!> \file GFS_GWD_generic.f
!! This file contains the CCPP-compliant orographic gravity wave
!! drag pre interstitial codes.

module GFS_GWD_generic_pre

contains

!> \section arg_table_GFS_GWD_generic_pre_init Argument Table
!!
      subroutine GFS_GWD_generic_pre_init()
      end subroutine GFS_GWD_generic_pre_init

!! \section arg_table_GFS_GWD_generic_pre_run Argument Table
!! | local_name     | standard_name                                                           | long_name                                                                                | units   | rank | type      | kind      | intent | optional |
!! |----------------|-------------------------------------------------------------------------|------------------------------------------------------------------------------------------|---------|------|-----------|-----------|--------|----------|
!! | im             | horizontal_loop_extent                                                  | horizontal dimension                                                                     | count   |    0 | integer   |           | in     | F        |
!! | levs           | vertical_dimension                                                      | vertical layer dimension                                                                 | count   |    0 | integer   |           | in     | F        |
!! | nmtvr          | number_of_statistical_measures_of_subgrid_orography                     | number of statistical measures of subgrid orography                                      | count   |    0 | integer   |           | in     | F        |
!! | mntvar         | statistical_measures_of_subgrid_orography                               | array of statistical measures of subgrid orography                                       | various |    2 | real      | kind_phys | in     | F        |
!! | hprime         | standard_deviation_of_subgrid_orography                                 | standard deviation of subgrid orography                                                  | m       |    1 | real      | kind_phys | out    | F        |
!! | oc             | convexity_of_subgrid_orography                                          | convexity of subgrid orography                                                           | none    |    1 | real      | kind_phys | out    | F        |
!! | oa4            | asymmetry_of_subgrid_orography                                          | asymmetry of subgrid orography                                                           | none    |    2 | real      | kind_phys | out    | F        |
!! | clx            | fraction_of_grid_box_with_subgrid_orography_higher_than_critical_height | horizontal fraction of grid box covered by subgrid orography higher than critical height | frac    |    2 | real      | kind_phys | out    | F        |
!! | theta          | angle_from_east_of_maximum_subgrid_orographic_variations                | angle with_respect to east of maximum subgrid orographic variations                      | degrees |    1 | real      | kind_phys | out    | F        |
!! | sigma          | slope_of_subgrid_orography                                              | slope of subgrid orography                                                               | none    |    1 | real      | kind_phys | out    | F        |
!! | gamma          | anisotropy_of_subgrid_orography                                         | anisotropy of subgrid orography                                                          | none    |    1 | real      | kind_phys | out    | F        |
!! | elvmax         | maximum_subgrid_orography                                               | maximum of subgrid orography                                                             | m       |    1 | real      | kind_phys | out    | F        |
!! | lssav          | flag_diagnostics                                                        | logical flag for storing diagnostics                                                     | flag    |    0 | logical   |           | in     | F        |
!! | ldiag3d        | flag_diagnostics_3D                                                     | flag for 3d diagnostic fields                                                            | flag    |    0 | logical   |           | in     | F        |
!! | dtdt           | tendency_of_air_temperature_due_to_model_physics                        | updated tendency of the temperature                                                      | K s-1   |    2 | real      | kind_phys | in     | F        |
!! | dt3dt          | cumulative_change_in_temperature_due_to_orographic_gravity_wave_drag    | cumulative change in temperature due to orographic gravity wave drag                     | K       |    2 | real      | kind_phys | inout  | F        |
!! | dtf            | time_step_for_dynamics                                                  | dynamics timestep                                                                        | s       |    0 | real      | kind_phys | in     | F        |
!! | errmsg         | ccpp_error_message                                                      | error message for error handling in CCPP                                                 | none    |    0 | character | len=*     | out    | F        |
!! | errflg         | ccpp_error_flag                                                         | error flag for error handling in CCPP                                                    | flag    |    0 | integer   |           | out    | F        |
!!
!!  \section general General Algorithm
!!  \section detailed Detailed Algorithm
!!  @{
      subroutine GFS_GWD_generic_pre_run(                               &
     &           im, levs, nmtvr, mntvar,                               &
     &           hprime, oc, oa4, clx, theta,                           &
     &           sigma, gamma, elvmax, lssav, ldiag3d,                  &
     &           dtdt, dt3dt, dtf, errmsg, errflg)

      use machine, only : kind_phys
      implicit none

      integer, intent(in) :: im, levs, nmtvr
      real(kind=kind_phys), intent(in) :: mntvar(im,nmtvr)

      real(kind=kind_phys), intent(out) ::                              &
     &  hprime(im), oc(im), oa4(im,4), clx(im,4),                       &
     &  theta(im), sigma(im), gamma(im), elvmax(im)

      logical, intent(in) :: lssav, ldiag3d
      real(kind=kind_phys), intent(in) :: dtdt(im,levs)
      ! dt3dt only allocated only if ldiag3d is .true.
      real(kind=kind_phys), intent(inout) :: dt3dt(:,:)
      real(kind=kind_phys), intent(in) :: dtf

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer :: i, k

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (nmtvr == 14) then  ! current operational - as of 2014
        hprime(:) = mntvar(:,1)
        oc(:)     = mntvar(:,2)
        oa4(:,1)  = mntvar(:,3)
        oa4(:,2)  = mntvar(:,4)
        oa4(:,3)  = mntvar(:,5)
        oa4(:,4)  = mntvar(:,6)
        clx(:,1)  = mntvar(:,7)
        clx(:,2)  = mntvar(:,8)
        clx(:,3)  = mntvar(:,9)
        clx(:,4)  = mntvar(:,10)
        theta(:)  = mntvar(:,11)
        gamma(:)  = mntvar(:,12)
        sigma(:)  = mntvar(:,13)
        elvmax(:) = mntvar(:,14)
      elseif (nmtvr == 10) then
        hprime(:) = mntvar(:,1)
        oc(:)     = mntvar(:,2)
        oa4(:,1)  = mntvar(:,3)
        oa4(:,2)  = mntvar(:,4)
        oa4(:,3)  = mntvar(:,5)
        oa4(:,4)  = mntvar(:,6)
        clx(:,1)  = mntvar(:,7)
        clx(:,2)  = mntvar(:,8)
        clx(:,3)  = mntvar(:,9)
        clx(:,4)  = mntvar(:,10)
      elseif (nmtvr == 6) then
        hprime(:) = mntvar(:,1)
        oc(:)     = mntvar(:,2)
        oa4(:,1)  = mntvar(:,3)
        oa4(:,2)  = mntvar(:,4)
        oa4(:,3)  = mntvar(:,5)
        oa4(:,4)  = mntvar(:,6)
        clx(:,1)  = 0.0
        clx(:,2)  = 0.0
        clx(:,3)  = 0.0
        clx(:,4)  = 0.0
      else
        hprime = 0
        oc     = 0
        oa4    = 0
        clx    = 0
        theta  = 0
        gamma  = 0
        sigma  = 0
        elvmax = 0
      endif   ! end if_nmtvr

      if (lssav) then
        if (ldiag3d) then
          do k=1,levs
            do i=1,im
              dt3dt(i,k) = dt3dt(i,k) - dtdt(i,k)*dtf
            enddo
          enddo
        endif
      endif

      end subroutine GFS_GWD_generic_pre_run
!> @}

! \ingroup GFS_ogwd
! \brief Brief description of the subroutine
!
!> \section arg_table_GFS_GWD_generic_pre_finalize Argument Table
!!
      subroutine GFS_GWD_generic_pre_finalize()
      end subroutine GFS_GWD_generic_pre_finalize

end module GFS_GWD_generic_pre
