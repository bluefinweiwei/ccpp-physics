!>\file rrtmg_sw_post
!! This file contains
      module rrtmg_sw_post
      contains

!>\defgroup rrtmg_sw_post GFS RRTMG scheme post
!! @{
!> \section arg_table_rrtmg_sw_post_init Argument Table
!!
      subroutine rrtmg_sw_post_init ()
      end subroutine rrtmg_sw_post_init
! PGI compiler does not accept lines longer than 264 characters, remove during pre-processing
#ifndef __PGI
!> \section arg_table_rrtmg_sw_post_run Argument Table
!! | local_name     | standard_name                                                                                  | long_name                                                                    | units    | rank |  type                 |   kind    | intent | optional |
!! |----------------|------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------|----------|------|-----------------------|-----------|--------|----------|
!! | Model          | GFS_control_type_instance                                                                      | Fortran DDT containing FV3-GFS model control parameters                      | DDT      |    0 | GFS_control_type      |           | in     | F        |
!! | Grid           | GFS_grid_type_instance                                                                         | Fortran DDT containing FV3-GFS grid and interpolation related data           | DDT      |    0 | GFS_grid_type         |           | in     | F        |
!! | Diag           | GFS_diag_type_instance                                                                         | Fortran DDT containing FV3-GFS diagnotics data                               | DDT      |    0 | GFS_diag_type         |           | inout  | F        |
!! | Radtend        | GFS_radtend_type_instance                                                                      | Fortran DDT containing FV3-GFS fields targetted for diagnostic output        | DDT      |    0 | GFS_radtend_type      |           | inout  | F        |
!! | Coupling       | GFS_coupling_type_instance                                                                     | Fortran DDT containing FV3-GFS fields to/from coupling with other components | DDT      |    0 | GFS_coupling_type     |           | inout  | F        |
!! | im             | horizontal_loop_extent                                                                         | horizontal loop extent                                                       | count    |    0 | integer               |           | in     | F        |
!! | ltp            | extra_top_layer                                                                                | extra top layers                                                             | none     |    0 | integer               |           | in     | F        |
!! | nday           | daytime_points_dimension                                                                       | daytime points dimension                                                     | count    |    0 | integer               |           | in     | F        |
!! | lm             | number_of_vertical_layers_for_radiation_calculations                                                        | number of vertical layers for radiation calculation                          | count    |    0 | integer               |           | in     | F        |
!! | kd             | vertical_index_difference_between_inout_and_local                                              | vertical index difference between in/out and local                           | index    |    0 | integer               |           | in     | F        |
!! | htswc          | tendency_of_air_temperature_due_to_shortwave_heating_on_radiation_time_step                    | total sky heating rate due to shortwave radiation                            | K s-1    |    2 | real                  | kind_phys | in     | F        |
!! | htsw0          | tendency_of_air_temperature_due_to_shortwave_heating_assuming_clear_sky_on_radiation_time_step | clear sky heating rates due to shortwave radiation                           | K s-1    |    2 | real                  | kind_phys | in     | F        |
!! | sfcalb1        | surface_albedo_due_to_near_IR_direct                                                           | surface albedo due to near IR direct beam                                    | frac     |    1 | real                  | kind_phys | in     | F        |
!! | sfcalb2        | surface_albedo_due_to_near_IR_diffused                                                         | surface albedo due to near IR diffused beam                                  | frac     |    1 | real                  | kind_phys | in     | F        |
!! | sfcalb3        | surface_albedo_due_to_UV_and_VIS_direct                                                        | surface albedo due to UV+VIS direct beam                                     | frac     |    1 | real                  | kind_phys | in     | F        |
!! | sfcalb4        | surface_albedo_due_to_UV_and_VIS_diffused                                                      | surface albedo due to UV+VIS diffused beam                                   | frac     |    1 | real                  | kind_phys | in     | F        |
!! | scmpsw         | components_of_surface_downward_shortwave_fluxes                                                | derived type for special components of surface downward shortwave fluxes     | W m-2    |    1 | cmpfsw_type           |           | inout  | F        |
!! | errmsg         | ccpp_error_message                                                                             | error message for error handling in CCPP                                     | none     |    0 | character             | len=*     | out    | F        |
!! | errflg         | ccpp_error_flag                                                                                | error flag for error handling in CCPP                                        | flag     |    0 | integer               |           | out    | F        |
!!
#endif
      subroutine rrtmg_sw_post_run (Model, Grid, Diag, Radtend, Coupling,  &
                 im, ltp, nday, lm, kd, htswc, htsw0,                      &
                 sfcalb1, sfcalb2, sfcalb3, sfcalb4, scmpsw, errmsg, errflg)

      use machine,                   only: kind_phys
      use module_radsw_parameters,   only: topfsw_type, sfcfsw_type,   &
                                           cmpfsw_type
      use GFS_typedefs,              only: GFS_coupling_type,          &
                                           GFS_control_type,           &
                                           GFS_grid_type,              &
                                           GFS_radtend_type,           &
                                           GFS_diag_type

      implicit none
      type(GFS_control_type),         intent(in)    :: Model
      type(GFS_coupling_type),        intent(inout) :: Coupling
      type(GFS_radtend_type),         intent(inout) :: Radtend
      type(GFS_grid_type),            intent(in)    :: Grid
      type(GFS_diag_type),            intent(inout) :: Diag
      integer, intent(in)                           :: im, lm, kd, nday, ltp
      type(cmpfsw_type), dimension(size(Grid%xlon,1)), intent(inout) :: scmpsw
      real(kind=kind_phys), dimension(Size(Grid%xlon,1), Model%levr+LTP), intent(in) ::  htswc, htsw0
      real(kind=kind_phys), dimension(size(Grid%xlon,1)), intent(in) :: sfcalb1, sfcalb2, sfcalb3, sfcalb4
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      ! Local variables
      integer :: i, k1, k

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (Model%lsswr) then
        if (nday > 0) then
          do k = 1, LM
            k1 = k + kd
            Radtend%htrsw(1:im,k) = htswc(1:im,k1)
          enddo
          ! We are assuming that radiative tendencies are from bottom to top 
          ! --- repopulate the points above levr i.e. LM
          if (lm < Model%levs) then
            do k = lm,Model%levs
              Radtend%htrsw (1:im,k) = Radtend%htrsw (1:im,LM)
            enddo
          endif

          if (Model%swhtr) then
            do k = 1, lm
               k1 = k + kd
               Radtend%swhc(1:im,k) = htsw0(1:im,k1)
             enddo
             ! --- repopulate the points above levr i.e. LM
             if (lm < Model%levs) then
               do k = lm,Model%levs
                 Radtend%swhc(1:im,k) = Radtend%swhc(1:im,LM)
               enddo
             endif
          endif

!  --- surface down and up spectral component fluxes
!>  - Save two spectral bands' surface downward and upward fluxes for
!!    output.

          do i=1,im
            Coupling%nirbmdi(i) = scmpsw(i)%nirbm
            Coupling%nirdfdi(i) = scmpsw(i)%nirdf
            Coupling%visbmdi(i) = scmpsw(i)%visbm
            Coupling%visdfdi(i) = scmpsw(i)%visdf

            Coupling%nirbmui(i) = scmpsw(i)%nirbm * sfcalb1(i)
            Coupling%nirdfui(i) = scmpsw(i)%nirdf * sfcalb2(i)
            Coupling%visbmui(i) = scmpsw(i)%visbm * sfcalb3(i)
            Coupling%visdfui(i) = scmpsw(i)%visdf * sfcalb4(i)
          enddo

        else                   ! if_nday_block

          Radtend%htrsw(:,:) = 0.0

          Radtend%sfcfsw = sfcfsw_type( 0.0, 0.0, 0.0, 0.0 )
          Diag%topfsw    = topfsw_type( 0.0, 0.0, 0.0 )
          scmpsw         = cmpfsw_type( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 )

          do i=1,im
            Coupling%nirbmdi(i) = 0.0
            Coupling%nirdfdi(i) = 0.0
            Coupling%visbmdi(i) = 0.0
            Coupling%visdfdi(i) = 0.0

            Coupling%nirbmui(i) = 0.0
            Coupling%nirdfui(i) = 0.0
            Coupling%visbmui(i) = 0.0
            Coupling%visdfui(i) = 0.0
          enddo

          if (Model%swhtr) then
            Radtend%swhc(:,:) = 0
          endif

        endif                  ! end_if_nday

! --- radiation fluxes for other physics processes
        do i=1,im
          Coupling%sfcnsw(i) = Radtend%sfcfsw(i)%dnfxc - Radtend%sfcfsw(i)%upfxc
          Coupling%sfcdsw(i) = Radtend%sfcfsw(i)%dnfxc
        enddo

      endif                                ! end_if_lsswr

      end subroutine rrtmg_sw_post_run
 
!> \section arg_table_rrtmg_sw_post_finalize Argument Table
!!
      subroutine rrtmg_sw_post_finalize ()
      end subroutine rrtmg_sw_post_finalize
!! @}
      end module rrtmg_sw_post
