!> \file cu_gf_driver_pre.F90
!!  Contains code related to GF convective schemes to be used within the GFS physics suite.

!> This module contains code related to GF convective schemes to be used within the GFS physics suite.
module cu_gf_driver_pre

   implicit none

   private

   public :: cu_gf_driver_pre_run

   contains

!>\ingroup cu_gf_group
!> \section arg_table_cu_gf_driver_pre_run Argument Table
!! \htmlinclude cu_gf_driver_pre_run.html
!!
   subroutine cu_gf_driver_pre_run (cactiv, cactiv_m, conv_act, conv_act_m,         &
                                    ntsmoke, ntdust, ntcoarsepm, chem3d, gq0, errmsg, errflg)

      use machine, only: kind_phys

      implicit none

      integer,          intent(out) :: cactiv(:)
      integer,          intent(out) :: cactiv_m(:)
      integer,          intent(in)  :: ntsmoke, ntdust, ntcoarsepm
!$acc declare copyout(cactiv,cactiv_m)
      real(kind_phys),  intent(in)  :: conv_act(:)
      real(kind_phys),  intent(in)  :: conv_act_m(:)
      real(kind_phys),  intent(inout), optional :: chem3d(:,:,:)
      real(kind_phys),  intent(inout) :: gq0(:,:,:)
!$acc declare copyin(conv_act,conv_act_m) copy(chem3d,gq0)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0


!$acc kernels
      cactiv(:)=nint(conv_act(:))
      cactiv_m(:)=nint(conv_act_m(:))

      if (present(chem3d)) then
       chem3d(:,:,1) = gq0(:,:,ntsmoke)
       chem3d(:,:,2) = gq0(:,:,ntdust)
       chem3d(:,:,3) = gq0(:,:,ntcoarsepm)
      endif
!$acc end kernels

   end subroutine cu_gf_driver_pre_run

end module cu_gf_driver_pre
