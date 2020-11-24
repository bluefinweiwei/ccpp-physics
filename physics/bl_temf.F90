!> \file bl_temf.F90
!! This file contains the CCPP-compliant Total Energy - Mass Flux (TEMF) PBL scheme.
!! About this scheme:
!! Initial implementation 2010 by Wayne Angevine, CIRES/NOAA ESRL.
!! References:
!!    Angevine et al., 2010, MWR
!!    Angevine, 2005, JAM
!!    Mauritsen et al., 2007, JAS
!----------------------------------------------------------------------
    module bl_temf

        use machine, only : kind_phys

        use physcons, only :   rcp    => con_rocp, &
                               cpv    => con_cvap, &
                               omega  => con_omega

        implicit none

        public bl_temf_init, bl_temf_run, bl_temf_finalize

        private

        real, parameter     :: hmax = 4000.0  ! Max hd,hct WA 11/20/09
        real, parameter     :: Ceps = 0.070
        real, parameter     :: Cgamma = Ceps
        real, parameter     :: Cphi = Ceps
        real, parameter     :: ftau0 = 0.17
        real, parameter     :: fth0 = 0.145
        real, parameter     :: Cf = 0.185
        real, parameter     :: CN = 2.0
        real, parameter     :: lasymp = 200.0 ! Asymptotic length scale WA 11/20/09
        real, parameter     :: Cc = 3.0       ! Prefactor for convective length scale
        real, parameter     :: PrT0 = Cphi/Ceps * ftau0**2 / 2. / fth0**2
        real, parameter     :: visc_temf = 1.57e-4   ! WA TEST bigger minimum K
        real, parameter     :: conduc_temf = 1.57e-4 / 0.733
        real, parameter     :: Cw = 0.5       ! Prefactor for surface wUPD
        real, parameter     :: Cdelt = 0.006  ! Prefactor for detrainment rate
        real, parameter     :: CM = 0.03      ! Proportionality constant for subcloud MF
        real, parameter     :: TEmin = 1e-3   ! minimum TE value

        contains

! -----------------------------------------------------------------------
! CCPP entry points for TEMF
! -----------------------------------------------------------------------
    subroutine bl_temf_init ()

    end subroutine bl_temf_init

    subroutine bl_temf_finalize()

    end subroutine bl_temf_finalize

!> \defgroup TEMF bl_temf Main
!! \brief This subroutine contains all of the logic for the
!! TEMF PBL scheme.
!!
!> \section arg_table_bl_temf_run Argument Table
!! \htmlinclude bl_temf_run.html
!!
! -------------------------------------------------------------------------------------------
    subroutine bl_temf_run (im,km,dt,nvdiff,ntqvx,ntcwx,&
                            r_d,r_v,cp,g,ep1,ep2,xlv,karman,&
                            tsk,hfx,qfx,zorl,zsrf,sinlat,&
                            exner,ugrs,vgrs,tgrs,qgrs,p2d,phil,prsi,&
                            kpbl1d,hpbl,t2,u10,v10,&
                            rubltenx, rvbltenx, rtbltenx, rqbltenx,     &
                            lssav, ldiag3d, qdiag3d,                                    &
                            flag_for_pbl_generic_tend, du3dt_PBL, dv3dt_PBL,             &
                            dt3dt_PBL, dq3dt_PBL, &
                            errmsg,errflg)

    implicit none

    !WL* or put this in CCPP physics namelist?
    logical, parameter  :: MFopt = .true.  ! Use mass flux or not
    !*WL

    ! input variables
    integer,                                        intent(in)      ::  im,km,nvdiff,ntqvx,ntcwx
    logical,                                        intent(in)      ::  lssav, ldiag3d, qdiag3d, flag_for_pbl_generic_tend
    real(kind=kind_phys),                           intent(in)      ::  dt,ep1,ep2,r_d,r_v,g,cp,xlv,karman
    real(kind=kind_phys),   dimension(im),          intent(in)      ::  tsk,hfx,qfx,zorl,zsrf,sinlat
    real(kind=kind_phys),   dimension(im,km),       intent(in)      ::  phil,p2d,exner,ugrs,vgrs,tgrs
    real(kind=kind_phys),   dimension(im,km+1),     intent(in)      ::  prsi
    real(kind=kind_phys),   dimension(im,km,nvdiff),intent(in)      ::  qgrs


    ! input/output variables
    real(kind=kind_phys),   dimension(im),          intent(inout)   ::  hpbl,u10,v10,t2
    real(kind=kind_phys),   dimension(im,km),       intent(inout)   ::  rubltenx,rvbltenx,rtbltenx
    real(kind=kind_phys),   dimension(im,km,nvdiff),intent(inout)   ::  rqbltenx

    real(kind=kind_phys),   dimension(im,km),       intent(inout)   ::  du3dt_PBL, dv3dt_PBL, dt3dt_PBL, dq3dt_PBL

    ! output variables
    integer,                dimension(im),          intent(out)     ::  kpbl1d

    ! error messages
    character(len=*),                               intent(out)     ::  errmsg
    integer,                                        intent(out)     ::  errflg

    ! local variables
    integer i,k,l,kt
    real    tvcon,Cepsmf,red_fact,sigq2,rst  ! red_fact for reducing MF components
    logical is_convective
    integer,               dimension(im)      :: hdidx,hctidx,lclidx,hmax_idx,h0idx,htidx,tval
    real(kind=kind_phys),  dimension(im)      :: fCor,znt,z0t,sfcTHVF,h0,wstr,hct,hd,lcl,ht,ang,wm,&
                                                 convection_TKE_surface_src,sfcFTE,ust
    real(kind=kind_phys),  dimension(im)      :: cfm_temfx,hd_temfx,lcl_temfx,hct_temfx
    real(kind=kind_phys),  dimension(im,km)   :: theta,thetal,thetav
    real(kind=kind_phys),  dimension(im,km)   :: qt,rv,rl,rt
    real(kind=kind_phys),  dimension(im,km)   :: chi_poisson,gam
    real(kind=kind_phys),  dimension(im,km)   :: u_temf,v_temf
    real(kind=kind_phys),  dimension(im,km)   :: z2d,zm,dzm,zt,dzt
    real(kind=kind_phys),  dimension(im,km)   :: dthdz,dqtdz,dudz,dvdz
    real(kind=kind_phys),  dimension(im,km)   :: rho,lepsmin,N2,S,Ri,beta,deltmf,ratio,ftau,fth
    real(kind=kind_phys),  dimension(im,km)   :: TE2,TKE
    real(kind=kind_phys),  dimension(im,km)   :: ustrtilde,linv,leps,kmm,kh,&
                                                 lconv,kh_conv,kmm_conv,Fz,epsmf,B,Bmoist
    real(kind=kind_phys),  dimension(im,km)   :: rUPD,wupd_dry,UUPD,VUPD,thetavUPD,&
                                                 thetavUPDmoist,TEUPD,qlUPD,TUPD,rstUPD,rlUPD,qstUPD,dthUPDdz,&
                                                 dwUPDdz,dqtup_temfxdz,dUUPDdz,dVUPDdz,dTEUPDdz,dwUPDmoistdz
    real(kind=kind_phys),  dimension(im,km)   :: MFCth,MFCq,MFCu,MFCv,MFCql,MFCthv,MFCTE,dMdz
    real(kind=kind_phys),  dimension(im,km)   :: au,sigq,qst,satdef,QFK,uwk,vwk,beta1,alpha2,beta2,THVF,&
                                                 buoy_src,srcs
    real(kind=kind_phys),  dimension(im,km)   :: te_temfx,qf_temfx,kh_temfx,kmm_temfx,&
                                                 thup_temfx,qtup_temfx,&
                                                 qlup_temfx,wupd_temfx,mf_temfx,&
                                                 cf3d_temfx,shf_temfx,&
                                                 uw_temfx,vw_temfx
    real(kind=kind_phys),  dimension(im,km)   :: u_new,v_new,theta_new,qvx_new,qcx_new

    !!!!!!!!!!!!!!!!!!!!!!! Actual Code !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !--- initialize CCPP error handling variables
    errmsg = ''
    errflg = 0
!----------------------------------------------------------------------
! Below is from the original model_bl_temf.F in WRF:
! Grid staggering:  Matlab version has mass and turbulence levels.
! WRF has full levels (with w) and half levels (u,v,theta,q*).  Both
! sets of levels use the same indices (1:km).  See pbl_driver or
! WRF Physics doc for (a few) details.
! So *mass levels correspond to half levels.*
! WRF full levels are ignored, we define our own turbulence levels
! in order to put the first one below the first half level.
! Another difference is that
! the Matlab version (and the Mauritsen et al. paper) consider the
! first mass level to be at z0 (effectively the surface).  WRF considers
! the first half level to be above the effective surface.  The first half
! level, at k=1, has nonzero values of u,v for example.  Here we convert
! all incoming variables to internal ones with the correct indexing
! in order to make the code consistent with the Matlab version.  We
! already had to do this for thetal and qt anyway, so the only additional
! overhead is for u and v.
! I use suffixes m for mass and t for turbulence as in Matlab for things
! like indices.
! Note that zsrf is the terrain height ASL, from Registry variable ht.
! Translations (Matlab to WRF):
!     dzt -> calculated below
!     dzm -> not supplied, calculated below
!     k -> karman
!     z0 -> znt
!     z0t -> not in WRF, calculated below
!     zt -> calculated below
!     zm -> (z2d - zsrf) but NOTE zm(1) is now z0 (znt) and zm(2) is
!                                 z2d(1) - zsrf
!
! WA I take the temperature at z0 to be
! TSK.  This isn't exactly robust.
! WA 2/16/11 removed calculation of flhc, flqc which are not needed here.
! These should be removed from the calling sequence someday.
!
! Other notes:
! - I have often used 1 instead of kts below, because the scheme demands
!   to know where the surface is.  It won't work if kts .NE. 1.
!


   do i = 1,im      ! Main loop

      do k=1,km
        theta(i,k) = tgrs(i,k)/exner(i,k)
        z2d(i,k) = phil(i,k)/g  !WL* geopotential heights at model full levels (z2d)
        te_temfx(i,k) = TEmin   !WL* initial value for te_temfx temfinit
      end do

      ! Get incoming surface theta from TSK (WA for now)
      thetal(i,1) = tsk(i) / exner(i,1)  ! WA really should use Exner func. at z0
      qt(i,1) = qgrs(i,1,ntqvx)
      rv(i,1) = qt(i,1) / (1.-qt(i,1))   ! Water vapor
      rl(i,1) = 0.
      rt(i,1) = rv(i,1) + rl(i,1)        ! Total water (without ice)
      chi_poisson(i,1) = rcp * (1.+rv(i,1)/ep2) / (1.+rv(i,1)*cpv/cp)
      gam(i,1) = rv(i,1) * r_v / (cp + rv(i,1)*cpv)
      thetav(i,1) = thetal(i,1) * (1. + 0.608*qt(i,1) - qgrs(i,1,ntcwx))  ! WA 4/6/10 allow environment liquid
      ! Convert incoming theta to thetal and qv,qc to qt
      ! NOTE this is where the indexing gets changed from WRF to TEMF basis
      do k = 2,km
         ! Convert specific humidities to mixing ratios
         rv(i,k) = qgrs(i,k-1,ntqvx) / (1.-qgrs(i,k-1,ntqvx))   ! Water vapor
         rl(i,k) = qgrs(i,k-1,ntcwx) / (1.-qgrs(i,k-1,ntcwx))   ! Liquid water
         rt(i,k) = rv(i,k) + rl(i,k)        ! Total water (without ice)
         chi_poisson(i,k) = rcp * (1.+rv(i,k)/ep2) / (1.+rv(i,k)*cpv/cp)
         gam(i,k) = rt(i,k) * r_v / (cp + rt(i,k)*cpv)
         thetal(i,k) = (tgrs(i,k-1) / exner(i,k-1)) * &
            ((ep2+rv(i,k))/(ep2+rt(i,k)))**chi_poisson(i,k) * &
            (rv(i,k)/rt(i,k))**(-gam(i,k)) * exp( -xlv*rl(i,k) / &
            ((cp + rt(i,k)*cpv) * tgrs(i,k)))
         qt(i,k) = qgrs(i,k-1,ntqvx) + qgrs(i,k-1,ntcwx)
         thetav(i,k) = thetal(i,k) * (1. + 0.608*qt(i,k) - qgrs(i,k-1,ntcwx))  ! WA 4/6/10 allow environment liquid
      end do

      ! Convert incoming u,v to internal u_temf, v_temf
      ! NOTE this is where the indexing gets changed from WRF to TEMF basis
      u_temf(i,1) = 0.   ! zero winds at z0
      v_temf(i,1) = 0.
      do k = 2,km
         u_temf(i,k) = ugrs(i,k-1)
         v_temf(i,k) = vgrs(i,k-1)
      end do


      !WL* convert roughness from cm to m
      znt(i) = zorl(i)/100.
      !*WL (don't see any reason of converting znt to z0t but keep it for now to be consistent with the original code)
      z0t(i) = znt(i)

      ! Get delta height at half (mass) levels (zm)
      zm(i,1) = znt(i) !WL* roughness length as 1st mass level
      dzt(i,1) = z2d(i,1) - zsrf(i) - zm(i,1)
      ! Get height and delta at turbulence levels (zt)
      zt(i,1) = (z2d(i,1) - zsrf(i) - znt(i)) / 2.
      do kt = 2,km
         zm(i,kt) = z2d(i,kt-1) - zsrf(i) ! Convert indexing from WRF to TEMF
         zt(i,kt) = (zm(i,kt) + z2d(i,kt) - zsrf(i)) / 2.
         dzm(i,kt) = zt(i,kt) - zt(i,kt-1)
         dzt(i,kt) = z2d(i,kt+1) - z2d(i,kt)
      end do
      dzm(i,1) = dzm(i,2)  ! WA why?
      dzt(i,km) = dzt(i,km-1)    ! WA 12/23/09

      ! Gradients at first level
      dthdz(i,1) = (thetal(i,2)-thetal(i,1)) / (zt(i,1) * log10(zm(i,2)/z0t(i)))
      dqtdz(i,1) = (qt(i,2)-qt(i,1)) / (zt(i,1) * log10(zm(i,2)/z0t(i)))
      dudz(i,1) = (u_temf(i,2)-u_temf(i,1)) / (zt(i,1) * log10(zm(i,2)/znt(i)))
      dvdz(i,1) = (v_temf(i,2)-v_temf(i,1)) / (zt(i,1) * log10(zm(i,2)/znt(i)))

      ! Surface thetaV flux from Stull p.147

      !WL* calculate rho and get qf_temfx(i,1) from upward latent heat flux qfx
      tvcon    = (1.+ep1 * qgrs(i,1,ntqvx))
      rho(i,1) = p2d(i,1)/(r_d * tgrs(i,1) * tvcon)
      qf_temfx(i,1) = qfx(i)/rho(i,1)
      sfcTHVF(i) = hfx(i)/(rho(i,1)*cp) * (1.+0.608*(qgrs(i,1,ntqvx)+qgrs(i,1,ntcwx))) + 0.608*thetav(i,1)*qf_temfx(i,1)
      !*WL

      ! WA use hd_temf to calculate w* instead of finding h0 here????
      ! Watch initialization!
      h0idx(i) = 1
      h0(i) = zm(i,1)

      lepsmin(i,1) = 0.

      ! WA 2/11/13 find index just above hmax for use below
      hmax_idx(i) = km-1

      do k = 2,km-1
         lepsmin(i,k) = 0.

         ! Mean gradients
         dthdz(i,k) = (thetal(i,k+1) - thetal(i,k)) / dzt(i,k)
         dqtdz(i,k) = (qt(i,k+1) - qt(i,k)) / dzt(i,k)
         dudz(i,k) = (u_temf(i,k+1) - u_temf(i,k)) / dzt(i,k)
         dvdz(i,k) = (v_temf(i,k+1) - v_temf(i,k)) / dzt(i,k)

         ! Find h0 (should eventually be interpolated for smoothness)
         if (thetav(i,k) > thetav(i,1) .AND. h0idx(i) .EQ. 1) then
         ! WA 9/28/11 limit h0 as for hd and hct
            if (zm(i,k) < hmax) then
               h0idx(i) = k
               h0(i) = zm(i,k)
            else
               h0idx(i) = k
               h0(i) = hmax
            end if
         end if
         ! WA 2/11/13 find index just above hmax for use below
         if (zm(i,k) > hmax) then
            hmax_idx(i) = min(hmax_idx(i),k)
         end if
      end do

      ! Gradients at top level

      dthdz(i,km) = dthdz(i,km-1)
      dqtdz(i,km) = dqtdz(i,km-1)
      dudz(i,km) = dudz(i,km-1)
      dvdz(i,km) = dvdz(i,km-1)

      if ( hfx(i) > 0.) then
         ! wstr(i) = (g * h0(i) / thetav(i,2) * shf_temfx(i,1) ) ** (1./3.)
         wstr(i) = (g * h0(i) / thetav(i,2) * hfx(i)/(rho(i,1)*cp) ) ** (1./3.)
      else
         wstr(i) = 0.
      end if


      ! Set flag convective or not for use below
      is_convective = wstr(i) > 0. .AND. MFopt .AND. dthdz(i,1)<0. .AND. dthdz(i,2)<0.  ! WA 12/16/09 require two levels of negative (unstable) gradient

      !WL*
      print*,'WL***',is_convective
      !*WL


      ! Find stability parameters and length scale (on turbulence levels)
      do kt = 1,km-1
         N2(i,kt) = 2. * g / (thetav(i,kt) + thetav(i,kt+1))*dthdz(i,kt)
         S(i,kt) = sqrt(dudz(i,kt)**2. + dvdz(i,kt)**2.)
         !WL* 
         if (S(i,kt) /= 0.) Ri(i,kt) = N2(i,kt) / S(i,kt)**2.
         !*WL
         if (S(i,kt) < 1e-15) then
            if (N2(i,kt) >= 0) then
               Ri(i,kt) = 10.
            else
               Ri(i,kt) = -1.
            end if
         end if
         beta(i,kt) = 2. * g / (thetav(i,kt)+thetav(i,kt+1))
         if (Ri(i,kt) > 0) then
            ratio(i,kt) = Ri(i,kt)/(Cphi**2.*ftau0**2./(2.*Ceps**2.*fth0**2.)+3.*Ri(i,kt))
            ftau(i,kt) = ftau0 * ((3./4.) / (1.+4.*Ri(i,kt)) + 1./4.)
            fth(i,kt) = fth0 / (1.+4.*Ri(i,kt))
            TE2(i,kt) = 2. * te_temfx(i,kt) * ratio(i,kt) * N2(i,kt) / beta(i,kt)**2.
         else
            ratio(i,kt) = Ri(i,kt)/(Cphi**2.*ftau0**2./(-2.*Ceps**2.*fth0**2.)+2.*Ri(i,kt))
            ftau(i,kt) = ftau0
            fth(i,kt) = fth0
            TE2(i,kt) = 0.
         end if
         TKE(i,kt) = te_temfx(i,kt) * (1. - ratio(i,kt))
         ustrtilde(i,kt) = sqrt(ftau(i,kt) * TKE(i,kt))
         if (N2(i,kt) > 0.) then
            !WL* Coriolis parameter
            fCor(i) = 2 * omega * sinlat(i)
            !*WL
            linv(i,kt) = 1./karman / zt(i,kt) + abs(fCor(i)) / &
               (Cf*ustrtilde(i,kt)) + &
               sqrt(N2(i,kt))/(CN*ustrtilde(i,kt)) + 1./lasymp
         else
            linv(i,kt) = 1./karman / zt(i,kt) + abs(fCor(i)) / &
               (Cf*ustrtilde(i,kt)) + 1./lasymp
         end if
         leps(i,kt) = 1./linv(i,kt)
         leps(i,kt) = max(leps(i,kt),lepsmin(i,kt))
      end do
      S(i,km) = 0.0
      N2(i,km) = 0.0
      TKE(i,km) = 0.0
      linv(i,km) = linv(i,km-1)
      leps(i,km) = leps(i,km-1)


      ! Find diffusion coefficients
      ! First use basic formulae for stable and neutral cases,
      ! then for convective conditions, and finally choose the larger
      ! WA 12/23/09 use convective form up to hd/2 always
      ! WA 12/28/09 after restructuring, this block is above MF block,
      ! so hd is not yet available for this timestep, must use h0,
      ! or use hd from previous timestep but be careful about initialization.
      do kt = 1,km-1    ! WA 12/22/09
         ! WA 4/8/10 remove beta term to avoid negative and huge values
         ! of km due to very small denominator.  This is an interim fix
         ! until we find something better (more theoretically sound).
         ! km(i,kt) = TKE(i,kt)**1.5 * ftau(i,kt)**2. / (-beta(i,kt) * fth(i,kt) * sqrt(TE2(i,kt)) + Ceps * sqrt(TKE(i,kt)*te_temfx(i,kt)) / leps(i,kt))
         !WL* replace origional "km" with "kmm" to avoid conflicting with vertical dimension "km"
         kmm(i,kt) = TKE(i,kt)**1.5 * ftau(i,kt)**2. / (Ceps * sqrt(TKE(i,kt)*te_temfx(i,kt)) / leps(i,kt))
         kh(i,kt) = 2. * leps(i,kt) * fth(i,kt)**2. * TKE(i,kt) / sqrt(te_temfx(i,kt)) / Cphi
         if ( is_convective) then
            ! WA 2/20/14 trap rare "equality" of h0 and zt (only when h0 = hmax)
            if (kt <= h0idx(i) .AND. h0(i)-zt(i,kt) > 1e-15) then
               lconv(i,kt) = 1. / (1. / (karman*zt(i,kt)) + Cc / (karman * (h0(i) - zt(i,kt))))
            else
               lconv(i,kt) = 0.
            end if
            ! WA 12/15/09 use appropriate coeffs to match kh_conv and kh at neutral
            kh_conv(i,kt) = ftau0**2. / Ceps / PrT0 * sqrt(TKE(i,kt)) * lconv(i,kt)
            if (kh_conv(i,kt) < 0.) then
               kh_conv(i,kt) = 0.
            end if
            kmm_conv(i,kt) = PrT0 * kh_conv(i,kt)
            if (zt(i,kt) <= h0(i)/2.) then
               kmm(i,kt) = kmm_conv(i,kt)
               kh(i,kt) = kh_conv(i,kt)
            end if
            if (zt(i,kt) > h0(i)/2. .AND. kt <= h0idx(i)) then
               kmm(i,kt) = max(kmm(i,kt),kmm_conv(i,kt),visc_temf)
               kh(i,kt) = max(kh(i,kt),kh_conv(i,kt),conduc_temf)
            end if
         end if  ! is_convective
         kmm(i,kt) = max(kmm(i,kt),visc_temf)
         kh(i,kt) = max(kh(i,kt),conduc_temf)
         Fz(i,kt) = -kh(i,kt) * dthdz(i,kt)  ! Diffusive heat flux
      end do
      kmm(i,km) = kmm(i,km-1)
      kh(i,km) = kh(i,km-1)
      Fz(i,km) = 0.0


      !*** Mass flux block starts here ***

      if ( is_convective) then

         Cepsmf = 2. / max(200.,h0(i))
         Cepsmf = max(Cepsmf,0.002)
         do k = 1,km
            ! Calculate lateral entrainment fraction for subcloud layer
            ! epsilon and delta are defined on mass grid (half levels)
            epsmf(i,k) = Cepsmf
         end do

         ! Initialize updraft
         thup_temfx(i,1) = thetal(i,1)    ! No excess
         qtup_temfx(i,1) = qt(i,1)            ! No excess
         rUPD(i,1) = qtup_temfx(i,1) / (1. - qtup_temfx(i,1))
         wupd_temfx(i,1) = Cw * wstr(i)
         wupd_dry(i,1) = Cw * wstr(i)
         UUPD(i,1) = u_temf(i,1)
         VUPD(i,1) = v_temf(i,1)
         thetavUPD(i,1) = thup_temfx(i,1) * (1. + 0.608*qtup_temfx(i,1))  ! WA Assumes no liquid
         thetavUPDmoist(i,1) = thup_temfx(i,1) * (1. + 0.608*qtup_temfx(i,1))  ! WA Assumes no liquid
         TEUPD(i,1) = te_temfx(i,1) + g / thetav(i,1) * sfcTHVF(i)
         qlUPD(i,1) = qgrs(i,1,ntcwx)  ! WA allow environment liquid
         TUPD(i,1) = thup_temfx(i,1) * exner(i,1)
         rstUPD(i,1) = rsat(p2d(i,1),TUPD(i,1),ep2)
         rlUPD(i,1) = 0.

         ! Calculate updraft parameters counting up
         do k = 2,km
            ! WA 2/11/13 use hmax index to prevent oddness high up
            if ( k < hmax_idx(i)) then
               dthUPDdz(i,k-1) = -epsmf(i,k) * (thup_temfx(i,k-1) - thetal(i,k-1))
               thup_temfx(i,k) = thup_temfx(i,k-1) + dthUPDdz(i,k-1) * dzm(i,k-1)
               dqtup_temfxdz(i,k-1) = -epsmf(i,k) * (qtup_temfx(i,k-1) - qt(i,k-1))
               qtup_temfx(i,k) = qtup_temfx(i,k-1) + dqtup_temfxdz(i,k-1) * dzm(i,k-1)
               thetavUPD(i,k) = thup_temfx(i,k) * (1. + 0.608*qtup_temfx(i,k))  ! WA Assumes no liquid
               B(i,k-1) = g * (thetavUPD(i,k) - thetav(i,k)) / thetav(i,k)
               if ( wupd_dry(i,k-1) < 1e-15 ) then
                  wupd_dry(i,k) = 0.
               else
                  dwUPDdz(i,k-1) = -2. *epsmf(i,k)*wupd_dry(i,k-1) + 0.33*B(i,k-1)/wupd_dry(i,k-1)
                  wupd_dry(i,k) = wupd_dry(i,k-1) + dwUPDdz(i,k-1) * dzm(i,k-1)
               end if
               dUUPDdz(i,k-1) = -epsmf(i,k) * (UUPD(i,k-1) - u_temf(i,k-1))
               UUPD(i,k) = UUPD(i,k-1) + dUUPDdz(i,k-1) * dzm(i,k-1)
               dVUPDdz(i,k-1) = -epsmf(i,k) * (VUPD(i,k-1) - v_temf(i,k-1))
               VUPD(i,k) = VUPD(i,k-1) + dVUPDdz(i,k-1) * dzm(i,k-1)
               dTEUPDdz(i,k-1) = -epsmf(i,k) * (TEUPD(i,k-1) - te_temfx(i,k-1))
               TEUPD(i,k) = TEUPD(i,k-1) + dTEUPDdz(i,k-1) * dzm(i,k-1)
               ! Alternative updraft velocity based on moist thetav
               ! Need thetavUPDmoist, qlUPD
               rUPD(i,k) = qtup_temfx(i,k) / (1. - qtup_temfx(i,k))
               ! WA Updraft temperature assuming no liquid
               TUPD(i,k) = thup_temfx(i,k) * exner(i,k)
               ! Updraft saturation mixing ratio
               rstUPD(i,k) = rsat(p2d(i,k-1),TUPD(i,k),ep2)
               ! Correct to actual temperature (Sommeria & Deardorff 1977)
               beta1(i,k) = 0.622 * (xlv/(r_d*TUPD(i,k))) * (xlv/(cp*TUPD(i,k)))
               rstUPD(i,k) = rstUPD(i,k) * (1.0+beta1(i,k)*rUPD(i,k)) / (1.0+beta1(i,k)*rstUPD(i,k))
               qstUPD(i,k) = rstUPD(i,k) / (1. + rstUPD(i,k))
               if (rUPD(i,k) > rstUPD(i,k)) then
                  rlUPD(i,k) = rUPD(i,k) - rstUPD(i,k)
                  qlUPD(i,k) = rlUPD(i,k) / (1. + rlUPD(i,k))
                  thetavUPDmoist(i,k) = (thup_temfx(i,k) + ((xlv/cp)*qlUPD(i,k)/exner(i,k))) * &
                                        (1. + 0.608*qstUPD(i,k) - qlUPD(i,k))
               else
                  rlUPD(i,k) = 0.
                  qlUPD(i,k) = qgrs(i,k-1,ntcwx)!qcx(i,k-1)   ! WA 4/6/10 allow environment liquid
                  thetavUPDmoist(i,k) = thup_temfx(i,k) * (1. + 0.608*qtup_temfx(i,k))
               end if
               Bmoist(i,k-1) = g * (thetavUPDmoist(i,k) - thetav(i,k)) / thetav(i,k)
               if ( wupd_temfx(i,k-1) < 1e-15 ) then
                  wupd_temfx(i,k) = 0.
               else
                  dwUPDmoistdz(i,k-1) = -2. *epsmf(i,k)*wupd_temfx(i,k-1) + 0.33*Bmoist(i,k-1)/wupd_temfx(i,k-1)
                  wupd_temfx(i,k) = wupd_temfx(i,k-1) + dwUPDmoistdz(i,k-1) * dzm(i,k-1)
               end if
            else
               thup_temfx(i,k) = thetal(i,k)
               qtup_temfx(i,k) = qt(i,k)
               wupd_dry(i,k) = 0.
               UUPD(i,k) = u_temf(i,k)
               VUPD(i,k) = v_temf(i,k)
               TEUPD(i,k) = te_temfx(i,k)
               qlUPD(i,k) = qgrs(i,k-1,ntcwx)!qcx(i,k-1)
               wupd_temfx(i,k) = 0.
            end if
         end do

         ! Find hd based on wUPD
         if (wupd_dry(i,1) == 0.) then
            hdidx(i) = 1
         else
            hdidx(i) = km  ! In case wUPD <= 0 not found
            do k = 2,km
               ! if (wupd_dry(i,k) <= 0.) then
               if (wupd_dry(i,k) <= 0. .OR. zm(i,k) > hmax) then
                  hdidx(i) = k
                  !WL--> GO TO not allowed in CCPP
                  !goto 100   ! FORTRAN made me do it!
                  exit
                  !<--WL
               end if
            end do
         end if
         hd(i) = zm(i,hdidx(i))
         kpbl1d(i) = hdidx(i)
         hpbl(i) = hd(i)       ! WA 5/11/10 hpbl is height.  Should still be replaced by something that works whether convective or not.

         ! Find LCL, hct, and ht
         lclidx(i) = km   ! In case LCL not found
         do k = 1,km
            if ( k < hmax_idx(i) .AND. rUPD(i,k) > rstUPD(i,k)) then
               lclidx(i) = k
               !WL--> no GO TO allowed in CCPP
               !goto 200
                exit
               !<--WL
            end if
         end do
         lcl(i) = zm(i,lclidx(i))

         if (hd(i) > lcl(i)) then   ! Forced cloud (at least) occurs
            ! Find hct based on wUPDmoist
            if (wupd_temfx(i,1) == 0.) then
               hctidx(i) = 1
            else
               hctidx(i) = km  ! In case wUPD <= 0 not found
               do k = 2,km
                  if (wupd_temfx(i,k) <= 0. .OR. zm(i,k) > hmax) then
                     hctidx(i) = k
                     !WL--> no GO TO allowed in CCPP
                     !goto 300   ! FORTRAN made me do it!
                     exit
                     !<--WL
                  end if
               end do
            end if
            hct(i) = zm(i,hctidx(i))
            if (hctidx(i) <= hdidx(i)+1) then   ! No active cloud
               hct(i) = hd(i)
               hctidx(i) = hdidx(i)
            else
            end if
         else   ! No cloud
            hct(i) = hd(i)
            hctidx(i) = hdidx(i)
         end if
         ht(i) = max(hd(i),hct(i))
         htidx(i) = max(hdidx(i),hctidx(i))

         ! Now truncate updraft at ht with taper
         do k = 1,km
            if (zm(i,k) < 0.9*ht(i)) then  ! Below taper region
               tval(i) = 1
            else if (zm(i,k) >= 0.9*ht(i) .AND. zm(i,k) <= 1.0*ht(i)) then
               ! Within taper region
               tval(i) = 1. - ((zm(i,k) - 0.9*ht(i)) / (1.0*ht(i) - 0.9*ht(i)))
            else  ! Above taper region
               tval(i) = 0.
            end if
            thup_temfx(i,k) = tval(i) * thup_temfx(i,k) + (1-tval(i))*thetal(i,k)
            thetavUPD(i,k) = tval(i) * thetavUPD(i,k) + (1-tval(i))*thetav(i,k)
            qtup_temfx(i,k) = tval(i) * qtup_temfx(i,k) + (1-tval(i)) * qt(i,k)
            ! WA 6/21/13 was a subscript error when k=1
            if (k > 1) then
               qlUPD(i,k) = tval(i) * qlUPD(i,k) + (1-tval(i)) * qgrs(i,k-1,ntcwx)
            end if
            UUPD(i,k) = tval(i) * UUPD(i,k) + (1-tval(i)) * u_temf(i,k)
            VUPD(i,k) = tval(i) * VUPD(i,k) + (1-tval(i)) * v_temf(i,k)
            TEUPD(i,k) = tval(i) * TEUPD(i,k) + (1-tval(i)) * te_temfx(i,k)
            if (zm(i,k) > ht(i)) then  ! WA this is just for cleanliness
               wupd_temfx(i,k) = 0.
               dwUPDmoistdz(i,k) = 0.
               wupd_dry(i,k) = 0.
               dwUPDdz(i,k) = 0.
            end if
         end do

         ! Calculate lateral detrainment rate for cloud layer
         deltmf(i,1) = Cepsmf
         do k = 2,km-1
            if (hctidx(i) > hdidx(i)+1) then      ! Some cloud
               deltmf(i,k) = 0.9 * Cepsmf + Cdelt * (atan((zm(i,k)-(lcl(i)+(hct(i)-lcl(i))/1.5))/ &
                                                          ((hct(i)-lcl(i))/8))+(3.1415926/2))/3.1415926
            else if (k < hdidx(i)) then   ! No cloud, below hd
               deltmf(i,k) = Cepsmf + 0.05 * 1. / (hd(i) - zm(i,k))
            else if (k >= hdidx(i)) then    ! No cloud, above hd
               deltmf(i,k) = deltmf(i,k-1)
            end if
         end do

         ! Calculate mass flux (defined on turbulence levels)
         mf_temfx(i,1) = CM * wstr(i)
         do kt = 2,km-1
            dMdz(i,kt) = (epsmf(i,kt) - deltmf(i,kt)) * mf_temfx(i,kt-1) * dzt(i,kt)
            mf_temfx(i,kt) = mf_temfx(i,kt-1) + dMdz(i,kt)
         end do

         ! WA 12/28/09 If mass flux component > diffusive
         ! component at second level,
         ! reduce M to prevent a stable layer
         MFCth(i,2) = mf_temfx(i,2) * (thup_temfx(i,2)-thetal(i,2) + thup_temfx(i,3)-thetal(i,3)) / 2.
         if (MFCth(i,2) > Fz(i,2)) then
            red_fact = Fz(i,2) / MFCth(i,2)
            do kt = 1,km
               mf_temfx(i,kt) = mf_temfx(i,kt) * red_fact
            end do
         end if  ! Reduce M to prevent stable layer at second level

         ! Calculate mass flux contributions to fluxes (defined on turb levels)
         ! Use log interpolation at first level
         MFCth(i,1) = mf_temfx(i,1) * (thup_temfx(i,1)-thetal(i,1) &
            + (thup_temfx(i,2)-thetal(i,2) - &
            (thup_temfx(i,1)-thetal(i,1))) * log(zt(i,1)/znt(i))/log(zm(i,2)/znt(i)))
         MFCq(i,1) = mf_temfx(i,1) * (qtup_temfx(i,1)-qt(i,1) &
            + (qtup_temfx(i,2)-qt(i,2) - &
            (qtup_temfx(i,1)-qt(i,1))) * log(zt(i,1)/znt(i))/log(zm(i,2)/znt(i)))
         MFCu(i,1) = mf_temfx(i,1) * (UUPD(i,1)-u_temf(i,1) &
            + (UUPD(i,2)-u_temf(i,2) - &
            (UUPD(i,1)-u_temf(i,1))) * log(zt(i,1)/znt(i))/log(zm(i,2)/znt(i)))
         MFCv(i,1) = mf_temfx(i,1) * (VUPD(i,1)-v_temf(i,1) &
            + (VUPD(i,2)-v_temf(i,2) - &
            (VUPD(i,1)-v_temf(i,1))) * log(zt(i,1)/znt(i))/log(zm(i,2)/znt(i)))
         MFCql(i,1) = mf_temfx(i,1) * (qlUPD(i,1)-qgrs(i,1,ntcwx) &
            + (qlUPD(i,2)-qgrs(i,2,ntcwx)) - &
            (qlUPD(i,1)-qgrs(i,1,ntcwx)) * log(zt(i,1)/znt(i))/log(zm(i,2)/znt(i)))
         MFCTE(i,1) = mf_temfx(i,1) * (TEUPD(i,1)-te_temfx(i,1) &
            + (TEUPD(i,2)-te_temfx(i,2) - &
            (TEUPD(i,1)-te_temfx(i,1))) * log(zt(i,1)/znt(i))/log(zm(i,2)/znt(i)))  ! WA Check this
         do kt = 2,km-1
            MFCth(i,kt) = mf_temfx(i,kt) * (thup_temfx(i,kt)-thetal(i,kt) + thup_temfx(i,kt+1)-thetal(i,kt+1)) / 2.
            MFCq(i,kt) = mf_temfx(i,kt) * (qtup_temfx(i,kt)-qt(i,kt) + qtup_temfx(i,kt+1)-qt(i,kt+1)) / 2.
            MFCu(i,kt) = mf_temfx(i,kt) * (UUPD(i,kt)-u_temf(i,kt) + UUPD(i,kt+1)-u_temf(i,kt+1)) / 2.
            MFCv(i,kt) = mf_temfx(i,kt) * (VUPD(i,kt)-v_temf(i,kt) + VUPD(i,kt+1)-v_temf(i,kt+1)) / 2.
            MFCql(i,kt) = mf_temfx(i,kt) * (qlUPD(i,kt)-qgrs(i,kt-1,ntcwx) + qlUPD(i,kt+1)-qgrs(i,kt,ntcwx)) / 2.
            MFCTE(i,kt) = mf_temfx(i,kt) * (TEUPD(i,kt)-te_temfx(i,kt)) ! TE is on turb levels
         end do
         MFCth(i,km) = 0
         MFCq(i,km) = 0
         MFCu(i,km) = 0
         MFCv(i,km) = 0
         MFCql(i,km) = 0
         MFCTE(i,km) = 0

         ! Calculate cloud fraction (on mass levels)
         cf3d_temfx(i,1) = 0.0
         cfm_temfx(i) = 0.0
         do k = 2,km
            ! if (wupd_temfx(i,k-1) >= 1.0e-15 .AND. wupd_temfx(i,k) >= 1.0e-15 .AND. .NOT. isnan(wupd_temfx(i,k-1)) .AND. .NOT. isnan(wupd_temfx(i,k))) then
            if (wupd_temfx(i,k-1) >= 1.0e-15 .AND. wupd_temfx(i,k) >= 1.0e-15) then
               au(i,k) = ((mf_temfx(i,k-1)+mf_temfx(i,k))/2.0) / ((wupd_temfx(i,k-1)+wupd_temfx(i,k))/2.0)  ! WA average before divide, is that best?
            else
               au(i,k) = 0.0
            end if
            sigq2 = au(i,k) * (qtup_temfx(i,k)-qt(i,k))
            if (sigq2 > 0.0) then
               sigq(i,k) = sqrt(sigq2)
            else
               sigq(i,k) = 0.0
            end if
            ! rst = rsat(p2d(i,k),theta(i,k)*exner(i,k),ep2)
            rst = rsat(p2d(i,k-1),theta(i,k-1)*exner(i,k-1),ep2)
            qst(i,k) = rst / (1. + rst)
            satdef(i,k) = qt(i,k) - qst(i,k)
            if (satdef(i,k) <= 0.0) then
               if (sigq(i,k) > 1.0e-15) then
                  cf3d_temfx(i,k) = max(0.5 + 0.36 * atan(1.55*(satdef(i,k)/sigq(i,k))),0.0)
               else
                  cf3d_temfx(i,k) = 0.0
               end if
            else
               cf3d_temfx(i,k) = 1.0
            end if
            if (zm(i,k) < lcl(i)) then
               cf3d_temfx(i,k) = 0.0
            end if
            ! Put max value so far into cfm
            if (zt(i,k) <= hmax) then
               cfm_temfx(i) = max(cf3d_temfx(i,k),cfm_temfx(i))
            end if
         end do

      else    ! not is_convective, no MF components
         do kt = 1,km
            MFCth(i,kt) = 0
            MFCq(i,kt) = 0
            MFCu(i,kt) = 0
            MFCv(i,kt) = 0
            MFCql(i,kt) = 0
            MFCTE(i,kt) = 0
         end do
         lcl(i) = zm(i,km-1)
         hct(i) = zm(i,1)
         hctidx(i) = 1
         hd(i) = zm(i,1)
         hdidx(i) = 1
         ht(i) = hd(i)
         ! Cloud fraction calculations
         cf3d_temfx(i,1) = 0.0
         cfm_temfx(i) = 0.0
         do k = 2,km
            if (qgrs(i,k-1,ntcwx) > 1.0e-15) then
               cf3d_temfx(i,k) = 1.0
            else
               cf3d_temfx(i,k) = 0.0
            end if
            ! Put max value so far into cfm
            if (zt(i,k) <= hmax) then
               cfm_temfx(i) = max(cf3d_temfx(i,k),cfm_temfx(i))
            end if
         end do

      end if   ! MF components or not
      cf3d_temfx(i,km) = 0.0
      ! Mass flux block ends here

      ! Flux profiles
      do kt = 2,km
         ! Fz(i,kt) = -kh(i,kt) * dthdz(i,kt)
         shf_temfx(i,kt) = Fz(i,kt) + MFCth(i,kt)
         QFK(i,kt) = -kh(i,kt) * dqtdz(i,kt)
         qf_temfx(i,kt) = QFK(i,kt) + MFCq(i,kt)
         uwk(i,kt) = -kmm(i,kt) * dudz(i,kt)
         uw_temfx(i,kt) = uwk(i,kt) + MFCu(i,kt)
         vwk(i,kt) = -kmm(i,kt) * dvdz(i,kt)
         vw_temfx(i,kt) = vwk(i,kt) + MFCv(i,kt)
      end do

      ! Surface momentum fluxes
      ! WA TEST 11/7/13 use w* as a component of the mean wind inside the
      ! u* calculation instead of in the velocity scale below (Felix)
      ! ust(i) = sqrt(ftau(i,1)/ftau0) * sqrt(u_temf(i,2)**2. + v_temf(i,2)**2.) * leps(i,1) / log(zm(i,2)/znt(i)) / zt(i,1)
      ust(i) = sqrt(ftau(i,1)/ftau0) * sqrt(u_temf(i,2)**2. + v_temf(i,2)**2. + (0.5*wstr(i))**2.) * leps(i,1) / log(zm(i,2)/znt(i)) / zt(i,1)

      ang(i) = atan2(v_temf(i,2),u_temf(i,2))
      uw_temfx(i,1) = -cos(ang(i)) * ust(i)**2.
      vw_temfx(i,1) = -sin(ang(i)) * ust(i)**2.

      ! Calculate mixed scaling velocity (Moeng & Sullivan 1994 JAS p.1021)
      ! Replaces ust everywhere
      ! WA TEST 11/7/13 back to wm = u* but with "whole" wind in u* above
      wm(i) = ust(i)
      ! WA 7/23/10 reduce velocity scale to fix excessive fluxes
      ! wm(i) = 0.5 * (1./5. * (wstr(i)**3. + 5. * ust(i)**3.)) ** (1./3.)

      ! Specified flux versions (flux is modified by land surface)
      ! WA 5/31/13 use whole surface flux to improve heat conservation
      shf_temfx(i,1) = hfx(i)/(rho(i,1)*cp)
      qf_temfx(i,1) = qfx(i)/rho(i,1)
      Fz(i,1) = shf_temfx(i,1) - MFCth(i,1)
      QFK(i,1) = qf_temfx(i,1) - MFCq(i,1)

      ! Calculate thetav and its flux
      ! From Lewellen 2004 eq. 3
      ! WA The constants and combinations below should probably be replaced
      ! by something more consistent with the WRF system, but for now
      ! I don't want to take the chance of breaking something.
      do kt = 2,km-1
         alpha2(i,kt) = 0.61 * (thetal(i,kt) + thetal(i,kt+1)) / 2.
         beta2(i,kt) = (100000. / prsi(i,kt))**0.286 * 2.45e-6 / 1004.67 - 1.61 * (thetal(i,kt) + thetal(i,kt+1)) / 2.
      end do
      alpha2(i,1) = 0.61 * (thetal(i,1) + (thetal(i,2)-thetal(i,1)) * (zt(i,2) - znt(i)) / (zm(i,2) - znt(i)))
      alpha2(i,km) = 0.61 * thetal(i,km)
      beta2(i,1) = (100000. / prsi(i,1))**0.286 * 2.45e-6 / &
         1004.67 - 1.61 * (thetal(i,1) + (thetal(i,2) - thetal(i,1)) &
         * (zt(i,2) - znt(i)) / (zm(i,2) - znt(i)))
      beta2(i,km) = (100000. / prsi(i,km))**0.286 * 2.45e-6 / 1004.67 - 1.61 * thetal(i,km)
      if ( is_convective ) then ! Activate EDMF components
         do kt = 1,km-1
            MFCthv(i,kt) = (1. + 0.61 * (qtup_temfx(i,kt)+qtup_temfx(i,kt+1))) / 2. * MFCth(i,kt) + &
                           alpha2(i,kt) * MFCq(i,kt) + beta2(i,kt) * MFCql(i,kt)
         end do
         MFCthv(i,km) = 0.
      else    ! No MF components
         do kt = 1,km
            MFCthv(i,kt) = 0.
         end do
      end if

      do kt = 1,km
         THVF(i,kt) = (1. + 0.61 * qt(i,kt)) * Fz(i,kt) + alpha2(i,kt) * QFK(i,kt) + MFCthv(i,kt)
      end do

      ! Update mean variables:
      ! This is done with implicit solver for diffusion part followed
      ! by explicit solution for MF terms.
      ! Note that Coriolis terms that were source terms for u and v
      ! in Matlab code are now handled by WRF outside this PBL context.

      u_new(i,:) = u_temf(i,:)
      call solve_implicit_temf(kmm(i,1:km-1),u_new(i,2:km), &
         uw_temfx(i,1),dzm(i,1:km-1),dzt(i,1:km-1),1,km-1,dt,.FALSE.)
      do k = 2,km-1
         u_new(i,k) = u_new(i,k) + dt * (-(MFCu(i,k)-MFCu(i,k-1))) / dzm(i,k)
      end do

      v_new(i,:) = v_temf(i,:)
      call solve_implicit_temf(kmm(i,1:km-1),v_new(i,2:km), &
         vw_temfx(i,1),dzm(i,1:km-1),dzt(i,1:km-1),1,km-1,dt,.FALSE.)
      do k = 2,km-1
         v_new(i,k) = v_new(i,k) + dt * (-(MFCv(i,k)-MFCv(i,k-1))) / dzm(i,k)
      end do

      call solve_implicit_temf(kh(i,1:km-1),thetal(i,2:km),Fz(i,1),dzm(i,1:km-1),&
                               dzt(i,1:km-1),1,km-1,dt,.FALSE.)
      do k = 2,km-1
         thetal(i,k) = thetal(i,k) + dt * (-(MFCth(i,k)-MFCth(i,k-1))) / dzm(i,k)
      end do

      call solve_implicit_temf(kh(i,1:km-1),qt(i,2:km),QFK(i,1),dzm(i,1:km-1),&
                               dzt(i,1:km-1),1,km-1,dt,.FALSE.)
      do k = 2,km-1
         qt(i,k) = qt(i,k) + dt * (-(MFCq(i,k)-MFCq(i,k-1))) / dzm(i,k)
      end do

      ! Stepping TE forward is a bit more complicated
      te_temfx(i,1) = ust(i)**2. / ftau(i,1) * (1. + ratio(i,1))
      if ( is_convective ) then
         ! WA currently disabled if MFopt=false, is that right?
         convection_TKE_surface_src(i) = 2. * beta(i,1) * shf_temfx(i,1)
      else
         convection_TKE_surface_src(i) = 0.
      end if
      te_temfx(i,1) = max(te_temfx(i,1), &
                          (leps(i,1) / Cgamma * (ust(i)**2. * S(i,1) + convection_TKE_surface_src(i)))**(2./3.))
      if (te_temfx(i,1) > 20.0) then
         te_temfx(i,1) = 20.0    ! WA 9/28/11 limit max TE
      end if
      sfcFTE(i) = -(kmm(i,1)+kmm(i,2)) / 2. * (te_temfx(i,2)-te_temfx(i,1)) / dzm(i,2)

      do kt = 1,km
         if (THVF(i,kt) >= 0) then
            buoy_src(i,kt) = 2. * g / thetav(i,kt) * THVF(i,kt)
         else
            buoy_src(i,kt) = 0.  ! Cancel buoyancy term when locally stable
         end if
         srcs(i,kt) = -uw_temfx(i,kt) * dudz(i,kt) - vw_temfx(i,kt) * dvdz(i,kt) - &
                      Cgamma * te_temfx(i,kt)**1.5 * linv(i,kt) + buoy_src(i,kt)
      end do
      call solve_implicit_temf((kmm(i,1:km-1)+kmm(i,2:km))/2.0, &
         te_temfx(i,2:km),sfcFTE(i),dzt(i,2:km),dzt(i,1:km-1),1,km-1,dt,.false.)
      do kt = 2,km-1
         te_temfx(i,kt) = te_temfx(i,kt) + dt * srcs(i,kt)
         te_temfx(i,kt) = te_temfx(i,kt) + dt * (-(MFCTE(i,kt)-MFCTE(i,kt-1))) / dzt(i,kt)
         if (te_temfx(i,kt) < TEmin) te_temfx(i,kt) = TEmin
      end do
      te_temfx(i,km) = 0.0
      do kt = 2,km-1
         if (te_temfx(i,kt) > 20.0) then
            te_temfx(i,kt) = 20.0    ! WA 9/29/11 reduce limit max TE from 30
         end if
      end do

      ! Done with updates, now convert internal variables back to WRF vars
      do k = 1,km
         ! Populate kh_temfx, kmm_temfx from kh, km
         kh_temfx(i,k) = kh(i,k)
         kmm_temfx(i,k) = kmm(i,k)
      end do

      ! Convert thetal, qt back to theta, qv, qc
      ! See opposite conversion at top of subroutine
      ! WA this accounts for offset of indexing between
      ! WRF and TEMF, see notes at top of this routine.
      call thlqt2thqvqc(thetal(i,2:km),qt(i,2:km), &
         theta_new(i,1:km-1),qvx_new(i,1:km-1),qcx_new(i,1:km-1), &
         p2d(i,1:km-1),exner(i,1:km-1),1,km-1,ep2,xlv,cp)

      do k = 1,km-1
         ! Calculate tendency terms
         ! WA this accounts for offset of indexing between
         ! WRF and TEMF, see notes at top of this routine.
         rubltenx(i,k) = rubltenx(i,k) + (u_new(i,k+1) - u_temf(i,k+1)) / dt
         rvbltenx(i,k) = rvbltenx(i,k) + (v_new(i,k+1) - v_temf(i,k+1)) / dt
         rtbltenx(i,k) = rtbltenx(i,k) + (theta_new(i,k) - theta(i,k)) * exner(i,k) / dt ! convert from pot temp to air temp

          !WL*do l=1,nvdiff-1
            rqbltenx(i,k,ntqvx) = rqbltenx(i,k,ntqvx) + (qvx_new(i,k) - qgrs(i,k,ntqvx)) / dt
            rqbltenx(i,k,ntcwx) = rqbltenx(i,k,ntcwx) + (qcx_new(i,k) - qgrs(i,k,ntcwx)) / dt
            !WL*if (l==ntqvx) then
            !    rqbltenx(i,k,ntqvx) = rqbltenx(i,k,ntqvx) + (qvx_new(i,k) - qgrs(i,k,ntqvx)) / dt
            !WL*else if (l==ntcwx) then
            !    rqbltenx(i,k,ntcwx) = rqbltenx(i,k,ntcwx) + (qcx_new(i,k) - qgrs(i,k,ntcwx)) / dt
            !WL*else
            !    rqbltenx(i,k,l) = rqbltenx(i,k,l) + (0. - qgrs(i,k,l)) / dt
            !WL*endif
          !WL*enddo
         !WL* Cumulative_change of u, v, t, q due to PBL parameterization
          if(lssav .and. ldiag3d .and. .not. flag_for_pbl_generic_tend) then
            dt3dt_pbl(i,k) = dt3dt_pbl(i,k) + (theta_new(i,k) - theta(i,k)) * exner(i,k)
            du3dt_pbl(i,k) = du3dt_pbl(i,k) + (u_new(i,k+1) - u_temf(i,k+1))
            dv3dt_pbl(i,k) = dv3dt_pbl(i,k) + (v_new(i,k+1) - v_temf(i,k+1))
            if (qdiag3d) then
              dq3dt_pbl(i,k) = dq3dt_pbl(i,k) + (qvx_new(i,k) - qgrs(i,k,ntqvx))
            end if
          end if
         !*WL
      end do
      rubltenx(i,km) = 0.
      rvbltenx(i,km) = 0.
      rtbltenx(i,km) = 0.
      rqbltenx(i,km,nvdiff) = 0.

      ! Populate surface exchange coefficient variables to go back out
      ! for next time step of surface scheme
      ! WA 2/16/11 removed, not needed in BL

      ! Populate 10 m winds and 2 m temp
      ! WA Note this only works if first mass level is above 10 m
      u10(i) = u_new(i,2) * log(10.0/znt(i)) / log(zm(i,2)/znt(i))
      v10(i) = v_new(i,2) * log(10.0/znt(i)) / log(zm(i,2)/znt(i))
      t2(i) = (tsk(i)/exner(i,1) + (theta_new(i,1) - tsk(i)/exner(i,1)) * log(2.0/z0t(i)) / log(zm(i,2)/z0t(i))) * exner(i,1)  ! WA this should also use pi at z0

      ! Populate diagnostic variables
      hd_temfx(i) = hd(i)
      lcl_temfx(i) = lcl(i)
      hct_temfx(i) = hct(i)

      ! Send updraft liquid water back
      if ( is_convective) then
         do k = 1,km-1
            qlup_temfx(i,k) = qlUPD(i,k)
         end do
      else
         qlup_temfx(i,1) = qgrs(i,1,ntcwx)!qcx(i,1)
         do k = 2,km-1
            qlup_temfx(i,k) = qgrs(i,k-1,ntcwx)!qcx(i,k-1)
         end do
      end if
      qlup_temfx(i,km) = qgrs(i,km,ntcwx)

   end do  ! Main (i) loop


   return

   end subroutine bl_temf_run


!--------------------------------------------------------------------
!
   subroutine thlqt2thqvqc(thetal,qt,theta,qv,qc,p,piex,kbot,ktop,ep2,Lv,Cp)

!  Calculates theta, qv, qc from thetal, qt.
!  Originally from RAMS (subroutine satadjst) by way of Hongli Jiang.

   implicit none
   integer, intent(in   ) :: kbot, ktop
   real,    dimension( kbot:ktop ), intent(in   ) :: thetal, qt
   real,    dimension( kbot:ktop ), intent(  out) :: theta, qv, qc
   real,    dimension( kbot:ktop ), intent(in   ) :: p, piex
   real,    intent(in   ) :: ep2, Lv, Cp
!
!  Local variables
   integer :: k, iterate
   real :: T1, Tt
   real, dimension( kbot:ktop) :: rst
   real, dimension( kbot:ktop) :: Tair, rc, rt, rv
!
   do k = kbot,ktop
      T1 = thetal(k) * piex(k)   ! First guess T is just thetal converted to T
      Tair(k) = T1
      rt(k) = qt(k) / (1. - qt(k))

      do iterate = 1,20
         rst(k) = rsat(p(k),Tair(k),ep2)
         rc(k) = max(rt(k) - rst(k), 0.)
         Tt = 0.7*Tair(k) + 0.3*T1 * (1.+Lv*rc(k) / (Cp*max(Tair(k),253.)))
         if ( abs(Tt - Tair(k)) < 0.001) exit !GOTO 100
         Tair(k) = Tt
      end do
!100   continue
      rv(k) = rt(k) - rc(k)
      qv(k) = rv(k) / (1. + rv(k))
      qc(k) = rc(k) / (1. + rc(k))
      theta(k) = Tair(k) / piex(k)
   end do ! k loop
   return
   end subroutine thlqt2thqvqc
!

!--------------------------------------------------------------------
!
   subroutine solve_implicit_temf(Khlf,psi_n,srf_flux,dzm,dzt,kbot,ktop,dt,print_flag)

!  Implicit solution of vertical diffusion for conserved variable
!  psi given diffusivity Khlf on turbulence levels,
!  and surface flux srf_flux.
!  dzm is delta_z of mass levels, dzt is delta_z of turbulence levels.
!  dt is timestep (s).

   implicit none
   integer  :: kbot, ktop
   logical  :: print_flag
   real :: srf_flux, dt
   real,    dimension( kbot:ktop ), intent(in   ) :: Khlf
   real,    dimension( kbot:ktop ), intent(in   ) :: dzm, dzt
   real,    dimension( kbot:ktop ), intent(inout) :: psi_n
!
!  Local variables
   integer :: k
   real,    dimension( kbot:ktop ) :: AU, BU, CU, YU
!
   AU(kbot) = Khlf(kbot) / (dzm(kbot)*dzt(kbot))
   BU(kbot) = -1.0/dt - Khlf(kbot+1)/(dzm(kbot+1)*dzt(kbot+1))
   CU(kbot) = Khlf(kbot+1)/(dzm(kbot)*dzt(kbot+1))
   YU(kbot) = -psi_n(kbot)/dt - srf_flux/dzm(kbot)

   do k = kbot+1,ktop-1
      ! Subdiagonal (A) vector
      AU(k) = Khlf(k) / (dzm(k) * dzt(k))
      ! Main diagonal (B) vector
      BU(k) = -1.0/dt - (Khlf(k)/dzt(k) + Khlf(k+1)/dzt(k+1)) / dzm(k)
      ! Superdiagonal (C) vector
      CU(k) = Khlf(k+1) / (dzm(k)*dzt(k+1))
      ! Result vector
      YU(k) = -psi_n(k)/dt
   end do ! k loop

   AU(ktop) = 0.
   BU(ktop) = -1.0 / dt
   YU(ktop) = -psi_n(ktop) / dt

   ! Compute result with tridiagonal routine
   psi_n = trid(AU,BU,CU,YU,kbot,ktop)

   return
   end subroutine solve_implicit_temf

   function trid(a,b,c,r,kbot,ktop)

!  Solves tridiagonal linear system.
!  Inputs are subdiagonal vector a, main diagonal b, superdiagonal c,
!  result r, column top and bottom indices kbot and ktop.
!  Originally from Numerical Recipes section 2.4.

   implicit none
   real,    dimension( kbot:ktop ) :: trid
   integer  :: kbot, ktop
   real,    dimension( kbot:ktop ), intent(in   ) :: a, b, c, r
!
!  Local variables
   integer :: k
   real    :: bet
   real,    dimension( kbot:ktop ) :: gam, u
!
   bet = b(kbot)
   u(kbot) = r(kbot) / bet

   do k = kbot+1,ktop
      gam(k) = c(k-1) / bet
      bet = b(k) - a(k)*gam(k)
      u(k) = (r(k) - a(k)*u(k-1)) / bet
   end do

   do k = ktop-1,kbot,-1
      u(k) = u(k) - gam(k+1)*u(k+1)
   end do

   trid = u

   return
   end function trid

   !--------------------------------------------------------------------
!
   real function rsat(p,T,ep2)

!  Calculates the saturation mixing ratio with respect to liquid water
!  Arguments are pressure (Pa) and absolute temperature (K)
!  Uses the formula from the ARM intercomparison setup.
!  Converted from Matlab by WA 7/28/08

implicit none
real p, T, ep2
real temp, x
real, parameter :: c0 = 0.6105851e+3
real, parameter :: c1 = 0.4440316e+2
real, parameter :: c2 = 0.1430341e+1
real, parameter :: c3 = 0.2641412e-1
real, parameter :: c4 = 0.2995057e-3
real, parameter :: c5 = 0.2031998e-5
real, parameter :: c6 = 0.6936113e-8
real, parameter :: c7 = 0.2564861e-11
real, parameter :: c8 = -0.3704404e-13

temp = T - 273.15

x =c0+temp*(c1+temp*(c2+temp*(c3+temp*(c4+temp*(c5+temp*(c6+temp*(c7+temp*c8)))))))
rsat = ep2*x/(p-x)

return
end function rsat

!-------------------------------------------------------------------------------
end module bl_temf
!-------------------------------------------------------------------------------
