      module Models
          implicit none 
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @author 
!> Martin Boeff, Anxin Ma
! DESCRIPTION: 
!> Flow and hardening laws
!> @param[in] nslip            Number of slip systems
!> @param[in] nprops            Number of material parameters
!> @param[in] props            Material parameters
!> @param[in] IVB                    Critical resolved shear stress \f$ \tau_c \f$ is passed into routine
!> @param[in] tau                    Resolved shear stress \f$ \tau \f$ is passed into routine
!> @param[in] tau_nsd(nslip)         Resolved shear stress \f$ \tau \f$ (non-Schmidt)
!> @param[in] crsF_mat(nslip,nslip)  Hardening matrix \f$ h_{\alpha\beta} \f$
!> @param[in] IIVB(nslip)            Hardening due to first gradient of FP
!> @param[in] BIVB(nslip)            Hardening due to second gradient of FP
!> @param[out] dgmdt            Shear rate \f$\dot{\gamma}\f$
!> @param[out] ddgmdt_dtau        \f$\frac{\dot{\gamma}}{\dot{\tau}}\f$
!> @param[out] ddgmdt_dIVB        \f$\frac{\dot{\gamma}}{\dot{\tau_C}}\f$
!> @param[out] dIVBdt        Critical resolved shear stress \f$\tau_C\f$
!> @param[out] ddIVBdt_ddgmdt    \f$\frac{\tau}{\dot{\gamma}}\f$
!> @param[out] ddIVBdt_dIVB        \f$\frac{\tau}{\dot{\tau_C}}\f$
!> @param[out] ctrl_warningflag      ctrl_warningflag is 112 if \f$\tau<10^{-10}\tau_c^0\f$ and \f$\tau>\tau_{c,s}\f$
!---------------------------------------------------------------------------
          contains
          subroutine cal_flow_harden_phenom(nslip, nprops, props, IVB, grainID, tau, tau_nsd, 
     &                                      crsF_mat, IIVB,backstressn, modeltype, dgmdt,
     &                                      ddgmdt_dtau,ddgmdt_dIVB, dIVBdt,ddIVBdt_ddgmdt,
     &                                      ddIVBdt_dIVB,dbackstressdt, ctrl_tempwarningflag, statev)
              use GlobalVariables, only: ctrl_isexplicit, ctrl_nonlocal, ctrl_infinity
              use MaterialData, only: materialproperties, grainDia
              integer, intent(in) ::  nslip                       !Number of slip systems
              integer(4), intent(in) ::  modeltype                   !Material model
              integer, intent(in) ::  nprops                      !Number of material parameters
              real(8), intent(in) ::  props(nprops)               !Material parameters
              real(8), intent(in) ::  IVB(nslip)                  !Critical resolved shear stress
              real(8), intent(in) ::  tau(nslip)                  !Resolved shear stress
              real(8), intent(in) ::  tau_nsd(nslip)              !Resolved shear stress (non-Schmidt)
              real(8), intent(in) ::  crsF_mat(nslip,nslip)       !hardening matrix h_alphabeta
              real(8), intent(in) ::  IIVB(nslip)                 !Hardening due to first gradient of FP
              !real(8), intent(in) ::  BIVB(nslip)                !Hardening due to second gradient of FP
              real(8), intent(in) ::  backstressn(nslip)          !Backstress calculated by phanomenological law
              integer, intent(in) ::  grainID                     !Grain ID
              !real(8), intent(in) ::  IVBdt(nslip)                !CHANGED added
              real(8), intent(in) ::  statev(115)
              !Out                                                                                  
              real(8), intent(out) :: dgmdt(nslip)                !Shear rate
              real(8), intent(out) :: ddgmdt_dtau(nslip)          !Derivatives of dgmdt w.r.t. dtau
              real(8), intent(out) :: ddgmdt_dIVB(nslip)          !Derivatives of dgmdt w.r.t. IVB
              real(8), intent(out) :: dIVBdt(nslip)               !Derivatives of IVB w.r.t. dt
              real(8), intent(out) :: ddIVBdt_ddgmdt(nslip,nslip) !Derivatives of divbdt w.r.t. dgmdt
              real(8), intent(out) :: ddIVBdt_dIVB(nslip,nslip)   !Derivatives of divbdt w.r.t. IVB 
              real(8), intent(out) :: dbackstressdt(nslip) 
              integer, intent(out) :: ctrl_tempwarningflag            !Numerical error flag
!-------------------------------------------------------------------------------
!-------------Local Parameters
              real(8) :: x1,x2                !Support variables
              integer :: i, j, k, is, js, ns  !Loop variables
              real(8) :: IVB_eff(nslip)       ! Effective critical resolved shear stress
              real(8) :: tau_eff(nslip)       ! Effective resolved shear stress
!-------------Material Properties
              real(8) :: pw_fl             !Power term in shear rate equation
              real(8) :: shrt_0            !Shear rate gamma 0
              real(8) :: hdrt_0            !Hardening rate h0 from Roters p.50 eq.6.10
              real(8) :: crss_0            !Initial critical resolved shear stress
              real(8) :: crss_s            !Saturated critical resolved shear stress
              real(8) :: pw_hd             !Power term in satuarion equation for slip resistance
              real(8) :: Adir, Adyn
              real(8) :: C1,C2,NEWC3,C_T1,C_T2,C_T3,T,m   !Harris
              real(8) :: PS_rate,P
!------------------------------------------------------------------------------
!-------------Initializing values
              ctrl_tempwarningflag = 0
              dgmdt = 0.d0
              ddgmdt_dtau =0.d0
              ddgmdt_dIVB = 0.d0
              dIVBdt = 0.d0
              ddIVBdt_ddgmdt = 0.d0
              ddIVBdt_dIVB = 0.d0
              pw_fl   = materialproperties(modeltype,5)        
              shrt_0  = materialproperties(modeltype,6)
              hdrt_0  = materialproperties(modeltype,7)
              if (modeltype ==1  .or. modeltype ==2)  then
                  crss_0  = materialproperties(modeltype,8)
                  crss_s  = materialproperties(modeltype,9)
                  pw_hd   = materialproperties(modeltype,10)
                  Adir    = materialproperties(modeltype,11)
                  Adyn    = materialproperties(modeltype,12)
              else
                  crss_0  = materialproperties(modeltype,8)+materialproperties(modeltype,9)/sqrt(grainDia(grainID))
                  crss_s  = materialproperties(modeltype,10)
                  pw_hd   = materialproperties(modeltype,11)
                  Adir    = materialproperties(modeltype,12)
                  Adyn    = materialproperties(modeltype,13)
              endif

              if (modeltype == 5) then
                  C1          = materialproperties(modeltype,16)
                  C2          = materialproperties(modeltype,17)
                  NEWC3       = materialproperties(modeltype,18)
                  C_T1        = materialproperties(modeltype,19)
                  C_T2        = materialproperties(modeltype,20)
                  C_T3        = materialproperties(modeltype,21)
                  T           = materialproperties(modeltype,22)
                  PS_rate     = materialproperties(modeltype,23)
                  m           = materialproperties(modeltype,24)
              endif

            !   if (modeltype == 2) then
            !     write(*,*), "shrt_0: ",shrt_0
            !     write(*,*), "hdrt_0: ",hdrt_0
            !     write(*,*), "crss_0: ",crss_0
            !     write(*,*), "crss_s: ",crss_s
            !     write(*,*), "pw_hd: ",pw_hd
            !     write(*,*), "Adir: ",Adir
            !     write(*,*), "Adyn: ",Adyn
            !     write(*,*), "C1: ",C1
            !     write(*,*), "C2: ",C2
            !     write(*,*), "NEWC3: ",NEWC3
            !     write(*,*), "C_T1: ",C_T1
            !     write(*,*), "C_T2: ",C_T2
            !     write(*,*), "C_T3: ",C_T3
            !     write(*,*), "T: ",T
            !     write(*,*), "PS_rate: ",PS_rate
            !     write(*,*), "m: ",m
              
            !     call XIT
            !   endif
!------------------------------------------------------------------------------
!-------------Check resolved shear stress and resistence
              do is=1,nslip
!------------------x1=crss_0*1.d-10
                  x1=0
                  x2=crss_s
                  if(IVB(is)<=x1 .or. IVB(is)>x2)then
!---------------------write(*,*) 'IVB(is): ', IVB(is)
                      ctrl_tempwarningflag=112
                      return
                  endif
              enddo
!-------------Define effective resolved shear stress and effective hardening parameter    
              if (ctrl_nonlocal .eqv. .true.) then
                  IVB_eff= IVB + IIVB
!-----------------tau_eff= tau + tau_nsd
                  tau_eff= tau 
              endif
              if (modeltype == 1) then
                  IVB_eff= IVB 
                  tau_eff= tau
              endif
              if (modeltype == 2) then
                  IVB_eff= IVB  
                  tau_eff= tau - backstressn
              endif
              if (modeltype == 3 .OR. modeltype == 4) then
                  IVB_eff= IVB 
                  tau_eff= tau + tau_nsd - backstressn
              endif
              if (modeltype == 5) then
                  pw_fl=1
                  IVB_eff= IVB 
                  tau_eff= tau - backstressn + tau_nsd
              endif
!-------------------------------------------------------------------------------
!-------------Shear rate, derivative of shear rate w.r.t. pk2i,IVB
              do is=1,nslip
                  dgmdt(is) = shrt_0*(abs(tau_eff(is)/IVB_eff(is)))**pw_fl*sign(1.d0,tau_eff(is))
                  if (abs(dgmdt(is)) > ctrl_infinity) then !!! infinity check
                      !write(*,*) 'dgmdt(is): ', dgmdt(is)
                      ctrl_tempwarningflag=113   
                      return
                  endif
                  if(ctrl_isexplicit/=1)then
                      ddgmdt_dtau(is) = pw_fl/IVB_eff(is)*shrt_0*(abs(tau_eff(is)/IVB_eff(is)))**(pw_fl-1)
                      if (abs(ddgmdt_dtau(is)) > ctrl_infinity) then !!! infinity check
                          ctrl_tempwarningflag=114   
                          return
                       endif
                      ddgmdt_dIVB(is) = -pw_fl*dgmdt(is)/IVB_eff(is)
                      if (abs(ddgmdt_dIVB(is)) > ctrl_infinity) then !!! infinity check
                          ctrl_tempwarningflag=115   
                          return
                       endif
                  endif
              enddo
!------------------------------------------------------------------------------
!-------------Evolution rate, derivative of evolution rate w.r.t. pk2i,IVB
              do is=1,nslip
                  do js=1,nslip
                      x1=1-IVB(js)/crss_s
                      divbdt(is) = divbdt(is) + crsF_mat(is,js) * hdrt_0 * abs(dgmdt(js))*x1**pw_hd
                      if (abs(divbdt(is)) > ctrl_infinity) then !!! infinity check
                          ctrl_tempwarningflag=116   
                          return
                      endif
                      if(ctrl_isexplicit/=1)then
                          ddIVBdt_ddgmdt(is,js) = crsF_mat(is,js)*hdrt_0*sign(1.d0,tau_eff(js))*x1**pw_hd 
                          if (abs(ddIVBdt_ddgmdt(is,js)) > ctrl_infinity) then !!! infinity check
                              ctrl_tempwarningflag=117   
                              return
                          endif
                          ddIVBdt_dIVB(is,js) = crsF_mat(is,js)*hdrt_0*ddgmdt_dIVB(js)*sign(1.d0,tau_eff(js))*
     &                                          x1**pw_hd -crsF_mat(is,js)*hdrt_0*abs(dgmdt(js)) *x1**(pw_hd-1)*
     &                                          pw_hd/crss_s
                          if (abs(ddIVBdt_dIVB(is,js)) > ctrl_infinity) then !!! infinity check
                              ctrl_tempwarningflag=118  
                              return
                          endif
                      endif
                  enddo
              enddo
!-------------------------------------------------------------------------------
!-------------Evolution of backstress
              do is=1,nslip
                  dbackstressdt(is) = Adir*dgmdt(is) - Adyn*backstressn(is)*abs(dgmdt(is))
!-----------------if (dbackstressdt(is).GT.0.000001) then
!---------------------print*, dbackstressdt(is)
!-----------------endif
              enddo
              return
          end subroutine cal_flow_harden_phenom

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @author 
!> Philipp Engels, Anxin Ma
! DESCRIPTION: 
!> Flow and hardening laws
!> @param[in] nslip            Number of slip systems
!> @param[in] nprops            Number of Material parameters
!> @param[in] props            Material parameters
!> @param[in] crsF_mat                 Hardening matrix q_alpha_beta from Roters Book p.50 eq.6.10
!> @param[in] IVB                      Critical resolved shear stress \f$ \tau_c \f$ is passed into routine
!> @param[in] tau                      Resolved shear stress \f$ \tau \f$ is passed into routine
!> @param[out] dgmdt            shear rate \f$\dot{\gamma}\f$
!> @param[out] ddgmdt_dtau        \f$\frac{\dot{\gamma}}{\dot{\tau}}\f$
!> @param[out] ddgmdt_dIVB        \f$\frac{\dot{\gamma}}{\dot{\tau_C}}\f$
!> @param[out] dIVBdt               critical resolved shear stress \f$\tau_C\f$
!> @param[out] ddIVBdt_ddgmdt           \f$\frac{\tau}{\dot{\gamma}}\f$
!> @param[out] ddIVBdt_dIVB        \f$\frac{\tau}{\dot{\tau_C}}\f$
!> @param[out] ctrl_warningflag        ctrl_warningflag is 112 if \f$\tau<10^{-10}\tau_c^0\f$ and \f$\tau>\tau_{c,s}\f$
!---------------------------------------------------------------------------
          subroutine cal_flow_harden_dislocation(nslip,nprops,props,IVB, tau,modeltype,
     &                                           dgmdt,ddgmdt_dtau,ddgmdt_dIVB,dIVBdt,ddIVBdt_ddgmdt,
     &                                           ddIVBdt_dIVB,ctrl_warningflag)
              use GlobalVariables, only: ctrl_isexplicit
              use MaterialData, only: materialproperties
              integer, intent(in) :: nslip
              integer, intent(in) :: nprops
              real(8), intent(in) :: props(nprops)
!-------------real(8), intent(in) :: crsF_mat(nslip,nslip) !Wird noch nicht benoetigt
              real(8), intent(in) :: IVB(nslip)
              real(8), intent(in) :: tau(nslip)
              integer, intent(in) :: modeltype
!-------------Out
              real(8), intent(out) :: dgmdt(nslip)
              real(8), intent(out) :: ddgmdt_dtau(nslip)
              real(8), intent(out) :: ddgmdt_dIVB(nslip)
              real(8), intent(out) :: dIVBdt(nslip)
              real(8), intent(out) :: ddIVBdt_ddgmdt(nslip,nslip)
              real(8), intent(out) :: ddIVBdt_dIVB(nslip,nslip)
              integer, intent(out) :: ctrl_warningflag         
!---------------------------------------------------------------------------
!-------------Local Parameters
!-------------integer :: i,j,is,js ns         !Loop variables
              integer :: is
              real(8) :: x1                   !Work variable
              real(8) :: IVB_eff(nslip)              ! Effective critical resolved shear stress
              real(8) :: tau_eff(nslip)              ! Effective resolved shear stress
!-------------real(8) :: IIVB(nslip), BIVB(nslip)    ! For gradient calculations
!-------------Material Properties
              real(8) :: C11,C12,C44          !Material elastic stiffness if c44=2*(c11-c12) -> iso
              real(8) :: pw_fl                !Power term in shear rate equation
              real(8) :: brgvec               !Burgers vector length
              real(8) :: v0                   !Reference dislocation velocity
              real(8) :: c3                   !Taylor hardening parameter
              real(8) :: km1, km2             !Kocks-Mecking law parameters
!-------------real(8) :: phi1,phi,phi2        !Initial euler angles
              real(8) :: Gmod                 !Gmod
!---------------------------------------------------------------------------
!-------------Initializing values
              ctrl_warningflag = 0
              dgmdt = 0.d0
              ddgmdt_dtau = 0.d0
              ddgmdt_dIVB = 0.d0
              dIVBdt = 0.d0
              ddIVBdt_ddgmdt = 0.d0
              ddIVBdt_dIVB = 0.d0
              C11   = materialproperties(modeltype,1)
              C12   = materialproperties(modeltype,2)
              C44   = materialproperties(modeltype,3)
              pw_fl = materialproperties(modeltype,5)
              brgvec= materialproperties(modeltype,6)
              v0    = materialproperties(modeltype,7) 
              km1   = materialproperties(modeltype,9) 
              km2   = materialproperties(modeltype,10) 
              c3    = materialproperties(modeltype,11)
              Gmod = 1./( 1./C44 +2./3.*(2./(C11-C12)) - 1./C44)
!-------------Define effective resolved shear stress and effective hardening parameter        
              if (modeltype == 1) then
                  IVB_eff= IVB 
                  tau_eff= tau
              endif
              if (modeltype == 2) then
                  IVB_eff= IVB 
                  tau_eff= tau
              endif
!---------------------------------------------------------------------------
!-------------Shear rate, derivative of shear rate w.r.t. pk2i,IVB
              do is=1,nslip
                  if(IVB(is)<0.0) then
                      ctrl_warningflag = 220
                      call inout_warningflaglocal(ctrl_warningflag)
                  end if
                  x1 = sqrt(IVB_eff(is))*c3*brgvec*Gmod ! -> Taylor-type tau_c
                  dgmdt(is) = IVB_eff(is)*brgvec*v0*(abs(tau_eff(is))/x1)**pw_fl*sign(1.d0,tau_eff(is))
                  divbdt(is) = (km1*sqrt(IVB_eff(is))-km2*IVB_eff(is))*abs(dgmdt(is))
                  if(ctrl_isexplicit/=1)then
!---------------------Derivatives w.r.t. shear rate
                      ddgmdt_dtau(is) = brgvec*v0*pw_fl*IVB_eff(is)*x1**(-pw_fl)*abs(tau_eff(is))**(pw_fl-1.d0)
                      ddgmdt_dIVB(is) = -0.5d0*brgvec*(pw_fl-2.d0)*v0*abs(tau_eff(is))**pw_fl*x1**(-pw_fl)*sign(1.d0,tau_eff(is))
!---------------------Derivatives w.r.t. dislocation evolution
!---------------------Simply Fortran says there is an Syntax error below !!!
                      ddIVBdt_ddgmdt(is,is) = (km1*sqrt(IVB_eff(is))-km2*IVB_eff(is))*sign(1.d0,dgmdt(is))
                      ddIVBdt_dIVB(is,is) = (1.d0/(2.d0*sqrt(IVB_eff(is))*km1-km2)*
     &                                      abs(dgmdt(is))-(km1*sqrt(IVB_eff(is))-km2*
     &                                      IVB_eff(is))*(brgvec*v0*(pw_fl-2.d0)/(2.d0*x1**pw_fl)*
     &                                      abs(tau_eff(is))**pw_fl))
                  endif
              enddo
              return
          end subroutine cal_flow_harden_dislocation
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      end module Models
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++