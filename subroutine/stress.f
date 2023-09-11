! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @author 
!> Martin Boeff
! DESCRIPTION: 
!> Assembles isotropic elastic stiffness matrix
!> @param[in] materialproperties   Array of material properties
!> @param[out] C                   Elastic Stiffness Matrix
!---------------------------------------------------------------------------
      subroutine AssembleStiffness(modeltype,C)
          use MaterialData, only: materialproperties
!---------------------------------------------------------------------------
          implicit none
          integer, intent(in) :: modeltype
          real(8), intent(out):: C(6,6)
!---------Local variable
          integer :: warningflag
          real(8) :: c11, c12, c44
!---------------------------------------------------------------------------
          c11 = materialproperties(modeltype,1)
          c12 = materialproperties(modeltype,2)
          c44 = materialproperties(modeltype,3)
!---------------------------------------------------------------------------
          C=0.0
          C(1,1)=c11; C(1,2)=c12; C(1,3)=c12
          C(2,1)=c12; C(2,2)=c11; C(2,3)=c12
          C(3,1)=c12; C(3,2)=c12; C(3,3)=c11
          C(4,4)=c44*2.0
          C(5,5)=c44*2.0
          C(6,6)=c44*2.0
      end subroutine AssembleStiffness
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @author 
!> Mohamed Sharaf, Martin Boeff, Philipp Engels, Aenne Köster
! DESCRIPTION: 
!> Assembles Schmid matrix for every slip system saved in the form: [nslip,3,3]
!> @param[in] nslip                Number of slip systems
!> @param[in] modeltype            Number of material model
!> @param[out] smdMi               List of Schmidt Matrices for every slip system
!> @param[out] smdVi2              List of Schmidt Matrices for every slip system saved as vector
!> @param[out] smnVi2              List of matrices of plane normals for every slip system saved as vector
!---------------------------------------------------------------------------
      subroutine AssembleLattice(nslip, modeltype, smdMi, smdVi2,smnVi2, 
     &                           n_smdVi2, crsF_mat, vn_mat_Norm, 
     &                           vd_mat_Norm, vl_mat_Norm)
          use MaterialData, only: materialproperties
          use GlobalVariables, only: ctrl_nonlocal
!-------------------------------------------------------------------------
          implicit none
          integer, intent(in) :: nslip
          integer, intent(in) :: modeltype
          real(8), intent(out):: smdMi(nslip,3,3)  !Schmid matrix
          real(8), intent(out):: smdVi2(nslip,6)   !Schmid matrix saved as vector, 4-6 is doubled
          real(8), intent(out):: smnVi2(nslip,6)   !Matrix of plane normals saved as a vector, 4-6 is doubled
          real(8), intent(out):: n_smdVi2(nslip,6) !Schmid matrix saved as vector, 4-6 is doubled
          real(8), intent(out):: crsF_mat(nslip,nslip)   !Hardening matrix
          real(8), intent(out):: vl_mat_Norm(nslip,3)    !Norm of slip plane normal
!---------Local parameters
          integer :: is, i, j                !Loop parameters
          integer :: warningflag
          real(8) :: vd_mat(nslip,3)         !Slip direction
          real(8) :: vn_mat(nslip,3)         !Slip plane normal
          real(8) :: vn1_mat(nslip,3)        !Additional slip plane normal for BCC
          real(8) :: vl_mat(nslip,3)         !Slip plane normal
          real(8) :: vd_mat_Norm(nslip,3)    !Norm of slip direction
          real(8) :: vn_mat_Norm(nslip,3)    !Norm of slip plane normal
          real(8) :: vn1_mat_Norm(nslip,3)   !Norm of additional BCC slip plane normals
          real(8) :: x3
          real(8) :: vn1d_mat(nslip,3)       !vector of outerproduct of 
          real(8) :: vn1d_mat_Norm(nslip,3)  !
          real(8) :: n_smdMi(nslip,3,3)      !
          real(8) :: smdSMi(nslip,3,3)       !Symmetric part of Schmidt Matrix
          real(8) :: smdVi1(nslip,6)         !Schmidt matrix saved as vector
          real(8) :: smnMi(nslip,3,3)        !Matrix of plane normals
          real(8) :: smnSMi(nslip,3,3)       !Symmetric part of matrix of plane normals
          real(8) :: smnVi1(nslip,6)         !Matrix of plane normals saved as vector
          real(8) :: n_smdVi1(nslip,6)       ! 
!---------------------------------------------------------------------------
!---------Start Definition of slip direction, slip plane normal and assembly of hardening matrix
! --------Slip direction and slip plane normal
          if (modeltype ==1  .or. modeltype ==2)  then
              vd_mat(1,:) =[ 0,  1, -1] ; vn_mat(1,:) =[ -1,  1,  1]
              vd_mat(2,:) =[ 1,  0,  1] ; vn_mat(2,:) =[ -1,  1,  1]
              vd_mat(3,:) =[ 1,  1,  0] ; vn_mat(3,:) =[ -1,  1,  1]
              vd_mat(4,:) =[ 0,  1, -1] ; vn_mat(4,:) =[  1,  1,  1]
              vd_mat(5,:) =[ 1,  0, -1] ; vn_mat(5,:) =[  1,  1,  1]
              vd_mat(6,:) =[ 1, -1,  0] ; vn_mat(6,:) =[  1,  1,  1]
              vd_mat(7,:) =[ 0,  1,  1] ; vn_mat(7,:) =[  1,  1, -1]
              vd_mat(8,:) =[ 1,  0,  1] ; vn_mat(8,:) =[  1,  1, -1]
              vd_mat(9,:) =[ 1, -1,  0] ; vn_mat(9,:) =[  1,  1, -1]
              vd_mat(10,:)=[ 0,  1,  1] ; vn_mat(10,:)=[  1, -1,  1]
              vd_mat(11,:)=[ 1,  0, -1] ; vn_mat(11,:)=[  1, -1,  1]
              vd_mat(12,:)=[ 1,  1,  0] ; vn_mat(12,:)=[  1, -1,  1]
          else
              vn_mat( 1,:)=[ 0,  1,  1] ; vd_mat( 1,:)=[ 1,  1, -1]
              vn_mat( 2,:)=[ 1,  0,  1] ; vd_mat( 2,:)=[ 1,  1, -1]
              vn_mat( 3,:)=[ 1, -1,  0] ; vd_mat( 3,:)=[ 1,  1, -1]
              vn_mat( 4,:)=[ 0,  1, -1] ; vd_mat( 4,:)=[ 1, -1, -1]
              vn_mat( 5,:)=[ 1,  0,  1] ; vd_mat( 5,:)=[ 1, -1, -1]
              vn_mat( 6,:)=[ 1,  1,  0] ; vd_mat( 6,:)=[ 1, -1, -1]
              vn_mat( 7,:)=[ 0,  1,  1] ; vd_mat( 7,:)=[ 1, -1,  1]
              vn_mat( 8,:)=[ 1,  0, -1] ; vd_mat( 8,:)=[ 1, -1,  1]
              vn_mat( 9,:)=[ 1,  1,  0] ; vd_mat( 9,:)=[ 1, -1,  1]
              vn_mat(10,:)=[ 0,  1, -1] ; vd_mat(10,:)=[ 1,  1,  1]
              vn_mat(11,:)=[ 1,  0, -1] ; vd_mat(11,:)=[ 1,  1,  1]
              vn_mat(12,:)=[ 1, -1,  0] ; vd_mat(12,:)=[ 1,  1,  1]
              vn1_mat(1,:)=[-1,  1,  0]
              vn1_mat(2,:)=[ 0, -1,  1]
              vn1_mat(3,:)=[ 1,  0, -1]
              vn1_mat(4,:)=[ 0,  1, -1]
              vn1_mat(5,:)=[-1, -1,  0]
              vn1_mat(6,:)=[ 1,  0,  1]
              vn1_mat(7,:)=[ 1, -1,  0]
              vn1_mat(8,:)=[ 0,  1,  1]
              vn1_mat(9,:)=[-1,  0, -1]
              vn1_mat(10,:)=[ 0, -1, -1]
              vn1_mat(11,:)=[ 1,  1,  0]
              vn1_mat(12,:)=[-1,  0,  1]
          endif
!---------Calculate line direction of something
          vl_mat=0.0
          vl_mat_Norm=0.0
          if (ctrl_nonlocal .eqv. .true.) then
              do is=1,nslip
                  vl_mat(is,1)=vn_mat(is,2)*vd_mat(is,3)-vn_mat(is,3)*vd_mat(is,2)
                  vl_mat(is,2)=vn_mat(is,3)*vd_mat(is,1)-vn_mat(is,1)*vd_mat(is,3)
                  vl_mat(is,3)=vn_mat(is,1)*vd_mat(is,2)-vn_mat(is,2)*vd_mat(is,1)
                  x3=sqrt(sum(vl_mat(is,:)**2))
                  vl_mat_Norm(is,:)=vl_mat(is,:)/x3
              enddo
          endif
!---------End Definition of slip direction and slip plane normal and assembly of hardening matrix
!---------------------------------------------------------------------------
          smdMi  = 0.0; smdSMi = 0.0; smnMi  = 0.0; smnSMi = 0.0; n_smdMi = 0.0 
          smdVi1 = 0.0; smdVi2 = 0.0; smnVi1 = 0.0; smnVi2 = 0.0
          vl_mat = 0.0; vl_mat_Norm = 0.0
          vn1d_mat = 0.0; vn1d_mat_Norm = 0.0 
          vn_mat_Norm = 0.0; vn1_mat_Norm = 0.0
          vd_mat_Norm = 0.0; 
          n_smdVi1 = 0.0; n_smdVi2 = 0.0 
          
          do is=1, nslip
!-------------Normalize vectors
              vd_mat_Norm(is,:)=vd_mat(is,:)/sqrt(sum(vd_mat(is,:)**2)) 
              vn_mat_Norm(is,:)=vn_mat(is,:)/sqrt(sum(vn_mat(is,:)**2))
              if (modeltype == 3 .or. modeltype == 4) then
!-----------------Calculate line direction needed for non-Schmidt definition
                  call crossproduct(vl_mat(is,:), vn_mat(is,:), vd_mat(is,:))
                  call crossproduct(vn1d_mat(is,:), vn1_mat(is,:),vd_mat(is,:))
                  vn1d_mat_Norm(is,:) = vn1d_mat(is,:)/sqrt(sum(vn1d_mat(is,:)**2))
                  vn1_mat_Norm(is,:) = vn1_mat(is,:)/sqrt(sum(vn1_mat(is,:)**2))
                  vl_mat_Norm(is,:)=vl_mat(is,:)/sqrt(sum(vl_mat(is,:)**2))
              end if
!             NEU MUSS GETESTET WERDEN
!             n_smdMi_mat(is,i,j)=
!    &            0.61d0*outerprod(vd_mat(is,:),vn1_mat(is,:))+
!    &            0.23d0*outerprod(vl_mat(is,:),vn_mat(is,:))+
!    &            0.55d0*outerprod(vn1d_mat(is,:),vn1_mat(is,:))+
!    &            0.11d0*outerprod(vn_mat(is,:),vn_mat(is,:))+
!    &            0.09d0*outerprod(vl_mat(is,:),vl_mat(is,:))+
!    &            (-0.2d0)*outerprod(vd_mat(is,:),vd_mat(is,:))
              do i=1,3
                  do j=1,3
                      smdMi(is,i,j)=vd_mat_Norm(is,i)*vn_mat_Norm(is,j)
                      smnMi(is,i,j)=vn_mat_Norm(is,i)*vn_mat_Norm(is,j)
                      if (modeltype == 3 .or. modeltype == 4) then
                          n_smdMi(is,i,j)=materialproperties(modeltype,16)*(vd_mat(is,i)*vn1_mat(is,j))
     &                    +materialproperties(modeltype,17)*(vl_mat(is,i)*vn_mat(is,j))
     &                    +materialproperties(modeltype,18)*(vn1d_mat(is,i)*vn1_mat(is,j))
     &                    +materialproperties(modeltype,19)*(vn_mat(is,i)*vn_mat(is,j))
     &                    +materialproperties(modeltype,20)*(vl_mat(is,i)*vl_mat(is,j))
     &                    +materialproperties(modeltype,21)*(vd_mat(is,i)*vd_mat(is,j))
                      end if
                  enddo
              end do
!-------------Symmetrize  for vectorization
              smdSMi(is,:,:) = (smdMi(is,:,:) + transpose(smdMi(is,:,:)))/2.0
              smnSMi(is,:,:) = (smnMi(is,:,:) + transpose(smnMi(is,:,:)))/2.0
              n_smdMi(is,:,:) = (n_smdMi(is,:,:) + transpose(n_smdMi(is,:,:)))/2.0
!-------------Vectorize
              call icams_conv33to6(smdSMi(is,:,:),smdVi1(is,:))
              call icams_conv33to6(smnSMi(is,:,:),smnVi1(is,:))
              call icams_conv33to6(n_smdMi(is,:,:),n_smdVi1(is,:))
!-------------Doubling the off-diagonal entries
              smdVi2(is,1:3)=smdVi1(is,1:3) 
              smdVi2(is,4:6)=smdVi1(is,4:6)*2.0
              smnVi2(is,1:3)=smnVi1(is,1:3) 
              smnVi2(is,4:6)=smnVi1(is,4:6)*2.0
              n_smdVi2(is,1:3) = n_smdVi1(is,1:3)
              n_smdVi2(is,4:6) = n_smdVi1(is,4:6)*2.0
          end do

!---------Assemble hardening matrix
          do i=1,nslip
              do j=1,nslip
                  if(sum(abs(vn_mat(i,:)-vn_mat(j,:)))<1.d-10) then
                      crsF_mat(i,j) = 1.0
                  else
                      crsF_mat(i,j) = 1.4
                  endif
              enddo
          enddo
      end subroutine AssembleLattice
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @author 
!> Mohamed Sharaf, Martin Boeff, Philipp Engels, Anxin Ma
! DESCRIPTION: 
!> Stress calculation
!> @param[in]  nslip                   Number of slip systems
!> @param[in]  nprops                  number of material parameters             
!> @param[in]  modeltype               Number of material model
!> @param[in]  props                   Material parameters
!> @param[in]  eang00(4)               Initial Euler angles
!> @param[in]  eang0(4)                Euler angles from time step \f$t_n\f$
!> @param[in]  crsF_mat                Hardening matrix q_alpha_beta (Roters p.50 eq.6.10)
!> @param[in]  IVB                     Critical resolved shear stress \f$t_{n}\f$
!> @param[in]  tau                     Shear stress at \f$t_{n}\f$
!> @param[in]  Fg                      Plastic defomation gradient at \f$t_{n}\f$
!> @param[in]  STFei26                 Elastic stiffness matrix
!> @param[in]  smdMi                   List of Schmidt Matrices for every glide system
!> @param[in]  smdVi2                  List of Schmidt Matrices for every glide system saved as Vector
!> @param[in]  smnVi2                  List of Matrices of normal to slip planes saved as Vector
!> @param[in]  pk2i_0                  2.Piola Kirchhoff stress at \f$t_{n}\f$
!> @param[in]  csM0                    Cauchy stress at \f$t_{n}\f$
!> @param[in]  modeltype               Material Model: 1-dislocation based ,2-pheno,3-phenoBCC
!> @param[in]  IVB0                    Critical resolved shear stress at \f$t_{n}\f$
!> @param[in]  IIVB                    Hardening due to first gradient of FP
!> @param[in]  BIVB                    Hardening due to second gradient of FP
!> @param[inout] p_slip                Accumulated plastic slip
!> @param[inout] p_fs                  Fatemi-Socie-inspired parameter at the 6 unparallel slip planes of BCC
!> @param[out] cs                      Cauchy stress vector at \f$t_{n+1}\f$
!> @param[out] pk2i                    2.Piola Kirchhoff stress at \f$t_{n+1}\f$
!> @param[out] IVB                     Critical resolved shear stress at \f$t_{n+1}\f$
!> @param[out] tau                     Resolved shear stress at \f$t_{n+1}\f$
!> @param[out] sigma_n                 Resolved normal stress at \f$t_{n+1}\f$
!> @param[out] Fe                      Elastic deformation gradient at \f$t_{n+1}\f$
!> @param[out] Fp                      Plastic deformation gradient at \f$t_{n+1}\f$
!> @param[out] Fp                      Plastic deformation gradient at \f$t_{n+1}\f$
!> @param[out] eang                    Euler angles at \f$t_{n+1}\f$
!> @param[out] MatJacb(6,6)            Material Jacobian at \f$t_{n+1}\f$
!> @param[out] ctrl_warningflag        Numerical error flag
!---------------------------------------------------------------------------
      subroutine cal_stress_mul(nslip, nprops, props,eang00,eang0, crsF_mat,dt1,
     &    Fp0, Fg, STFei26, smdMi, smdVi2, smnVi2, n_smdVi2,pk2i_0, csM0, IVB0, grainID, 
     &    IIVB,backstressn,cs,pk2i,IVB,tau,sigma_n,Fe,Fp,eang,MatJacb, modeltype,ctrl_warningflag,
     &    p_slip,p_fs,dgmdt,backstress,gam,p_irr,statev)
          use GlobalVariables, only: ctrl_isexplicit, ctrl_toler_NRloop, ctrl_NRmaxIterations, ctrl_infinity
          use MaterialData, only: materialproperties, grainDia
          use Models, only: cal_flow_harden_phenom, cal_flow_harden_dislocation
!--------------------------------------------------------------------------
          implicit none
!---------Start Material properties
          integer, intent(in) ::  nslip                   !Number of slip systems
          integer, intent(in) ::  modeltype                 !Material Model: 1-dislocation based ,2-pheno,3-phenoBCC
          integer, intent(in) ::  nprops                    !Number of material parameters             
          real(8), intent(in) ::  props(nprops)             !Material parameters
          real(8), intent(in) ::  eang00(4)                 !Initial Euler angles
          real(8), intent(in) ::  eang0(4)                  !Euler angles from time step tn
          real(8), intent(in) ::  statev(115)
!--------------------------------------------------------------------------
!---------End Material properties
          real(8), intent(in) ::  dt1                       !Time step    
          real(8), intent(in) ::  crsF_mat(nslip,nslip)     !Hardening matrix
          real(8), intent(in) ::  Fp0(3,3)                  !Plastic deformation gradient at tn
          real(8), intent(in) ::  Fg(3,3)                   !Total deformation gradient
          real(8), intent(in) ::  STFei26(6,6)              !Elastic stiffness matrix
          real(8), intent(in) ::  smdMi(nslip,3,3)          !Schmidt matrix
          real(8), intent(in) ::  smdVi2(nslip,6)           !Schmidt matrix saved as vector
          real(8), intent(in) ::  smnVi2(nslip,6)           !Matrix of normal to slip planes saved as a vector
          real(8), intent(in) ::  n_smdVi2(nslip,6)         !Non-Schmidt matrix saved as vector
          real(8), intent(in) ::  pk2i_0(6)                 !2.Piola Kirchhoff at tn
          real(8), intent(in) ::  csM0(3,3)                 !Cauchy stress matrix
          real(8), intent(in) ::  IVB0(nslip)               !Critical resolved shear stress at tn
          integer, intent(in) ::  grainID                   !Grain ID
          real(8), intent(in) ::  IIVB(nslip)               !Hardening due to first gradient of FP
          real(8), intent(in) ::  backstressn(nslip)        !backstress at time tn
!---------real(8), intent(in) ::  BIVB(nslip)               !Hardening due to second gradient of FP
          real(8), intent(inout) :: p_slip                  !Accumulated plastic slip
          real(8), intent(inout) :: gam(nslip)              !Accumulated shear
          real(8), intent(inout) :: p_fs(nslip/2)           !Fatemi-Socie-inspired parameter at the 6 unparallel slip planes of BCC
!---------Out
          real(8), intent(out) :: cs(6)                     !Cauchy stress vector at tn+1
          real(8), intent(out) :: pk2i(6)                   !2.Piola Kirchhoff stress
          real(8), intent(out) :: IVB(nslip)                !Critical resolved shear stress at tn+1
          real(8), intent(out) :: tau(nslip)                !Resolved shear stress at tn+1    
          real(8), intent(out) :: sigma_n(nslip)            !Resolved normal stress at tn+1
          real(8), intent(out) :: Fe(3,3)                   !Elastic deformation gradient at tn+1
          real(8), intent(out) :: Fp(3,3)                   !Plastic deformation gradient at tn+1
          real(8), intent(out) :: eang(4)                   !Euler angles at tn+1
          real(8), intent(out) :: MatJacb(6,6)              !Material Jacobian at tn+1
          integer, intent(out) :: ctrl_warningflag          !Numerical error flag
          real(8), intent(out) :: dgmdt(nslip)              !Shear rate \f$\dot{\gamma}^\alpha\f$
          real(8), intent(out) :: backstress(nslip)         !backstress at time tn+1
          real(8), intent(out) :: p_irr(nslip)              !Shear irreversibility
!------------------------------------------------------------------------------
!---------Local defined variables
          integer :: i,j,k,l,m,n,i1,j1,k1,l1,m1,n1,is,js      !Loop variables
          integer :: iNRloop                                  !Loop integer newton raphson loop
          integer :: Iexp_abq                                 !Flag for Umat (=0) or Vumat (=1)
          integer :: ctrl_tempwarningflag                     !Flag for inversion
          real(8) :: det_Fg, det_Fe                           !Determinats of deformation gradient
          real(8) :: x1,x2,x3,x4,y1,z1                        !Help functions
          real(8) :: TFg(3,3)                                 !Transposed deformation gradient
          real(8) :: CGE(3,3)                                 !Cauchy-Green strain tensor
          real(8) :: CGEe_max(3,3)                            !Trial Cauchy-Green strain tensor
          real(8) :: IFp0(3,3), TIFp0(3,3), IFg(3,3)          !Plastic deformation gradient at tn (I-inversion, T-transposed)
          real(8) :: MX1(3,3),MX2(3,3),MX3(3,3),MX4(3,3)      !Help functions for shorter formulas
          real(8) :: STFrlx6(nslip,6)                         !
          real(8) :: pk2i_max(6)                              !Piola Kirchhoff stress with plastic deformation gradient from tn!?
          real(8) :: tau_nsd(nslip)                           !non-Schmidt resolved shear stress at tn+1
!---------Possible Module values
          real(8) :: XI33(3,3)                                !Unit Matrices
          integer :: ib1(9),ib2(9)                            !Abaqus loop variables
!---------Start cal flow harden parameters
          real(8) :: dIVBdt(nslip)                            !Critical resolved shear stress \f$\tau_C\f$ with respect to dt
!---------real(8) :: dgmdt(nslip)                             !Shear rate \f$\dot{\gamma}^\alpha\f$ 
          real(8) :: ddgmdt_dtau(nslip)                      !Derivatives of dgmdt w.r.t. dtau
          real(8) :: ddgmdt_dIVB(nslip)                  !Derivatives of dgmdt w.r.t. dIVBdt
          real(8) :: ddIVBdt_ddgmdt(nslip,nslip)              !Derivatives of dIVBdt w.r.t. dgmdt
          real(8) :: ddIVBdt_dIVB(nslip,nslip)                !Derivatives of dIVBdt w.r.t. dIVBdt
!---------Start Newton iteration variables
          real(8) :: dpk2i(6)                                 !Increment in 2.Piola Kirchhoff stress
          real(8) :: dIVB(nslip)                              !Increment in Critical resolved shear stress
          real(8) :: ddgmdt_dpk2i(nslip,6)                    !Derivative of gammaDot w.r.t. 2.PK-Stress in Intermediate
          real(8) :: ddIVBdt_dpk2i(nslip,6)                   !Derivative of c.r.s.s w.r.t. 2.PK-Stress in Intermediate
          real(8) :: GV1(6)                                   !1. Equation in Newton iteration
          real(8) :: GV2(nslip)                               !2. Equation in Newton iteration
!---------Values from NR-Subroutine
          real(8) :: refv_pk2i                                !
          real(8) :: refv_IVB                                 !
          real(8) :: eqM6nGv1(6,nslip)
          real(8) :: eqM66Gv1(6,6)
          real(8) :: IeqM66Gv1(6,6)
          real(8) :: eqMn6Gv2(nslip,6)
          real(8) :: IeqMnnGv2(nslip,nslip)
!---------Updated values after Newton Iteration
          real(8) :: Lp(3,3)                                  !Plastic velocity gradient
          real(8) :: DPv(6)                                   !
          real(8) :: det_Fp, IFp(3,3),TIFp(3,3)               !Plastic deformation gradient tn+1
          real(8) :: pk2i_M(3,3)                              !Piola Kirchhoff stress saved as Matrix
          real(8) :: csM(3,3)                                 !Cauchy Stress at tn+1
          real(8) :: pk2r_M(3,3), pk2r(6)                     !2.Piola Kirchhoff stress in Referenze config
          real(8) :: dgam(nslip)                              !shear increment
!---------Material Stiffness
          real(8) :: idGv2_dIVB(nslip,nslip)                  !
          real(8) :: dGv2_dpk2i(nslip,6)                      !
          real(8) :: Gmod                                     !Gmod
          real(8) :: dbackstressdt(nslip)                     !
          real(8) :: gam_c                                    !gam_c
          real(8) :: pw_irr                                   !pw_irr
!---------------------------------------------------------------------------  
!---------Get material properties
          Gmod = 1./( 1./materialproperties(modeltype,3) +2./3.*(2./(materialproperties(modeltype,1)
     &           - materialproperties(modeltype,2))) - 1./materialproperties(modeltype,3))
            gam_c = materialproperties(modeltype,14)
            pw_irr = materialproperties(modeltype,15)
!---------write(*,*)'in stress statev(116)==',statev(116)
!---------write(*,*)'statev(114)==',statev(114)
          
!---------------------------------------------------------------------------
!---------Assemble Unit Matrices
          call IdentityMatrix(3,XI33)
!---------Abaqus stress vector sequence
          ib1=[1,2,3,1,1,2,2,3,3]
          ib2=[1,2,3,2,3,3,1,1,2]
          backstress = 0.0
!---------------------------------------------------------------------------
          if ((modeltype .NE. 3) .AND. (modeltype .NE. 4)) then !non-Schmidt behaviour
              tau_nsd=0.0
          endif
!---------------------------------------------------------------------------
!---------Check distortion of Fg
          ctrl_warningflag=0
          call get_determinant(Fg,det_Fg)
          if(det_Fg < 1.d-10) then
              ctrl_warningflag=20 
              return
          endif

!---------Determinant of elastic def.grad. equals det of total def.grad
          det_Fe=det_Fg
!---------Transpose Deformation gradient
          TFg=transpose(Fg)
!---------Calculate Cauchy-Green strain
          CGE=matmul(TFg,Fg)
!---------------------------------------------------------------------------
!---------Check if Fp0 is invertable
          ctrl_tempwarningflag=0
!---------Inverse of deformation gradient
          call gaussj(Fp0,3,IFp0,ctrl_tempwarningflag)
          if(ctrl_tempwarningflag/=0)then
              ctrl_warningflag=1  !Fp0 is non-invertible
              return
          endif
!---------Transpose Inverse
          TIFp0=transpose(IFp0)
!---------------------------------------------------------------------------
!---------Calculate trial Piola Kirchhoff stress by pk2i_max(6)=C*(CGE_max(6)-I(6))/2 

!---------Maximum Cauchy strain
          CGEe_max = matmul( matmul(TIFp0,CGE),IFp0 )
          MX1 = (CGEe_max-XI33)/2

!---------Assemble trial piola Kirchhoff stress
          do i=1,6
              pk2i_max(i)=0
              do j=1,6
                  pk2i_max(i) = pk2i_max(i) + STFei26(i,j) * MX1(ib1(j),ib2(j))
              enddo
          enddo
!---------------------------------------------------------------------------
!---------Calculate STFrlx6=C*( CGE_max*smdMi + (CGE_max*smdMi)^T )/2    
          do is=1,nslip
              MX1=matmul(CGEe_max,smdMi(is,:,:))
              MX2=(MX1+transpose(MX1))/2
              do i=1,6
                  STFrlx6(is,i)=0
                  do j=1,6
                      STFrlx6(is,i) = STFrlx6(is,i) + STFei26(i,j)*MX2(ib1(j),ib2(j))
                  enddo
              enddo
          enddo
!---------------------------------------------------------------------------
!---------Calculate shear irreversibility
          p_irr = (gam/gam_c)**pw_irr
          p_irr = 1.-p_irr
!---------------------------------------------------------------------------
!---------Begin newton raphson method to solve pk2i, IVB
!---------Initialite pk2i and IVB from tn
          pk2i = pk2i_0
          IVB  = IVB0
          do iNRloop=1,ctrl_NRmaxIterations
!-------------Calculate the actual resolved shear and normal stresses (intital guess with pk2i from tn?!)
              do is=1,nslip
                  tau(is) = dot_product(pk2i,smdVi2(is,:))
!-----------------tau(is) = dot_product(pk2i,smdVi2(is,:))/p_irr(is)
                  sigma_n(is)= dot_product(pk2i,smnVi2(is,:))
              enddo
              if (modeltype==3 .or. modeltype == 4) then !non-Schmidt behaviour
                  do is=1,nslip
                      tau_nsd(is) = dot_product(pk2i,n_smdVi2(is,:))
                  enddo
              end if
!-------------Set initial error flag for Function cal_flow_harden to zero        
              ctrl_tempwarningflag=0
!-------------Call derivatives used for newton iteration from cal_flow_harden
!-------------Distinguish between modeltypes 
              if (modeltype /= 1) then                
                  call cal_flow_harden_phenom(nslip,nprops,props,IVB, grainID, tau, tau_nsd, crsF_mat,
     &                                       IIVB,backstressn,modeltype, dgmdt, ddgmdt_dtau,ddgmdt_dIVB,
     &                                        dIVBdt,ddIVBdt_ddgmdt, ddIVBdt_dIVB,dbackstressdt,
     &                                        ctrl_tempwarningflag,statev)
              refv_IVB  = Gmod * 1.d-6 ! residuum nominator for GV2
              else
                  call cal_flow_harden_dislocation(nslip,nprops,props,IVB, tau,modeltype, dgmdt,
     &                                             ddgmdt_dtau,ddgmdt_dIVB,dIVBdt,
     &                                             ddIVBdt_ddgmdt,ddIVBdt_dIVB,ctrl_tempwarningflag)
!-----------------residuum nominator for GV2
                  refv_IVB = materialproperties(modeltype,8)+materialproperties(modeltype,9)/sqrt(grainDia(grainID))
              endif
!-------------Transfer cal_flow_harden error flag to global error flag
              if(ctrl_tempwarningflag/=0)then
                  ctrl_warningflag=ctrl_tempwarningflag  !!! shearrate non-cvg,stop
                  return
              endif
!-------------Calculate derivative of gammaDot w.r.t. 2.Piola-Kirchhoff stress in intermediate config
              do is=1,nslip
                  ddgmdt_dpk2i(is,:) = ddgmdt_dtau(is) * smdVi2(is,:)
              enddo
!-------------Calculate derivative of c.r.s.s. w.r.t. 2.Piola-Kirchhoff stress in intermediate config
              ddIVBdt_dpk2i(1:nslip,:) = matmul(ddIVBdt_ddgmdt(1:nslip,1:nslip),ddgmdt_dpk2i(1:nslip,:))
!-------------Assemble GV1
              Gv1 = +pk2i - pk2i_max + matmul(transpose(STFrlx6(1:nslip,:)),dgmdt(1:nslip)) * dt1
!-------------Assemble GV2
              Gv2(1:nslip) = +IVB(1:nslip) - IVB0(1:nslip) - dIVBdt(1:nslip)*dt1
!-------------Toleranz criteriums for GV1 
              refv_pk2i = Gmod * 1.d-6
              x1 = sum(abs(Gv1(1:6))) / refv_pk2i
              x2 = sum(abs(Gv2(1:nslip))) / refv_IVB
!-------------Convergence check
              if( (x1+x2) < ctrl_toler_NRloop .or. ctrl_isexplicit==1 ) then
!-----------------Update Piola Kirchhoff stress
                  pk2i = +pk2i_max - matmul(transpose(STFrlx6(1:nslip,:)),dgmdt(1:nslip))*dt1
!-----------------Update Critical resolved shear stress
                  IVB(1:nslip) = IVB0(1:nslip) + dIVBdt(1:nslip) * dt1
                  do is=1,nslip
                      if (IVB(is) > 1.d4 .or. IVB(is)<=0) then
                          ctrl_warningflag=112  !!! shearrate non-cvg,stop
                          return
                      endif
                  enddo
                  goto 101 !Converge jump out NR loop and continue
              endif
!-------------Call Newton Raphson Matrix
              ctrl_tempwarningflag=0
              call cal_NRmatrix(nslip, dt1, STFrlx6, ddgmdt_dpk2i, ddgmdt_dIVB,ddIVBdt_dpk2i, ddIVBdt_dIVB,
     &                          eqM6nGv1, eqM66Gv1, IeqM66Gv1, eqMn6Gv2, IeqMnnGv2, idGv2_dIVB,dGv2_dpk2i, 
     &                          ctrl_tempwarningflag)
              if(ctrl_tempwarningflag/=0)then
                  ctrl_warningflag=ctrl_tempwarningflag  !NR matrix cannot be calculated,stop
                  return
              endif
!-------------Calculate increment in 2.Piola Kirchhoff stress and c.r.s.s.
              dpk2i = matmul(IeqM66Gv1,-(Gv1-matmul(eqM6nGv1(:,1:nslip),Gv2(1:nslip))))
              dIVB(1:nslip) = matmul(IeqMnnGv2(1:nslip,1:nslip),-(Gv2(1:nslip)-matmul(eqMn6Gv2(1:nslip,:),Gv1)))
!-------------Update 2.Piola Kirchhoff stress and critical resolved shear stress during Newton Iteration 
              pk2i           = pk2i + dpk2i 
              IVB(1:nslip)   = IVB(1:nslip) + dIVB(1:nslip)
!-------------do is=1,nslip
!---------if (dbackstressdt(is).GT.1) then
!---------print*, dbackstressdt(is)*dt1
!---------endif
!---------enddo
          enddo
          ctrl_warningflag=5  !Non-covergence for ctrl_NRmaxIterations
          return
101       continue
!---------------------------------------------------------------------------
!---------Restore shear irreversibility value
          p_irr = 1.-p_irr 
!---------Solve cauchy stress    
!---------Convert 2.Pki-Vector to matrix
          call icams_conv6to33(pk2i,pk2i_M)
!---------Calculate plastic velocity gradient
          do i=1,3
              do j=1,3
                  Lp(i,j)=dot_product(dgmdt(:),smdMi(:,i,j))
              enddo
          enddo

          MX1=matmul(Lp,Fp0)
          MX2=(MX1+transpose(MX1))/2

          do i=1,6
!----------------DPv(i)=dot_product(dgmdt(:),smdVi1(:,i))
                 DPv(i)=MX2(ib1(i),ib2(i))
          enddo
!---------Calculate plastic deformation gradient for tn+1
          Fp=matmul(XI33+Lp*dt1,Fp0)
          call get_determinant(Fp,det_Fp)
          if(det_Fp==0)then
              ctrl_warningflag=22  !Fp is non-invertible,stop current time step
              return
          endif
          Fp=Fp/det_Fp**(1.0/3.0)
!---------Set Determinante of Fp to one
          det_Fp=1.0
          ctrl_tempwarningflag=0
          call gaussj(Fp,3,IFp,ctrl_tempwarningflag)
          if(ctrl_tempwarningflag/=0)then
              ctrl_warningflag=2  !Fp is non-invertible,stop current time step
              return
          endif
          TIFp=transpose(IFp)
          Fe=matmul(Fg,IFp)
          call get_determinant(Fe,det_Fe)
!---------Calculate Cauchy stress in matrix notation
          csM=matmul(matmul(Fe,pk2i_M),transpose(Fe))/det_Fe
!---------Convert Cauchy Stress-Matrix to Vector
          call icams_conv33to6(csM,cs)
          ctrl_tempwarningflag=0
          call gaussj(Fg,3,IFg,ctrl_tempwarningflag)
          if(ctrl_tempwarningflag/=0)then
              ctrl_warningflag=3  !Fg is non-invertible,stop current time step
              return
          endif
!---------Calculate Piola Kirchhoff stress in referenze configuration
          pk2r_M=matmul(matmul(IFg,csM),transpose(IFg))*det_Fe
          call icams_conv33to6(pk2r_M,pk2r)
!---------shear rate at tn+1
          dgam = dgmdt * dt1
!---------------------------------------------------------------------------
!---------Calculate Euler angles and misorientation
!---------Euler angles at tn
          eang(1:3) = eang0(1:3)
!---------Calculate Euler angles based on elastic deformation gradient
          call caleulang(Fe,eang(1:3),ctrl_tempwarningflag)
          if(ctrl_tempwarningflag/=0) then
              eang(1:3) = eang0(1:3)
          end if
!---------Calculate missorientation based on euler angles
          call icams_misori(eang00(1:3),eang(1:3),eang(4))
!---------------------------------------------------------------------------
!---------Calculate accumulated plastic slip       
          call cal_pl_slip(Lp,dt1,p_slip)
          gam = gam + abs(dgam)
          do is=1,nslip
              if (abs(gam(is)) > ctrl_infinity) then !! infinity check
                  ctrl_warningflag = 119
                  return
              endif
          enddo
          backstress = backstressn + dbackstressdt * dt1
!---------------------------------------------------------------------------
!---------Calculate Fatemi-Socie-inspired parameter)
          p_fs(1)=p_fs(1)+(1.+0.4*sigma_n(1)/0.000355)*(abs(dgam(1))+abs(dgam( 7)))
          p_fs(2)=p_fs(2)+(1.+0.4*sigma_n(2)/0.000355)*(abs(dgam(2))+abs(dgam( 5)))
          p_fs(3)=p_fs(3)+(1.+0.4*sigma_n(3)/0.000355)*(abs(dgam(3))+abs(dgam(12)))
          p_fs(4)=p_fs(4)+(1.+0.4*sigma_n(4)/0.000355)*(abs(dgam(4))+abs(dgam(10)))
          p_fs(5)=p_fs(5)+(1.+0.4*sigma_n(6)/0.000355)*(abs(dgam(6))+abs(dgam( 9)))
          p_fs(6)=p_fs(6)+(1.+0.4*sigma_n(8)/0.000355)*(abs(dgam(8))+abs(dgam(11)))
!---------------------------------------------------------------------------
!---------calculate material stiffness
!---------ctrl_isexplicit=0                         !VUmat

          if(ctrl_isexplicit==0)then
              ctrl_tempwarningflag=0
!-------------Call NRmatrix again to get values at tn+1
              call cal_NRmatrix(nslip, dt1, STFrlx6, ddgmdt_dpk2i, ddgmdt_dIVB,ddIVBdt_dpk2i, 
     &                          ddIVBdt_dIVB, eqM6nGv1, eqM66Gv1, IeqM66Gv1, eqMn6Gv2, 
     &                          IeqMnnGv2, idGv2_dIVB, dGv2_dpk2i, ctrl_tempwarningflag)
              if(ctrl_tempwarningflag/=0)then
                  ctrl_warningflag=420+ctrl_tempwarningflag  !NR matrix cannot be got,stop
                  return
              endif
!-------------Calculate Material Jacobian
              call cal_MatStiffness(nslip,dt1,IdGv2_dIVB,dGv2_dpk2i,Fp0,TIFp0,IFp0,TIFp,IFp,Fg,
     &                              det_Fg, STFei26,smdMi,IeqM66Gv1, ddgmdt_dpk2i,ddgmdt_dIVB,
     &                              pk2i_M,csM0,dgmdt,MatJacb)
          endif
          return
      end subroutine cal_stress_mul
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @author 
!> Mohamed Sharaf, Martin Boeff, Anxin Ma
! DESCRIPTION: 
!> Newton Raphson Iteration
!> @param[in]  nslip                          Number of slip systems
!> @param[in]  dt1                            Time step
!> @param[in]  STFrlx6                        
!> @param[in]  ddgmdt_dpk2i                   Derivatives of dgmdt w.r.t. PK2
!> @param[in]  ddgmdt_dIVB                    Derivatives of dgmdt w.r.t. IVB
!> @param[in]  ddIVBdt_dpk2i                  Derivatives of dIVBdt w.r.t. PK2
!> @param[in]  ddIVBdt_dIVB                   Derivatives of dIVBdt w.r.t. IVB
!> @param[out] eqM6nGv1(6,nslip)              
!> @param[out] eqM66Gv1(6,6)                  
!> @param[out] IeqM66Gv1(6,6)                 
!> @param[out] eqMn6Gv2(nslip,6)              
!> @param[out] IeqMnnGv2(nslip,nslip)         
!> @param[out] idGv2_dIVB(nslip,nslip)        
!> @param[out] ctrl_warningflag               Numerical error flag
!---------------------------------------------------------------------------  
      subroutine cal_NRmatrix(nslip, dt1, STFrlx6, ddgmdt_dpk2i, ddgmdt_dIVB,ddIVBdt_dpk2i, 
     &                        ddIVBdt_dIVB, eqM6nGv1, eqM66Gv1, IeqM66Gv1, eqMn6Gv2, IeqMnnGv2, 
     &                        idGv2_dIVB,dGv2_dpk2i, ctrl_warningflag)
          implicit none
!---------Start Material properties
          integer, intent(in)::  nslip                       !Number of slip systems     
          real(8), intent(in)::  dt1                         !Time step
          real(8) ::  STFrlx6(nslip,6)            !
          real(8), intent(in)::  ddgmdt_dpk2i(nslip,6)       !Derivatives of dgmdt w.r.t. PK2
          real(8), intent(in)::  ddgmdt_dIVB(nslip)          !Derivatives of dgmdt w.r.t. IVB
          real(8), intent(in)::  ddIVBdt_dpk2i(nslip,6)      !Derivatives of dIVBdt w.r.t. PK2
          real(8), intent(in)::  ddIVBdt_dIVB(nslip,nslip)   !Derivatives of dIVBdt w.r.t. IVB
          !Out
          real(8), intent(out):: eqM6nGv1(6,nslip)           !
          real(8), intent(out):: eqM66Gv1(6,6)               !
          real(8), intent(out):: IeqM66Gv1(6,6)              !
          real(8), intent(out):: eqMn6Gv2(nslip,6)           !
          real(8), intent(out):: IeqMnnGv2(nslip,nslip)      !
          real(8), intent(out):: idGv2_dIVB(nslip,nslip)     !Inverse of derivative of second residuum w.r.t. to IVB
          real(8), intent(out):: dGv2_dpk2i(nslip,6)         !Derivative of second residuum w.r.t. to 2PK
          integer, intent(out):: ctrl_warningflag            !Numerical error flag
!------------------------------------------------------------------------------  
!---------Local parameters
          integer i,is                      !Loop variables
          integer ctrl_tempwarningflag      !Numerical error flag
          real(8) XI66(6,6)                 !Identity matrix
          real(8) XInn(nslip,nslip)         !
          real(8) dGv1_dpk2i(6,6)           !Derivative of first residuum w.r.t. to IVB
          real(8) dGv1_dIVB(6,nslip)        !Derivative of first residuum w.r.t. to IVB
          real(8) idGv1_dpk2i(6,6)          !Inverse of derivative of first residuum w.r.t. to IVB
          real(8) dGv2_dIVB(nslip,nslip)    !Derivative of second residuum w.r.t. to IVB
          real(8) eqMnnGv2(nslip,nslip)     !
!------------------------------------------------------------------------------  
          call IdentityMatrix(6,XI66)
          call IdentityMatrix(nslip,XInn)
!------------------------------------------------------------------------------  
!---------calculate dG1_dpk2i(6,6),dG1_dIVB(6,:),IdG1_dpk2i(6,6)
          ctrl_warningflag=0
          dGv1_dpk2i = XI66 + matmul(transpose(STFrlx6(1:nslip,:)),ddgmdt_dpk2i(1:nslip,:))*dt1
          do is=1,nslip
              dGv1_dIVB(:,is) = STFrlx6(is,:)*ddgmdt_dIVB(is)*dt1
          enddo
!---------Inversion
          ctrl_tempwarningflag=0
          call gaussj(dGv1_dpk2i,6,idGv1_dpk2i,ctrl_tempwarningflag)
          if(ctrl_tempwarningflag/=0)then
              ctrl_warningflag=41  !dGv1_dpk2i non-invertable,stop
              return
          endif
!------------------------------------------------------------------------------  
!---------Calculate dG2_dpk2i(:,6),dG2_dIVB(:,:),IdG2_dIVB(:,:)
          
          dGv2_dpk2i(1:nslip,:) = -ddIVBdt_dpk2i(1:nslip,:) * dt1
          dGv2_dIVB(1:nslip,1:nslip) = +XInn(1:nslip,1:nslip) - ddIVBdt_dIVB(1:nslip,1:nslip) * dt1
          
          ctrl_tempwarningflag=0
          call gaussj(dGv2_dIVB(1:nslip,1:nslip), nslip, idGv2_dIVB(1:nslip,1:nslip), ctrl_tempwarningflag)
          if(ctrl_tempwarningflag/=0)then
              ctrl_warningflag=42  !dGv2_dIVB non-invertable,stop
              return
          endif
!------------------------------------------------------------------------------  
!---------Calculate eqM66Gv1,eqMnnGv2,IeqM66Gv1,IeqMnnGv2
          eqM6nGv1 = matmul(dGv1_dIVB(:,1:nslip),IdGv2_dIVB(1:nslip,1:nslip))
          eqM66Gv1 = -matmul(matmul(dGv1_dIVB(:,1:nslip),IdGv2_dIVB(1:nslip,1:nslip)),dGv2_dpk2i(1:nslip,:)) + dGv1_dpk2i
          ctrl_tempwarningflag=0
          call gaussj(eqM66Gv1,6,IeqM66Gv1,ctrl_tempwarningflag)
          if(ctrl_tempwarningflag/=0)then
              ctrl_warningflag=43  !eqM66Gv1 non-invertable,stop
              return
          endif
          eqMn6Gv2(1:nslip,:)=matmul(dGv2_dpk2i(1:nslip,:),IdGv1_dpk2i)
          eqMnnGv2(1:nslip,1:nslip)=-matmul(matmul(dGv2_dpk2i(1:nslip,:),IdGv1_dpk2i),
     &                              dGv1_dIVB(:,1:nslip))+dGv2_dIVB(1:nslip,1:nslip)
          ctrl_tempwarningflag=0
          call gaussj(eqMnnGv2(1:nslip,1:nslip),nslip,IeqMnnGv2(1:nslip,1:nslip),ctrl_tempwarningflag)
          if(ctrl_tempwarningflag/=0)then
              ctrl_warningflag=44  !eqMnnGv2 non-invertable,stop
              return
          endif
          return
      endsubroutine cal_NRmatrix

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @author 
!> Martin Boeff, Anxin Ma
! DESCRIPTION: 
!> THIS ROUTINE calculates material jacb for umat STFjc_66 \n
!> THIS ROUTINE also can calculates material jacb for uel STFtk_66
!> @param[in]  nslip                          Number of slip systems
!> @param[in]  dt1                            !Time increment
!> @param[in]  IdGv2_dIVB                     Inverse of derivative of second residuum w.r.t. to IVB
!> @param[in]  dGv2_dpk2i                     Derivative of second residuum w.r.t. to 2PK
!> @param[in]  Fp0(3,3)                       Plastic deformation gradient \f$F^p_n\f$
!> @param[in]  TIFp0(3,3)                     Tranpose of inverse of pl. deformation gradient \f$F_n^{-T}\f$
!> @param[in]  IFp0(3,3)                      Inverse of pl. deformation gradient \f$F_{p,n}^{-1}\f$
!> @param[in]  TIFp                           Tranpose of inverse of pl. deformation gradient \f$F_{p,n+1}^{-T}\f$
!> @param[in]  IFp(3,3)                       Inverse of  pl. deformation gradient \f$F_{p,n+1}^{-1}\f$
!> @param[in]  Fg(3,3)                        Deformation gradient \f$F_{n+1}\f$
!> @param[in]  STFei26(6,6)                   Elastic stiffness matrix
!> @param[in]  smdMi(nslip,3,3)               Schmidt matrix
!> @param[in]  IeqM66Gv1(6,6)                 
!> @param[in]  ddgmdt_dpk2i(nslip,6)          Derivative of dgmdt w.r.t. to 2PK
!> @param[in]  ddgmdt_dIVB(nslip)             Derivative of dgmdt w.r.t. to IVB
!> @param[in]  pk2i_M(3,3)                    2.Piola Kirchhoff stress intermediate at tn+1 Matrix
!> @param[in]  csM0(3,3)                      Cauchy stress saved as Matrix at \f$t_n\f$
!> @param[in]  dgmdt(nslip)                   Shear increment
!> @param[out] MatJacb(6,6)                   Material Jacobian matrix at \f$t_{n+1}\f$
!---------------------------------------------------------------------------  
      subroutine cal_MatStiffness(nslip, dt1, IdGv2_dIVB, dGv2_dpk2i, Fp0, TIFp0, IFp0, TIFp, IFp, 
     &                            Fg, det_Fg, STFei26, smdMi, IeqM66Gv1, ddgmdt_dpk2i, ddgmdt_dIVB, 
     &                            pk2i_M, csM0, dgmdt, MatJacb)
          implicit none
!---------Input
          integer, intent(in) :: nslip                            !Slip systems
          real(8), intent(in) :: dt1                              !Time increment
          real(8), intent(in) :: IdGv2_dIVB(nslip,nslip)          !Inverse of derivative of second residuum w.r.t. IVB
          real(8), intent(in) :: dGv2_dpk2i(nslip,6)              !Derivative of second residuum w.r.t. to 2PK
          real(8), intent(in) :: Fp0(3,3)                         !Plastic deformation gradient at t_n
          real(8), intent(in) :: TIFp0(3,3)                       !Tranpose of inverse of pl. deformation gradient at t_n
          real(8), intent(in) :: IFp0(3,3)                        !Inverse of pl. deformation gradient at t_n
          real(8), intent(in) :: TIFp(3,3)                        !Tranpose Inverse of deformation gradient at tn+1
          real(8), intent(in) :: IFp(3,3)                         !Inverse of deformation gradient at tn+1
          real(8), intent(in) :: Fg(3,3)                          !Deformation gradient tn+1
          real(8), intent(in) :: det_Fg                           !Determinante of deformation gradient
          real(8), intent(in) :: STFei26(6,6)                     !Elastic stiffness matrix
          real(8), intent(in) :: smdMi(nslip,3,3)                 !Schmidt matrix
          real(8), intent(in) :: IeqM66Gv1(6,6)
          real(8), intent(in) :: ddgmdt_dpk2i(nslip,6)            !Derivative of dgmdt w.r.t. to 2PK
          real(8), intent(in) :: ddgmdt_dIVB(nslip)               !Derivative of dgmdt w.r.t. to IVB
          real(8), intent(in) :: pk2i_M(3,3)                      !2. Piola Kirchhoff stress intermediate at tn+1 Matrix
          real(8), intent(in) :: csM0(3,3)                        !Cauchy stress saved as Matrix at tn
          real(8), intent(in) :: dgmdt(nslip)                     !Shear increment
!---------Out
          real(8), intent(out)::MatJacb(6,6)                      !Material Jacobian at tn+1
!------------------------------------------------------------------------------  
!---------Local variables
          integer  i,j,k,l,m,n,i1,j1,k1,l1,m1,n1,is,js           !Loop variables
          integer  ib1(9),ib2(9)                                 !Abaqus loop variables
          real(8)  x1,y1,z1
          real(8)  M1_66(6,6),M1_96(9,6),M1_99(9,9)
          real(8)  M2_66(6,6),M2_96(9,6),M2_99(9,9)
          real(8)  M3_66(6,6),M3_96(9,6),M3_99(9,9)
          real(8)  M1_3333(3,3,3,3),M2_3333(3,3,3,3),M3_3333(3,3,3,3)
          real(8)  MX1(3,3),MX2(3,3),MX3(3,3),MX4(3,3)
          real(8)  XI33(3,3)
          real(8)  dIVB_dpk2i(nslip,6)              !
          real(8)  dCGEe_mx_dE(6,6)                 !
          real(8)  dpk2i_mx_dE(6,6)                 !
          real(8)  dSTFrlx6_dE(nslip,6,6)           !
          real(8)  dGv1_dE(6,6)                     !
          real(8)  dTdgmdt_dpk2i(nslip,9)
          real(8)  dpk2i_dE(6,6)                    !
          real(8)  dIFp_dE(9,6)                     !
          real(8)  dTIFp_dE(9,6)                    !
          real(8)  dFp_dE(9,6)
          real(8)  ddgmdt_dE(nslip,6)               !
          real(8)  dpk2r_dE(6,6)                    !
          real(8)  STFjc_66(6,6)                    !
          real(8)  STFtk_66(6,6)                    !
!-----------------------------------------------------------------------------
!---------Abaqus stress vector sequence
          ib1=[1,2,3,1,1,2,2,3,3]
          ib2=[1,2,3,2,3,3,1,1,2]
!-----------------------------------------------------------------------------
          dIVB_dpk2i=-matmul(IdGv2_dIVB,dGv2_dpk2i)
          call IdentityMatrix(3,XI33)
!-----------------------------------------------------------------------------
!         jacb#1: dpk2i_dE
!         #1.2:   dGv1_dE
!         #1.2.1:   dCGEe_mx_dE(6,6) = d(Fe_mx^T*Fe_mx)_dE
!-----------------------------------------------------------------------------
          do i=1,6
              do j=1,6
                  if(j<=3)then
                      dCGEe_mx_dE(i,j) = 2*TIFp0(ib1(i),ib1(j)) * IFp0(ib2(j),ib2(i))
                  else
                      dCGEe_mx_dE(i,j) = (+2*TIFp0(ib1(i),ib1(j  ))*
     &                                   IFp0(ib2(j  ), ib2(i))+2*
     &                                   TIFp0(ib1(i),ib1(j+3))*IFp0(ib2(j+3),ib2(i)))/2
                  endif
              enddo
          enddo
!-----------------------------------------------------------------------------
!         jacb#1: dpk2i_dE
!         #1.2:   dGv1_dE
!         #1.2.2:   dpk2i_mx_dE(6,6)
!-----------------------------------------------------------------------------
          dpk2i_mx_dE=matmul( STFei26 , dCGEe_mx_dE )/2
!-----------------------------------------------------------------------------
!         jacb#1: dpk2i_dE
!         #1.2:   dGv1_dE
!         #1.2.3:   dSTFrlx6_dE(nslip,6,6)
!-----------------------------------------------------------------------------
          do is=1,nslip
              MX1=matmul(IFp0,smdMi(is,:,:))
              MX2=transpose(MX1)
              do i=1,6
                  do j=1,6
                      M1_66(i,j) = +2*TIFp0(ib1(i),ib1(j))*MX1 (ib2(j),ib2(j))+2*
     &                             MX2  (ib1(i),ib1(j))*IFp0(ib2(j),ib2(j))
                  enddo
              enddo
              dSTFrlx6_dE(is,:,:) = matmul(STFei26, M1_66)/2
          enddo
!-----------------------------------------------------------------------------
!         jacb#1: dpk2i_dE
!         #1.2:   dGv1_dE
!         #1.2.4:   dGv1_dE(6,6)
!-----------------------------------------------------------------------------
          do i=1,6
              do j=1,6
                  dGv1_dE(i,j)= -dpk2i_mx_dE(i,j) + dot_product(dgmdt,dSTFrlx6_dE(:,i,j))*dt1
              enddo
          enddo
!-----------------------------------------------------------------------------
!         jacb#1: dpk2i_dE
!         #1.3:   dpk2i_dE(6,6)
!-----------------------------------------------------------------------------
          dpk2i_dE = -matmul(IeqM66Gv1,dGv1_dE)
!-----------------------------------------------------------------------------
!         jacb#2: dpk2r_dE
!         #2.1:  dIFp_dE(9,6),dTIFp_dE(9,6)
!-----------------------------------------------------------------------------
          M1_99=0
          M2_99=0
          M3_99=0
          do is=1,nslip
              MX1=matmul(IFp0,smdMi(is,:,:))
              MX2=transpose(MX1)
              MX4=matmul(smdMi(is,:,:),Fp0)
              do i=1,9
                  if(i<=6) i1=i
                  if(i >6) i1=i-3
                  MX3(ib1(i),ib2(i)) = +ddgmdt_dpk2i(is,i1)+ddgmdt_dIVB(is)*dIVB_dpk2i(is,i1)
                  dTdgmdt_dpk2i(is,i) = +ddgmdt_dpk2i(is,i1)+ddgmdt_dIVB(is)*dIVB_dpk2i(is,i1)
              enddo
              do i=1,9
                  do j=1,9
                      M1_99(i,j) = M1_99(i,j)+MX1(ib1(i),ib2(i))*MX3(ib1(j),ib2(j))
                      M2_99(i,j) = M2_99(i,j)+MX2(ib1(i),ib2(i))*MX3(ib1(j),ib2(j))
                      M3_99(i,j) = M3_99(i,j)+MX4(ib1(i),ib2(i))*MX3(ib1(j),ib2(j))
                  enddo
              enddo
          enddo
          dIFp_dE  = 0
          dTIFp_dE = 0
          M1_96(1:6,:) = dpk2i_dE
          M1_96(7,:)=M1_96(4,:)
          M1_96(8,:)=M1_96(5,:)
          M1_96(9,:)=M1_96(6,:)
          dIFp_dE  = -matmul(M1_99,M1_96)*dt1
          dTIFp_dE = -matmul(M2_99,M1_96)*dt1
          dFp_dE   = -matmul(M3_99,M1_96)*dt1
          ddgmdt_dE=  matmul(dTdgmdt_dpk2i,M1_96)*dt1
          dFp_dE(:,4:6)    = dFp_dE(:,4:6)*2
          ddgmdt_dE(:,4:6) = ddgmdt_dE(:,4:6)*2
!-----------------------------------------------------------------------------
!         jacb#2: dpk2r_dE
!         #2.2: dpk2r_dE(6,6)
!-----------------------------------------------------------------------------
          dpk2r_dE=0
          MX1=matmul(pk2i_M,TIFp)
          MX2=matmul(IFp,pk2i_M)
          do i=1,6
              do k=1,6
                  do m=1,9
                      if(m<=6)m1=m
                      if(m >6)m1=m-3
                      dpk2r_dE(i,k)= dpk2r_dE(i,k) + XI33(ib1(i),ib1(m))*MX1(ib2(m),ib2(i))*
     &                               dIFp_dE(m ,k) + IFp(ib1(i), ib1(m))*TIFp(ib2(m),ib2(i))*
     &                               dpk2i_dE(m1,k) + MX2(ib1(i),ib1(m))*XI33(ib2(i),ib2(m))* dTIFp_dE(m,k)
                        enddo
                 enddo
          enddo
!-------------------------------------------------------------------------------------------
!         Calculate stiffness for Jaummann rate of Cauchy stress STFjc_66 for user material
!-------------------------------------------------------------------------------------------
          do i=1,9
                 do j=1,9
                        if(i<=6) i1=i
                        if(i >6) i1=i-3
                        if(j<=6) j1=j
                        if(j >6) j1=j-3
                        M1_3333(ib1(i),ib2(i),ib1(j),ib2(j))=dpk2r_dE(i1,j1)
                 enddo
          enddo
         
          M2_3333=0
          do i=1,3
             do j=1,3
                 do k=1,3
                     do l=1,3
                         x1=0
                         do i1=1,3
                             do j1=1,3
                                 do k1=1,3
                                     do l1=1,3
                                         x1=x1+M1_3333(i1,j1,k1,l1)*Fg(i,i1)*Fg(j,j1)*Fg(k,k1)*Fg(l,l1)
                                     enddo
                                  enddo
                              enddo
                          enddo
                      M2_3333(i,j,k,l)=x1/det_Fg+XI33(i,k)*csM0(l,j)+csM0(i,k)*XI33(l,j)+csM0(i,j)*XI33(k,l)
                      enddo
                  enddo
              enddo
          enddo
          do i=1,6
              do j=1,6
                  STFjc_66(I,J)=M2_3333(ib1(i),ib2(i),ib1(j),ib2(j))
              enddo
          enddo
!-----------------------------------------------------------------------------------------
!         Calculate stiffness for Trusdell rate of Kirchhoff stress STFtk_66 for user element
!-----------------------------------------------------------------------------------------
          do i=1,9
             do j=1,9
                 if(i<=6) i1=i
                 if(i >6) i1=i-3
                 if(j<=6) j1=j
                 if(j >6) j1=j-3
                 M1_3333(ib1(i),ib2(i),ib1(j),ib2(j)) = dpk2r_dE(i1,j1)
             enddo
          enddo
          M2_3333=0
          do i=1,3
              do j=1,3
                  do k=1,3
                      do l=1,3
                          do k1=1,3
                              do l1=1,3
                                  M2_3333(i,j,k,l)=M2_3333(i,j,k,l)+M1_3333(i,j,k1,l1)*Fg(k,k1)*Fg(l,l1)
                              enddo
                          enddo
                      enddo
                  enddo
              enddo
          enddo
          
          do i=1,6
              do j=1,6
                  STFtk_66(I,J)=M2_3333(ib1(i),ib2(i),ib1(j),ib2(j))
              enddo
          enddo
!-----------------------------------------------------------------------------
!         Save Material tangent matrix
!-----------------------------------------------------------------------------
          MatJacb=STFjc_66
          return
      endsubroutine cal_MatStiffness
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------------
!     Equivalent strain   
!-----------------------------------------------------------------------------  
      subroutine cal_eqpl_strain(Fp,peeq)
          implicit none
          real(8), intent(in)::  Fp(3,3)                         !Plastic deformation gradient
          real(8), intent(inout):: peeq                          !Plastic equivalent strain at tn1
          
!---------Local parameters
          real(8) Cp(3,3)                                        !Plastic Cauchy strain
          real(8) XI33(3,3)                                      !Identity matrix
          real(8) Ep(3,3) ,Ep_norm                               !Plastic Green-Lagrange strain
!---------Calculate plastic Cauchy strain
          Cp = matmul(transpose(Fp),Fp)
!---------Calculate plastic Green-Langrange strain
          call IdentityMatrix(3,XI33)
          Ep = 0.5 * ( Cp - XI33)
!---------Calculate amount of EP for use for peeq
          call icams_matrixnorm(Ep,3,3,Ep_norm)
!---------Calculate peeq
          peeq = peeq + sqrt(2./3.)*Ep_norm
!---------write(*,*)'peeq=',peeq
          return
      endsubroutine cal_eqpl_strain
!-----------------------------------------------------------------------------
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------------
!      Accumulated plastic slip 
!-----------------------------------------------------------------------------  
      subroutine cal_pl_slip(Lp,dt1,p_slip)
          implicit none 
          real(8), intent(in)::  Lp(3,3)                         !Plastic deformation gradient 
          real(8), intent(in)::  dt1                             !Time increment
          real(8), intent(inout):: p_slip                        !Plastic equivalent strain at tn1
          
!---------Local parameters
          real(8) Lp_norm                                        !Plastic velocity gradient after Matrix verjüngung
!---------Calculate amount of EP for use for peeq
          call icams_matrixnorm(Lp,3,3,Lp_norm)
!---------Calculate peeq
          p_slip = p_slip + sqrt(2./3.)*Lp_norm*dt1
          return
      end subroutine cal_pl_slip
!-----------------------------------------------------------------------------
