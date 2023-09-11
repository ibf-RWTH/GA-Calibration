      include "mod_global.f"
      include "precmod.f"
      include "stringmod.f"
      include "mod_materialdata.f"
      include "models.f"
      include "inout.f"
      include "stress.f"
      include "other_code.f"
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @author 
!> Martin Boeff, Philipp Engels, Anxin Ma, Philipp Schwittek
!
! DESCRIPTION: 
!> Uexternaldb - Calls external routines at:
!> lop = 0 beginning of analysis,                               
!> 1 start of increment,                                        
!> 2 end of increment,                                          
!> 3 end of analysis,
!> 4 beginning of restart 
!> 5 start of step 
!> 6 end of step
! REVISION HISTORY:
!
!> @param[in] lop              
!> @param[in] lrestart
!> @param[in] time      
!> @param[in] dtime
!> @param[in] kstep
!> @param[in] dtime 
!> @param[in] kinc  
!> @return Starts a new job if wanted
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine uexternaldb(lop,lrestart,time,dtime,kstep,kinc)

          use GlobalVariables, only: ctrl_printinfostoterminal, 
     &                               ctrl_nonlocal,
     &                               ctrl_averageNodal, ctrl_abaqus, 
     &                               ctrl_jobdir, ctrl_jobname, 
     &                               F_ProcessSleep, F_GradientRead
          use MaterialData
          use GlobalVariables2
!------------------------------------------------------------------------------
          implicit none
          integer :: lop,k
          integer :: lrestart
          integer :: kstep 
          integer :: kinc
          integer :: ierr
          integer :: ctrl_warningflag   !flag shows if numerical problems occur
          integer :: lctrl_jobname, lctrl_jobdir
          integer :: irank 
          integer :: clock_max
          integer :: clock_rate
          integer :: clock_reading
          integer :: ip, io_er, read_er
          character(len=8) :: crank
          character(len=8) :: ckstep
          character(len=8) :: ckinc
          character(len=8) :: cmaxNumNodesInElement
          character(len=8) :: cmaxNode
          character(len=1000) :: command
          character(len=100) :: Fname, Fname2
          character(len=512) :: cwd
          character(len=256) :: data_file_name
          character(len=256) :: data_file_path
          real(8) :: time(2)
          real(8) :: dtime
          logical :: exists
          real(8) :: ctime, cctime, ccctime
          real(8) :: wtime
!
!---------Writes current step/increment to string
          write(ckstep,'(I5)') kstep
          write(ckinc,'(I5)') kinc
          ckstep=adjustl(ckstep)
          ckinc=adjustl(ckinc)
          ckstep=trim(ckstep)
          ckinc=trim(ckinc)
!---------Get current directory and jobname
          call getjobname(ctrl_jobname, lctrl_jobname)
          call getoutdir(ctrl_jobdir, lctrl_jobdir)
!---------Get Current cpus rank
          call getrank(irank) 
          write(crank,'(I5)') irank
!---------Start of the calculation
          if(lop==0)then
              call inout_startscreen()
!-------------read material parameters from input-files
              call inout_readmaterialparameters() !read material properties
              call inout_readgrainorientations()  !Call orientions of grain ID
!             write(*,*) orientations
          endif
!---------Start of the increment
          if(lop==1) then
              call system_clock ( clock_reading, clock_rate, clock_max )
              wtime = real ( clock_reading, kind = 8 )/ real ( clock_rate, kind = 8 )
              a0time = wtime
              call inout_checkwarnings()
              call inout_printincrementinfos(dtime, time, kinc)
          endif
!---------End of the increment
          !if(lop==2) then
          !    if (RF .ne. 0) then
          !        open(105,file= data_file_path, status="old", position="append", action="write")
          !        write(105,*) time(2), RF 
          !        close(105)
          !    endif
          !endif
!---------End of the calculation
          if(lop==3) then
              call inout_endscreen()
          endif
          return
      end subroutine uexternaldb

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @author 
!> Mohamed Sharaf, Martin Boeff, Philipp Engels, Anxin Ma
! DESCRIPTION: 
!> ctrl_abaqus UMAT
!> @param[in] ndi                  number of stress components
!> @param[in] nshr                 number of engineering shear stress components
!> @param[in] ntens                size of the stress array (ndi + nshr)
!> @param[in] nstatv               number state variables
!> @param[in] nprops               number of material constants
!> @param[in] layer                layer number
!> @param[in] kspt                 section point number within the current layer
!> @param[in] kstep                step number
!> @param[in] noel                 element number
!> @param[in] npt                  integration point number
!> @param[in] kinc                 increment number
!---------------------------------------------------------------------------
!> @param[in] drpldt               jacobian drpl_dt
!> @param[in] dtime                time increment dt
!> @param[in] temp                 temperature at t0
!> @param[in] dtemp                increment of temperature.
!> @param[in] celent               characteristic element length
!> @param[in] sse                  specific elastic strain energy
!> @param[in] spd                  specific plastic dissipation
!> @param[in] scd                  specific creep dissipation
!> @param[in] rpl                  volumetric heat generation per unit time
!> @param[in] pnewdt               dt_next/dt_now
!---------------------------------------------------------------------------
!> @param[in] ddsdde(ntens,ntens)  jacobian ds_de
!> @param[in] statev(nstatv)       state variables
!> @param[in] props (nprops)       material constants 
!> @param[in] ddsddt(ntens)        jacobian ds_dt
!> @param[in] drplde(ntens)        jacobian drpl_de
!> @param[in] stress(ntens)        stress tensor
!> @param[in] stran (ntens)        strains at t0
!> @param[in] dstran(ntens)        strain increments
!> @param[in] dfgrd0(3,3)          deformation gradient at t0
!> @param[in] dfgrd1(3,3)          deformation gradient at t0+dt
!> @param[in] drot  (3,3)          rotation increment matrix
!> @param[in] coords(3)            coordinates of this point
!> @param[in] time  (2)            1:step time; 2:total time, At t0
!> @param[in] predef(1)            predefined field variables at t0
!> @param[in] dpred (1)            incr of predefined field vrbs
!---------------------------------------------------------------------------
!> @return Stress update \f$\sigma_{n+1}\f$, Material jacobian update \f$ \frac{\delta\sigma}{\delta\varepsilon} \f$
!---------------------------------------------------------------------------
      subroutine umat(stress,statev,ddsdde,
     &                sse,spd,scd,rpl,
     &                ddsddt,drplde,drpldt,
     &                stran,dstran,
     &                time,dtime,
     &                temp,dtemp,
     &                predef,dpred,
     &                cmname,
     &                ndi,nshr,ntens,nstatv,
     &                props,nprops,coords,drot,
     &                pnewdt,
     &                celent,
     &                dfgrd0,dfgrd1,
     &                noel,npt,layer,kspt,kstep,kinc)
!---------------------------------------------------------------------------
!         include 'aba_param.inc'
          use MaterialData, only: materialproperties, orientations, grainDia
          use GlobalVariables, only: ctrl_nonlocal, F_GradientRead
!---------Umat provided variables
          implicit none
          character*80 :: cmname         !user defined material name
          integer :: ndi                  !number of stress components
          integer :: nshr                 !number of engineering shear stress components
          integer :: ntens                !size of the stress array (ndi + nshr)
          integer :: nstatv               !number state variables
          integer :: nprops               !number of material constants
          integer :: layer                !layer number
          integer :: kspt                 !section point number within the current layer
          integer :: kstep                !step number
          integer :: noel                 !element number
          integer :: npt                  !integration point number
          integer :: kinc                 !increment number
          real(8) :: drpldt               !jacobian drpl_dt
          real(8) :: dtime                !time increment dt
          real(8) :: temp                 !temperature at t0
          real(8) :: dtemp                !increment of temperature.
          real(8) :: celent               !characteristic element length
          real(8) :: sse                  !specific elastic strain energy
          real(8) :: spd                  !specific plastic dissipation
          real(8) :: scd                  !specific creep dissipation
          real(8) :: rpl                  !volumetric heat generation per unit time
          real(8) :: pnewdt               !dt_next/dt_now
          real(8) :: ddsdde(ntens,ntens)  !jacobian ds_de
          real(8) :: statev(nstatv)       !state variables
          real(8) :: props (nprops)       !material constants 
          real(8) :: ddsddt(ntens)        !jacobian ds_dt
          real(8) :: drplde(ntens)        !jacobian drpl_de
          real(8) :: stress(ntens)        !stress tensor
          real(8) :: stran (ntens)        !strains at t0
          real(8) :: dstran(ntens)        !strain increments
          real(8) :: dfgrd0(3,3)          !deformation gradient at t0
          real(8) :: dfgrd1(3,3)          !deformation gradient at t0+dt
          real(8) :: drot(3,3)            !rotation increment matrix
          real(8) :: coords(3)            !coordinates of this point
          real(8) :: time  (2)            !1:step time; 2:total time, At t0
          real(8) :: predef(1)            !predefined field variables at t0
          real(8) :: dpred (1)            !incr of predefined field vrbs
          real(8) :: crss_0
          real(8) :: PS_rate
!------------------------------------------------------------------------------------
!---------Material Properties
          integer :: modeltype           !Flag Materialmodel (see documentation)
          integer :: grainID
          integer :: nslip                !Material slip systems
          real(8) :: C11,C12,C44          !Material elastic stiffness if c44=2*(c11-c12) -> iso
          real(8) :: phi1,phi,phi2        !Initial euler angles
          real(8) :: eang00(4)            !Initial euler angles
          real(8) :: Gmod                 !Gmod
!------------------------------------------------------------------------------------
!---------Local defined state variables (0: for tn)
          real(8) :: eang0(4), eang(4)    !Euler angles
          real(8) :: pk2i(6), pk2i_0(6)   !Piola Kirchhoff stress in intermediate configuration
          real(8) :: cs(6), cs0(6)        !Cauchy stress in current configuration Voigt notation
          real(8) :: csM(3,3), csM0(3,3)  !Cauchy stress in current configuration Matrix notation
          real(8) :: Fg(3,3)              !Total deformation gradient
          real(8) :: Fg0(3,3)             !Total deformation gradient at t=n
          real(8) :: Fe(3,3),TFe(3,3)     !Elastic part of deformation grad and transposed
          real(8) :: Fe0(3,3),TFe0(3,3)   !Elastic part of deformation grad and transposed at t=n
          real(8) :: Fp(3,3),TFp(3,3)     !Plasitc part of deformation grad and transposed
          real(8) :: Fp0(3,3),TFp0(3,3)   !Plasitc part of deformation grad and transposed at t=n
          real(8) :: mx33_1(3,3)          !Orientation Matrix assambled with Euler angles
          real(8) :: MatJacb0(6,6),MatJacb(6,6)  !Material jacobian
          real(8) :: STFei26(6,6)         !Material elastic stiffness matrix
!------------------------------------------------------------------------------------
          integer :: i,j                  !Loop integer
          integer :: is,js                !Loop integer typically used for slip systems
          integer :: ib1(9),ib2(9)        !ctrl_abaqus loop variables for conversion from tensor to voigt vector
          integer :: ctrl_isexplicit      !Type of stress calculation (0--implicit, 1--explicit)
          integer :: irank                 !CPU Type
          integer :: noasl                !Number of active slip systems
          real(8) :: dt1                  !Time increment
          real(8) :: peeq                 !plastic equivalent strain
          real(8) :: p_slip               !Plastic equivalent slip
          real(8) :: p_fs(6)              !Fatemi-Socie inspired parameter
!-------------------------------------------------------------------------------------
!---------Values on slip systems
          integer :: ctrl_warningflag                !Warning flag
          real(8), allocatable, dimension(:) :: IVB              !Critical resolved shear stress
          real(8), allocatable, dimension(:) :: IVB0             !Critical resolved shear stress 
          real(8), allocatable, dimension(:) :: tau              !Resolved shear stress
          real(8), allocatable, dimension(:) :: sigma_n          !Resolved normal stress
          real(8), allocatable, dimension(:) :: tau_nsd          !Non-Schmidt resolved shear stress
          real(8), allocatable, dimension(:,:,:) :: smdMi        !List of Schmid matrix
          real(8), allocatable, dimension(:,:) :: smdVi2         !List of Schmid matrix saved as Vector
          real(8), allocatable, dimension(:,:) :: smnVi2         !List of matrix of plane normals saved as Vector
          real(8), allocatable, dimension(:,:) :: n_smdVi2   
          real(8), allocatable, dimension(:,:) :: crsF_mat
          real(8), allocatable, dimension(:) :: IIVB             !Isotropic hardening due to first gradient
          real(8), allocatable, dimension(:) :: BIVB             !Isotropic hardening due to first gradient
          real(8), allocatable, dimension(:,:) :: vl_mat_Norm    
          real(8), allocatable, dimension(:,:) :: vn_mat_Norm    
          real(8), allocatable, dimension(:,:) :: vd_mat_Norm
          real(8), allocatable, dimension(:) :: dgmdt
          real(8), allocatable, dimension(:) :: backstress
          real(8), allocatable, dimension(:) :: backstressn
          real(8), allocatable, dimension(:) :: gam              !Accumulated shear
          real(8), allocatable, dimension(:) :: p_irr            !Shear irreversibility
!---------character(len=100) :: Fname2
!------------------------------------------------------------------------------------
!---------Implementation of Damage by Karl Gillner
          real(8):: PEEQN1
          real(8):: DEPL
          real(8):: D1
          real(8):: D2
          real(8):: GF
          real(8):: TOLER
          real(8):: SOD
          real(8):: DCS1
          real(8):: DCS2
          real(8):: DCS3
          real(8):: DCS4
          real(8):: DCS5
          real(8):: DCS6
          real(8):: SEFNM
          real(8):: ETA
          real(8):: DAMINI
          real(8):: DDAM
          real(8):: DAMAGE
          real(8):: SIGM
          real(8):: old_peeq=0
          real(8):: new_peeq=0
!------------------------------------------------------------------------------------
          call getrank(irank)
!---------Initializations
          dt1=dtime
          ib1=[1,2,3,1,1,2,2,3,3] 
          ib2=[1,2,3,2,3,3,1,1,2]
!---------Initialize type of calculation
          grainID = int(props(1))
          modeltype = int(props(2)) ! where 1=dislo and 2=FCC/phenom and 3,4 and 5 = BCC/phenom
          ctrl_warningflag = 0
!---------Start calling required constant material property for allocation
          nslip = int(materialproperties(modeltype,4)) !int(props(4))
!---------Allocate slip system variables 
          allocate(crsF_mat(nslip,nslip)); crsF_mat = 0.0
          allocate(IVB(nslip)); IVB = 0.0
          allocate(IVB0(nslip)); IVB0 = 0.0
          allocate(tau(nslip)); tau = 0.0
          allocate(sigma_n(nslip)); sigma_n = 0.0
          allocate(smdMi(nslip,3,3)); smdMi = 0.0
          allocate(smdVi2(nslip,6)); smdVi2 = 0.0
          allocate(smnVi2(nslip,6)); smnVi2 = 0.0
          allocate(n_smdVi2(nslip,6)); n_smdVi2 = 0.0
          allocate(vl_mat_Norm(nslip,3)); vl_mat_Norm = 0.0
          allocate(vn_mat_Norm(nslip,3)); vn_mat_Norm = 0.0
          allocate(vd_mat_Norm(nslip,3)); vd_mat_Norm = 0.0
          allocate(dgmdt(nslip)); dgmdt = 0.0
          allocate(backstress(nslip)); backstress = 0.0
          allocate(backstressn(nslip)); backstressn = 0.0
          allocate(gam(nslip)); gam = 0.0
          allocate(p_irr(nslip)); p_irr = 0.0
!---------Start Initializing the state variables at the beginning
!---------Assemble schmid matrix
          call AssembleLattice(nslip,modeltype,smdMi,smdVi2,smnVi2,n_smdVi2,crsF_mat,vn_mat_Norm,vd_mat_Norm,vl_mat_Norm)
!---------Assemble elastic stiffness matrix
          call AssembleStiffness(modeltype, STFei26)
!---------Initialize undeformed euler angles
          eang00(1) = orientations(grainID,1)  !phi1
          eang00(2) = orientations(grainID,2)  !phi
          eang00(3) = orientations(grainID,3)  !phi2   
          if(time(2)==0) then
            
!-------------For t=0 set all state variables to zero
              statev(1:nstatv)=0.0
              
!-------------Calculate Orientation Matrix mx33_1 from initial euler angles
              call icams_eang2Q(eang00(1),eang00(2),eang00(3),mx33_1)
!-------------Save initial Euler angles
!-------------statev(1)=phi1; statev(2)=phi; statev(3)=phi2
              statev(1)=eang00(1); statev(2)=eang00(2); statev(3)=eang00(3)
!-------------Save mx33_1 as initial Fe and Fp
              do i=1,9
                  statev(10+i)=mx33_1(ib1(i),ib2(i))        !Fe
                  statev(19+i)=mx33_1(ib2(i),ib1(i))        !Fp
              enddo
!-------------Save initial statevs for IVB as crss_0
              if (modeltype ==1  .or. modeltype ==2)  then
                statev(28+0*nslip+1 : 28+1*nslip)=materialproperties(modeltype,8)
              else
                statev(28+0*nslip+1 : 28+1*nslip)=materialproperties(
     &            modeltype,8)+materialproperties(modeltype,9)/
     &            sqrt(grainDia(grainID)) !crss_0/rho_0
              endif
          endif
!---------End Initializing the state variables
!------------------------------------------------------------------------------------
!---------Start calling the State Variables from previous time step tn
          eang0  = statev( 1: 4)
          pk2i_0 = statev( 5:10)
          Fg0    = dfgrd0
          Fg     = dfgrd1
          cs0    = stress
          backstressn = statev(28+2*nslip+4+1:28+2*nslip+4+nslip)
!---------Save Cauchy stress in matrix notation
          do i=1,9
              j=i; if(i>6) j=i-3
              csM0(ib1(i),ib2(i))=cs0(j)
          enddo
!---------Call deformation gradient from statevs
          do i=1,9
              Fe0(ib1(i),ib2(i)) = statev(10+i)
              Fe (ib1(i),ib2(i)) = statev(10+i)
              Fp0(ib1(i),ib2(i)) = statev(19+i)
              Fp (ib1(i),ib2(i)) = statev(19+i)
          enddo
!---------Retrieve resolved stresses
          do is=1,nslip
              IVB0(is) = statev(28+0*nslip+is)
              IVB (is) = statev(28+0*nslip+is)
              tau (is) = statev(28+1*nslip+is)
          enddo
!---------Call equivalent plastic strain
          peeq = statev(28+2*nslip+2)
!---------------------------------------------------------------------------
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          PEEQN1 = peeq
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!---------Call equivalent plastic slip
          p_slip = statev(28+2*nslip+3)
!---------Retrieve accumulated shear
          do is=1,nslip
              gam(is) = statev(28+4*nslip+ 4+is)
          enddo
!---------Retrieve Fatemi-Socie inspired parameter
          do i=1,6
              p_fs(i) = statev(28+6*nslip+ 4+i)        !! SDV105:SDV110
          enddo
!---------End calling the State Variables
!---------------------------------------------------------------------------
!---------Start stress calculation 
!---------write(*,*),c11,c22,c44
          call cal_stress_mul(nslip,nprops,props,eang00,eang0,
     &                        crsF_mat,dt1,Fp0, Fg, STFei26, smdMi,
     &                        smdVi2,smnVi2,n_smdVi2,pk2i_0,csM0,
     &                        IVB0,grainID,IIVB,backstressn,cs,pk2i,IVB,
     &                        tau,sigma_n,Fe,Fp,eang,MatJacb,modeltype,
     &                        ctrl_warningflag, p_slip,p_fs,
     &                        dgmdt,backstress,gam,p_irr,statev)
!---------End stress calculation
!----------------------------------------------------------------------------
!---------Save Stress, Material-Jacobian and SDVs
          if(ctrl_warningflag .ne. 0)then
              call inout_warningflaglocal(ctrl_warningflag, noel, npt)
              pnewdt=0.5d0
              return
          else
!-------------Calculate plastic equivalent strain
              call cal_eqpl_strain(Fp,peeq)
              
!-------------start of calculation of plastic strain rate
              new_peeq = peeq
              old_peeq = statev(114)
              PS_rate = (new_peeq-old_peeq)/dtime
              statev(114) = new_peeq
              statev(115) = PS_rate 
!-------------end of calculation of plastic strain rate
!-------------Save Cauchy stress
              stress  = cs
!-------------Save Material-Jacobian
              ddsdde  = MatJacb
!-------------Save SDVs
              statev( 1: 4)   = eang                    !!--> 1- 4    Euler angles
              statev( 5:10)   = pk2i                    !!--> 5-10    2. Piola Kirchhoff stress
!-------------Save Deformation gradients
              do i=1,9
                  statev(10+i)=Fe(ib1(i),ib2(i))     !!-->11-19    Fe
                  statev(19+i)=Fp(ib1(i),ib2(i))     !!-->20-28    Fp
              enddo
!---------------------------------------------------------------------------
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
              DEPL = peeq - PEEQN1                  !! Calculating equivalent plastic strain increment
              D1 = 0.000000003
              D2 = 0.3
              GF = 200.0
              TOLER = 10.0**(-12.0)
              SOD = statev(28+6*nslip+12)
              DAMAGE = statev(28+6*nslip+13)
!-------------Calculating stress triaxiality
              SIGM = (cs(1)+cs(2)+cs(3))/3.0
              DCS1 = cs(1) - SIGM
              DCS2 = cs(2) - SIGM
              DCS3 = cs(3) - SIGM
              DCS4 = cs(4)
              DCS5 = cs(5)
              DCS6 = cs(6)
              SEFNM = DCS1**2.0 + DCS2**2.0 + DCS3**2.0 + 2.0*
     &                DCS4**2.0 + 2.0*DCS5**2.0 + 2.0*DCS6**2.0
              SEFNM = SQRT((3.0/2.0)*SEFNM)
              ETA = SIGM/SEFNM
              DAMINI = D1 * EXP(-D2*ETA)
              IF((peeq.GT.DAMINI).AND.(SOD.LE.0.0)) THEN
                  SOD = SEFNM
              DDAM = (SOD*CELENT)/(2.0*GF)*SQRT(2.0/3.0)
              ELSE IF((peeq.GT.DAMINI).AND.(SOD.GT.0.0)) THEN
                  DDAM = (SOD*CELENT)/(2.0*GF)*SQRT(2.0/3.0)
              END IF
!-------------Calculating Damage
              DAMAGE = DAMAGE + DDAM * DEPL
!-------------Rapid Stiffness degradation due to Fracture Initiation
!-------------IF (DAMAGE.GT.TOLER) THEN
!-------------    ddsdde  = 0.000000001
!-------------     WRITE(6,*) "D_Element=", NOEL
!-------------     WRITE(6,*) "D_INCREMENT=", KINC
!-------------END IF
!-------------Save informations on slip systems
              do is=1,nslip
                  statev(28+0*nslip+is)=IVB(is)!!*(1.0 - DAMAGE)      !!-->29-40    IVB- Critical resolved shear stress with damage
                  statev(28+1*nslip+is)=tau(is)      !!-->41-52    tau- Resolved shear stress
              enddo
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!----------------------------------------------------------------------------------------
!------------- The next 5 lines containing the do-loop are the original once without any change
!-------------Save informations on slip systems
!-------------         do is=1,nslip
!-------------                statev(28+0*nslip+is)=IVB(is)      !!-->29-40    IVB- Critical resolved shear stress
!-------------                statev(28+1*nslip+is)=tau(is)      !!-->41-52    tau- Resolved shear stress
!-------------         enddo
              statev(28+2*nslip+1) = sum(abs(IVB(1:nslip)))      !! total dislocation density
              statev(28+2*nslip+2) = peeq                         !! plastic equivalent strain
              statev(28+2*nslip+3) = p_slip                       !! SDV55 Equivalent plastic slip
!-------------Save Fatemi-Socie-inspired parameter
              do i=1,6
                  statev(28+6*nslip+ 4+i) = p_fs(i)                 !! SDV105:SDV110
              enddo
              statev(28+6*nslip+11)=max(p_fs(1),p_fs(2),p_fs(3),p_fs(4),p_fs(5),p_fs(6)) 
              noasl=0
              do is=1,nslip
                  if (dgmdt(is).GT. 0.00001) then
                      noasl = noasl +1
                  endif
              enddo
              statev(28+2*nslip+ 4) = noasl
              do is=1,nslip
                  statev(28+2*nslip+4+is) = backstress(is)
              enddo
              do is=1,nslip
                  statev(28+3*nslip+4+is)=dgmdt(is)
              enddo
              do is=1,nslip
                  statev(28+4*nslip+ 4+is) = gam(is)                !! Accumulated shear
              enddo
              do is=1,nslip
                  statev(28+5*nslip+ 4+is) = p_irr(is)              !! Shear irreversibility
              enddo
              !
              !---------------------------------------------------------------------------
              !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
              statev(28+6*nslip+12)= SOD
              statev(28+6*nslip+13)= DAMAGE   !!###LAST### SDV113
              !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          endif
!---------End Save Stress, Material-Jacobian and SDVs
          return
      end subroutine umat