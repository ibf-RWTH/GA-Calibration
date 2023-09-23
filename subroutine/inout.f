! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @author 
!> Philipp Engels
! DESCRIPTION: 
!> Displays starting screen to .dat-file and terminal (if allowed)
!--------------------------------------------------------------------------
      subroutine inout_startscreen()
          use GlobalVariables, only: ctrl_glob_startingtime, 
     &    ctrl_printinfostoterminal,ctrl_versiondate, ctrl_versioninfos
          implicit none
          integer :: clock_max, clock_reading, clock_rate
          integer :: ctrl_numprocesses
          integer :: ctrl_warningflag
          integer :: irank
          real(8) :: clock_reading_r
          real(8) :: clock_rate_r
!---------Calculate initial time
          call system_clock(clock_reading, clock_rate, clock_max)
          clock_rate_r = real(clock_reading)
          clock_rate_r = real(clock_rate)
          ctrl_glob_startingtime = clock_reading_r/clock_rate_r
!---------Get CPU data (only for MPI)
          call getnumcpus(ctrl_numprocesses) 
          call getrank(irank) 
          if (irank==0) then
              write(6,*) '                                                                  '
              write(6,*) '                                                                  '
              write(6,*) '    --------------------------------------------------------------'
              write(6,*) '    --------------------Crystal Plasticity Code-------------------'
              write(6,*) '    -------------------------MMM/ICAMS/RUB------------------------'
              write(6,*) '    --------------------------------------------------------------'
              write(6,*) '    --------------------------------------------------------------'
              write(6,*) '      Version: ', trim(ctrl_versiondate)
              write(6,*) '      Infos: ', trim(ctrl_versioninfos)
              write(6,*) '      Number of used CPUs (MPI only): ', ctrl_numprocesses
              write(6,*) '                                                                  '
              write(6,*) '      <Start>                                                     '
              write(6,*) '                                                                  '
              if (ctrl_printinfostoterminal .eqv. .true.) then
                  write(*,*) '                                                                  '
                  write(*,*) '                                                                  '
                  write(*,*) '    --------------------------------------------------------------'
                  write(*,*) '    --------------------Crystal Plasticity Code-------------------'
                  write(*,*) '    -------------------------MMM/ICAMS/RUB------------------------'
                  write(*,*) '    --------------------------------------------------------------'
                  write(*,*) '    --------------------------------------------------------------'
                  write(*,*) '      Version: ', trim(ctrl_versiondate)
                  write(*,*) '      Infos: ', trim(ctrl_versioninfos)
                  write(*,*) '      Number of used CPUs (MPI only): ', ctrl_numprocesses
                  write(*,*) '                                                                  '
                  write(*,*) '      <Start>                                                     '
                  write(*,*) '                                                                  '
              end if 
              if (ctrl_numprocesses>1) then 
                  ctrl_warningflag = 31
                  call inout_warningflagglobal(ctrl_warningflag)
              end if
          end if
      end subroutine inout_startscreen

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @author 
!> Philipp Engels
! DESCRIPTION: 
!> Prints increment informations to .dat-file and terminal (if allowed)
!> @param[in] dtime         Time step duration
!> @param[in] time          Time information vector of ABAQUS
!> @param[in] kinc          Increment number
!--------------------------------------------------------------------------
      subroutine inout_printincrementinfos(dtime, time, kinc)
          use GlobalVariables, only: ctrl_printinfostoterminal
!--------------------------------------------------------------------------
          implicit none
          integer, intent(in) :: kinc         !Time step duration
          real(8), intent(in) :: dtime         !Time step duration
          real(8), intent(in), dimension(2) :: time !Time information vector of ABAQUS
!---------Local variables
          integer :: i,j         !Loop integers
          integer :: irank
          
          call getrank(irank) 
          
          write(6,62) kinc, dtime, time(2)
          if (ctrl_printinfostoterminal .eqv. .true. .and. irank==0) then
              write(*,62) kinc, dtime, time(2)
!-------------write(*,*) ' Increment duration: ', dtime,' Total time: ', time(2)
          end if
 62       format (I4,' / Increment duration: ', F10.6, 's / Total time: ', F10.4, 's')
          return
      end subroutine inout_printincrementinfos

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @author 
!> Philipp Engels
! DESCRIPTION: 
!> Displays ending screen to .dat-file and terminal (if allowed)
!--------------------------------------------------------------------------
      subroutine inout_endscreen()
          use GlobalVariables, only: ctrl_glob_startingtime, ctrl_displayedwarnings,ctrl_printinfostoterminal
!--------------------------------------------------------------------------
          implicit none
          integer :: clock_max, clock_reading, clock_rate
          integer :: irank
          real(8) :: ctrl_finishingtime
          call system_clock(clock_reading, clock_rate, clock_max)
          ctrl_finishingtime = real (clock_reading, kind=8)/real(clock_rate, kind=8) - ctrl_glob_startingtime
          call getrank(irank)
          if (irank == 0) then
              if (ctrl_printinfostoterminal .eqv. .true.) then
                  write(*,*) '                                                                  '
                  write(*,92), ctrl_finishingtime
                  write(*,93), ctrl_displayedwarnings
                  write(*,*) '                                                                  '
                  write(*,*) '    <End>                                                         '
                  write(*,*) '    ------------------Crystal Plasticity Code---------------------'
                  write(*,*) '    --------------------------------------------------------------'
                  write(*,*) '    --------------------------------------------------------------'
              end if
          write(6,*) '                                                                  '
          write(6,92) ctrl_finishingtime
          write(6,93) ctrl_displayedwarnings
          write(*,*) '                                                                  '
          write(6,*) '    <End>                                                         '
          write(6,*) '    --------------------------------------------------------------'
          write(6,*) '    ------------------Crystal Plasticity Code---------------------'
          write(6,*) '    --------------------------------------------------------------'
          write(6,*) '    --------------------------------------------------------------'
 92       format ('    Calculation time: ', F10.3, ' seconds')
 93       format ( '    Collected warnings: ', I5)
          end if
      end subroutine inout_endscreen

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @author 
!> Philipp Engels
! DESCRIPTION: 
!> Surpress warnings if ctrl_displayedwarnings > ctrl_maxwarningstodisplay
!--------------------------------------------------------------------------
      subroutine inout_checkwarnings()
          use GlobalVariables, only: ctrl_displayedwarnings, ctrl_maxwarningstodisplay, ctrl_printinfostoterminal        
!--------------------------------------------------------------------------
          implicit none
          integer :: irank
          call getrank(irank)
          if (ctrl_displayedwarnings>=ctrl_maxwarningstodisplay .and. ctrl_displayedwarnings < 666 .and. irank == 0) then
              if (ctrl_printinfostoterminal .eqv. .true.) then
                  write(*,*),'To many warnings! Further warnings surpressed!'
              end if
              write(6,*),'To many warnings! Further warnings surpressed!'
              ctrl_displayedwarnings = 666 
          end if
      end subroutine inout_checkwarnings

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @author 
!> Philipp Engels
! DESCRIPTION: 
!> Displays warning flag returned from stress calculation to .dat-file and terminal (if allowed)
!> @param[in] ctrl_warningflag         Control flag
!--------------------------------------------------------------------------
      subroutine inout_warningflagglobal(ctrl_warningflag)
          use GlobalVariables, only: ctrl_displayedwarnings, ctrl_printinfostoterminal
!--------------------------------------------------------------------------
          implicit none
          integer, intent(in) :: ctrl_warningflag    ! Control flag
          select case(ctrl_warningflag)    
              case(31)
                  write(6,*) 'Warning: Multiple CPUs only supported with iterative solver!'
                  if (ctrl_displayedwarnings /=666 .and. ctrl_printinfostoterminal .eqv. .true.) then
                      write(*,*) 'Warning: Multiple CPUs only supported with iterative solver!'
                  end if
          end select
          ctrl_displayedwarnings = ctrl_displayedwarnings + 1
      end subroutine inout_warningflagglobal

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @author 
!> Philipp Engels
! DESCRIPTION: 
!> Displays warning flag returned from stress calculation to .dat-file and terminal (if allowed)
!> @param[in] ctrl_warningflag         Control flag
!> @param[in] noel               Element Number
!> @param[in] npt                Integration Point
!--------------------------------------------------------------------------
      subroutine inout_warningflaglocal(ctrl_warningflag, noel, npt)
          use GlobalVariables, only: ctrl_displayedwarnings, ctrl_printinfostoterminal
!--------------------------------------------------------------------------
          implicit none
          integer, intent(in) :: ctrl_warningflag    ! Control flag
          integer, intent(in) :: noel                    ! Element Number
          integer, intent(in) :: npt                     ! Integration Point
          select case(ctrl_warningflag)
              case(112)
                write(6,"('Warning: Abnormal IVB magnitude at EL ', I5,' / GP ', I2)") noel, npt  
                if (ctrl_displayedwarnings /=666 .and. ctrl_printinfostoterminal .eqv. .true.) then
                     write(*,"('Warning: Abnormal IVB magnitude at EL', I5,' / GP ', I2)") noel, npt
                end if
              case(113)
                write(6,"('Warning: Abnormal dgmdt magnitude at EL ', I5,' / GP ', I2)") noel, npt  
                if (ctrl_displayedwarnings /=666 .and. ctrl_printinfostoterminal .eqv. .true.) then
                     write(*,"('Warning: Abnormal dgmdt magnitude at EL', I5,' / GP ', I2)") noel, npt
                end if  
              case(114)
                write(6,"('Warning: Abnormal ddgmdt_dtau magnitude at EL ', I5,' / GP ', I2)") noel, npt  
                if (ctrl_displayedwarnings /=666 .and. ctrl_printinfostoterminal .eqv. .true.) then
                     write(*,"('Warning: Abnormal ddgmdt_dtau magnitude at EL', I5,' / GP ', I2)") noel, npt
                end if 
              case(115)
                write(6,"('Warning: Abnormal ddgmdt_dIVB magnitude at EL ', I5,' / GP ', I2)") noel, npt  
                if (ctrl_displayedwarnings /=666 .and. ctrl_printinfostoterminal .eqv. .true.) then
                     write(*,"('Warning: Abnormal ddgmdt_dIVB magnitude at EL', I5,' / GP ', I2)") noel, npt
                end if
              case(116)
                write(6,"('Warning: Abnormal divbdt magnitude at EL ', I5,' / GP ', I2)") noel, npt  
                if (ctrl_displayedwarnings /=666 .and. ctrl_printinfostoterminal .eqv. .true.) then
                     write(*,"('Warning: Abnormal divbdt magnitude at EL', I5,' / GP ', I2)") noel, npt
                end if
              case(117)
                write(6,"('Warning: Abnormal ddIVBdt_ddgmdt magnitude at EL ', I5,' / GP ', I2)") noel, npt  
                if (ctrl_displayedwarnings /=666 .and. ctrl_printinfostoterminal .eqv. .true.) then
                     write(*,"('Warning: Abnormal ddIVBdt_ddgmdt magnitude at EL', I5,' / GP ', I2)") noel, npt
                end if
              case(118)
                write(6,"('Warning: Abnormal ddIVBdt_dIVB magnitude at EL ', I5,' / GP ', I2)") noel, npt  
                if (ctrl_displayedwarnings /=666 .and. ctrl_printinfostoterminal .eqv. .true.) then
                     write(*,"('Warning: Abnormal ddIVBdt_dIVB magnitude at EL', I5,' / GP ', I2)") noel, npt
                end if
              case(119)
                write(6,"('Warning: Abnormal gam magnitude at EL ', I5,' / GP ', I2)") noel, npt  
                if (ctrl_displayedwarnings /=666 .and. ctrl_printinfostoterminal .eqv. .true.) then
                     write(*,"('Warning: Abnormal gam magnitude at EL', I5,' / GP ', I2)") noel, npt
                end if
              case(5)
                if (ctrl_displayedwarnings /=666 .and. ctrl_printinfostoterminal .eqv. .true.) then
                  write(*,"('Warning: Non-covergence for ctrl_NRmaxIterations at EL ', I5,' / GP ', I2)") noel, npt  
                end if 
                write(6,"('Warning: Non-covergence for ctrl_NRmaxIterations at EL ', I5,' / GP ', I2)") noel, npt
              case(1)
                 if (ctrl_displayedwarnings /=666 .and. ctrl_printinfostoterminal .eqv. .true.) then
                   write(*,"('Warning: IFp0 not invertable at EL ', I5,' / GP ', I2)") noel, npt  
                 end if
                 write(6,"('Warning: IFp0 not invertable at EL ', I5,' / GP ', I2)") noel, npt
              case(20)
                if (ctrl_displayedwarnings /=666 .and. ctrl_printinfostoterminal .eqv. .true.) then
                  write(*,"('Warning: at EL ', I5,' / GP ', I2)") noel, npt  
                end if
                 write(6,"('Warning: at EL ', I5,' / GP ', I2)") noel, npt
              case(2)
                 if (ctrl_displayedwarnings /=666 .and. ctrl_printinfostoterminal .eqv. .true.) then
               write(*,"('Warning: Fp is non-invertible at EL ', I5,' / GP ', I2)") noel, npt
                 end if
                 write(6,"('Warning: Fp is non-invertible at EL ', I5,' / GP ', I2)") noel, npt    
              case(22)
                 if (ctrl_displayedwarnings /=666 .and. ctrl_printinfostoterminal .eqv. .true.) then
                   write(*,"('Warning: Fp is non-invertible at EL ', I5,' / GP ', I2)") noel, npt 
                 end if
                 write(6,"('Warning: Fp is non-invertible at EL ', I5,' / GP ', I2)") noel, npt
              case(41)
                if (ctrl_displayedwarnings /=666 .and. ctrl_printinfostoterminal .eqv. .true.) then
                  write(*,"('Warning: dGv1_dpk2i non-invertible at EL ', I5,' / GP ', I2)") noel, npt 
                end if
                write(6,"('Warning: dGv1_dpk2i non-invertible at EL ', I5,' / GP ', I2)") noel, npt
              case(42)
                if (ctrl_displayedwarnings /=666 .and. ctrl_printinfostoterminal .eqv. .true.) then
                  write(*,"('Warning: dGv2_dIVB non-invertible at EL ', I5,' / GP ', I2)") noel, npt  
                end if
                write(6,"('Warning: dGv2_dIVB non-invertible at EL ', I5,' / GP ', I2)") noel, npt 
              case(43)
                if (ctrl_displayedwarnings /=666 .and. ctrl_printinfostoterminal .eqv. .true.) then
                  write(*,"('Warning: eqM66Gv1 non-invertible at EL ', I5,' / GP ', I2)") noel, npt  
                end if
                write(6,"('Warning: eqM66Gv1 non-invertible at EL ', I5,' / GP ', I2)") noel, npt 
              case(44)
                if (ctrl_displayedwarnings /=666 .and. ctrl_printinfostoterminal .eqv. .true.) then
                  write(*,"('Warning: eqMnnGv2 non-invertible at EL ', I5,' / GP ', I2)") noel, npt  
                end if
                write(6,"('Warning: eqMnnGv2 non-invertible at EL ', I5,' / GP ', I2)") noel, npt 
              case(3)
                if (ctrl_displayedwarnings /=666 .and. ctrl_printinfostoterminal .eqv. .true.) then  
                  write(*,"('Warning: Fg non-invertible at EL ', I5,' / GP ', I2)") noel, npt  
                end if
                write(6,"('Warning: Fg non-invertible at EL ', I5,' / GP ', I2)") noel, npt 
              case(46)
                if (ctrl_displayedwarnings /=666 .and. ctrl_printinfostoterminal .eqv. .true.) then
                  write(*,"('Warning:  at EL ', I5,'-', I2)") noel, npt
                end if
                write(6,"('Warning:  at EL ', I5,'-', I2)") noel, npt   
              case(220)
                if (ctrl_printinfostoterminal .eqv. .true.) then
                  write(*,*) 'Error: negative dislocation density detected. EXIT!'
                end if
                write(6,*) 'Error: negative dislocation density detected. EXIT!'
          end select 
          ctrl_displayedwarnings = ctrl_displayedwarnings + 1
      end subroutine inout_warningflaglocal
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @author 
!> Philipp Engels
! DESCRIPTION: 
!> Displays error flag, terminates calculation!
!> @param[in] ctrl_warningflag         Control flag
!> @param[in] noel               Element Number
!> @param[in] npt                Integration Point
!--------------------------------------------------------------------------
      subroutine inout_errorflag(ctrl_warningflag)
          use GlobalVariables, only: ctrl_printinfostoterminal
          !--------------------------------------------------------------------------
          implicit none
          integer, intent(in) :: ctrl_warningflag    !Control flag
          select case(ctrl_warningflag)
              case(17)
                if (ctrl_printinfostoterminal .eqv. .true.) then
                  write(*,*) 'Error: Unspecified allocation problem! EXIT!'
                end if
                write(6,*) 'Error: Unspecified allocation problem! EXIT!'
              case(55)
                if (ctrl_printinfostoterminal .eqv. .true.) then
                  write(*,*) 'Error: matdata.inp not found! EXIT!'
                end if
                write(6,*) 'Error: matdata.inp not found! EXIT!'
              case(56)
                if (ctrl_printinfostoterminal .eqv. .true.) then
                  write(*,*) 'Error: could not parse value from inputfile! EXIT!'
                end if
                write(6,*) 'Error: could not parse value from inputfile! EXIT!'
              case(57)
                if (ctrl_printinfostoterminal .eqv. .true.) then
                  write(*,*) 'Error: no material descriptor found! EXIT!'
                end if
                write(6,*) 'Error: no material descriptor found! EXIT!'
              case(58)
                if (ctrl_printinfostoterminal .eqv. .true.) then
                  write(*,*) 'Error: graindata.inp not found! EXIT!'
                end if
                write(6,*) 'Error: graindata.inp not found! EXIT!'
          end select
         CALL XIT
      end subroutine inout_errorflag

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @author 
!> Philipp Engels
! DESCRIPTION: 
!> Displays success!
!> @param[in] ctrl_warningflag   Control flag
!> @param[in] noel               Element Number
!> @param[in] npt                Integration Point
!--------------------------------------------------------------------------
      subroutine inout_success(ctrl_warningflag)
          use GlobalVariables, only: ctrl_printinfostoterminal
!--------------------------------------------------------------------------
          implicit none
          integer, intent(in) :: ctrl_warningflag    !Control flag
          select case(ctrl_warningflag)
              case(9999)
                if (ctrl_printinfostoterminal .eqv. .true.) then
                    write(*,*) 'Ok: .'
                end if
                write(6,*) 'Ok: .'
          end select
      end subroutine inout_success

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @author 
!> Philipp Engels
! DESCRIPTION: 
!> Prints Matrix to .dat-file and terminal
!> @param[in] M33         3x3 Matrix to be displayed
!--------------------------------------------------------------------------
      subroutine inout_print_33matrix(M33)
!--------------------------------------------------------------------------
          implicit none
          integer :: i,j
          real(8), intent(in), dimension(3,3) :: M33
          do i=1,3
              write(*,60,advance = 'NO') (M33(i,j),j=1,3)
              write(6,60,advance = 'NO') (M33(i,j),j=1,3)
          end do
 60       format ('|',3(ES15.4,' '),'|',/)    
          return
      end subroutine inout_print_33matrix

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @author 
!> Philipp Engels
! DESCRIPTION: 
!> Prints Matrix to .dat-file and terminal
!> @param[in] M66         6x6 Matrix to be displayed
!--------------------------------------------------------------------------
      subroutine inout_print_66matrix(M66)
!--------------------------------------------------------------------------
          implicit none
          integer :: i,j
          real(8), intent(in), dimension(6,6) :: M66
          
          do i=1,6
              write(*,61,advance = 'NO') (M66(i,j),j=1,6)
              write(6,61,advance = 'NO') (M66(i,j),j=1,6)
          end do
 61       format (6(ES15.4,' , '),/)
          return
      end subroutine inout_print_66matrix
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++