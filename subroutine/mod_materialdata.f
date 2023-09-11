! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @author 
!> Martin Boeff, Philipp Engels, Philipp Schwittek
!
! DESCRIPTION: 
!> Global module \n
!> pay attention: ONLY VALUES SAVED IN MODULE WICH ARE GLOBALY DEFINED (for all materials)\n
!> and are NOT changed during one time step
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      module MaterialData
        implicit none
        character(len=50) :: materialpropnames(30)
        real(8) :: materialproperties(5,30) = 0.0
        real(8) :: orientations(1000,4) = 0.0
        real(8) :: grainDia(1000) = 6.0
        contains
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @author 
!> Philipp Engels
! DESCRIPTION: 
!> Reads material properties from matdata.inp in the working directory. 
!> To parse the input file to the fortran routine, the library "Fortran String Utilities" of
!> Dr. George Benthien are used. 
!> http://www.gbenthien.net/strings/index.html
!--------------------------------------------------------------------------
          subroutine inout_readmaterialparameters()
              use GlobalVariables, only: ctrl_jobdir
              use Strings, only: removesp, readline, match, parse, value
!----------------------------------------------------------------------------
              implicit none
!-------------Local variables
              integer :: i,k             !Loop variable
              integer :: delpos            !Delimiter position
              integer :: io_error          !Inout error flag
              integer :: ios               !Return flag
              integer :: ntokens           !Number of individualized strings
              integer :: matindex          !Index of red material
              integer :: ctrl_warningflag  !Warning flag
!-------------real(8) :: parameters(cpmaterialprops)
              character(len=350) :: path   !
!-------------character(len=255) :: pathtemp   !
              character(len=200) :: line
              character(len=50) :: tokens(2) ! An array of strings to hold the individual fields on the line
              character(len=50) :: materialtokens(5) ! An array of strings to hold the individual fields on the line
              path =  ctrl_jobdir // '/matdata.inp'    !windows-linux
              call removesp(path) !Remove blanks       !windows-linux
              matindex = 0
              open(unit=20,file=path,status='old',action='read',iostat=io_error)
              if (io_error == 0) then
                  k = 0
                  do i=1,250
                      delpos = 0
                      !read line
                      call readline(20, line, ios)  
                      call match(line, 1, delpos) !detect delimiter <>
                      if (delpos > 1) then
                          call parse(line, ':', materialtokens, ntokens)
                          call value(materialtokens(3), matindex, ios)
                        !   write(*,*), "materialstoken", materialtokens(3)
                          k = 1
                          cycle
                      end if  
                      if (ios < 0) then !End of file
                          close(unit=20)
                          exit
                      elseif (ios > 0) then
                          ctrl_warningflag = 56   !Value could not be parsed
                          call inout_errorflag(ctrl_warningflag)
                      elseif(k/=0) then
                          call parse(line, ':', tokens, ntokens)
                          call value(tokens(2), materialproperties(matindex,k), ios)
                        !   write(*,*), "matindex: ", matindex
                        !   write(*,*), "k: ", k
                        !   write(*,*), "tokens(2): ", tokens(2)
                        !   write(*,*), "materialproperties(matindex,k): ", materialproperties(matindex,k)
                          materialpropnames(k) = tokens(1)
                          k = k + 1
                          
                      else
                          ctrl_warningflag = 57   !No material delimiter found!
                          call inout_errorflag(ctrl_warningflag)
                      end if
                  end do
                ! write(*,*), materialproperties
                !   k = 1
                !   do i=1,30
                !     write(*,*),"i: ", i
                !     write(*,*), "matindex: ", 1
                !     write(*,*), "materialpropnames: ", materialpropnames(i)
                !     write(*,*), "matproperties: ", materialproperties(1,i)
                !     k = k+1
                !   end do
                !   k = 1
                !   do i=1,30
                !     write(*,*), "matindex: ", 2
                !     write(*,*), "materialpropnames: ", materialpropnames(i)
                !     write(*,*), "matproperties: ", materialproperties(2,i)
                !     k = k+1
                !   end do
                !   k = 1
                !   do i=1,30
                !     write(*,*), "matindex: ", 3
                !     write(*,*), "materialpropnames: ", materialpropnames(i)
                !     write(*,*), "matproperties: ", materialproperties(3,i)
                !     k = k+1
                !   end do
                !   k = 1
                !   do i=1,30
                !     write(*,*), "matindex: ", 4
                !     write(*,*), "materialpropnames: ", materialpropnames(i)
                !     write(*,*), "matproperties: ", materialproperties(4,i)
                !     k = k+1
                !   end do
                !   k = 1
                !   do i=1,30
                !     write(*,*), "matindex: ", 5
                !     write(*,*), "materialpropnames: ", materialpropnames(i)
                !     write(*,*), "matproperties: ", materialproperties(5,i)
                !     k = k+1
                !   end do
                !   call XIT
              else
                  ctrl_warningflag = 55  !File could not be opened
                  write(*,*) path
                  call inout_errorflag(ctrl_warningflag)
              end if
              
          end subroutine inout_readmaterialparameters

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @author 
!> Philipp Engels
! DESCRIPTION: 
!> Reads material properties from graindata.inp in the working directory. 
!> To parse the input file to the fortran routine, the library "Fortran String Utilities" of
!> Dr. George Benthien is used. 
!> http://www.gbenthien.net/strings/index.html
!--------------------------------------------------------------------------
          subroutine inout_readgrainorientations() !grainID, orientations
!-------------use GlobalVariables, only: ctrl_printinfostoterminal
              use Strings, only: readline, parse, value
!------------------------------------------------------------------------------
              implicit none
!-------------Local variables
              integer :: i
              integer :: io_error
              integer :: ios
              integer :: iosparse(4)
              integer :: ntokens
              integer :: grainIDtemp
              integer :: ctrl_warningflag
!-------------real(8) :: parameters(cpmaterialprops)
              character*256 :: path
              character*256 :: line
              character(len=50) ::  tokens(6) ! An array of strings to hold the individual fields on the line
              call getoutdir(path, ios) ! Returns working dir of calculation 
              write(6,*) path
              path = trim(path) // '/graindata.inp'
              orientations = 0.0
              grainDia = 6.0
              open(unit=21,file=path,status='old',action='read',iostat=io_error)
              if (io_error == 0) then
                  do i=1,1000
                      call readline(21, line, ios)
                      
                      if (ios < 0) then !End of file
                          close(unit=21)
                          exit
                      end if
                      call parse(line, ':', tokens, ntokens)
                      call value(tokens(2), grainIDtemp, ios)
                      call value(tokens(3), orientations(grainIDtemp,1), iosparse(1))
                      call value(tokens(4), orientations(grainIDtemp,2), iosparse(2))
                      call value(tokens(5), orientations(grainIDtemp,3), iosparse(3))
                      call value(tokens(6), grainDia(grainIDtemp)      , iosparse(4))
                      if (sum(iosparse)+ios /= 0.0) then
                          ctrl_warningflag = 56   !Value could not be parsed
                          call inout_errorflag(ctrl_warningflag)
                      end if
                  end do
              else
                  ctrl_warningflag = 58  !File could not be opened
                 call inout_errorflag(ctrl_warningflag)
              end if
          end subroutine inout_readgrainorientations
      end module MaterialData
