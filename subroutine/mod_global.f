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
      module GlobalVariables
          implicit none
          character(len=50) :: ctrl_versiondate = '17.09.2020'
          character(len=150) :: ctrl_versioninfos = "Gradient calculation not fully implemented yet."
          character(len=11) :: ctrl_abaqus = 'abq6111.exe'
          character(len=12) :: ctrl_averageNodal= '1'    ! 2 - volume average, 1 - nodal average, 0 - no average
          character(len=256) :: ctrl_jobname
          character(len=256) :: ctrl_jobdir
          character(len=100) :: F_ProcessSleep
          character(len=100) :: F_GradientRead
          logical :: ctrl_printinfostoterminal = .true.
          logical :: ctrl_nonlocal = .false.
          integer :: ctrl_displayedwarnings = 0
          integer :: ctrl_maxwarningstodisplay = 50
          integer :: ctrl_isexplicit= 0                  !
          integer :: ctrl_NRmaxIterations = 800          ! Maximum number of local Newton iterations
          real(8) :: ctrl_toler_NRloop = 1.d-5           ! Toleranz for NI in calculation of stress increment
          real(8) :: ctrl_glob_startingtime
          real(8) :: ctrl_infinity = 1.d20               ! Toleranz for large number to catch infinity-error
      end module GlobalVariables
        
      module GlobalVariables2
          implicit none
          real (8) :: a0time
          real (8) :: a1time
          real (8) :: a2time
          real (8) :: a3time
      end module GlobalVariables2
