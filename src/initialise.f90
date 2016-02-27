Module initialise

! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! Module: initialise
! Updated: 18th May 2015
! --------------------------------------------------------------!
! Description:
! Makes directories used when program runs
! Creates the output files with headers
! Defines the ProgramTime() function
! --------------------------------------------------------------!

! Setup Modules
  Use kinds          ! from libBP
  Use strings        ! from libBP
  Use general        ! from libBP
  Use printMod       ! from libBP
  Use globals        ! declare all globals
! Force declaration of all variables
  Implicit None
! Include MPI header
  Include 'mpif.h'
! Privacy of functions/subroutines/variables
  Private
! Subroutines
  Public :: runInitialise        !Subroutine

! ------------------------------------------------------------------------!
!                                                                         !
! MODULE SUBROUTINES                                                      !
!                                                                         !
!                                                                         !
! ------------------------------------------------------------------------!

  contains

! Run all the input subroutines
  Subroutine runInitialise()
! Force declaration of all variables
    Implicit None
! Declare private variables
!
    Call initProgram()
    Call programHeader()
  End Subroutine runInitialise

!=======================================================

  Subroutine initProgram()
  ! Force declaration of all variables
    Implicit None
  ! Declare private variables
    Integer(kind=StandardInteger) :: error
  !
  ! MPI variables (public)
    Call MPI_Comm_size(MPI_COMM_WORLD,mpiProcessCount,error)
    Call MPI_Comm_rank(MPI_COMM_WORLD,mpiProcessID,error)
  ! Initialise variables
    inputFileName = BlankString(inputFileName)
  ! Check for input file, exit if no file
    Call get_command_argument(1,inputFileName)
    If(inputFileName(1:4).eq."    ")Then
      If(mpiProcessID.eq.0)Then
        Print *,"No input file, exiting."
      End If
      Call Exit(0)
    End If
  ! Read in input file
      Call readFile(inputFileName, inputFileData, inputFileRows)
  ! get the working directory
      Call getcwd(currentWorkingDirectory)
  ! Set output and temp/scratch directories
      outputDirectory = Trim(currentWorkingDirectory)//"/output"
      tempDirectory = Trim(currentWorkingDirectory)//"/temp"
  End Subroutine initProgram


  Subroutine programHeader()
  ! Force declaration of all variables
    Implicit None
  ! Declare private variables
  !
    Call printBR(80, "=") ! ====
    print *,"                            Fitting Program                             "
    Call printBR(80, "=") ! ====
    print *,"MpiProcesses: ",mpiProcessCount
    print *,"Input File:       ",trim(inputFileName)
    print *,"Temp Dir:         ",trim(outputDirectory)
    print *,"Output Dir:       ",trim(tempDirectory)
    Call printBR(80, "=") ! ====
  End Subroutine programHeader

End Module initialise
