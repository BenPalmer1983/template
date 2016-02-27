Module globals

! --------------------------------------------------------------!
! Global Variables
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!

! Declare all global variables

! ----------------------------------------
! Updated: 21st February 2016
! ----------------------------------------

! Setup Modules
  Use kinds
  Use strings

! Force declaration of all variables
  Implicit None

!================================================================
! Declare Variables
!================================================================

! Standard Variables
!--------------------------------------------------------------
  Integer(kind=StandardInteger) :: mpiProcessCount, mpiProcessID
  Character(len=64) :: compileLine, inputFileName
  Character(len=255) :: currentWorkingDirectory, outputDirectory, tempDirectory
  Integer(kind=StandardInteger) :: inputFileRows
  Character(len=255), Dimension(1:4096) :: inputFileData
  Real(kind=DoubleReal) :: programStartTime, programEndTime, programRunTime

! Standard Parameters
!--------------------------------------------------------------


!================================================================
! Set scope: Public
!================================================================

  Public :: mpiProcessCount, mpiProcessID
  Public :: compileLine, inputFileName
  Public :: currentWorkingDirectory, outputDirectory, tempDirectory
  Public :: inputFileRows, inputFileData
  Public :: programStartTime, programEndTime, programRunTime


!================================================================
! Subroutines
!================================================================

  Contains

! Init global variables
  Subroutine initGlobals()
    Implicit None
! Initialise Subroutine Variable
    Call cpu_time(programStartTime)
    compileLine = "01:00:08  22/12/2015"  !--compile-line-replace
    inputFileData = BlankStringArray(inputFileData)
  End Subroutine initGlobals

End Module globals
