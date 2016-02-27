Module types
! --------------------------------------------------------------!
! Data types
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
! Define data types
! ----------------------------------------
! Updated: 12th Aug 2014
! ----------------------------------------
! Setup Modules
  Use kinds
! Force declaration of all variables
  Implicit None
!
  Type :: chart
    Character(len=128) ::    tempDirectory
    Character(len=128) ::    outputDirectory
    Character(len=64) ::     outputName
    Character(len=32) ::     title
    Character(len=32) ::     xAxis
    Character(len=32) ::     yAxis
    Real(kind=DoubleReal) :: xMin=1.1D99
    Real(kind=DoubleReal) :: xMax=-1.1D99
    Real(kind=DoubleReal) :: yMin=1.1D99
    Real(kind=DoubleReal) :: yMax=-1.1D99
    Logical ::               cleanPyFile=.true.
    Integer(kind=StandardInteger), Dimension(1:10) :: rowStart = -1
    Integer(kind=StandardInteger), Dimension(1:10) :: rowEnd = -1
    Integer(kind=StandardInteger), Dimension(1:10) :: colX = 1
    Integer(kind=StandardInteger), Dimension(1:10) :: colY = 2
  End Type

  
End Module types
