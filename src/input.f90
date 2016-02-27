Module input
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! Module: loadData
! Updated: 18th February 2016
! --------------------------------------------------------------!
! Description:
! Output data
! --------------------------------------------------------------!

! Setup Modules
  Use kinds          ! from libBP
  Use strings        ! from libBP
  Use general        ! from libBP
  Use globals
! Force declaration of all variables
  Implicit None
! Privacy of variables/functions/subroutines
  Private
! Public Subroutines
  Public :: readInput

  Contains

! ---------------------------------------------------------------------------------------------------
! Save to specific file
! ---------------------------------------------------------------------------------------------------

  Subroutine readInput()
! Saves the eam file to the output directory
    Implicit None   ! Force declaration of all variables
! Private variables
    Character(len=255) :: fileRow, fileRowU
    Integer(kind=StandardInteger) :: i, j
! Read through input file
    Do i=1,inputFileRows
      fileRow = inputFileData(i)
      fileRowU = StrToUpper(fileRow)
      !If(fileRowU(1:1).eq." ")Then
      !End If
    End Do
  End Subroutine readInput

End Module input
