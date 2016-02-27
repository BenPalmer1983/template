Module output
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! Module: loadData
! Updated: 18th February 2016
! --------------------------------------------------------------!
! Description:
! Output data
! --------------------------------------------------------------!

! Setup Modules
  Use kinds
! Force declaration of all variables
  Implicit None
! Privacy of variables/functions/subroutines
  Private
! Public Subroutines
  Public :: outputSubroutine

  Contains

! ---------------------------------------------------------------------------------------------------
! Save to specific file
! ---------------------------------------------------------------------------------------------------

  Subroutine outputSubroutine(fileName)
! Saves the eam file to the output directory
    Implicit None   ! Force declaration of all variables
! Private variables
    Character(*) :: fileName
  End Subroutine outputSubroutine

End Module output
