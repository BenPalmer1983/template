Module finalise
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! Module: finalise
! Updated: 18th May 2015
! --------------------------------------------------------------!
! Description:
! Loads isotope data into memory
! --------------------------------------------------------------!

! Setup Modules
  Use kinds          ! from libBP
  Use printMod       ! from libBP
  Use globals        ! declare all globals

! Force declaration of all variables
  Implicit None
! Include MPI header
  Include 'mpif.h'
! Privacy of functions/subroutines/variables
  Private
! Subroutines
  Public :: runFinalise        !Subroutine

! ------------------------------------------------------------------------!
!                                                                         !
! MODULE SUBROUTINES                                                      !
!                                                                         !
!                                                                         !
! ------------------------------------------------------------------------!

  contains

! Run all the input subroutines
  Subroutine runFinalise()
! Force declaration of all variables
    Implicit None
! Declare private variables
! Time
    Call cpu_time(programEndTime)
    programRunTime = programEndTime - programStartTime
! Internal subroutine variables
    If(mpiProcessID.eq.0)Then
      print *,""
      print *,""
      print *,""
      print *,""
      print *,""
      Call printBR(80, "=") ! ====
      print *,"                           Run Time"
      print *,"Total time:                ",programRunTime
      Call printBR(80, "=") ! ====
    End If
  End Subroutine runFinalise

End Module finalise
