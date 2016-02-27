PROGRAM template
! University of Birmingham
! Ben Palmer
!
! This package has been designed with several aims in mind.
! 1.
!
! Internal units:
! - energy, eV
! - length, Angstrom
! - forces, ev/Angstrom
! - stress, ev/ang^3
! - mass, atomic mass A
! - acceleration, ang/ns^2
! - time, ns  (nano seconds)
!
! Printed units:
! - energy, eV
! - length, Angstrom
! - forces, ev/Ang
! - stress, GPa
!
! Calculations:
! F = m a
! m = F (1/a)
! a = F/m        a [ang/ns^2] = 9.6485E9 * (F[ev/ang]/m[amu])
!
! Setup Modules
  Use kinds          ! data kinds
  Use globals        ! declare all globals
  Use output         ! output data
  Use initialise     ! initialise program
  Use input         ! output data
  Use finalise       ! finalise program

! Force declaration of all variables
  Implicit None
! Include MPI header
  Include 'mpif.h'
! Variables
  Integer(kind=StandardInteger) :: error
!======================================
! Init MPI
  Call MPI_Init(error)
!======================================
  Call runInitialise()
  Call initGlobals()
!--------------------------------------------------------------------------------------------------------
! Program Starts !
!----------------!

! Read in input/data files
! --------------------------------------------------------------
  Call readInput()


! --------------------------------------------------------------
  Call runFinalise()

!======================================
! Finalise MPI
  Call MPI_Finalize(error)
!======================================

End Program template
