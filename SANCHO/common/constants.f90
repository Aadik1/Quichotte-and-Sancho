!==============================================================================
! common/constants.f90  --  module Constants
!
! Working precision and physical constants. Natural units: hbar = e = 1,
! energies in eV, voltages in V, temperature in K.
!
! The files in common/ are IDENTICAL in QUICHOTTE and SANCHO
! (checked by `make check-common`).
!==============================================================================
module Constants
   use, intrinsic :: iso_fortran_env, only: int64
   implicit none
   private
   public :: dp, int64, hbar, kb, pi, im

   integer, parameter :: dp = kind(1.d0)              !...double precision

   real(dp), parameter :: hbar = 1.d0                 !...reduced Planck constant (natural units)
   real(dp), parameter :: kb = 8.6173303d-5           !...Boltzmann constant [eV/K]
   real(dp), parameter :: pi = 4.d0*atan(1.d0)
   complex(dp), parameter :: im = (0.d0, 1.d0)        !...imaginary unit
end module Constants
