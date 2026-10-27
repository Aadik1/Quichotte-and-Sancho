!==============================================================================
! common/scf_checks.f90  --  module ScfChecks
!
! Physical sanity checks on the output Green's functions of one SCF step.
! For every basis state i (spin-orbital in QUICHOTTE, site in SANCHO):
!     occupation       n_i = -i hbar (delta/2pi) sum_w G^<_ii(w)
!     spectral weight  W_i = -(hbar delta/pi) sum_w Im G^r_ii(w)   (= 1 if the grid is wide enough)
!
!   runaway     NaN/Inf, or n_i outside [-1, 2]  -> the SCF should be abandoned
!   unphysical  n_i outside [0, 1] (+-0.05) or W_i < 0.5 for some i
!               -> a converged solution of this kind is not trusted
! The worst state (index, n, W) is returned for the error message.
!
! The files in common/ are IDENTICAL in QUICHOTTE and SANCHO
! (checked by `make check-common`).
!==============================================================================
module ScfChecks
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use Constants, only: dp, hbar, pi, im
   implicit none
   private
   public :: check_gf_sanity

contains

   subroutine check_gf_sanity(Gr, GL, nb, nw, delta, runaway, unphysical, i_worst, n_worst, w_worst)
      integer, intent(in) :: nb, nw
      complex(dp), intent(in) :: Gr(nb, nb, nw), GL(nb, nb, nw)
      real(dp), intent(in) :: delta
      logical, intent(out) :: runaway, unphysical
      integer, intent(out) :: i_worst
      real(dp), intent(out) :: n_worst, w_worst
      integer :: i
      real(dp) :: nocc, wsum

      runaway = .false.; unphysical = .false.
      i_worst = 1; n_worst = 0.5d0; w_worst = 1.d0
      do i = 1, nb
         nocc = real(-im*hbar*(delta/(2.d0*pi))*sum(GL(i, i, :)))
         wsum = -hbar*(delta/pi)*sum(aimag(Gr(i, i, :)))
         if (.not. ieee_is_finite(nocc) .or. .not. ieee_is_finite(wsum)) then
            runaway = .true.; i_worst = i; n_worst = nocc; w_worst = wsum
            return
         end if
         if (nocc .lt. -1.d0 .or. nocc .gt. 2.d0) runaway = .true.
         if (nocc .lt. -0.05d0 .or. nocc .gt. 1.05d0 .or. wsum .lt. 0.5d0) unphysical = .true.
         if (abs(nocc - 0.5d0) .gt. abs(n_worst - 0.5d0) .or. &
             (wsum .lt. w_worst .and. abs(n_worst - 0.5d0) .le. 0.55d0)) then
            i_worst = i; n_worst = nocc; w_worst = wsum
         end if
      end do
      if (runaway) unphysical = .true.
   end subroutine check_gf_sanity

end module ScfChecks
