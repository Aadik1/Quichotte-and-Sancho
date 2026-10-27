!==============================================================================
! common/broyden_mixing.f90  --  module BroydenMixing
!
! Modified Broyden mixing (D.D. Johnson, Phys. Rev. B 38, 12807 (1988)) for
! the fixed point x = x_out(x), i.e. F(x) = x_out(x) - x = 0.
!
!   call broyden_solve(residual, x, n, alpha_mix, nhist, maxscf, tolf, &
!                      log_unit, converged, n_iter, fmax)
!
! residual(n, x, F, abort) evaluates F at x (one SCF step). Setting
! abort = .true. (e.g. NaN, runaway occupations) ends the iteration with
! converged = .false.
!
! Converged when max|F| < tolf, tested BEFORE any mixing, so on success x is
! the input whose own output differs from it by less than tolf everywhere.
! No Jacobian is stored: each step costs ONE residual call plus BLAS work on
! an (n x nhist) history kept in a ring buffer. The Gram matrix dF^T dF is
! updated one column per step. nhist = 0 gives plain linear mixing.
! Memory: 2*nhist*n + 3*n doubles (+ the caller's x).
!
! One line per step is written to log_unit: step, max|F|, rms|F|, history used.
!
! The files in common/ are IDENTICAL in QUICHOTTE and SANCHO
! (checked by `make check-common`).
!==============================================================================
module BroydenMixing
   use Constants, only: dp
   implicit none
   private
   public :: broyden_solve

   abstract interface
      subroutine residual_function(n, x, F, abort)
         import :: dp
         integer, intent(in) :: n
         real(dp), intent(in) :: x(n)
         real(dp), intent(out) :: F(n)
         logical, intent(out) :: abort
      end subroutine residual_function
   end interface

contains

   subroutine broyden_solve(residual, x, n, alpha_mix, nhist, maxscf, tolf, log_unit, converged, n_iter, fmax)
      procedure(residual_function) :: residual
      integer, intent(in) :: n, nhist, maxscf, log_unit
      real(dp), intent(inout) :: x(n)
      real(dp), intent(in) :: alpha_mix, tolf
      logical, intent(out) :: converged
      integer, intent(out) :: n_iter
      real(dp), intent(out) :: fmax

      real(dp), parameter :: w0 = 0.01d0          !...Johnson's regulariser
      real(dp), allocatable :: F(:), Fold(:), xold(:), dF(:, :), U(:, :), a(:, :), cvec(:), gram(:, :)
      real(dp) :: rms, fnorm
      integer :: its, i, slot, nstored, nh, info_l, lda
      logical :: abort
      real(dp), external :: dnrm2

      lda = max(nhist, 1)
      allocate (F(n), Fold(n), xold(n), dF(n, nhist), U(n, nhist), a(lda, lda), cvec(lda), gram(lda, lda))

      converged = .false.
      nstored = 0
      fmax = huge(1.d0)
      write (log_unit, '(a)') '#  step        max|F|          rms|F|   nhist'

      do its = 1, maxscf
         n_iter = its

         !...residual at the current input
         call residual(n, x, F, abort)
         if (abort) exit

         fmax = maxval(abs(F))
         rms = dnrm2(n, F, 1)/sqrt(dble(n))
         write (log_unit, '(i6,2es16.6,i6)') its, fmax, rms, min(nstored, nhist)
         flush (log_unit)

         if (fmax .lt. tolf) then
            converged = .true.
            exit
         end if

         !...store the new (dF, u) pair in the ring buffer
         if (its .gt. 1 .and. nhist .gt. 0) then
            slot = mod(nstored, nhist) + 1
            dF(:, slot) = F - Fold
            fnorm = dnrm2(n, dF(:, slot), 1)
            if (fnorm .gt. 0.d0) then
               dF(:, slot) = dF(:, slot)/fnorm
               U(:, slot) = alpha_mix*dF(:, slot) + (x - xold)/fnorm
               nstored = nstored + 1
               !...only the new row/column of the Gram matrix dF^T dF changes
               nh = min(nstored, nhist)
               call dgemv('T', n, nh, 1.d0, dF, n, dF(:, slot), 1, 0.d0, cvec, 1)
               gram(1:nh, slot) = cvec(1:nh)
               gram(slot, 1:nh) = cvec(1:nh)
            end if
         end if
         nh = min(nstored, nhist)

         xold = x
         Fold = F

         !...gamma = (w0^2 I + dF^T dF)^{-1} dF^T F
         if (nh .gt. 0) then
            a(1:nh, 1:nh) = gram(1:nh, 1:nh)
            do i = 1, nh
               a(i, i) = a(i, i) + w0*w0
            end do
            call dgemv('T', n, nh, 1.d0, dF, n, F, 1, 0.d0, cvec, 1)
            call dposv('U', nh, 1, a, lda, cvec, lda, info_l)
            if (info_l .ne. 0) then
               write (*, *) 'WARNING: Broyden history ill-conditioned, resetting'
               nstored = 0
               nh = 0
            end if
         end if

         !...x_new = x + alpha F - sum_i gamma_i u_i
         x = x + alpha_mix*F
         if (nh .gt. 0) call dgemv('N', n, nh, -1.d0, U, n, cvec, 1, 1.d0, x, 1)
      end do

      deallocate (F, Fold, xold, dF, U, a, cvec, gram)
   end subroutine broyden_solve

end module BroydenMixing
