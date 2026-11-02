!==============================================================================
! common/broyden_mixing.f90  --  module BroydenMixing
!
! Modified Broyden mixing (D.D. Johnson, Phys. Rev. B 38, 12807 (1988)) for
! the fixed point x = x_out(x), i.e. F(x) = x_out(x) - x = 0.
!
!   call broyden_solve(residual, x, n, alpha_mix, nhist, maxscf, tolf, &
!                      log_unit, converged, n_iter, fmax [, single_history])
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
! single_history = .true. (optional, default .false.) stores the history columns
! dF, u in single precision: memory 3*n doubles + 2*nhist*n singles, i.e. the
! history costs half. They only steer the update direction, so rounding them
! changes the path to the fixed point, not the fixed point; x, F and every dot
! product stay in double precision. With the default the iteration is
! bit-for-bit the original one.
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

   subroutine broyden_solve(residual, x, n, alpha_mix, nhist, maxscf, tolf, log_unit, converged, n_iter, fmax, &
                            single_history)
      procedure(residual_function) :: residual
      integer, intent(in) :: n, nhist, maxscf, log_unit
      real(dp), intent(inout) :: x(n)
      real(dp), intent(in) :: alpha_mix, tolf
      logical, intent(out) :: converged
      integer, intent(out) :: n_iter
      real(dp), intent(out) :: fmax
      logical, intent(in), optional :: single_history

      real(dp), parameter :: w0 = 0.01d0          !...Johnson's regulariser
      real(dp), allocatable :: F(:), Fold(:), xold(:), dF(:, :), U(:, :), a(:, :), cvec(:), gram(:, :)
      real, allocatable :: dFs(:, :), Us(:, :)    !...single-precision history (single_history)
      real(dp) :: rms, fnorm, rf
      integer :: its, i, k, slot, nstored, nh, info_l, lda
      logical :: abort, sp
      real(dp), external :: dnrm2

      sp = .false.
      if (present(single_history)) sp = single_history
      lda = max(nhist, 1)
      allocate (F(n), Fold(n), xold(n), a(lda, lda), cvec(lda), gram(lda, lda))
      if (sp) then
         allocate (dFs(n, lda), Us(n, lda))
      else
         allocate (dF(n, nhist), U(n, nhist))
      end if

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
            if (sp) then
               Fold = F - Fold                    !...dF in double, before rounding
               fnorm = dnrm2(n, Fold, 1)
               if (fnorm .gt. 0.d0) then
                  rf = 1.d0/fnorm
                  !$OMP PARALLEL DO
                  do i = 1, n
                     dFs(i, slot) = real(Fold(i)*rf)
                     Us(i, slot) = real(alpha_mix*Fold(i)*rf + (x(i) - xold(i))*rf)
                  end do
                  !$OMP END PARALLEL DO
                  nstored = nstored + 1
                  nh = min(nstored, nhist)
                  do k = 1, nh
                     cvec(k) = dot_ss(n, dFs(:, k), dFs(:, slot))
                  end do
                  gram(1:nh, slot) = cvec(1:nh)
                  gram(slot, 1:nh) = cvec(1:nh)
               end if
            else
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
            if (sp) then
               do k = 1, nh
                  cvec(k) = dot_sd(n, dFs(:, k), F)
               end do
            else
               call dgemv('T', n, nh, 1.d0, dF, n, F, 1, 0.d0, cvec, 1)
            end if
            call dposv('U', nh, 1, a, lda, cvec, lda, info_l)
            if (info_l .ne. 0) then
               write (*, *) 'WARNING: Broyden history ill-conditioned, resetting'
               nstored = 0
               nh = 0
            end if
         end if

         !...x_new = x + alpha F - sum_i gamma_i u_i
         if (sp) then
            !$OMP PARALLEL DO PRIVATE(k)
            do i = 1, n
               x(i) = x(i) + alpha_mix*F(i)
               do k = 1, nh
                  x(i) = x(i) - dble(Us(i, k))*cvec(k)
               end do
            end do
            !$OMP END PARALLEL DO
         else
            x = x + alpha_mix*F
            if (nh .gt. 0) call dgemv('N', n, nh, -1.d0, U, n, cvec, 1, 1.d0, x, 1)
         end if
      end do

      deallocate (F, Fold, xold, a, cvec, gram)
      if (allocated(dF)) deallocate (dF, U)
      if (allocated(dFs)) deallocate (dFs, Us)
   end subroutine broyden_solve

   real(dp) function dot_ss(n, u, v)
      !...double-precision sum of a product of two single-precision columns
      integer, intent(in) :: n
      real, intent(in) :: u(n), v(n)
      integer :: i
      real(dp) :: s
      s = 0.d0
      !$OMP PARALLEL DO REDUCTION(+:s)
      do i = 1, n
         s = s + dble(u(i))*dble(v(i))
      end do
      !$OMP END PARALLEL DO
      dot_ss = s
   end function dot_ss

   real(dp) function dot_sd(n, u, v)
      !...double-precision sum of a single-precision column times a double vector
      integer, intent(in) :: n
      real, intent(in) :: u(n)
      real(dp), intent(in) :: v(n)
      integer :: i
      real(dp) :: s
      s = 0.d0
      !$OMP PARALLEL DO REDUCTION(+:s)
      do i = 1, n
         s = s + dble(u(i))*v(i)
      end do
      !$OMP END PARALLEL DO
      dot_sd = s
   end function dot_sd

end module BroydenMixing
