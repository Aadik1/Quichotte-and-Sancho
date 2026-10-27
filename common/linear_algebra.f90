!==============================================================================
! common/linear_algebra.f90  --  module LinearAlgebra
!
! Small linear-algebra helpers shared by QUICHOTTE and SANCHO:
!   * LAPACK wrappers   : hermitian_eigensystem, invert_complex
!   * matrix utilities  : hermitian_conjugate, matrix_trace
!   * 3-vector utilities: cross_product, norm3
!   * Broyden packing   : pack_gf, unpack_gf
!
! Everything lives in one module so the compiler checks every call's
! argument list.
!
! The files in common/ are IDENTICAL in QUICHOTTE and SANCHO
! (checked by `make check-common`).
!==============================================================================
module LinearAlgebra
   use, intrinsic :: iso_fortran_env, only: error_unit
   use Constants, only: dp
   implicit none
   private
   public :: hermitian_eigensystem, invert_complex
   public :: hermitian_conjugate, matrix_trace
   public :: cross_product, norm3
   public :: pack_gf, unpack_gf

contains

   !===========================================================================
   ! LAPACK wrappers
   !===========================================================================

   subroutine hermitian_eigensystem(A, N, E)
      !...Hermitian eigenproblem A*X = E*X (LAPACK zheev).
      !   On exit the columns of A hold the eigenvectors (A is overwritten)
      !   and E holds the eigenvalues in ascending order.
      implicit none
      integer, intent(in) :: N
      complex(dp), intent(inout) :: A(N, N)
      real(dp), intent(out) :: E(N)
      integer :: INFO, LWORK
      complex(dp), allocatable :: WORK(:)
      real(dp), allocatable :: RWORK(:)

      !...workspace query first, then the real call
      allocate (WORK(1), RWORK(max(1, 3*N - 2)))
      call ZHEEV('V', 'U', N, A, N, E, WORK, -1, RWORK, INFO)
      LWORK = int(WORK(1))
      deallocate (WORK); allocate (WORK(LWORK))
      call ZHEEV('V', 'U', N, A, N, E, WORK, LWORK, RWORK, INFO)
      deallocate (WORK, RWORK)

      if (INFO .ne. 0) then
         write (error_unit, *) 'ERROR: ZHEEV eigenproblem failed, INFO =', INFO
         stop 1
      end if
   end subroutine hermitian_eigensystem

   subroutine invert_complex(N, A, INFO)
      !...Inverse of a general complex N x N matrix (LAPACK zgetrf + zgetri).
      !   A is overwritten by its inverse. Stops on a singular matrix.
      !   Thread-safe: all work arrays are local (called inside OpenMP loops).
      implicit none
      integer, intent(in) :: N
      complex(dp), intent(inout) :: A(N, N)
      integer, intent(out) :: INFO
      complex(dp), allocatable :: WORK(:)
      integer, allocatable :: IPIV(:)
      integer :: LWORK

      allocate (IPIV(N), WORK(1))

      call ZGETRF(N, N, A, N, IPIV, INFO)          !...LU factorisation
      if (INFO .ne. 0) then
         write (error_unit, *) 'ERROR: ZGETRF failed (singular matrix?), INFO =', INFO
         stop 1
      end if

      call ZGETRI(N, A, N, IPIV, WORK, -1, INFO)   !...workspace query
      LWORK = int(real(WORK(1)))
      deallocate (WORK); allocate (WORK(LWORK))
      call ZGETRI(N, A, N, IPIV, WORK, LWORK, INFO)
      deallocate (WORK, IPIV)

      if (INFO .ne. 0) then
         write (error_unit, *) 'ERROR: ZGETRI failed, INFO =', INFO
         stop 1
      end if
   end subroutine invert_complex

   !===========================================================================
   ! Matrix and vector utilities
   !===========================================================================

   subroutine hermitian_conjugate(A, N, B)
      !...B = A^dagger
      implicit none
      integer, intent(in) :: N
      complex(dp), intent(in) :: A(N, N)
      complex(dp), intent(out) :: B(N, N)
      integer :: i, j

      do i = 1, N
         do j = 1, N
            B(i, j) = conjg(A(j, i))
         end do
      end do
   end subroutine hermitian_conjugate

   subroutine matrix_trace(A, n, trace)
      !...trace = Tr(A)
      implicit none
      integer, intent(in) :: n
      complex(dp), intent(in) :: A(n, n)
      complex(dp), intent(out) :: trace
      integer :: i

      trace = (0.d0, 0.d0)
      do i = 1, n
         trace = trace + A(i, i)
      end do
   end subroutine matrix_trace

   subroutine cross_product(a, b, c)
      !...c = a x b
      implicit none
      real(dp), intent(in) :: a(3), b(3)
      real(dp), intent(out) :: c(3)

      c(1) = a(2)*b(3) - a(3)*b(2)
      c(2) = a(3)*b(1) - a(1)*b(3)
      c(3) = a(1)*b(2) - b(1)*a(2)
   end subroutine cross_product

   real(dp) function norm3(r)
      !...Euclidean length |r|
      implicit none
      real(dp), intent(in) :: r(3)

      norm3 = sqrt(r(1)**2 + r(2)**2 + r(3)**2)
   end function norm3

   !===========================================================================
   ! Packing of the Green's functions into one real vector for Broyden mixing
   !
   ! Per frequency point iw the block has length nblock = 2*n_full + n_strict_upper + n_upper:
   !
   !   [ Re G^r(i,j)  (n_full)  | Im G^r(i,j)  (n_full)  |
   !     Re G^<(i<j)  (n_strict_upper) | Im G^<(i<=j) (n_upper) ]
   !
   ! with n_full = N^2, n_strict_upper = N(N-1)/2, n_upper = n_strict_upper + N.
   ! Only the independent part of G^< is stored, because G^< is
   ! anti-Hermitian: G^<(j,i) = -conjg(G^<(i,j)) (so Re G^<(i,i) = 0).
   !===========================================================================

   subroutine pack_gf(Input1, Input2, output, vec_len, n_full, n_strict_upper, nb, Nw)
      !...(G^r, G^<) on the omega grid  ->  real vector 'output'
      implicit none
      integer, intent(in) :: nb, Nw, vec_len, n_full, n_strict_upper
      complex(dp), intent(in) :: Input1(nb, nb, Nw)   !...G^r
      complex(dp), intent(in) :: Input2(nb, nb, Nw)   !...G^<
      real(dp), intent(out) :: output(vec_len)
      integer :: i, j, iw, n, n1, nblock, offset

      nblock = 2*n_full + n_strict_upper + (n_strict_upper + nb)

      do iw = 1, Nw
         offset = (iw - 1)*nblock

         n = 0
         do i = 1, nb
            do j = 1, nb
               n = n + 1
               output(offset + n) = real(Input1(i, j, iw))
               output(offset + n + n_full) = aimag(Input1(i, j, iw))
            end do
         end do

         n = 2*n_full + n_strict_upper
         n1 = 2*n_full
         do i = 1, nb
            do j = i, nb
               n = n + 1
               if (i .ne. j) then
                  n1 = n1 + 1
                  output(offset + n1) = real(Input2(i, j, iw))
               end if
               output(offset + n) = aimag(Input2(i, j, iw))
            end do
         end do
      end do
   end subroutine pack_gf

   subroutine unpack_gf(Output1, Output2, input, vec_len, n_full, n_strict_upper, nb, Nw)
      !...real vector 'input'  ->  (G^r, G^<); inverse of pack_gf.
      !   The lower triangle of G^< is rebuilt from anti-Hermiticity.
      implicit none
      integer, intent(in) :: vec_len, nb, Nw, n_full, n_strict_upper
      real(dp), intent(in) :: input(vec_len)
      complex(dp), intent(out) :: Output1(nb, nb, Nw)  !...G^r
      complex(dp), intent(out) :: Output2(nb, nb, Nw)  !...G^<
      complex(dp), parameter :: im = (0.0d0, 1.0d0)
      integer :: i, j, iw, n, n1, nblock, offset
      real(dp) :: re

      nblock = 2*n_full + n_strict_upper + (n_strict_upper + nb)

      do iw = 1, Nw
         offset = (iw - 1)*nblock

         n = 0
         do i = 1, nb
            do j = 1, nb
               n = n + 1
               Output1(i, j, iw) = input(offset + n) + im*input(offset + n + n_full)
            end do
         end do

         n = 2*n_full + n_strict_upper
         n1 = 2*n_full
         do i = 1, nb
            do j = i, nb
               n = n + 1
               re = 0.0d0
               if (i .ne. j) then
                  n1 = n1 + 1
                  re = input(offset + n1)
               end if
               Output2(i, j, iw) = re + im*input(offset + n)
               Output2(j, i, iw) = -conjg(Output2(i, j, iw))
            end do
         end do
      end do
   end subroutine unpack_gf

end module LinearAlgebra
