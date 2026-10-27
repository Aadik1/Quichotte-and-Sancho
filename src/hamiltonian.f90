!==============================================================================
! src/hamiltonian.f90  --  module Hamiltonian   (QUICHOTTE)
!
! Central-region Hamiltonian of a helical molecule with spin-orbit coupling:
!
!   H = sum_i E_CC n_i
!     + t_hop * sum_<i,i+1>   (spin-diagonal nearest-neighbour hopping)
!     + i*lamb * sum_<i,i+2>  (d_{i,i+1} x d_{i,i+2}) . sigma
!
! with d the unit bond vectors between sites. Basis: see parameters.f90.
!==============================================================================
module Hamiltonian
   use Constants, only: dp, pi, im
   use LinearAlgebra, only: cross_product, norm3
   use Parameters
   implicit none

   complex(dp), allocatable, dimension(:, :) :: H   !...central Hamiltonian (Norb x Norb)

contains

   subroutine helix_coordinate(r, i)
      !...Cartesian position of site i = 1 .. N_ions*N_turns on the helix:
      !     x = R cos(hand * 2 pi (i-1)/N_ions),  y = R sin(hand * 2 pi (i-1)/N_ions)
      !     z = (i-1) * dz,   dz = hel_length / (N_turns*(N_ions-1))
      !   i.e. hel_length/N_turns is the rise over the first N_ions-1 bonds of a
      !   turn, so for one turn hel_length is the end-to-end length exactly; for
      !   N_turns turns the end-to-end length is hel_length*(N_sites-1)/(N_turns*(N_ions-1)).
      implicit none
      integer, intent(in) :: i
      real(dp), intent(out) :: r(3)

      r(1) = hel_radius*cos(hand*(i - 1)*2*pi/N_ions)
      r(2) = hel_radius*sin(hand*(i - 1)*2*pi/N_ions)
      r(3) = hel_length*(i - 1)/(N_turns*(N_ions - 1))
   end subroutine helix_coordinate

   subroutine soc_matrix_elements(r, r1, r2, ii, jj)
      !...SOC block between site i (spin-orbitals ii-1, ii) and site i+2
      !   (jj-1, jj), from the positions r, r1, r2 of sites i, i+1, i+2:
      !   H_ij = i*lamb * (d1 x d2) . sigma
      implicit none
      real(dp), intent(in) :: r(3), r1(3), r2(3)
      integer, intent(in) :: ii, jj
      real(dp) :: d1(3), d2(3), v(3)

      d1 = (r - r1)/norm3(r - r1)
      d2 = (r - r2)/norm3(r - r2)
      call cross_product(d1, d2, v)
      H(ii - 1, jj - 1) = im*lamb*v(3)
      H(ii - 1, jj) = im*lamb*(v(1) - im*v(2))
      H(ii, jj - 1) = im*lamb*(v(1) + im*v(2))
      H(ii, jj) = im*lamb*(-v(3))
   end subroutine soc_matrix_elements

   subroutine make_hermitian()
      !...fill the lower triangle of H from the upper triangle
      implicit none
      integer :: i, j

      do i = 1, Norb
         do j = i, Norb
            H(j, i) = conjg(H(i, j))
         end do
      end do
   end subroutine make_hermitian

   subroutine build_helix_hamiltonian()
      !...Hamiltonian of the analytic helix (positions from helix_coordinate)
      implicit none
      integer :: i, ii, jj, Nat
      real(dp) :: r(3), r1(3), r2(3)

      H = (0.0d0, 0.0d0)
      Nat = Norb/2

      do i = 1, Nat
         ii = 2*i

         !...on-site
         H(ii - 1, ii - 1) = E_CC
         H(ii, ii) = E_CC

         !...nearest-neighbour hopping
         if (i .le. Nat - 1) then
            jj = 2*(i + 1)
            H(ii - 1, jj - 1) = t_hop
            H(ii, jj) = t_hop
         end if

         !...next-nearest-neighbour spin-orbit coupling
         if (i .le. Nat - 2) then
            jj = 2*(i + 2)
            call helix_coordinate(r, i)
            call helix_coordinate(r1, i + 1)
            call helix_coordinate(r2, i + 2)
            call soc_matrix_elements(r, r1, r2, ii, jj)
         end if
      end do

      call make_hermitian()
   end subroutine build_helix_hamiltonian

   subroutine build_hamiltonian_from_xyz(xyz_file)
      !...Same model as build_helix_hamiltonian, but site positions are read from an
      !   .xyz file (e.g. one written by tools/helix_xyz.f90).
      !   The file must contain Norb/2 atoms; H must already be allocated.
      !   NOT USED by default -- see the commented call in quichotte.f90.
      implicit none
      character(len=*), intent(in) :: xyz_file
      integer :: i, ii, jj, Nat, ierr, unit_xyz
      real(dp), allocatable :: coords(:, :)
      character(len=2) :: atom_type
      character(len=200) :: comment_line

      open (newunit=unit_xyz, file=trim(xyz_file), status='old', action='read', iostat=ierr)
      if (ierr /= 0) then
         print *, 'ERROR: cannot open XYZ file: ', trim(xyz_file)
         stop 1
      end if

      read (unit_xyz, *) Nat
      if (2*Nat .ne. Norb) then
         print *, 'ERROR: ', trim(xyz_file), ' has', Nat, 'atoms but the input gives', Norb/2, 'sites'
         stop 1
      end if
      allocate (coords(3, Nat))
      read (unit_xyz, '(A)') comment_line
      do i = 1, Nat
         read (unit_xyz, *) atom_type, coords(1:3, i)
      end do
      close (unit_xyz)

      H = (0.0d0, 0.0d0)
      do i = 1, Nat
         ii = 2*i

         H(ii - 1, ii - 1) = E_CC
         H(ii, ii) = E_CC

         if (i .le. Nat - 1) then
            jj = 2*(i + 1)
            H(ii - 1, jj - 1) = t_hop
            H(ii, jj) = t_hop
         end if

         if (i .le. Nat - 2) then
            jj = 2*(i + 2)
            call soc_matrix_elements(coords(:, i), coords(:, i + 1), coords(:, i + 2), ii, jj)
         end if
      end do

      call make_hermitian()
      deallocate (coords)
   end subroutine build_hamiltonian_from_xyz

end module Hamiltonian
