!==============================================================================
! src/hamiltonian.f90  --  module Hamiltonian   (SANCHO)
!
! Central-region Hamiltonian: tight-binding chain of Nsites sites,
!   H = sum_i E_CC n_i + t_hop sum_<i,i+1> (c_i^dag c_{i+1} + h.c.)
! plus an EXPERIMENTAL builder from atomic positions (not used yet).
!==============================================================================
module Hamiltonian
   use Constants, only: dp
   implicit none

   complex(dp), allocatable, dimension(:, :) :: H   !...central Hamiltonian (Nsites x Nsites)

contains

   subroutine build_chain_hamiltonian(HCC, n, ECC, thop)
      !...nearest-neighbour chain: HCC(i,i) = ECC, HCC(i,i+-1) = thop
      implicit none
      integer, intent(in) :: n
      real(dp), intent(in) :: ECC, thop
      complex(dp), intent(out) :: HCC(n, n)
      integer :: i, j

      HCC = 0.d0
      do i = 1, n
         HCC(i, i) = ECC
         do j = 1, n
            if (i .eq. j + 1 .or. j .eq. i + 1) then
               HCC(i, j) = thop
            end if
         end do
      end do
   end subroutine build_chain_hamiltonian

   !============================================================================
   ! EXPERIMENTAL (not called): general Hamiltonian from atomic positions.
   ! On-site energies by element, hopping thop between atoms closer than
   ! 'cutoff'. The element table and cutoff are placeholders.
   !============================================================================

   real(dp) function onsite_energy(symbol)
      implicit none
      character(len=*), intent(in) :: symbol
      select case (trim(symbol))
      case ('C')
         onsite_energy = 1.d0
      case ('H')
         onsite_energy = 2.d0
      case ('O')
         onsite_energy = 3.d0
      case default
         onsite_energy = 0.0d0
      end select
   end function onsite_energy

   subroutine build_hamiltonian_from_xyz(Hxyz, N, coords, elements, thop)
      implicit none
      integer, intent(in) :: N
      real(dp), intent(in) :: coords(3, N), thop
      character(len=2), intent(in) :: elements(N)
      complex(dp), intent(out) :: Hxyz(N, N)
      real(dp), parameter :: cutoff = 0.5d0
      integer :: i, j
      real(dp) :: dist

      Hxyz = (0.d0, 0.d0)
      do i = 1, N
         Hxyz(i, i) = onsite_energy(elements(i))
         do j = i + 1, N
            dist = norm2(coords(:, i) - coords(:, j))
            if (dist .le. cutoff) then
               Hxyz(i, j) = thop
               Hxyz(j, i) = thop
            end if
         end do
      end do
   end subroutine build_hamiltonian_from_xyz

end module Hamiltonian
