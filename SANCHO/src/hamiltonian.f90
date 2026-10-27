!==============================================================================
! src/hamiltonian.f90  --  module Hamiltonian   (SANCHO)
!
! Central-region Hamiltonian H (Nsites x Nsites). The input key 'hamiltonian'
! chooses the model; build_hamiltonian (below) calls the matching routine:
!
!   chain   build_chain_hamiltonian    nearest-neighbour chain (default)
!           H = sum_i E_CC n_i + t_hop sum_<i,i+1> (c_i^dag c_{i+1} + h.c.)
!   custom  build_custom_hamiltonian   YOUR model: edit the routine at the end
!
! Rules every model must follow (details: manual, "Using your own Hamiltonian"):
!   * one orbital per site, index i = site i; the same H for both spins
!   * Hermitian with a real diagonal: fill BOTH triangles
!   * energies in eV; no Hubbard/Hartree terms (U enters via the self-energy)
!   * left lead on site 1, right lead on site Nsites
!
! Also here: an EXPERIMENTAL builder from atomic positions (not selectable yet).
!==============================================================================
module Hamiltonian
   use Constants, only: dp
   use Parameters
   implicit none

   complex(dp), allocatable, dimension(:, :) :: H   !...central Hamiltonian (Nsites x Nsites)

contains

   !===========================================================================
   ! Dispatcher: builds the model chosen by the input key 'hamiltonian'.
   ! To add a model of your own: write its routine below, then add ONE case
   ! here, e.g.   case ('ssh');  call build_ssh_hamiltonian()
   ! and select it with  hamiltonian = ssh  (plus N_sites) in input.dat.
   ! H must already be allocated as H(Nsites, Nsites).
   !===========================================================================
   subroutine build_hamiltonian()
      select case (trim(hamiltonian_model))
      case ('chain')
         call build_chain_hamiltonian(H, Nsites, E_CC, t_hop)
      case ('custom')
         call build_custom_hamiltonian()
      case default
         print *, 'ERROR: unknown hamiltonian = ', trim(hamiltonian_model)
         print *, '       built in: chain, custom (see src/hamiltonian.f90)'
         stop 1
      end select
   end subroutine build_hamiltonian

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


   !===========================================================================
   ! YOUR MODEL  (input: hamiltonian = custom)
   !
   ! Replace the body with your own Hamiltonian. As shipped it is the plain
   ! nearest-neighbour chain (E_CC on site, t_hop between neighbours), written
   ! out explicitly, so 'hamiltonian = custom' runs out of the box.
   !
   ! Rules: one orbital per site, index i = 1..Nsites; the same H for both spins;
   ! Hermitian with a real diagonal -- SANCHO has no make_hermitian, so set
   ! H(j,i) = conjg(H(i,j)) yourself; energies in eV; no Hubbard/Hartree terms;
   ! leads attach to site 1 and site Nsites.
   ! New parameters: declare them in parameters.f90 and read them in read_input
   ! (input_output.f90) with  call input_get('my_key', my_var[, default=...]).
   !===========================================================================
   subroutine build_custom_hamiltonian()
      integer :: i

      H = (0.d0, 0.d0)
      do i = 1, Nsites
         H(i, i) = E_CC                        !...on-site energy
      end do
      do i = 1, Nsites - 1
         H(i, i + 1) = t_hop                   !...bond i -> i+1
         H(i + 1, i) = conjg(H(i, i + 1))      !...keep H Hermitian
      end do
   end subroutine build_custom_hamiltonian

end module Hamiltonian
