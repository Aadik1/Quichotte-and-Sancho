!==============================================================================
! src/observables.f90  --  module Observables   (QUICHOTTE)
!
! All observables come from the FINAL sweep of solve_scf (scf.f90), i.e. from
! the Green's functions built with the converged self-energies:
!   Current              spin-resolved current through the LEFT lead (from trU, trD)
!   spin_texture         <S_alpha> and spin-resolved populations of every site (from G_less_obs)
!   write_spectral       A_ii(omega) = -2 Im G^r_ii(omega) of every spin-orbital (from dosd)
!   atom_texture         per ATOM: <S_alpha>, populations and orbital moment <L_alpha>
!                        (multi-orbital Hamiltonians: dftb, dftb_hs)
!==============================================================================
module Observables
   use Constants, only: dp, hbar, pi, im
   use Parameters
   use GreenFunctions
   implicit none

   complex(dp), allocatable, dimension(:, :) :: S_a   !...<S_alpha> per site (Norb, 3)

contains

   subroutine Current(Volt, J_up, J_down)
      !...J_s = (delta/hbar) sum_w trans_s(w)   (spin up and spin down), where
      !   trans_s(w) = Tr_s { i GammaL [ (f_L - 1) G^< - f_L G^> ] } / 2pi was stored
      !   by the final sweep (G_at_energy, mode M_FINAL)
      implicit none
      real(dp), intent(in) :: Volt
      real(dp), intent(out) :: J_up, J_down
      integer :: iw

      J_up = 0.d0; J_down = 0.d0
      do iw = 1, n_omega
         J_up = J_up + trU(iw)
         J_down = J_down + trD(iw)
      end do

      J_up = J_up*(delta/hbar)
      J_down = J_down*(delta/hbar)
   end subroutine Current

   subroutine spin_texture(Volt, S_alpha, N, unit_num)
      !...<S_alpha>_i = -i (hbar^2/2) Tr[sigma_alpha G^<_ii(t=0)] and the
      !...spin-resolved populations n_is = -i hbar G^<_is,is(t=0) of every site
      !...(G^<(t=0) = G_less_obs, from the final sweep);
      !...one line per site is appended to unit_num (populations.dat)
      implicit none
      integer, intent(in) :: N, unit_num
      real(dp), intent(in) :: Volt
      integer :: i, s1, s2
      complex(dp) :: S_alpha(N, 3), traces(3)
      complex(dp) :: pop_up, pop_down, pref, pop_tot

      S_alpha = (0.d0, 0.d0)
      pref = -im*0.5d0*hbar*hbar

      do i = 1, Norb, 2
         traces = (0.d0, 0.d0)

         do s1 = 0, 1
            do s2 = 0, 1
               traces(1) = traces(1) + pauli_x(s1 + 1, s2 + 1)*G_less_obs(i + s2, i + s1)
               traces(2) = traces(2) + pauli_y(s1 + 1, s2 + 1)*G_less_obs(i + s2, i + s1)
               traces(3) = traces(3) + pauli_z(s1 + 1, s2 + 1)*G_less_obs(i + s2, i + s1)
            end do
         end do

         S_alpha(i, 1) = traces(1)*pref
         S_alpha(i, 2) = traces(2)*pref
         S_alpha(i, 3) = traces(3)*pref

      end do

      pop_down = (0.d0, 0.d0); pop_up = (0.d0, 0.d0)
      !..write out the spin operator information
      do i = 1, Norb, 2

         pop_up = -im*G_less_obs(i, i)*hbar !...spin up at site i
         pop_down = -im*G_less_obs(i + 1, i + 1)*hbar !...spin down at site i
         pop_tot = pop_up + pop_down

         !...columns: V, site, Sx, Sy, Sz, n_down, n_up, n_total
         write (unit_num, *) Volt, (i + 1)/2, real(S_alpha(i, 1)), real(S_alpha(i, 2)), real(S_alpha(i, 3)), &
            real(pop_down), real(pop_up), real(pop_tot)
      end do
      flush (unit_num)
   end subroutine spin_texture


   subroutine write_spectral(Volt, unit_num)
      !...spectral function A_ii(w) = -2 Im G^r_ii(w) of every spin-orbital,
      !...one block per voltage (gnuplot: index k selects voltage k)
      real(dp), intent(in) :: Volt
      integer, intent(in) :: unit_num
      integer :: iw, i

      write (unit_num, '(a,f10.5,a)') '# V = ', Volt, '   columns: V, omega, A_ii for i = 1..Norb'
      do iw = 1, n_omega
         write (unit_num, *) Volt, omega(iw), (-2.d0*dosd(i, iw), i=1, Norb)
      end do
      write (unit_num, '(/)')
      flush (unit_num)
   end subroutine write_spectral

   subroutine atom_texture(Volt, unit_num)
      !...Per-atom observables from the equal-time density matrix
      !      rho_(a s),(b s') = <c+_(b s') c_(a s)> = -i hbar G^<_(a s),(b s')(t=0)
      !   summed over the orbitals of each atom:
      !      n_s     = sum_a rho_(a s),(a s)
      !      S_alpha = (1/2) sum_a Tr_spin[ sigma_alpha rho_aa ]                 [hbar]
      !      L_alpha = sum_(a,b in p shell) (L_alpha)_ab sum_s rho_(b s),(a s)   [hbar]
      !   with (L_k)_ij = -i eps_kij in the (p_x, p_y, p_z) basis and the DFTB+
      !   orbital order s, p_y, p_z, p_x. L is evaluated in the (Loewdin-)orthogonal
      !   orbitals the Hamiltonian is written in, so it is the orbital moment of
      !   those orbitals (close to, not identical with, the bare atomic p orbitals).
      !   One line per atom is appended to unit_num (populations_atoms.dat):
      !      V, atom, Sx, Sy, Sz, n_down, n_up, n_total, Lx, Ly, Lz
      implicit none
      real(dp), intent(in) :: Volt
      integer, intent(in) :: unit_num
      integer :: at, i, s1, s2, k, p, q, ip, iq
      integer, parameter :: cart(3) = [2, 3, 1]   !...p_y, p_z, p_x -> y, z, x
      real(dp) :: Sv(3), Lv(3), n_up, n_dn
      complex(dp) :: rho2(2, 2), rhop, Lk(3, 3, 3)

      Lk = (0.d0, 0.d0)
      do k = 1, 3
         Lk(k, mod(k, 3) + 1, mod(k + 1, 3) + 1) = -im     !...(L_x)_yz = -i, ...
         Lk(k, mod(k + 1, 3) + 1, mod(k, 3) + 1) = im
      end do

      do at = 1, N_atoms
         Sv = 0.d0; Lv = 0.d0; n_up = 0.d0; n_dn = 0.d0
         do i = atom_first(at), atom_first(at) + atom_norb(at) - 1
            do s1 = 1, 2
               do s2 = 1, 2
                  rho2(s1, s2) = -im*hbar*G_less_obs(2*i - 2 + s1, 2*i - 2 + s2)
               end do
            end do
            n_up = n_up + real(rho2(1, 1)); n_dn = n_dn + real(rho2(2, 2))
            Sv(1) = Sv(1) + 0.5d0*real(sum(pauli_x*transpose(rho2)))
            Sv(2) = Sv(2) + 0.5d0*real(sum(pauli_y*transpose(rho2)))
            Sv(3) = Sv(3) + 0.5d0*real(sum(pauli_z*transpose(rho2)))
         end do
         if (atom_norb(at) .eq. 4) then
            do p = 1, 3
               do q = 1, 3
                  ip = atom_first(at) + p; iq = atom_first(at) + q    !...sites of p_y, p_z, p_x
                  rhop = -im*hbar*(G_less_obs(2*iq - 1, 2*ip - 1) + G_less_obs(2*iq, 2*ip))
                  do k = 1, 3
                     Lv(k) = Lv(k) + real(Lk(k, cart(p), cart(q))*rhop)
                  end do
               end do
            end do
         end if
         write (unit_num, *) Volt, at, Sv, n_dn, n_up, n_dn + n_up, Lv
      end do
      flush (unit_num)
   end subroutine atom_texture

end module Observables
