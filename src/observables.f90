!==============================================================================
! src/observables.f90  --  module Observables   (QUICHOTTE)
!
!   Current              spin-resolved current through the LEFT lead
!   current_integrand    its integrand at one energy
!   spin_texture         <S_alpha> and spin-resolved populations of every site
!   write_spectral       A_ii(omega) = -2 Im G^r_ii(omega) of every spin-orbital
!==============================================================================
module Observables
   use Constants, only: dp, hbar, pi, im
   use Parameters
   use GreenFunctions
   implicit none

   complex(dp), allocatable, dimension(:, :) :: S_a   !...<S_alpha> per site (Norb, 3)

contains

   subroutine current_integrand(iw, Volt, trans_up, trans_down)
      !...integrand of the current through the left lead at omega(iw),
      !   Tr_s { i GammaL [ (f_L - 1) G^< - f_L G^> ] } / 2pi, split by spin
      implicit none
      integer, intent(in) :: iw
      real(dp), intent(in) :: Volt
      real(dp), intent(out) :: trans_up, trans_down
      integer :: i
      complex(dp) :: trace1, trace2
      real(dp) :: w

      w = omega(iw)

      work1 = G_in%L(:, :, iw)
      work2 = G_in%G(:, :, iw)

      work3 = im*matmul(GammaL, (fermi(w, Volt) - 1.d0)*work1 - fermi(w, Volt)*work2)

      trace1 = (0.d0, 0.d0); trace2 = (0.d0, 0.d0)
      do i = 1, Norb, 2
         trace1 = trace1 + work3(i, i)          !...spin up
         trace2 = trace2 + work3(i + 1, i + 1)  !...spin down
      end do

      trans_up = real(trace1)/(2.d0*pi)
      trans_down = real(trace2)/(2.d0*pi)
   end subroutine current_integrand

   subroutine Current(Volt, J_up, J_down)
      !...J_s = (delta/hbar) sum_w trans_s(w)   (spin up and spin down)
      implicit none
      real(dp), intent(in) :: Volt
      real(dp), intent(out) :: J_up, J_down
      real(dp) :: trans_up, trans_down
      integer :: iw

      J_up = 0.d0; J_down = 0.d0
      do iw = 1, n_omega
         call current_integrand(iw, Volt, trans_up, trans_down)
         J_up = J_up + trans_up
         J_down = J_down + trans_down
      end do

      J_up = J_up*(delta/hbar)
      J_down = J_down*(delta/hbar)
   end subroutine Current

   subroutine spin_texture(Volt, S_alpha, N, unit_num)
      !...<S_alpha>_i = -i (hbar^2/2) Tr[sigma_alpha G^<_ii(t=0)] and the
      !...spin-resolved populations n_is = -i hbar G^<_is,is(t=0) of every site;
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
               traces(1) = traces(1) + pauli_x(s1 + 1, s2 + 1)*G_less_t0(i + s2, i + s1)
               traces(2) = traces(2) + pauli_y(s1 + 1, s2 + 1)*G_less_t0(i + s2, i + s1)
               traces(3) = traces(3) + pauli_z(s1 + 1, s2 + 1)*G_less_t0(i + s2, i + s1)
            end do
         end do

         S_alpha(i, 1) = traces(1)*pref
         S_alpha(i, 2) = traces(2)*pref
         S_alpha(i, 3) = traces(3)*pref

      end do

      pop_down = (0.d0, 0.d0); pop_up = (0.d0, 0.d0)
      !..write out the spin operator information
      do i = 1, Norb, 2

         pop_up = -im*G_less_t0(i, i)*hbar !...spin up at site i
         pop_down = -im*G_less_t0(i + 1, i + 1)*hbar !...spin down at site i
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
         write (unit_num, *) Volt, omega(iw), (-2.d0*aimag(G_in%R(i, i, iw)), i=1, Norb)
      end do
      write (unit_num, '(/)')
      flush (unit_num)
   end subroutine write_spectral

end module Observables
