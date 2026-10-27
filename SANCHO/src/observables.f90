!==============================================================================
! src/observables.f90  --  module Observables   (SANCHO)
!
!   Current              current through the LEFT lead (x2 for spin)
!   current_integrand    its integrand at one energy
!   write_populations    occupation per spin of every site
!   write_spectral       A_ii(omega) = -2 Im G^r_ii(omega) of every site
!==============================================================================
module Observables
   use Constants, only: dp, hbar, pi, im
   use LinearAlgebra, only: matrix_trace
   use Parameters
   use GreenFunctions
   implicit none

contains

   subroutine current_integrand(iw, Volt, trns)
      !...integrand of the current through the left lead at omega(iw):
      !   Tr{ i GammaL [ (f_L - 1) G^< - f_L G^> ] } / 2pi
      implicit none
      integer, intent(in) :: iw
      real(dp), intent(in) :: Volt
      real(dp), intent(out) :: trns
      real(dp) :: w
      complex(dp) :: trace1

      w = omega(iw)
      work3 = (0.d0, 0.d0)
      work1 = G_in%L(:, :, iw); work2 = G_in%G(:, :, iw)

      work3 = im*matmul(GammaL, (fermi(w, Volt) - 1.d0)*work1 - fermi(w, Volt)*work2)

      call matrix_trace(work3, Nsites, trace1)
      trns = real(trace1)/(2.d0*pi)

   end subroutine current_integrand

   real(dp) function Current(Volt)
      !...J = 2 (delta/hbar) sum_w current_integrand(w);  the 2 is the spin degeneracy
      implicit none
      real(dp), intent(in) :: Volt
      real(dp) :: J_L, trns
      integer :: iw

      J_L = 0.d0
      do iw = 1, n_omega
         call current_integrand(iw, Volt, trns)
         J_L = J_L + trns
      end do
      Current = 2.d0*J_L*(delta/hbar)
   end function Current

   subroutine write_populations(Volt, unit_num)
      !...occupation per spin n_i = -i hbar G^<_ii(t=0) of every site;
      !...columns: V, site, n_per_spin, n_total (= 2 n_per_spin)
      real(dp), intent(in) :: Volt
      integer, intent(in) :: unit_num
      integer :: i
      real(dp) :: n_i

      do i = 1, Nsites
         n_i = real(-im*hbar*G_less_t0(i))
         write (unit_num, *) Volt, i, n_i, 2.d0*n_i
      end do
      flush (unit_num)
   end subroutine write_populations

   subroutine write_spectral(Volt, unit_num)
      !...spectral function A_ii(w) = -2 Im G^r_ii(w) of every site,
      !...one block per voltage (gnuplot: index k selects voltage k)
      real(dp), intent(in) :: Volt
      integer, intent(in) :: unit_num
      integer :: iw, i

      write (unit_num, '(a,f10.5,a)') '# V = ', Volt, '   columns: V, omega, A_ii for i = 1..Nsites'
      do iw = 1, n_omega
         write (unit_num, *) Volt, omega(iw), (-2.d0*aimag(G_in%R(i, i, iw)), i=1, Nsites)
      end do
      write (unit_num, '(/)')
      flush (unit_num)
   end subroutine write_spectral

end module Observables
