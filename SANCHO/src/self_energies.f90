!==============================================================================
! src/self_energies.f90  --  module SelfEnergies   (SANCHO)
!
! Interaction (Hubbard) self-energies of the input Green's functions G_in,
! for ALL energies at once, in factorised (bubble) form:
!   build_self_energies   driver, called once per SCF step
!   compute_bubbles       particle-hole bubbles chi (order >= 2) and
!                         particle-particle bubbles phi (order 3)
!   second_order_sigma    Sig2R, Sig2L, Sig2G
!   third_order_sigma     Sig3R, Sig3L, Sig3G (corrected Appendix B)
!   alloc_sigma_work      work arrays (3rd-order ones only when order = 3)
! The first-order (Hartree-Fock) term is local and static; it is added
! directly in interacting_gf (scf.f90).
!
! Work arrays (frequency index first, so all inner loops are contiguous):
!   gR,gA,gL,gG (k,a,b)       = G_in%r/a/L/G (a,b,k)
!   chiX(m,a,b), m = -(n_omega-1):n_omega-1   particle-hole bubbles, nu = m*delta
!   phiX(n,a,b), n = 2:2*n_omega              particle-particle bubbles, Omega = 2*w_init + n*delta
!   SX(m,i,j)  = sum_k U_k D^X_ijk(m)   (X = G,L,R)
!   PX(n,i,j)  = sum_k U_k E^X_ijk(n)
! Sig2G / Sig3G are only needed for checks (Sigma^r - Sigma^a = Sigma^> - Sigma^<).
!==============================================================================
module SelfEnergies
   use Constants, only: dp, hbar, pi, im
   use Parameters
   use GreenFunctions
   implicit none

   complex(dp), allocatable, dimension(:, :, :) :: Sig2R, Sig2L, Sig2G
   complex(dp), allocatable, dimension(:, :, :) :: Sig3R, Sig3L, Sig3G
   complex(dp), allocatable, dimension(:, :, :) :: gR, gA, gL, gG
   complex(dp), allocatable, dimension(:, :, :) :: chiG, chiL, chiR, chiA
   complex(dp), allocatable, dimension(:, :, :) :: phiG, phiL, phiR, phiA
   complex(dp), allocatable, dimension(:, :, :) :: SG3, SL3, SR3, PG3, PL3, PR3

contains

!....Interaction self-energies for all omega at once (they need the whole grid).
!    Called after equal_time_lesser and before every interacting_gf sweep.
   subroutine build_self_energies()
      implicit none
      if (order .ge. 2) then
         call compute_bubbles(order) !...G copies, chi (and phi if order 3)
         call second_order_sigma()
      end if
      if (order .eq. 3) call third_order_sigma()
   end subroutine build_self_energies

   subroutine compute_bubbles(ord)
      !...Steps 0 and 1 of the factorised self-energies (shared by second and third order):
      !   frequency-first copies of G_in, the particle-hole bubbles chi (order >= 2) and the
      !   particle-particle bubbles phi (order 3 only).
      implicit none
      integer, intent(in) :: ord
      integer :: Na, Nw, a, b, m, n, k1, k2
      real(dp) :: pp
      complex(dp) :: s1, s2, s3, s4

      Na = Nsites; Nw = n_omega
      pp = delta/(2.d0*pi)
      call alloc_sigma_work(Na, Nw, ord)

      !...0. copies of the GFs with the frequency index first (contiguous inner loops)
      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(a,b)
      do b = 1, Na
         do a = 1, Na
            gR(:, a, b) = G_in%r(a, b, :); gA(:, a, b) = G_in%a(a, b, :)
            gL(:, a, b) = G_in%L(a, b, :); gG(:, a, b) = G_in%G(a, b, :)
         end do
      end do
      !$OMP END PARALLEL DO

      !...1a. particle-hole bubbles, chi(m,a,b), nu = m*delta; k2 = k1-m must stay on the grid
      !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(a,b,m,k1,k2,s1,s2,s3,s4)
      do b = 1, Na
         do a = 1, Na
            do m = -(Nw - 1), Nw - 1
               s1 = 0; s2 = 0; s3 = 0; s4 = 0
               do k1 = max(1, 1 + m), min(Nw, Nw + m)
                  k2 = k1 - m
                  s1 = s1 + gG(k1, a, b)*gL(k2, b, a)
                  s2 = s2 + gL(k1, a, b)*gG(k2, b, a)
                  s3 = s3 + gR(k1, a, b)*gL(k2, b, a) + gL(k1, a, b)*gA(k2, b, a)
                  s4 = s4 + gA(k1, a, b)*gL(k2, b, a) + gL(k1, a, b)*gR(k2, b, a)
               end do
               chiG(m, a, b) = s1*pp; chiL(m, a, b) = s2*pp
               chiR(m, a, b) = s3*pp; chiA(m, a, b) = s4*pp
            end do
         end do
      end do
      !$OMP END PARALLEL DO

      if (ord .lt. 3) return

      !...1b. particle-particle bubbles, phi(n,a,b), Omega = 2 w_init + n*delta; k2 = n-k1
      !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(a,b,n,k1,k2,s1,s2,s3,s4)
      do b = 1, Na
         do a = 1, Na
            do n = 2, 2*Nw
               s1 = 0; s2 = 0; s3 = 0; s4 = 0
               do k1 = max(1, n - Nw), min(Nw, n - 1)
                  k2 = n - k1
                  s1 = s1 + gG(k1, a, b)*gG(k2, a, b)
                  s2 = s2 + gL(k1, a, b)*gL(k2, a, b)
                  s3 = s3 + gR(k1, a, b)*gG(k2, a, b) + gL(k1, a, b)*gR(k2, a, b)
                  s4 = s4 + gA(k1, a, b)*gG(k2, a, b) + gL(k1, a, b)*gA(k2, a, b)
               end do
               phiG(n, a, b) = s1*pp; phiL(n, a, b) = s2*pp
               phiR(n, a, b) = s3*pp; phiA(n, a, b) = s4*pp
            end do
         end do
      end do
      !$OMP END PARALLEL DO

   end subroutine compute_bubbles

   subroutine second_order_sigma()

      implicit none
      integer :: Na, Nw, i, j, iw, m, kk
      real(dp) :: pp
      complex(dp) :: c2, accR, accL, accG

      Na = Nsites; Nw = n_omega
      pp = delta/(2.d0*pi)
      c2 = -(im*hbar)**2

      !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(i,j,iw,m,kk,accR,accL,accG)
      do j = 1, Na
         do i = 1, Na
            do iw = 1, Nw
               accR = 0; accL = 0; accG = 0
               do m = max(iw - Nw, -(Nw - 1)), min(iw - 1, Nw - 1)
                  kk = iw - m
                  accL = accL + gL(kk, i, j)*chiL(m, i, j)
                  accG = accG + gG(kk, i, j)*chiG(m, i, j)
                  accR = accR + gR(kk, i, j)*chiG(m, i, j) + gL(kk, i, j)*chiR(m, i, j)
               end do
               Sig2R(i, j, iw) = c2*U_site(i)*U_site(j)*pp*accR
               Sig2L(i, j, iw) = c2*U_site(i)*U_site(j)*pp*accL
               Sig2G(i, j, iw) = c2*U_site(i)*U_site(j)*pp*accG
            end do
         end do
      end do
      !$OMP END PARALLEL DO

   end subroutine second_order_sigma

   subroutine third_order_sigma()
      implicit none
      integer :: Na, Nw, i, j, k, m, n, iw, kk
      real(dp) :: pp
      complex(dp) :: c3, accR, accL, accG

      Na = Nsites; Nw = n_omega
      pp = delta/(2.d0*pi)
      c3 = -(im*hbar)**3

      !...steps 0 and 1 (G copies, chi and phi bubbles) are done in compute_bubbles(3)

      !...2. D and E (products in frequency), summed over the internal site k with U_k
      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i,j,k)
      do j = 1, Na
         do i = 1, Na
            SG3(:, i, j) = 0; SL3(:, i, j) = 0; SR3(:, i, j) = 0
            PG3(:, i, j) = 0; PL3(:, i, j) = 0; PR3(:, i, j) = 0
            do k = 1, Na
               SG3(:, i, j) = SG3(:, i, j) + U_site(k)*(chiR(:, i, k)*chiG(:, k, j) + chiG(:, i, k)*chiA(:, k, j))
               SL3(:, i, j) = SL3(:, i, j) + U_site(k)*(chiR(:, i, k)*chiL(:, k, j) + chiL(:, i, k)*chiA(:, k, j))
               SR3(:, i, j) = SR3(:, i, j) + U_site(k)*chiR(:, i, k)*chiR(:, k, j)
               PG3(:, i, j) = PG3(:, i, j) + U_site(k)*(phiR(:, i, k)*phiG(:, k, j) + phiG(:, i, k)*phiA(:, k, j))
               PL3(:, i, j) = PL3(:, i, j) + U_site(k)*(phiR(:, i, k)*phiL(:, k, j) + phiL(:, i, k)*phiA(:, k, j))
               PR3(:, i, j) = PR3(:, i, j) + U_site(k)*phiR(:, i, k)*phiR(:, k, j)
            end do
         end do
      end do
      !$OMP END PARALLEL DO

      !...3. convolution with the remaining GF (Gamma is a product in time)
      !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(i,j,iw,m,n,kk,accR,accL,accG)
      do j = 1, Na
         do i = 1, Na
            do iw = 1, Nw
               accR = 0; accL = 0; accG = 0
               !...particle-hole part: G_ij(omega - nu), kk = iw - m on the grid
               do m = max(iw - Nw, -(Nw - 1)), min(iw - 1, Nw - 1)
                  kk = iw - m
                  accL = accL + gL(kk, i, j)*SL3(m, i, j)
                  accG = accG + gG(kk, i, j)*SG3(m, i, j)
                  accR = accR + gR(kk, i, j)*SG3(m, i, j) + gL(kk, i, j)*SR3(m, i, j)
               end do
               !...particle-particle part: G_ji(Omega - omega), kk = n - iw on the grid
               do n = max(2, iw + 1), min(2*Nw, iw + Nw)
                  kk = n - iw
                  accL = accL + gG(kk, j, i)*PL3(n, i, j)
                  accG = accG + gL(kk, j, i)*PG3(n, i, j)
                  accR = accR + gA(kk, j, i)*PG3(n, i, j) + gG(kk, j, i)*PR3(n, i, j)
               end do
               Sig3R(i, j, iw) = c3*U_site(i)*U_site(j)*pp*accR
               Sig3L(i, j, iw) = c3*U_site(i)*U_site(j)*pp*accL
               Sig3G(i, j, iw) = c3*U_site(i)*U_site(j)*pp*accG
            end do
         end do
      end do
      !$OMP END PARALLEL DO

   end subroutine third_order_sigma

   subroutine alloc_sigma_work(Na, Nw, ord)
      !...allocate the 2nd-order work arrays (and the 3rd-order ones only if
      !...ord = 3) once; re-allocate if the grid changes
      implicit none
      integer, intent(in) :: Na, Nw, ord
      if (allocated(Sig2R)) then
         if (size(Sig2R, 1) == Na .and. size(Sig2R, 3) == Nw .and. (ord < 3 .or. allocated(Sig3R))) return
         deallocate (Sig2R, Sig2L, Sig2G, gR, gA, gL, gG, chiG, chiL, chiR, chiA)
         if (allocated(Sig3R)) deallocate (Sig3R, Sig3L, Sig3G, phiG, phiL, phiR, phiA, SG3, SL3, SR3, PG3, PL3, PR3)
      end if
      !...second order: Sigma, frequency-first GF copies, particle-hole bubbles
      allocate (Sig2R(Na, Na, Nw), Sig2L(Na, Na, Nw), Sig2G(Na, Na, Nw))
      allocate (gR(Nw, Na, Na), gA(Nw, Na, Na), gL(Nw, Na, Na), gG(Nw, Na, Na))
      allocate (chiG(-(Nw - 1):Nw - 1, Na, Na), chiL(-(Nw - 1):Nw - 1, Na, Na), &
                chiR(-(Nw - 1):Nw - 1, Na, Na), chiA(-(Nw - 1):Nw - 1, Na, Na))
      if (ord < 3) return
      !...third order: Sigma, particle-particle bubbles, D and E sums
      allocate (Sig3R(Na, Na, Nw), Sig3L(Na, Na, Nw), Sig3G(Na, Na, Nw))
      allocate (SG3(-(Nw - 1):Nw - 1, Na, Na), SL3(-(Nw - 1):Nw - 1, Na, Na), SR3(-(Nw - 1):Nw - 1, Na, Na))
      allocate (phiG(2:2*Nw, Na, Na), phiL(2:2*Nw, Na, Na), phiR(2:2*Nw, Na, Na), phiA(2:2*Nw, Na, Na))
      allocate (PG3(2:2*Nw, Na, Na), PL3(2:2*Nw, Na, Na), PR3(2:2*Nw, Na, Na))
   end subroutine alloc_sigma_work

end module SelfEnergies
