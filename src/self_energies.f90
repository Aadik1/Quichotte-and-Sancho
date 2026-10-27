!==============================================================================
! src/self_energies.f90  --  module SelfEnergies   (QUICHOTTE)
!
! Interaction (Hubbard) self-energies of the input Green's functions G_in:
!   sigma_first_order    first order (Hartree-Fock), static, on-site 2x2 blocks
!   sigma2_compute       second order, Sigma^r and Sigma^< for ALL energies at
!                        once by zero-padded FFT convolutions (FFTW)
!   sigma_retarded / sigma_lesser   total Sigma^r / Sigma^< at one energy
!   sigma2_init / sigma2_finalize   FFT plans and storage
!==============================================================================
module SelfEnergies
   use Constants, only: dp, int64, hbar, pi, im
   use Parameters
   use GreenFunctions
   implicit none

   !.....second-order self-energies for all omegas, built by FFT once per SCF iteration
   include 'fftw3.f'
   integer   :: Nf                     !...padded FFT length >= 2*n_omega - 1
   integer(int64) :: plan_fwd, plan_bwd
   complex(dp), allocatable, dimension(:, :, :) :: Fr, FL        !...(Nf,Norb,Norb): FFT of G^r, G^<
   complex(dp), allocatable, dimension(:, :, :) :: SigR2, SigL2  !...(Norb,Norb,n_omega)

! Every term of the old direct sums (Omega_r / int_SigLnG) is sum_{k1,k2} A(k1) B(k2) C(iw-k1+k2),
! i.e. A(t) B(-t) C(t) after a DFT along omega. Zero-padding to
! Nf >= 2*n_omega-1 makes this exact on the discrete grid (same k3 bounds).
! Spin enters only as 2x2 blocks per orbital pair (i,j):
!   chi1 = L_ji(-l).G_ij(l),  chi2 = L_ji(-l).r_ij(l) + a_ji(-l).L_ij(l),
!   chiL = G_ji(-l).L_ij(l)
!   Sigma^R_ij = c [ r_ij.adj(chi1) + L_ij.adj(chi2) ],  Sigma^<_ij = c L_ij.adj(chiL)
! adj(M) = Tr(M)*1 - M  = direct minus exchange.
!
! Only G^r and G^< are transformed. The other two follow EXACTLY, because the
! SCF builds G_in%a = (G_in%r)^dagger and G_in%G = G_in%L + G_in%r - G_in%a:
!   a_ij(l)  = [ r_ji(-l) ]^dagger            (from  Fa(l,i,j) = conj(Fr(-l,j,i)) )
!   G_ij(l)  = L_ij(l) + r_ij(l) - a_ij(l)    (the FFT is linear)


contains

   subroutine sigma_first_order(Sigma1)
      !...Sigma1 on the on-site 2x2 spin blocks, from G_less_t0 = G^<(t=0):
      !   same spin:      i*hbar*U_i * (G^<_ss - Tr G^<)
      !   opposite spin:  i*hbar*U_i *  G^<_ss'
      implicit none
      complex(dp), intent(out) :: Sigma1(:, :)
      integer :: i, s, s1
      complex(dp) :: Hartree

      Sigma1 = (0.d0, 0.d0); Hartree = (0.d0, 0.d0)

      do i = 1, Norb, 2 !..orbitals

         Hartree = G_less_t0(i, i) + G_less_t0(i + 1, i + 1)
         !.. spins
         do s = 0, 1
            do s1 = 0, 1
               if (s .eq. s1) then
                  Sigma1(i + s, i + s1) = im*hbar*U_site(i)*(G_less_t0(i + s, i + s1) - Hartree)
               else
                  Sigma1(i + s, i + s1) = im*hbar*U_site(i)*G_less_t0(i + s, i + s1)
               end if
            end do
         end do

      end do

   end subroutine sigma_first_order

   subroutine sigma_retarded(iw, SigmaR)
      !...retarded interaction self-energy at omega(iw):
      !   order 0: 0,  order 1: Sigma1 (Hartree-Fock),  order 2: Sigma1 + Sigma2^r
      implicit none
      integer, intent(in) :: iw
      complex(dp), intent(out) :: SigmaR(Norb, Norb)

      SigmaR = (0.d0, 0.d0)
      if (order .ge. 1) call sigma_first_order(SigmaR)
      if (order .eq. 2) SigmaR = SigmaR + SigR2(:, :, iw)
   end subroutine sigma_retarded

   subroutine sigma_lesser(iw, SigmaL)
      !...lesser interaction self-energy at omega(iw): non-zero only at order 2
      !   (the first-order Sigma is static and real-time local, so Sigma1^< = 0)
      implicit none
      integer, intent(in) :: iw
      complex(dp), intent(out) :: SigmaL(Norb, Norb)

      SigmaL = (0.d0, 0.d0)
      if (order .eq. 2) SigmaL = SigL2(:, :, iw)
   end subroutine sigma_lesser

   !....smallest n >= nmin whose prime factors are only 2,3,5,7 (FFT-friendly)
   integer function good_fft_size(nmin)
      implicit none
      integer :: nmin, n, m
      n = nmin
      do
         m = n
         do while (mod(m, 2) == 0); m = m/2; end do
         do while (mod(m, 3) == 0); m = m/3; end do
         do while (mod(m, 5) == 0); m = m/5; end do
         do while (mod(m, 7) == 0); m = m/7; end do
         if (m == 1) exit
         n = n + 1
      end do
      good_fft_size = n
   end function good_fft_size

   subroutine sigma2_init()
      implicit none
      complex(dp), allocatable :: tmp(:)

      Nf = good_fft_size(2*n_omega - 1)
      allocate (Fr(Nf, Norb, Norb), FL(Nf, Norb, Norb))
      allocate (SigR2(Norb, Norb, n_omega), SigL2(Norb, Norb, n_omega))

      !....plans are made once, outside any parallel region (planning is not
      !....thread-safe; dfftw_execute_dft on other arrays is)
      allocate (tmp(Nf))
      call dfftw_plan_dft_1d(plan_fwd, Nf, tmp, tmp, FFTW_FORWARD, FFTW_ESTIMATE + FFTW_UNALIGNED)
      call dfftw_plan_dft_1d(plan_bwd, Nf, tmp, tmp, FFTW_BACKWARD, FFTW_ESTIMATE + FFTW_UNALIGNED)
      deallocate (tmp)
      write (*, *) 'Sigma2FFT: n_omega =', n_omega, '  padded FFT length Nf =', Nf
   end subroutine sigma2_init

   subroutine sigma2_finalize()
      implicit none
      call dfftw_destroy_plan(plan_fwd)
      call dfftw_destroy_plan(plan_bwd)
      deallocate (Fr, FL, SigR2, SigL2)
   end subroutine sigma2_finalize

   !....zero-pad along omega and forward-transform, in place
   subroutine load_fft(G, ii, jj, dst)
      !....element (ii,jj) of G(:,:,1:n_omega), zero-padded to length Nf, then
      !....forward-transformed in place. Reads G directly (no array temporary
      !....for the strided slice G(ii,jj,:)).
      implicit none
      complex(dp), intent(in)    :: G(:, :, :)
      integer, intent(in)        :: ii, jj
      complex(dp), intent(inout) :: dst(Nf)
      integer :: iw
      do iw = 1, n_omega
         dst(iw) = G(ii, jj, iw)
      end do
      dst(n_omega + 1:Nf) = (0.d0, 0.d0)
      call dfftw_execute_dft(plan_fwd, dst, dst)
   end subroutine load_fft

   !....Call once per SCF iteration (after G_in is updated), BEFORE the loop
   !....over interacting_gf. Fills SigR2(:,:,iw) and SigL2(:,:,iw) for all iw
   !....(second-order parts only; Sigma1 is added in sigma_retarded).
   subroutine sigma2_compute(do_lesser)
      implicit none
      logical, intent(in) :: do_lesser
      integer :: ii, jj, io, jo, l, lm, sp, sp1
      real(dp)  :: pp
      complex(dp) :: c
      complex(dp), dimension(2, 2) :: r_ij, L_ij, a_ij, G_ij ! blocks at  +l
      complex(dp), dimension(2, 2) :: r_ji, L_ji, a_ji, G_ji ! blocks at  -l
      complex(dp), dimension(2, 2) :: chi1, chi2, chiL, sR, sL
      complex(dp), allocatable :: bufR(:, :, :), bufL(:, :, :)

      pp = delta/(2.d0*pi)

      !....1. forward FFT of every matrix element of G^r and G^< only
      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(ii,jj) SCHEDULE(static)
      do jj = 1, Norb
         do ii = 1, Norb
            call load_fft(G_in%r, ii, jj, Fr(:, ii, jj))
            call load_fft(G_in%L, ii, jj, FL(:, ii, jj))
         end do
      end do
      !$OMP END PARALLEL DO

      !....2. per orbital pair: 2x2 spin algebra at every l, then inverse FFT
      !$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(dynamic) &
      !$OMP& PRIVATE(io,jo,l,lm,sp,sp1,c,r_ij,L_ij,a_ij,G_ij,r_ji,L_ji,a_ji,G_ji, &
      !$OMP&         chi1,chi2,chiL,sR,sL,bufR,bufL)
      do jo = 1, Norb, 2
         do io = 1, Norb, 2

            if (U_site(io)*U_site(jo) == 0.d0) then
               SigR2(io:io + 1, jo:jo + 1, :) = (0.d0, 0.d0)
               SigL2(io:io + 1, jo:jo + 1, :) = (0.d0, 0.d0)
               cycle
            end if

            allocate (bufR(Nf, 2, 2), bufL(Nf, 2, 2))

            do l = 1, Nf
               lm = mod(Nf - (l - 1), Nf) + 1 ! index of -l  (time reversal)

               !....the only four loads: G^r and G^< blocks, (i,j) at +l and (j,i) at -l
               r_ij = Fr(l, io:io + 1, jo:jo + 1); L_ij = FL(l, io:io + 1, jo:jo + 1)
               r_ji = Fr(lm, jo:jo + 1, io:io + 1); L_ji = FL(lm, jo:jo + 1, io:io + 1)

               !....advanced and greater blocks, derived (no FFT, no storage)
               a_ij = conjg(transpose(r_ji)) ! a_ij(+l) = [r_ji(-l)]^dagger
               a_ji = conjg(transpose(r_ij)) ! a_ji(-l) = [r_ij(+l)]^dagger
               G_ij = L_ij + r_ij - a_ij ! G^> = G^< + G^r - G^a
               G_ji = L_ji + r_ji - a_ji

               !....Retarded Self-energy
               chi1 = matmul(L_ji, G_ij)
               chi2 = matmul(L_ji, r_ij) + matmul(a_ji, L_ij)
               sR = matmul(r_ij, adj2(chi1)) + matmul(L_ij, adj2(chi2))
               bufR(l, :, :) = sR

               if (do_lesser) then
                  chiL = matmul(G_ji, L_ij)
                  sL = matmul(L_ij, adj2(chiL))
                  bufL(l, :, :) = sL
               end if
            end do

            c = U_site(io)*U_site(jo)*hbar**2*pp*pp/dble(Nf) !...prefactors added here

            do sp1 = 1, 2
               do sp = 1, 2
                  call dfftw_execute_dft(plan_bwd, bufR(:, sp, sp1), bufR(:, sp, sp1))
                  SigR2(io + sp - 1, jo + sp1 - 1, :) = c*bufR(1:n_omega, sp, sp1)
                  if (do_lesser) then
                     call dfftw_execute_dft(plan_bwd, bufL(:, sp, sp1), bufL(:, sp, sp1))
                     SigL2(io + sp - 1, jo + sp1 - 1, :) = c*bufL(1:n_omega, sp, sp1)
                  end if
               end do
            end do

            deallocate (bufR, bufL)
         end do
      end do
      !$OMP END PARALLEL DO

   end subroutine sigma2_compute

   !....adjugate of a 2x2:  Tr(M)*1 - M   (direct minus exchange)
   pure function adj2(M) result(A)
      implicit none
      complex(dp), intent(in) :: M(2, 2)
      complex(dp) :: A(2, 2)
      A(1, 1) = M(2, 2); A(1, 2) = -M(1, 2)
      A(2, 1) = -M(2, 1); A(2, 2) = M(1, 1)
   end function adj2

end module SelfEnergies
