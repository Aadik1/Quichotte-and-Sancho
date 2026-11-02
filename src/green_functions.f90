!==============================================================================
! src/green_functions.f90  --  module GreenFunctions   (QUICHOTTE)
!
! Storage for the Green's functions and the quantities derived from them.
!
! MEMORY LAYOUT (see the manual, "Memory and the active subspace")
!   Only the ACTIVE states are stored on the energy grid: both spin states of
!   every site (orbital) that carries a Hubbard U, plus those of the contact
!   sites (non-zero GammaL/GammaR). The interaction self-energy lives on the
!   interacting sites, so the self-consistency closes on the active block:
!     order 2: G^r, G^< of the active block on the whole grid (G_in)
!     order 1: only G^<(t=0) on the interacting sites (G_less_t0)
!     order 0: nothing (no self-consistency)
!   The full Norb x Norb G^r exists only as a per-energy, per-thread temporary
!   (scf.f90, G_at_energy). G^a = (G^r)^dagger and G^> = G^< + G^r - G^a are
!   formed where needed and never stored. The observables (current, spin
!   texture, spectral function) come from one FINAL SWEEP with the converged
!   self-energies, which covers every site.
!
!   act(k)  full index of active state k (ascending); pos(i) = k, 0 if inactive
!   sup(:)  support of the lesser self-energy: the active states at order 2,
!           the contact states at orders 0 and 1
!   cl(:)   left-contact states (non-zero rows/columns of GammaL), for the current
!
! Routines:
!   fermi                Fermi-Dirac distribution
!   setup_storage        active subspace and allocations (after the grid is known)
!   equal_time_lesser    G^<(t=0) on the active sites, from G_in   (order 2)
!   gather_onsite        G^<(t=0) on chosen sites, from the per-energy blocks onsL
!==============================================================================
module GreenFunctions
   use Constants, only: dp, hbar, pi, im
   use Parameters
   use Hamiltonian, only: H
   implicit none

   type :: GreenFunctionSet
      complex(dp), allocatable, dimension(:, :, :) :: r, L   !...(nA, nA, n_omega), order 2 only
   end type GreenFunctionSet
   type(GreenFunctionSet) :: G_in        !...INPUT Green's functions of the current SCF step

   real(dp), allocatable, dimension(:) :: omega    !...energy grid (n_omega)
   real(dp), allocatable, dimension(:) :: U_site   !...Hubbard U per spin-orbital (set on the spin-up index)
   real(dp), allocatable, dimension(:) :: Ev       !...eigenvalues of H
   complex(dp), allocatable, dimension(:, :) :: Eigenvec
   complex(dp), allocatable, dimension(:, :) :: GammaL, GammaR   !...lead coupling matrices (Norb, Norb)

   !...equal-time lesser GF, on-site 2x2 spin blocks (Norb, Norb):
   complex(dp), allocatable, dimension(:, :) :: G_less_t0   !...INPUT of the SCF step (enters Sigma1)
   complex(dp), allocatable, dimension(:, :) :: G_less_obs  !...from the final sweep, every site (observables)

   !...active subspace
   integer :: nA, nS, nCL, nAs, n_int_sites
   integer, allocatable :: act(:), pos(:), sup(:), cl(:), asite(:)
   logical, allocatable :: site_int(:), site_act(:)

   !...fixed matrices: Hemb = -H + i(GammaL+GammaR)/2, first-order Sigma, lead couplings on sup / cl
   complex(dp), allocatable, dimension(:, :) :: Hemb, Sig1, GamL_s, GamR_s, GamL_cl

   !...per-energy results of a sweep (memory linear in the number of states)
   complex(dp), allocatable :: gRd(:, :), gLd(:, :)   !...(nA, n_omega) diagonals of the output G^r, G^<
   complex(dp), allocatable :: onsL(:, :)             !...(4*N_sites, n_omega) on-site 2x2 blocks of G^<
   real(dp), allocatable :: dosd(:, :)                !...(Norb, n_omega) Im G^r_ii (final sweep)
   real(dp), allocatable :: trU(:), trD(:)            !...(n_omega) current integrand, spin up / down

   integer :: INFO

contains

   !.....fermi-dirac distribution
   real(dp) function fermi(w, V)
      implicit none
      real(dp) :: arg, w, V

      arg = (w - mu + V/hbar)*beta
      !...exponent capped at 700 (exp(700) ~ 1e304): with wide spectra (e.g. DFTB, ~50 eV)
      !...exp would overflow to +Inf. The result is unchanged (f = 1e-304 instead of 0)
      !...but the debug build, which traps overflow, can run.
      fermi = 1.d0/(exp(min(hbar*arg, 700.d0)) + 1.d0)

   end function fermi

   !===========================================================================
   ! setup_storage: active subspace, allocations and a memory estimate.
   ! Call once, after H, GammaL/R, U_site and the energy grid are set.
   !===========================================================================
   subroutine setup_storage(n_vec_out, memory_gb)
      integer, intent(out) :: n_vec_out
      real(dp), intent(out) :: memory_gb
      logical, allocatable :: conL(:), conR(:)
      integer :: i, k, s, nfull, nsu
      real(dp) :: bytes

      allocate (site_int(N_sites), site_act(N_sites), pos(Norb), conL(Norb), conR(Norb))

      !.....interacting sites: U /= 0 and order >= 1
      do s = 1, N_sites
         site_int(s) = (order .ge. 1 .and. U_site(2*s - 1) .ne. 0.d0)
      end do
      n_int_sites = count(site_int)

      !.....contact states (non-zero row or column of Gamma); their sites are active too
      site_act = site_int
      do i = 1, Norb
         conL(i) = any(GammaL(i, :) .ne. (0.d0, 0.d0)) .or. any(GammaL(:, i) .ne. (0.d0, 0.d0))
         conR(i) = any(GammaR(i, :) .ne. (0.d0, 0.d0)) .or. any(GammaR(:, i) .ne. (0.d0, 0.d0))
         if (conL(i) .or. conR(i)) site_act((i + 1)/2) = .true.
      end do

      !.....active states and sites, ascending
      nA = 2*count(site_act)
      allocate (act(nA)); pos = 0; k = 0
      do i = 1, Norb
         if (site_act((i + 1)/2)) then
            k = k + 1; act(k) = i; pos(i) = k
         end if
      end do
      nAs = count(site_act)
      allocate (asite(nAs)); k = 0
      do s = 1, N_sites
         if (site_act(s)) then
            k = k + 1; asite(k) = s
         end if
      end do

      !.....left-contact states, and the support of the lesser self-energy
      nCL = count(conL)
      allocate (cl(nCL)); k = 0
      do i = 1, Norb
         if (conL(i)) then
            k = k + 1; cl(k) = i
         end if
      end do
      if (order .eq. 2) then
         nS = nA
         allocate (sup(nS)); sup = act
      else
         nS = count(conL .or. conR)
         allocate (sup(nS)); k = 0
         do i = 1, Norb
            if (conL(i) .or. conR(i)) then
               k = k + 1; sup(k) = i
            end if
         end do
      end if
      allocate (GamL_s(nS, nS), GamR_s(nS, nS), GamL_cl(nCL, nCL))
      GamL_s = GammaL(sup, sup); GamR_s = GammaR(sup, sup); GamL_cl = GammaL(cl, cl)

      allocate (Hemb(Norb, Norb), Sig1(Norb, Norb), G_less_t0(Norb, Norb), G_less_obs(Norb, Norb))
      Sig1 = (0.d0, 0.d0); G_less_t0 = (0.d0, 0.d0); G_less_obs = (0.d0, 0.d0)
      allocate (onsL(4*N_sites, n_omega), dosd(Norb, n_omega), trU(n_omega), trD(n_omega))
      onsL = (0.d0, 0.d0)
      if (order .ge. 1) allocate (gRd(nA, n_omega), gLd(nA, n_omega))

      !.....length of the SCF vector
      n_vec_out = 0
      if (order .eq. 1) then
         n_vec_out = 8*n_int_sites                 !...Re, Im of the 2x2 block of each interacting site
      else if (order .eq. 2) then
         nfull = nA*nA; nsu = (nfull - nA)/2
         n_vec_out = (2*nfull + nsu + nsu + nA)*n_omega
         allocate (G_in%r(nA, nA, n_omega), G_in%L(nA, nA, n_omega))
      end if

      !.....memory estimate (stored arrays; see the manual)
      bytes = 16.d0*4*N_sites*n_omega + 8.d0*Norb*n_omega + 16.d0*n_omega
      if (order .ge. 1) bytes = bytes + 2*16.d0*nA*n_omega
      if (order .eq. 2) bytes = bytes + 2*16.d0*nA*nA*n_omega &        !...G_in
                              + 2*16.d0*nA*nA*n_omega &                !...SigR2, SigL2
                              + 2*16.d0*nA*nA*(2*n_omega)              !...FFT arrays (padded ~2x)
      memory_gb = bytes/1024.d0**3

      write (*, '(a,i0,a,i0,a,i0,a,i0,a)') ' Active states: ', nA, ' of ', Norb, '  (', nAs, ' of ', &
         N_sites, ' sites)'
      write (*, '(a,i0,a,i0)') ' Interacting sites: ', n_int_sites, '   lesser self-energy support: ', nS
      deallocate (conL, conR)
   end subroutine setup_storage

   subroutine equal_time_lesser()
      !.....G_less_t0 = G^<(t=0) = (delta/2pi) sum_w G^<(w) from the input G_in,
      !.....on-site 2x2 spin blocks of the active sites (all that Sigma1 needs)
      implicit none
      integer :: is, a, b, k1
      complex(dp) :: s
      real(dp) :: pp

      pp = delta/(2.d0*pi)
      G_less_t0 = (0.d0, 0.d0)
      do is = 1, nAs
         do a = 2*asite(is) - 1, 2*asite(is)
            do b = 2*asite(is) - 1, 2*asite(is)
               s = (0.d0, 0.d0)
               do k1 = 1, n_omega
                  s = s + G_in%L(pos(a), pos(b), k1)
               end do
               G_less_t0(a, b) = s*pp
            end do
         end do
      end do
   end subroutine equal_time_lesser

   subroutine gather_onsite(mask, G)
      !.....G(a,b) = (delta/2pi) sum_w G^<_ab(w) on the 2x2 block of every masked site,
      !.....from the per-energy blocks onsL written by a sweep
      logical, intent(in) :: mask(N_sites)
      complex(dp), intent(inout) :: G(Norb, Norb)
      integer :: s, a, b, k, iw
      complex(dp) :: acc
      real(dp) :: pp

      pp = delta/(2.d0*pi)
      do s = 1, N_sites
         if (.not. mask(s)) cycle
         do b = 1, 2
            do a = 1, 2
               k = 4*(s - 1) + 2*(b - 1) + a
               acc = (0.d0, 0.d0)
               do iw = 1, n_omega
                  acc = acc + onsL(k, iw)
               end do
               G(2*s - 2 + a, 2*s - 2 + b) = acc*pp
            end do
         end do
      end do
   end subroutine gather_onsite

end module GreenFunctions
