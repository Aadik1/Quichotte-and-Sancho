!==============================================================================
! src/scf.f90  --  module SCF   (QUICHOTTE)
!
! Self-consistent solution at one voltage:
!   solve_scf          driver: starting point, iterate to convergence, final sweep
!   scf_residual2      one SCF map at order 2:  x = active-block G^r, G^< on the grid
!   scf_residual1      one SCF map at order 1:  x = G^<(t=0) on the interacting sites
!   sweep              one pass over the energy grid (OpenMP over energies)
!   G_at_energy        Dyson + Keldysh equations at one energy, keeping only what
!                      the current mode needs (see green_functions.f90 for the layout)
!
! The iteration itself is common/broyden_mixing.f90. Order 1 mixes G^<(t=0)
! because Sigma1 is static and depends on nothing else: the SCF vector is
! 8 numbers per interacting site instead of ~3 Norb^2 per energy.
! After convergence one FINAL sweep with the converged self-energies gives the
! current, the spin texture and the spectral function of every site.
!==============================================================================
module SCF
   use Constants, only: dp, hbar, pi, im
   use LinearAlgebra, only: invert_complex, pack_gf, unpack_gf
   use BroydenMixing, only: broyden_solve
   use ScfChecks, only: check_diag_sanity
   use Parameters
   use Hamiltonian, only: H
   use GreenFunctions
   use SelfEnergies
   use OMP_LIB
   implicit none

   !.....SCF settings (all settable in input.dat)
   !.....  alpha_mix    - linear-mixing weight / initial inverse Jacobian (-alpha_mix*I)
   !.....  nhist        - Broyden history length. nhist = 0 is plain linear mixing
   !.....  maxscf       - maximum SCF steps per voltage
   !.....  tolf         - order 2: converged when max|x_out - x| < tolf  (x = G^r, G^< on the grid)
   !.....  tolf_order1  - order 1: converged when max|x_out - x| < tolf_order1  (x = G^<(t=0))
   !.....  continuation - start each voltage from the previous converged solution
   !.....  broyden_single_precision - keep the Broyden history in single precision (half the memory)
   real(dp) :: alpha_mix = 0.5d0
   integer  :: nhist = 16, maxscf = 300
   real(dp) :: tolf = 1.d-6
   real(dp) :: tolf_order1 = 1.d-8
   logical  :: continuation = .true.
   logical  :: broyden_single_precision = .true.

   !.....order 2: packed vector per energy = 2*n_full + n_strict_upper + n_upper  (n_full = nA^2)
   integer :: n_full, n_strict_upper, n_upper, nblk, n_vec

   !.....modes of G_at_energy
   integer, parameter :: M_INIT2 = 1, M_SCF2 = 2, M_SCF1 = 3, M_FINAL = 4
   logical :: use_sigma = .false.

   !.....state shared between solve_scf and the residuals
   real(dp) :: V_scf
   integer  :: scf_iter

contains

!-----------------------------------------------------------------------------
! scf_setup_sizes: lengths of the packed order-2 vector (call after setup_storage)
!-----------------------------------------------------------------------------
   subroutine scf_setup_sizes(n_total)
      integer, intent(in) :: n_total
      n_full = nA*nA
      n_strict_upper = (n_full - nA)/2
      n_upper = n_strict_upper + nA
      nblk = 2*n_full + n_strict_upper + n_upper
      n_vec = n_total
   end subroutine scf_setup_sizes

!-----------------------------------------------------------------------------
! solve_scf: converge at voltage Volt, then run the final sweep.
!   fresh_start = .true.  : start from the NON-interacting solution
!   fresh_start = .false. : start from the previous voltage's converged input
!   converged   = .false. if the SCF did not converge, ran away (NaN, occupations
!                 far outside [0,1]) or converged to an unphysical solution.
!   On exit G_less_obs, dosd, trU, trD hold the observables of every site.
!-----------------------------------------------------------------------------
   subroutine solve_scf(Volt, fresh_start, converged)
      real(dp), intent(in) :: Volt
      logical, intent(in) :: fresh_start
      logical, intent(out) :: converged
      real(dp), allocatable :: x(:)
      real(dp) :: dummy(1)
      character(len=30) :: fn1
      logical :: runaway, unphysical
      integer :: n_iter, i_worst, s
      real(dp) :: fmax, n_worst, w_worst

      write (u_runtime, *) '........SCF Calculations at Voltage:', Volt, '..........'
      print *, '>>>>>>>>>>VOLTAGE:', Volt, merge(' (fresh start)          ', ' (from previous voltage)', fresh_start)

      Hemb = -H + 0.5d0*(im/hbar)*(GammaL + GammaR) !LK <========= must be +
      V_scf = Volt
      scf_iter = 0
      converged = .true.; n_iter = 0; fmax = 0.d0

      !.....one convergence log per voltage, e.g. scf_V_-0.950.dat
      write (fn1, '(f8.3)') Volt
      open (u_scflog, file='scf_V_'//trim(adjustl(fn1))//'.dat', status='replace')

      select case (order)
      case (1)
         if (fresh_start) then
            !.....G^<(t=0) of the non-interacting system
            G_less_t0 = (0.d0, 0.d0); use_sigma = .false.
            call sweep(M_SCF1, dummy)
            call gather_onsite(site_int, G_less_t0)
         end if
         allocate (x(n_vec))
         call pack_nil(G_less_t0, x)
         call broyden_solve(scf_residual1, x, n_vec, alpha_mix, nhist, maxscf, tolf_order1, u_scflog, &
                            converged, n_iter, fmax, single_history=broyden_single_precision)
         deallocate (x)
      case (2)
         if (.not. allocated(SigR2)) call sigma2_init()
         if (fresh_start) then
            !.....non-interacting active-block G^r, G^< on the whole grid
            use_sigma = .false.
            call sweep(M_INIT2, dummy)
         end if
         allocate (x(n_vec))
         call pack_gf(G_in%r, G_in%L, x, n_vec, n_full, n_strict_upper, nA, n_omega)
         call broyden_solve(scf_residual2, x, n_vec, alpha_mix, nhist, maxscf, tolf, u_scflog, &
                            converged, n_iter, fmax, single_history=broyden_single_precision)
         deallocate (x)
      case default
         !.....order 0: no self-consistency
         write (u_scflog, '(a)') '#  step        max|F|          rms|F|   nhist'
         write (u_scflog, '(i6,2es16.6,i6)') 1, 0.d0, 0.d0, 0
         n_iter = 1
      end select
      close (u_scflog)

      !.....final sweep with the converged self-energies (Sig1, SigR2, SigL2 of the last input)
      use_sigma = (order .ge. 1)
      call sweep(M_FINAL, dummy)
      G_less_obs = (0.d0, 0.d0)
      call gather_onsite([(.true., s=1, N_sites)], G_less_obs)
      call final_sanity(runaway, unphysical, i_worst, n_worst, w_worst)

      if (.not. converged) then
         write (*, '(a,i5,a,es11.3)') ' NOT CONVERGED after', n_iter, ' steps, max|F| =', fmax
      else if (unphysical) then
         converged = .false.
         write (*, '(a,i4,a,es11.3,a,es11.3)') ' CONVERGED TO AN UNPHYSICAL SOLUTION: state', i_worst, &
            ' n =', n_worst, ', spectral weight =', w_worst
      else
         write (*, '(a,i5,a)') ' ... CONVERGED in', n_iter, ' steps ...'
      end if
   end subroutine solve_scf

!-----------------------------------------------------------------------------
! scf_residual2 (order 2): input G_in = unpack(x); build its self-energies,
! solve Dyson/Keldysh at every energy, write pack(G_out) straight into F, F = F - x.
! abort = .true. if G_out ran away (NaN, occupations far outside [0,1]).
!-----------------------------------------------------------------------------
   subroutine scf_residual2(n, x, F, abort)
      integer, intent(in) :: n
      real(dp), intent(in) :: x(n)
      real(dp), intent(out) :: F(n)
      logical, intent(out) :: abort
      integer :: i_worst
      logical :: unphysical
      real(dp) :: n_worst, w_worst

      scf_iter = scf_iter + 1
      write (*, *) '.... ITERATION = ', scf_iter, ' ....'

      !.....input GFs from x
      call unpack_gf(G_in%R, G_in%L, x, n, n_full, n_strict_upper, nA, n_omega)

      !.....self-energies and output GFs
      call equal_time_lesser()
      call sigma_first_order(Sig1)
      call sigma2_compute(.true.)
      use_sigma = .true.
      call sweep(M_SCF2, F)

      call check_diag_sanity(gRd, gLd, nA, n_omega, delta, abort, unphysical, i_worst, n_worst, w_worst, act)
      if (abort) then
         write (*, '(a,i5,a,i4,a,es11.3,a,es11.3)') ' SCF RUNAWAY at iteration', scf_iter, &
            ': state', i_worst, ' n =', n_worst, ', spectral weight =', w_worst
         F = 0.d0
         return
      end if
      F = F - x
   end subroutine scf_residual2

!-----------------------------------------------------------------------------
! scf_residual1 (order 1): input G^<(t=0) = unpack(x) on the interacting sites;
! Sigma1 from it; one sweep gives the output G^<(t=0); F = pack(output) - x.
!-----------------------------------------------------------------------------
   subroutine scf_residual1(n, x, F, abort)
      integer, intent(in) :: n
      real(dp), intent(in) :: x(n)
      real(dp), intent(out) :: F(n)
      logical, intent(out) :: abort
      integer :: i_worst
      logical :: unphysical
      real(dp) :: n_worst, w_worst
      complex(dp), allocatable :: Gout(:, :)

      scf_iter = scf_iter + 1
      write (*, *) '.... ITERATION = ', scf_iter, ' ....'

      call unpack_nil(x, G_less_t0)
      call sigma_first_order(Sig1)
      use_sigma = .true.
      call sweep(M_SCF1, F)

      call check_diag_sanity(gRd, gLd, nA, n_omega, delta, abort, unphysical, i_worst, n_worst, w_worst, act)
      if (abort) then
         write (*, '(a,i5,a,i4,a,es11.3,a,es11.3)') ' SCF RUNAWAY at iteration', scf_iter, &
            ': state', i_worst, ' n =', n_worst, ', spectral weight =', w_worst
         F = 0.d0
         return
      end if
      allocate (Gout(Norb, Norb)); Gout = (0.d0, 0.d0)
      call gather_onsite(site_int, Gout)
      call pack_nil(Gout, F)
      F = F - x
      deallocate (Gout)
   end subroutine scf_residual1

   subroutine pack_nil(G, x)
      !...order-1 SCF vector: Re, Im of the 2x2 block of G^<(t=0) of each interacting site
      complex(dp), intent(in) :: G(Norb, Norb)
      real(dp), intent(out) :: x(*)
      integer :: s, a, b, k
      k = 0
      do s = 1, N_sites
         if (.not. site_int(s)) cycle
         do b = 2*s - 1, 2*s
            do a = 2*s - 1, 2*s
               x(k + 1) = real(G(a, b)); x(k + 2) = aimag(G(a, b)); k = k + 2
            end do
         end do
      end do
   end subroutine pack_nil

   subroutine unpack_nil(x, G)
      real(dp), intent(in) :: x(*)
      complex(dp), intent(out) :: G(Norb, Norb)
      integer :: s, a, b, k
      G = (0.d0, 0.d0); k = 0
      do s = 1, N_sites
         if (.not. site_int(s)) cycle
         do b = 2*s - 1, 2*s
            do a = 2*s - 1, 2*s
               G(a, b) = cmplx(x(k + 1), x(k + 2), kind=dp); k = k + 2
            end do
         end do
      end do
   end subroutine unpack_nil

!-----------------------------------------------------------------------------
! final_sanity: occupation and spectral weight of EVERY state after the final sweep
!-----------------------------------------------------------------------------
   subroutine final_sanity(runaway, unphysical, i_worst, n_worst, w_worst)
      logical, intent(out) :: runaway, unphysical
      integer, intent(out) :: i_worst
      real(dp), intent(out) :: n_worst, w_worst
      complex(dp), allocatable :: grd_all(:, :), gld_all(:, :)
      integer :: i, s, a

      allocate (grd_all(Norb, n_omega), gld_all(Norb, n_omega))
      grd_all = cmplx(0.d0, dosd, kind=dp)
      do i = 1, Norb
         s = (i + 1)/2; a = i - 2*(s - 1)
         gld_all(i, :) = onsL(4*(s - 1) + 2*(a - 1) + a, :)
      end do
      call check_diag_sanity(grd_all, gld_all, Norb, n_omega, delta, runaway, unphysical, i_worst, n_worst, w_worst)
      deallocate (grd_all, gld_all)
   end subroutine final_sanity

!-----------------------------------------------------------------------------
! sweep: G_at_energy at every energy (OpenMP over energies)
!-----------------------------------------------------------------------------
   subroutine sweep(mode, F)
      integer, intent(in) :: mode
      real(dp) :: F(*)
      integer :: iw

      !$OMP PARALLEL DO PRIVATE(iw)
      do iw = 1, n_omega
         call G_at_energy(iw, V_scf, mode, F)
      end do
      !$OMP END PARALLEL DO
   end subroutine sweep

!-----------------------------------------------------------------------------
! G_at_energy: Dyson and Keldysh equations at omega(iw),
!   G^r = [ hbar*w - H + (i/2hbar)(GammaL+GammaR) - Sigma^r ]^-1
!   G^< = G^r [ (i/hbar)(f_L GammaL + f_R GammaR) + Sigma^< ] G^a
! with the full Norb x Norb G^r as a per-thread temporary. Sigma^r lives on the
! active block, Sigma^< on sup, so G^<(X,Y) = G^r(X,sup) Sigma^<(sup,sup) [G^r(Y,sup)]^dagger.
! Kept, by mode:
!   M_INIT2  non-interacting active-block G^r, G^< -> G_in            (order 2 start)
!   M_SCF2   output active-block G^r, G^< packed straight into F, + diagonals
!   M_SCF1   output on-site G^< of the interacting sites, + diagonals of the active states
!   M_FINAL  on-site G^< of every site, Im G^r_ii of every state, current integrand
! Called inside an OpenMP loop: uses only local work arrays.
!-----------------------------------------------------------------------------
   subroutine G_at_energy(iw, Volt, mode, F)
      integer, intent(in) :: iw, mode
      real(dp), intent(in) :: Volt
      real(dp) :: F(*)
      integer :: i, k, s, a, lo, p1, off, info_loc
      real(dp) :: w, fL, fR
      complex(dp) :: trc_u, trc_d
      complex(dp), allocatable, dimension(:, :) :: w1, SL, T, blk, GLc, X, Mc
      complex(dp), allocatable, dimension(:, :, :) :: grA, glA

      w = omega(iw)
      allocate (w1(Norb, Norb))

      !............retarded
      w1 = Hemb
      if (use_sigma) then
         if (order .eq. 2) then
            w1(act, act) = w1(act, act) - (Sig1(act, act) + SigR2(:, :, iw))
         else
            w1(act, act) = w1(act, act) - Sig1(act, act)
         end if
      end if
      do i = 1, Norb
         w1(i, i) = w1(i, i) + hbar*w
      end do
      call invert_complex(Norb, w1, info_loc)

      !............lesser self-energy on sup: lead embedding (+ second order)
      fL = fermi(w, Volt); fR = fermi(w, 0.d0)
      allocate (SL(nS, nS))
      SL = im*(fL*GamL_s + fR*GamR_s)/hbar
      if (use_sigma .and. order .eq. 2) SL = SL + SigL2(:, :, iw)

      select case (mode)

      case (M_INIT2, M_SCF2)                  !...order 2: sup = act
         allocate (grA(nA, nA, 1), glA(nA, nA, 1), T(nA, nA))
         grA(:, :, 1) = w1(act, act)
         T = matmul(w1(act, act), SL)
         glA(:, :, 1) = matmul(T, conjg(transpose(grA(:, :, 1))))   !.. GL = Gr * SigmaL * Ga
         if (mode .eq. M_INIT2) then
            G_in%r(:, :, iw) = grA(:, :, 1)
            G_in%L(:, :, iw) = glA(:, :, 1)
         else
            off = (iw - 1)*nblk
            call pack_gf(grA, glA, F(off + 1:off + nblk), nblk, n_full, n_strict_upper, nA, 1)
            do k = 1, nA
               gRd(k, iw) = grA(k, k, 1); gLd(k, iw) = glA(k, k, 1)
            end do
         end if
         deallocate (grA, glA, T)

      case (M_SCF1)
         allocate (T(nA, nS), blk(2, 2))
         T = matmul(w1(act, sup), SL)
         do k = 1, nA
            gRd(k, iw) = w1(act(k), act(k))
            gLd(k, iw) = sum(T(k, :)*conjg(w1(act(k), sup)))
         end do
         do s = 1, N_sites
            if (.not. site_int(s)) cycle
            lo = 2*s - 1; p1 = pos(lo)
            blk = matmul(T(p1:p1 + 1, :), conjg(transpose(w1(lo:lo + 1, sup))))
            onsL(4*(s - 1) + 1:4*s, iw) = reshape(blk, [4])
         end do
         deallocate (T, blk)

      case (M_FINAL)
         allocate (T(Norb, nS), blk(2, 2))
         T = matmul(w1(:, sup), SL)
         do i = 1, Norb
            dosd(i, iw) = aimag(w1(i, i))
         end do
         do s = 1, N_sites
            lo = 2*s - 1
            blk = matmul(T(lo:lo + 1, :), conjg(transpose(w1(lo:lo + 1, sup))))
            onsL(4*(s - 1) + 1:4*s, iw) = reshape(blk, [4])
         end do
         !.....current integrand through the left lead (left-contact block):
         !.....Tr_s { i GammaL [ (f_L - 1) G^< - f_L G^> ] } / 2pi,  G^> = G^< + G^r - G^a
         allocate (GLc(nCL, nCL), X(nCL, nCL), Mc(nCL, nCL))
         GLc = matmul(T(cl, :), conjg(transpose(w1(cl, sup))))
         do a = 1, nCL
            do k = 1, nCL
               X(k, a) = (fL - 1.d0)*GLc(k, a) - fL*(GLc(k, a) + w1(cl(k), cl(a)) - conjg(w1(cl(a), cl(k))))
            end do
         end do
         Mc = im*matmul(GamL_cl, X)
         trc_u = (0.d0, 0.d0); trc_d = (0.d0, 0.d0)
         do a = 1, nCL
            if (mod(cl(a), 2) .eq. 1) then
               trc_u = trc_u + Mc(a, a)          !...spin up
            else
               trc_d = trc_d + Mc(a, a)          !...spin down
            end if
         end do
         trU(iw) = real(trc_u)/(2.d0*pi)
         trD(iw) = real(trc_d)/(2.d0*pi)
         deallocate (T, blk, GLc, X, Mc)

      end select

      deallocate (w1, SL)
   end subroutine G_at_energy

end module SCF
