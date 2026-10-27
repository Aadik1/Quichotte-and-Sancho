!==============================================================================
! src/scf.f90  --  module SCF   (QUICHOTTE)
!
! Self-consistent solution at one voltage:
!   solve_scf         driver: choose the starting point, iterate to convergence
!   scf_residual      one SCF map  x -> F(x) = x_out(x) - x
!   interacting_gf    Dyson + Keldysh equations at one energy (gives G_out)
!
! x is G_in packed into one real vector by pack_gf (common/linear_algebra.f90);
! the iteration itself is common/broyden_mixing.f90.
!==============================================================================
module SCF
   use Constants, only: dp, hbar, pi, im
   use LinearAlgebra, only: invert_complex, hermitian_conjugate, pack_gf, unpack_gf
   use BroydenMixing, only: broyden_solve
   use ScfChecks, only: check_gf_sanity
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
   !.....  tolf         - converged when max|x_out - x| < tolf
   !.....  continuation - start each voltage from the previous converged GFs
   real(dp) :: alpha_mix = 0.5d0
   integer  :: nhist = 16, maxscf = 300
   real(dp) :: tolf = 1.d-6
   logical  :: continuation = .true.

   !.....length of the packed vector, per energy: n_full + n_full + n_strict_upper + n_upper
   integer :: n_full, n_strict_upper, n_upper, n_vec

   !.....state shared between solve_scf and scf_residual
   real(dp) :: V_scf
   integer  :: scf_iter

contains

!-----------------------------------------------------------------------------
! solve_scf: converge the Green's functions at voltage Volt.
!   fresh_start = .true.  : start from the NON-interacting GFs
!   fresh_start = .false. : start from the current G_in (previous voltage)
!   converged   = .false. if the SCF did not converge, ran away (NaN, occupations
!                 far outside [0,1]) or converged to an unphysical solution.
!   On exit G_in holds the final GFs and G_less_t0 = G^<(t=0).
!-----------------------------------------------------------------------------
   subroutine solve_scf(Volt, fresh_start, converged)
      real(dp), intent(in) :: Volt
      logical, intent(in) :: fresh_start
      logical, intent(out) :: converged
      real(dp), allocatable :: x(:)
      character(len=30) :: fn1
      logical :: runaway, unphysical
      integer :: n_iter, i_worst
      real(dp) :: fmax, n_worst, w_worst

      write (u_runtime, *) '........SCF Calculations at Voltage:', Volt, '..........'
      print *, '>>>>>>>>>>VOLTAGE:', Volt, merge(' (fresh start)          ', ' (from previous voltage)', fresh_start)

      if (fresh_start) then
         call noninteracting_gr_ga()
         call noninteracting_gl_gg(Volt)
      end if
      if (order .eq. 2 .and. .not. allocated(SigR2)) call sigma2_init()

      !.....one convergence log per voltage, e.g. scf_V_-0.950.dat
      write (fn1, '(f8.3)') Volt
      open (u_scflog, file='scf_V_'//trim(adjustl(fn1))//'.dat', status='replace')

      allocate (x(n_vec))
      call pack_gf(G_in%R, G_in%L, x, n_vec, n_full, n_strict_upper, Norb, n_omega)
      V_scf = Volt
      scf_iter = 0
      call broyden_solve(scf_residual, x, n_vec, alpha_mix, nhist, maxscf, tolf, u_scflog, &
                         converged, n_iter, fmax)
      deallocate (x)
      close (u_scflog)

      !.....G_in = input of the last scf_residual call, G_out = its output
      call check_gf_sanity(G_out%R, G_out%L, Norb, n_omega, delta, runaway, unphysical, i_worst, n_worst, w_worst)
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
! scf_residual: one SCF step. Input G_in = unpack(x); build its self-energies,
! solve Dyson/Keldysh for G_out at every energy, return F = pack(G_out) - x.
! abort = .true. if G_out ran away (NaN, occupations far outside [0,1]).
!-----------------------------------------------------------------------------
   subroutine scf_residual(n, x, F, abort)
      integer, intent(in) :: n
      real(dp), intent(in) :: x(n)
      real(dp), intent(out) :: F(n)
      logical, intent(out) :: abort
      integer :: iw, i_worst
      logical :: unphysical
      real(dp) :: n_worst, w_worst

      scf_iter = scf_iter + 1
      write (*, *) '.... ITERATION = ', scf_iter, ' ....'

      !.....input GFs from x
      call unpack_gf(G_in%R, G_in%L, x, n, n_full, n_strict_upper, Norb, n_omega)
      call complete_input_gf()

      !.....self-energies and output GFs
      call equal_time_lesser()
      if (order .eq. 2) call sigma2_compute(.true.)

      !$OMP PARALLEL DO PRIVATE(iw)
      do iw = 1, n_omega
         call interacting_gf(iw, V_scf)
      end do
      !$OMP END PARALLEL DO

      call check_gf_sanity(G_out%R, G_out%L, Norb, n_omega, delta, abort, unphysical, i_worst, n_worst, w_worst)
      if (abort) then
         write (*, '(a,i5,a,i4,a,es11.3,a,es11.3)') ' SCF RUNAWAY at iteration', scf_iter, &
            ': state', i_worst, ' n =', n_worst, ', spectral weight =', w_worst
         F = 0.d0
         return
      end if

      call pack_gf(G_out%R, G_out%L, F, n, n_full, n_strict_upper, Norb, n_omega)
      F = F - x
   end subroutine scf_residual

!-----------------------------------------------------------------------------
! interacting_gf: interacting GFs at one energy omega(iw), from the self-energies of
! the current input G_in. Dyson and Keldysh equations:
!   G^r = [ hbar*w - H + (i/2hbar)(GammaL+GammaR) - Sigma^r ]^-1,  G^a = (G^r)^dagger
!   G^< = G^r [ (i/hbar)(f_L GammaL + f_R GammaR) + Sigma^< ] G^a
!   G^> = G^< + G^r - G^a
! Called inside an OpenMP loop over iw: uses only local work arrays.
!-----------------------------------------------------------------------------
   subroutine interacting_gf(iw, Volt)
      implicit none
      integer, intent(in) :: iw
      real(dp), intent(in) :: Volt
      integer :: i
      integer :: info_loc !...thread-local (module INFO is shared across OpenMP threads)
      real(dp) :: w
      complex(dp), allocatable, dimension(:, :) :: SigmaL, SigmaR, w1, w2

      allocate (SigmaL(Norb, Norb), SigmaR(Norb, Norb))
      allocate (w1(Norb, Norb), w2(Norb, Norb))

      !............interaction self-energies
      call sigma_retarded(iw, SigmaR)
      call sigma_lesser(iw, SigmaL)

      !............retarded and advanced
      w = omega(iw)
      w1 = -H + 0.5d0*(im/hbar)*(GammaL + GammaR) - SigmaR ! LK <== must be + for emb and minus for interaction sigma
      do i = 1, Norb
         w1(i, i) = w1(i, i) + hbar*w
      end do

      call invert_complex(Norb, w1, info_loc)
      call hermitian_conjugate(w1, Norb, w2)

      G_out%R(:, :, iw) = w1; G_out%A(:, :, iw) = w2

      !............lesser: lead embedding + interaction Sigma^<
      SigmaL = im*(fermi(w, Volt)*GammaL + fermi(w, 0.d0)*GammaR)/hbar + SigmaL

      G_out%L(:, :, iw) = matmul(matmul(w1, SigmaL), w2) !.. GL = Gr * SigmaL * Ga
      G_out%G(:, :, iw) = G_out%L(:, :, iw) + G_out%R(:, :, iw) - G_out%A(:, :, iw)

      deallocate (SigmaL, SigmaR, w1, w2)
   end subroutine interacting_gf

end module SCF
