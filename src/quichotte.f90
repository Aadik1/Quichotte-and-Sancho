!==============================================================================
! src/quichotte.f90  --  main program of QUICHOTTE
!   QUantum Interacting CHains: Out-of-equilibrium TransporT (spin-resolved)
!
! Spin-resolved I-V curve of a helical molecular junction with spin-orbit
! coupling and Hubbard interaction (0th, 1st or 2nd order in U), from
! self-consistent non-equilibrium Green's functions.
!
! Flow:
!   1. read input.dat, build H (model chosen by the key 'hamiltonian':
!      helix, xyz or custom -- see src/hamiltonian.f90) and the lead couplings
!   2. diagonalise H to choose the energy grid [Ev_min - dw, Ev_max + up]
!   3. for every voltage V = V, V+delv, ..., Vf:
!        solve_scf     -> converged G^r, G^<   (scf.f90)
!        Current       -> J_up, J_down         (observables.f90)
!        spin_texture  -> <S>, populations     (observables.f90)
!
! Output files (in the run directory):
!   Volt_Current_<order>.dat   V, J_up, J_down, J_total, converged (1/0)
!   populations.dat            V, site, Sx, Sy, Sz, n_down, n_up, n_total
!   spectral_function.dat      V, omega, A_ii(omega) for every spin-orbital
!   scf_V_<V>.dat              SCF convergence history at each voltage
!   info_Hamiltonian.dat       H, eigenvalues, eigenvectors
!   runtime_datasheet.dat      code version, voltages, wall time
!
! Exit codes: 0 = all voltages converged, 1 = bad input,
!             3 = finished, but some voltages did not converge.
!==============================================================================
program QUICHOTTE
   use Constants, only: dp
   use LinearAlgebra, only: hermitian_eigensystem
   use Parameters
   use Hamiltonian
   use GreenFunctions
   use SelfEnergies, only: SigR2, sigma2_finalize
   use SCF
   use Observables
   use InputOutput
   implicit none
   include 'version.inc'      !...code_version, written by the Makefile from `git describe`
   real(dp) :: V1, J_up, J_down, total_time
   integer :: k, i, start_tick, end_tick, rate, max_count, n_unconv
   character(len=30) :: vfn
   logical :: fresh, converged

   !...wall-clock timer
   open (u_runtime, file='runtime_datasheet.dat', status='replace')
   call system_clock(count_rate=rate, count_max=max_count)
   if (rate .eq. 0) then
      write (*, *) "Error: SYSTEM_CLOCK not supported or rate is 0."
      stop 1
   end if
   call system_clock(count=start_tick)
   write (u_runtime, *) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>', 'Runtime Data Sheet', '>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
   write (u_runtime, *) 'QUICHOTTE version: ', code_version
   write (*, *) 'QUICHOTTE version: ', code_version

   !=====================================================================
   ! 1. Input, Hamiltonian and lead couplings
   !=====================================================================
   call read_input()

   allocate (H(Norb, Norb))
   allocate (GammaL(Norb, Norb), GammaR(Norb, Norb))

   call build_hamiltonian()          !...model chosen by the input key 'hamiltonian'

   !...left lead couples to the first site (spin-dependent), right lead to the last site
   GammaL = (0.d0, 0.d0); GammaL(1, 1) = GammaL_up; GammaL(2, 2) = GammaL_down
   GammaR = (0.d0, 0.d0); GammaR(Norb - 1, Norb - 1) = GammaR_both; GammaR(Norb, Norb) = GammaR_both

   call setup_hubbard()

   !=====================================================================
   ! 2. Energy grid from the spectrum of H
   !=====================================================================
   allocate (Eigenvec(Norb, Norb), Ev(Norb))
   Eigenvec = H
   call hermitian_eigensystem(Eigenvec, Norb, Ev)

   w_init = Ev(1) - dw; w_fin = Ev(Norb) + up
   n_omega = int((w_fin - w_init)/delta)
   write (*, *) 'n_omega:', n_omega, 'w_fin:', w_fin, 'w_init:', w_init
   write (*, *) 'Temp:', T, 'delta:', delta, 'Hamiltonian Dimension:', Norb, '# of omegas:', n_omega

   allocate (omega(n_omega))
   do i = 1, n_omega
      omega(i) = w_init + i*delta
   end do

   call write_hamiltonian_info()
   deallocate (Eigenvec, Ev)

   !=====================================================================
   ! Storage: GFs on the whole grid (G_in = SCF input, G_out = SCF output)
   !=====================================================================
   allocate (G_in%r(Norb, Norb, n_omega), G_in%a(Norb, Norb, n_omega))
   allocate (G_in%L(Norb, Norb, n_omega), G_in%G(Norb, Norb, n_omega))
   allocate (G_out%R(Norb, Norb, n_omega), G_out%A(Norb, Norb, n_omega))
   allocate (G_out%L(Norb, Norb, n_omega), G_out%G(Norb, Norb, n_omega))
   allocate (G_less_t0(Norb, Norb))

   !...length of the flattened Broyden vector (layout: see pack_gf in linear_algebra.f90)
   n_full = Norb*Norb
   n_strict_upper = (n_full - Norb)/2
   n_upper = n_strict_upper + Norb
   n_vec = (n_full + n_full + n_strict_upper + n_upper)*n_omega
   write (*, *) 'n_vec', n_vec

   allocate (work1(Norb, Norb), work2(Norb, Norb), work3(Norb, Norb))
   allocate (S_a(Norb, 3))

   !=====================================================================
   ! 3. Voltage sweep
   !=====================================================================
   write (vfn, '(i0)') order
   open (u_iv, file='Volt_Current_'//trim(vfn)//'.dat', status='replace')
   open (u_pop, file='populations.dat', status='replace')
   open (u_spec, file='spectral_function.dat', status='replace')
   !...no header lines: the HPC sweep scripts count data lines in these files

   !...the first voltage starts from the non-interacting GFs; later ones from the
   !...previous voltage if continuation = T and that voltage converged
   fresh = .true.
   n_unconv = 0
   do k = 0, n_volt_steps
      V1 = V + k*delv
      !...H is built once above. If you sweep a Hamiltonian parameter with k
      !...(e.g. t_hop) instead of V, set it here and call build_hamiltonian() again.

      call solve_scf(V1, fresh, converged)
      if (.not. converged) n_unconv = n_unconv + 1
      fresh = (.not. continuation) .or. (.not. converged)

      call Current(V1, J_up, J_down)
      write (u_iv, *) V1, J_up, J_down, J_up + J_down, merge(1, 0, converged)
      flush (u_iv)

      call spin_texture(V1, S_a, Norb, u_pop)
      call write_spectral(V1, u_spec)

      print *, 'Progress:', 100.d0*k/max(n_volt_steps, 1), '%', J_up, J_down
   end do

   close (u_iv); close (u_pop); close (u_spec)

   !...wall time
   call system_clock(count=end_tick)
   total_time = real(end_tick - start_tick)/real(rate)/60.d0
   write (u_runtime, *) '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
   write (u_runtime, *) 'Total Runtime:', total_time, 'mins'
   write (u_runtime, *) '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
   write (*, *) 'Total Runtime:', total_time, 'mins'
   close (u_runtime)

   if (allocated(SigR2)) call sigma2_finalize()
   deallocate (G_out%L, G_out%G, G_out%R, G_out%A)
   deallocate (G_in%r, G_in%a, G_in%L, G_in%G)
   deallocate (work1, work2, work3)
   deallocate (H, U_site, omega, S_a)
   deallocate (GammaL, GammaR, G_less_t0)

   !...non-zero exit code so job scripts can tell a run with unconverged voltages
   if (n_unconv .gt. 0) then
      write (*, *) 'WARNING:', n_unconv, 'voltage(s) did not converge (column 5 = 0 in Volt_Current)'
      error stop 3
   end if
end program QUICHOTTE
