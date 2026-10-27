!==============================================================================
! src/parameters.f90  --  module Parameters   (QUICHOTTE)
!
! Run parameters (read from input.dat by input_output.f90), sizes derived
! from them, the Pauli matrices and the unit numbers of the output files.
!
! Basis: two spin-orbitals per site,
!     index 2*i-1 = site i, spin up
!     index 2*i   = site i, spin down
! so Norb = 2 * (number of sites) = 2 * N_ions * N_turns.
!==============================================================================
module Parameters
   use Constants, only: dp
   implicit none
   public

   !...sizes
   integer :: Norb           !...number of spin-orbitals = 2 * number of sites
   integer :: n_omega        !...number of points on the energy grid
   integer :: n_volt_steps   !...number of voltage steps (n_volt_steps + 1 voltages)
   integer :: N_ions         !...sites per turn of the helix
   integer :: N_turns        !...turns of the helix
   integer :: order          !...perturbation order in U: 0, 1 or 2

   !...voltage sweep, temperature, chemical potential
   real(dp) :: V, Vf, delv   !...first voltage, last voltage, voltage step [V]
   real(dp) :: T, mu, beta   !...temperature [K], chemical potential [eV], 1/(kB T)

   !...energy grid omega(i) = w_init + i*delta, i = 1..n_omega
   real(dp) :: delta, w_init, w_fin
   real(dp) :: dw, up        !...grid margins below/above the spectrum of H [eV]

   !...model
   real(dp) :: E_CC, t_hop   !...on-site energy, nearest-neighbour hopping [eV]
   real(dp) :: lamb          !...spin-orbit coupling strength [eV]
   real(dp) :: U_int         !...uniform Hubbard U [eV]
   real(dp) :: hel_radius    !...helix radius
   real(dp) :: hel_length    !...helix length scale (see helix_coordinate in hamiltonian.f90)
   real(dp) :: hand          !...helicity: +1 right-handed, -1 left-handed
   real(dp) :: GammaL_up, GammaL_down   !...left lead coupling to site 1, per spin [eV]
   real(dp) :: GammaR_both              !...right lead coupling to the last site, both spins [eV]

   !...Pauli matrices
   complex(dp), parameter :: pauli_x(2, 2) = reshape([(0.d0, 0.d0), (1.d0, 0.d0), &
                                                      (1.d0, 0.d0), (0.d0, 0.d0)], [2, 2])
   complex(dp), parameter :: pauli_y(2, 2) = reshape([(0.d0, 0.d0), (0.d0, 1.d0), &
                                                      (0.d0, -1.d0), (0.d0, 0.d0)], [2, 2])
   complex(dp), parameter :: pauli_z(2, 2) = reshape([(1.d0, 0.d0), (0.d0, 0.d0), &
                                                      (0.d0, 0.d0), (-1.d0, 0.d0)], [2, 2])

   !...Fortran unit numbers of the output files
   integer, parameter :: u_runtime = 212   !...runtime_datasheet.dat
   integer, parameter :: u_scflog = 18     !...scf_V_<V>.dat (one per voltage)
   integer, parameter :: u_iv = 30         !...Volt_Current_<order>.dat
   integer, parameter :: u_pop = 230       !...populations.dat
   integer, parameter :: u_spec = 71       !...spectral_function.dat
end module Parameters
