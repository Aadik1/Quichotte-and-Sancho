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
   integer :: N_sites        !...number of sites
   integer :: Norb           !...number of spin-orbitals = 2 * N_sites
   integer :: n_omega        !...number of points on the energy grid
   integer :: n_volt_steps   !...number of voltage steps (n_volt_steps + 1 voltages)
   integer :: N_ions         !...sites per turn of the helix   (hamiltonian = helix)
   integer :: N_turns        !...turns of the helix            (hamiltonian = helix)
   integer :: order          !...perturbation order in U: 0, 1 or 2

   !...which Hamiltonian to build (see build_hamiltonian in hamiltonian.f90)
   character(len=32)  :: hamiltonian_model = 'helix'   !...'helix', 'xyz', 'dftb', 'dftb_hs' or 'custom'
   character(len=256) :: xyz_file = 'molecule.xyz'     !...positions (xyz); elements (dftb_hs)
   character(len=256) :: ham_file = 'hamiltonian_blocks_fortran.txt'   !...dftb: blocks; dftb_hs: hamsqr1.dat
   character(len=256) :: overlap_file = 'oversqr.dat'  !...(hamiltonian = dftb_hs)
   logical :: write_hamiltonian = .false.              !...write H in full precision to hamiltonian_matrix.dat

   !...dftb_hs: elements and their shells / spin-orbit constants (from the input keys
   !   max_angular_momentum = 'C:p H:s' and soc_xi = 'C:0.00849', xi in eV)
   character(len=16) :: soc_form = 'dual'              !...'dual' (as DFTB+) or 'onsite'
   character(len=2), allocatable :: atom_elem(:)       !...element of atom a (N_atoms)
   character(len=2), allocatable :: elem_name(:)       !...elements given in max_angular_momentum
   integer, allocatable :: elem_lmax(:)                !...0 = s, 1 = s + p
   real(dp), allocatable :: elem_xi(:)                 !...p-shell spin-orbit constant [eV]

   !...atoms and orbitals. A 'site' is one spatial orbital = one spin pair (2i-1, 2i).
   !   For helix, xyz and custom every atom has one orbital (atom = site). For dftb an
   !   atom has several (C: s, p_y, p_z, p_x in DFTB+ order; H: s).
   integer :: N_atoms                               !...number of atoms
   integer, allocatable :: atom_norb(:)             !...spatial orbitals of atom a   (N_atoms)
   integer, allocatable :: atom_first(:)            !...first site (orbital) of atom a (N_atoms)
   integer, allocatable :: orb_atom(:)              !...atom of site i               (N_sites)

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
   integer, parameter :: u_popat = 231     !...populations_atoms.dat
end module Parameters
