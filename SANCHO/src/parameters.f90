!==============================================================================
! src/parameters.f90  --  module Parameters   (SANCHO)
!
! Run parameters (read from input.dat by input_output.f90), sizes derived
! from them and the unit numbers of the output files.
! One orbital per site; spin degeneracy enters only as a factor 2 in the current.
!==============================================================================
module Parameters
   use Constants, only: dp
   implicit none
   public

   !...sizes
   integer :: Nsites         !...number of sites in the chain
   integer :: n_omega        !...number of points on the energy grid
   integer :: n_volt_steps   !...number of voltage steps (n_volt_steps + 1 voltages)
   integer :: order          !...perturbation order in U: 0, 1, 2 or 3

   !...voltage sweep, temperature, chemical potential
   real(dp) :: V, Vf, delv   !...first voltage, last voltage, voltage step [V]
   real(dp) :: T, mu, beta   !...temperature [K], chemical potential [eV], 1/(kB T)

   !...energy grid omega(i) = w_init + i*delta, i = 1..n_omega
   real(dp) :: delta, w_init, w_fin
   real(dp) :: dw, up        !...grid margins below/above the spectrum of H [eV]

   !...model
   real(dp) :: E_CC, t_hop   !...on-site energy, nearest-neighbour hopping [eV]
   real(dp) :: U_int         !...uniform Hubbard U [eV]
   real(dp) :: GammaL_coupling   !...left lead coupling to site 1 [eV]
   real(dp) :: GammaR_coupling   !...right lead coupling to the last site [eV]

   !...Fortran unit numbers of the output files
   integer, parameter :: u_runtime = 22    !...runtime_datasheet.dat
   integer, parameter :: u_scflog = 17     !...scf_V_<V>.dat (one per voltage)
   integer, parameter :: u_iv = 30         !...Volt_Current_<order>.dat
   integer, parameter :: u_pop = 230       !...populations.dat
   integer, parameter :: u_spec = 71       !...spectral_function.dat
end module Parameters
