!==============================================================================
! src/input_output.f90  --  module InputOutput   (QUICHOTTE)
!
!   read_input               read input.dat (common/input_parser.f90), check values
!   setup_hubbard            per-site Hubbard U (uniform, or from a file)
!   write_hamiltonian_info   H, its eigenvalues/vectors -> info_Hamiltonian.dat
!
! Every input key is listed in read_input, with its default if it has one.
!==============================================================================
module InputOutput
   use, intrinsic :: iso_fortran_env, only: error_unit
   use Constants, only: dp, kb
   use InputParser
   use Parameters
   use Hamiltonian, only: H, xyz_site_count
   use GreenFunctions, only: U_site, Ev, Eigenvec
   use SCF, only: alpha_mix, nhist, maxscf, tolf, continuation
   implicit none
   private
   public :: read_input, setup_hubbard, write_hamiltonian_info

   character(len=256) :: hubbard_file = ''   !...optional: per-site U, one value per line

contains

   subroutine read_input()
      logical :: restart_unused

      call input_load([character(len=12) :: 'input.dat', 'inputSOC.dat'])

      !...voltage sweep, temperature, chemical potential
      call input_get('V', V)
      call input_get('Vf', Vf)
      call input_get('delv', delv)
      call input_get('T', T)
      call input_get('mu', mu)
      !...which Hamiltonian (see build_hamiltonian in hamiltonian.f90)
      call input_get('hamiltonian', hamiltonian_model, default='helix')
      call input_get('xyz_file', xyz_file, default='molecule.xyz')
      call input_get('N_sites', N_sites, default=0)
      !...helix geometry: required for hamiltonian = helix, optional otherwise
      if (trim(hamiltonian_model) .eq. 'helix') then
         call input_get('N_ions', N_ions)
         call input_get('N_turns', N_turns)
         call input_get('hel_length', hel_length)
         call input_get('hel_radius', hel_radius)
         call input_get('hand', hand)
      else
         call input_get('N_ions', N_ions, default=0)
         call input_get('N_turns', N_turns, default=0)
         call input_get('hel_length', hel_length, default=0.d0)
         call input_get('hel_radius', hel_radius, default=0.d0)
         call input_get('hand', hand, default=1.d0)
      end if
      !...Hamiltonian and interaction
      call input_get('E_CC', E_CC)
      call input_get('t_hop', t_hop)
      if (trim(hamiltonian_model) .eq. 'helix' .or. trim(hamiltonian_model) .eq. 'xyz') then
         call input_get('lamb', lamb)                   !...SOC: helix and xyz models
      else
         call input_get('lamb', lamb, default=0.d0)     !...available to your model if you want it
      end if
      call input_get('U_int', U_int, alias='Hubbard')
      call input_get('hubbard_file', hubbard_file, default='')
      call input_get('order', order)
      !...lead couplings
      call input_get('GammaL_up', GammaL_up)
      call input_get('GammaL_down', GammaL_down)
      call input_get('GammaR', GammaR_both, alias='GammaR_both')
      !...energy grid
      call input_get('dw', dw)
      call input_get('up', up)
      call input_get('delta', delta)
      !...SCF
      call input_get('alpha_mix', alpha_mix, default=0.5d0, alias='pulay')
      call input_get('nhist', nhist, default=16)
      call input_get('maxscf', maxscf, default=300)
      call input_get('tolf', tolf, default=1.d-6)
      call input_get('continuation', continuation, default=.true.)
      call input_get('restart', restart_unused, default=.false.)   !...old key, not used
      call input_finish()

      !...number of sites: from the helix, from the xyz file, or the key N_sites
      call set_number_of_sites()

      !...derived quantities
      beta = 1.d0/(kb*T)
      Norb = 2*N_sites          !...sites x 2 spins
      n_volt_steps = nint((Vf - V)/delv)

      call check_input()

      write (*, *) 'Input read from ', trim(input_file_name)
      write (*, *) 'Hamiltonian: ', trim(hamiltonian_model), '   N_sites:', N_sites
      write (*, *) 'T:', T, 'V:', V, 'Vf:', Vf, 'delv:', delv, 'mu:', mu
      write (*, *) 'order:', order, 'U_int:', U_int, 'Norb:', Norb
      write (*, *) 'dw:', dw, 'up:', up, 'delta:', delta
      write (*, *) 'SCF: alpha_mix =', alpha_mix, ' nhist =', nhist, ' maxscf =', maxscf, &
         ' tolf =', tolf, ' continuation = ', continuation
   end subroutine read_input

   subroutine set_number_of_sites()
      !...N_sites for each model:
      !     helix   N_ions*N_turns   (N_sites optional; must agree if given)
      !     xyz     atoms in xyz_file (N_sites optional; must agree if given)
      !     custom  the key N_sites  (required; also for any model you add)
      integer :: n_model

      select case (trim(hamiltonian_model))
      case ('helix')
         n_model = N_ions*N_turns
      case ('xyz')
         n_model = xyz_site_count(xyz_file)
      case default
         !...any other model (custom, or one you add to build_hamiltonian):
         !...its size comes from the key N_sites. Unknown names are rejected
         !...by build_hamiltonian in hamiltonian.f90.
         if (N_sites .lt. 1) then
            write (error_unit, *) 'ERROR: hamiltonian = ', trim(hamiltonian_model), ' needs the key N_sites'
            write (error_unit, *) '       (built-in models: helix, xyz, custom)'
            stop 1
         end if
         return
      end select

      if (N_sites .gt. 0 .and. N_sites .ne. n_model) then
         write (error_unit, '(a,i0,a,a,a,i0)') ' ERROR: N_sites = ', N_sites, ' but the ', &
            trim(hamiltonian_model), ' model has ', n_model
         stop 1
      end if
      N_sites = n_model
   end subroutine set_number_of_sites

   subroutine check_input()
      !...stop early, with a clear message, on inputs that cannot work
      logical :: bad

      bad = .false.
      if (trim(hamiltonian_model) .eq. 'helix' .and. (N_ions .lt. 2 .or. N_turns .lt. 1)) &
         call complain('need N_ions >= 2 and N_turns >= 1')
      if (N_sites .lt. 1) call complain('need at least one site')
      if (order .lt. 0 .or. order .gt. 2) call complain('order must be 0, 1 or 2')
      if (delta .le. 0.d0) call complain('delta (energy grid step) must be > 0')
      if (delv .le. 0.d0) call complain('delv (voltage step) must be > 0')
      if (Vf .lt. V) call complain('Vf must be >= V')
      if (T .le. 0.d0) call complain('T must be > 0')
      if (nhist .lt. 0) call complain('nhist must be >= 0')
      if (maxscf .lt. 1) call complain('maxscf must be >= 1')
      if (tolf .le. 0.d0) call complain('tolf must be > 0')
      if (bad) stop 1

   contains
      subroutine complain(msg)
         character(len=*), intent(in) :: msg
         write (error_unit, *) 'ERROR in input: ', msg
         bad = .true.
      end subroutine complain
   end subroutine check_input

   subroutine setup_hubbard()
      !...U_site(i) = U on the spin-up index of every site (i = 1, 3, 5, ...).
      !   Uniform U = U_int unless 'hubbard_file' is given in input.dat,
      !   in which case U of site k is read from line k of that file.
      implicit none
      integer :: i, u, ios

      allocate (U_site(Norb))
      U_site = 0.d0

      if (len_trim(hubbard_file) .eq. 0) then
         do i = 1, Norb, 2
            U_site(i) = U_int
         end do
      else
         open (newunit=u, file=trim(hubbard_file), status='old', action='read', iostat=ios)
         if (ios .ne. 0) then
            write (error_unit, *) 'ERROR: cannot open hubbard_file ', trim(hubbard_file)
            stop 1
         end if
         do i = 1, Norb, 2
            read (u, *, iostat=ios) U_site(i)
            if (ios .ne. 0) then
               write (error_unit, *) 'ERROR: ', trim(hubbard_file), ' needs one U per site (', Norb/2, ' lines)'
               stop 1
            end if
         end do
         close (u)
         write (*, *) 'Site-dependent U read from ', trim(hubbard_file)
      end if
   end subroutine setup_hubbard

   subroutine write_hamiltonian_info()
      !...H, a Hermiticity check, eigenvalues and eigenvectors -> info_Hamiltonian.dat
      implicit none
      integer :: i, j, u
      complex(dp) :: diff

      open (newunit=u, file='info_Hamiltonian.dat', status='replace')

      write (u, '(/a)') '... Hamiltonian:'
      do i = 1, Norb
         write (u, '(10(a,2f10.5,a))') (' [', H(i, j), '] ', j=1, Norb)
      end do

      write (u, '(/a)') '... Checking the Hamiltonian is Hermitian:'
      do i = 1, Norb
         do j = i, Norb
            diff = H(i, j) - conjg(H(j, i))
            write (u, *) i, j, diff
         end do
      end do

      write (u, '(/a)') '... Eigenvalues:'
      do i = 1, Norb
         write (u, *) i, '.', Ev(i)
      end do

      write (u, '(/a)') '... Eigenvectors:'
      do j = 1, Norb
         write (u, '(i3,10(a,2f10.5,a))') j, (' [', Eigenvec(i, j), '] ', i=1, Norb)
      end do

      write (u, *) 'Norb:', Norb
      close (u)
   end subroutine write_hamiltonian_info

end module InputOutput
