!==============================================================================
! src/input_output.f90  --  module InputOutput   (QUICHOTTE)
!
!   read_input               read input.dat (common/input_parser.f90), check values
!   setup_hubbard            per-site Hubbard U (uniform, or from a file)
!   setup_leads              GammaL, GammaR on all orbitals of the contact atoms
!   write_hamiltonian_info   H, its eigenvalues/vectors -> info_Hamiltonian.dat
!
! Every input key is listed in read_input, with its default if it has one.
!==============================================================================
module InputOutput
   use, intrinsic :: iso_fortran_env, only: error_unit
   use Constants, only: dp, kb
   use InputParser
   use Parameters
   use Hamiltonian, only: H, xyz_site_count, dftb_blocks_scan
   use GreenFunctions, only: U_site, Ev, Eigenvec, GammaL, GammaR
   use SCF, only: alpha_mix, nhist, maxscf, tolf, continuation
   implicit none
   private
   public :: read_input, setup_hubbard, setup_leads, write_hamiltonian_info

   character(len=256) :: hubbard_file = ''   !...optional: per-site U, one value per line
   character(len=256) :: lead_left_atoms = ''    !...atoms the left lead couples to  (default: 1)
   character(len=256) :: lead_right_atoms = ''   !...atoms the right lead couples to (default: N_atoms)

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
      call input_get('ham_file', ham_file, default='hamiltonian_blocks_fortran.txt')
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
      if (trim(hamiltonian_model) .eq. 'dftb') then
         call input_get('E_CC', E_CC, default=0.d0)     !...not used: H comes from ham_file
         call input_get('t_hop', t_hop, default=0.d0)
      else
         call input_get('E_CC', E_CC)
         call input_get('t_hop', t_hop)
      end if
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
      call input_get('lead_left_atoms', lead_left_atoms, default='')
      call input_get('lead_right_atoms', lead_right_atoms, default='')
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
      Norb = 2*N_sites          !...sites (spatial orbitals) x 2 spins
      call set_atom_map()
      n_volt_steps = nint((Vf - V)/delv)

      call check_input()

      write (*, *) 'Input read from ', trim(input_file_name)
      write (*, *) 'Hamiltonian: ', trim(hamiltonian_model), '   N_sites (orbitals):', N_sites, &
         '   N_atoms:', N_atoms
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
      integer :: n_model, n_at
      integer, allocatable :: norb(:)

      select case (trim(hamiltonian_model))
      case ('helix')
         n_model = N_ions*N_turns
      case ('xyz')
         n_model = xyz_site_count(xyz_file)
      case ('dftb')
         !...atoms and orbitals per atom from the block headers of ham_file
         call dftb_blocks_scan(ham_file, n_at, norb)
         n_model = sum(norb)
         N_atoms = n_at
         allocate (atom_norb(n_at)); atom_norb = norb
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

   subroutine set_atom_map()
      !...atom <-> site (spatial orbital) map. One orbital per atom for helix, xyz
      !   and custom; for dftb, atom_norb was read from ham_file.
      integer :: a, i

      if (.not. allocated(atom_norb)) then
         N_atoms = N_sites
         allocate (atom_norb(N_atoms)); atom_norb = 1
      end if
      allocate (atom_first(N_atoms), orb_atom(N_sites))
      i = 0
      do a = 1, N_atoms
         atom_first(a) = i + 1
         orb_atom(i + 1:i + atom_norb(a)) = a
         i = i + atom_norb(a)
      end do
   end subroutine set_atom_map

   subroutine setup_leads()
      !...Wide-band lead couplings on EVERY orbital of the contact atoms:
      !     left:  GammaL_up on the spin-up, GammaL_down on the spin-down state of each orbital
      !     right: GammaR on both spins of each orbital
      !   Contact atoms from lead_left_atoms / lead_right_atoms (lists such as '1 2 3'
      !   or '1:3'); defaults are atom 1 and atom N_atoms. With one orbital per atom and
      !   the defaults this is exactly the original coupling (states 1, 2 and Norb-1, Norb).
      integer, allocatable :: left(:), right(:)
      integer :: k, a, i

      call parse_atom_list(lead_left_atoms, 1, 'lead_left_atoms', left)
      call parse_atom_list(lead_right_atoms, N_atoms, 'lead_right_atoms', right)
      do k = 1, size(left)
         if (any(right .eq. left(k))) then
            write (error_unit, *) 'ERROR: atom', left(k), ' is in both lead_left_atoms and lead_right_atoms'
            stop 1
         end if
      end do

      GammaL = (0.d0, 0.d0); GammaR = (0.d0, 0.d0)
      do k = 1, size(left)
         a = left(k)
         do i = atom_first(a), atom_first(a) + atom_norb(a) - 1
            GammaL(2*i - 1, 2*i - 1) = GammaL_up
            GammaL(2*i, 2*i) = GammaL_down
         end do
      end do
      do k = 1, size(right)
         a = right(k)
         do i = atom_first(a), atom_first(a) + atom_norb(a) - 1
            GammaR(2*i - 1, 2*i - 1) = GammaR_both
            GammaR(2*i, 2*i) = GammaR_both
         end do
      end do
      write (*, '(a,*(1x,i0))') ' Left lead atoms: ', left
      write (*, '(a,*(1x,i0))') ' Right lead atoms:', right
   end subroutine setup_leads

   subroutine parse_atom_list(str, default_atom, key, list)
      !...'1 2 5', '1,2,5' or '1:4 7' -> integer list; empty -> [default_atom]
      character(len=*), intent(in) :: str, key
      integer, intent(in) :: default_atom
      integer, allocatable, intent(out) :: list(:)
      character(len=len(str)) :: s
      character(len=64) :: tok
      integer :: p, q, k, a, b, ios, n
      integer :: buf(10000)

      s = str
      do k = 1, len(s)
         if (s(k:k) .eq. ',') s(k:k) = ' '
      end do
      n = 0; p = 1
      do
         do while (p .le. len_trim(s))
            if (s(p:p) .ne. ' ') exit
            p = p + 1
         end do
         if (p .gt. len_trim(s)) exit
         q = index(s(p:), ' ') + p - 2
         if (q .lt. p) q = len_trim(s)
         tok = s(p:q); p = q + 1
         k = index(tok, ':')
         if (k .gt. 0) then
            read (tok(:k - 1), *, iostat=ios) a
            if (ios .eq. 0) read (tok(k + 1:), *, iostat=ios) b
         else
            read (tok, *, iostat=ios) a
            b = a
         end if
         if (ios .ne. 0 .or. a .lt. 1 .or. b .gt. N_atoms .or. b .lt. a) then
            write (error_unit, *) 'ERROR: ', key, ': bad atom or range "', trim(tok), '" (atoms are 1..', N_atoms, ')'
            stop 1
         end if
         do k = a, b
            if (n .lt. size(buf)) then
               n = n + 1; buf(n) = k
            end if
         end do
      end do
      if (n .eq. 0) then
         n = 1; buf(1) = default_atom
      end if
      allocate (list(n)); list = buf(1:n)
   end subroutine parse_atom_list

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
               write (error_unit, *) 'ERROR: ', trim(hubbard_file), ' needs one U per site/orbital (', Norb/2, ' lines)'
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
