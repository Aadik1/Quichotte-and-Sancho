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
   use Hamiltonian, only: H, xyz_site_count, dftb_blocks_scan, read_xyz_elements
   use GreenFunctions, only: U_site, Ev, Eigenvec, GammaL, GammaR
   use SCF, only: alpha_mix, nhist, maxscf, tolf, tolf_order1, continuation, broyden_single_precision
   implicit none
   private
   public :: read_input, setup_hubbard, setup_leads, write_hamiltonian_info

   character(len=256) :: hubbard_file = ''   !...optional: per-site U, one value per line
   character(len=256) :: lead_left_atoms = ''    !...atoms the left lead couples to  (default: 1)
   character(len=256) :: lead_right_atoms = ''   !...atoms the right lead couples to (default: N_atoms)
   character(len=256) :: max_l_string = ''       !...dftb_hs: shells per element, e.g. 'C:p H:s'
   character(len=256) :: soc_xi_string = ''      !...dftb_hs: p-shell xi per element [eV], e.g. 'C:0.00849'

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
      if (trim(hamiltonian_model) .eq. 'dftb_hs') then
         call input_get('ham_file', ham_file, default='hamsqr1.dat')
         call input_get('overlap_file', overlap_file, default='oversqr.dat')
         call input_get('max_angular_momentum', max_l_string)
         call input_get('soc_xi', soc_xi_string, default='')
         call input_get('soc_form', soc_form, default='dual')
      else
         call input_get('ham_file', ham_file, default='hamiltonian_blocks_fortran.txt')
         call input_get('overlap_file', overlap_file, default='oversqr.dat')
         call input_get('max_angular_momentum', max_l_string, default='')
         call input_get('soc_xi', soc_xi_string, default='')
         call input_get('soc_form', soc_form, default='dual')
      end if
      call input_get('write_hamiltonian', write_hamiltonian, default=.false.)
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
      !...peptide model (alpha-helix defaults; hel_radius and hand are shared with helix)
      call input_get('pep_twist', pep_twist, default=100.d0)
      call input_get('pep_rise', pep_rise, default=1.5d0)
      call input_get('hop_range', hop_range, default=1)
      call input_get('hop_decay', hop_decay, default=1.d0)
      call input_get('soc_model', soc_model, default='radial')
      !...Hamiltonian and interaction
      if (trim(hamiltonian_model) .eq. 'dftb' .or. trim(hamiltonian_model) .eq. 'dftb_hs') then
         call input_get('E_CC', E_CC, default=0.d0)     !...not used: H comes from the DFTB files
         call input_get('t_hop', t_hop, default=0.d0)
      else
         call input_get('E_CC', E_CC)
         call input_get('t_hop', t_hop)
      end if
      if (trim(hamiltonian_model) .eq. 'helix' .or. trim(hamiltonian_model) .eq. 'xyz' .or. &
          trim(hamiltonian_model) .eq. 'peptide') then
         call input_get('lamb', lamb)                   !...SOC: helix, xyz and peptide models
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
      call input_get('tolf_order1', tolf_order1, default=1.d-8)
      call input_get('broyden_single_precision', broyden_single_precision, default=.true.)
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
         ' tolf =', tolf, ' tolf_order1 =', tolf_order1, ' continuation = ', continuation, &
         ' single-precision history = ', broyden_single_precision
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
      case ('dftb_hs')
         !...atoms from xyz_file, orbitals per atom from max_angular_momentum
         call setup_dftb_elements()
         n_model = sum(atom_norb)
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

   subroutine setup_dftb_elements()
      !...hamiltonian = dftb_hs: elements of xyz_file, shells from max_angular_momentum
      !   ('C:p H:s'; s -> 1 orbital, p -> 4: s, p_y, p_z, p_x), p-shell xi from soc_xi
      !   ('C:0.00849', eV; elements not listed have xi = 0)
      character(len=8), allocatable :: names(:)
      character(len=32), allocatable :: vals(:)
      integer :: n, a, e, k, ios
      real(dp) :: xv

      call read_xyz_elements(xyz_file, N_atoms, atom_elem)
      do a = 1, N_atoms
         atom_elem(a) = element_case(atom_elem(a))
      end do

      call split_pairs(max_l_string, 'max_angular_momentum', names, vals, n)
      allocate (elem_name(n), elem_lmax(n), elem_xi(n)); elem_xi = 0.d0
      do e = 1, n
         elem_name(e) = element_case(names(e))
         select case (trim(adjustl(vals(e))))
         case ('s', 'S'); elem_lmax(e) = 0
         case ('p', 'P'); elem_lmax(e) = 1
         case default
            write (error_unit, *) 'ERROR: max_angular_momentum: only s and p shells are supported, not "', &
               trim(vals(e)), '" for ', trim(names(e))
            stop 1
         end select
      end do
      deallocate (names, vals)

      if (len_trim(soc_xi_string) .gt. 0) then
         call split_pairs(soc_xi_string, 'soc_xi', names, vals, n)
         do k = 1, n
            e = findloc(elem_name, element_case(names(k)), dim=1)
            read (vals(k), *, iostat=ios) xv
            if (e .eq. 0 .or. ios .ne. 0) then
               write (error_unit, *) 'ERROR: soc_xi: bad entry ', trim(names(k)), ':', trim(vals(k)), &
                  ' (element must also appear in max_angular_momentum)'
               stop 1
            end if
            elem_xi(e) = xv
         end do
      end if

      allocate (atom_norb(N_atoms))
      do a = 1, N_atoms
         e = findloc(elem_name, atom_elem(a), dim=1)
         if (e .eq. 0) then
            write (error_unit, *) 'ERROR: element ', trim(atom_elem(a)), ' of ', trim(xyz_file), &
               ' is missing from max_angular_momentum'
            stop 1
         end if
         atom_norb(a) = merge(4, 1, elem_lmax(e) .eq. 1)
      end do
   end subroutine setup_dftb_elements

   subroutine split_pairs(str, key, names, vals, n)
      !...'C:p H:s' or 'C:0.0085, N:0.012' -> names(:), vals(:)
      character(len=*), intent(in) :: str, key
      character(len=8), allocatable, intent(out) :: names(:)
      character(len=32), allocatable, intent(out) :: vals(:)
      integer, intent(out) :: n
      character(len=len(str)) :: s
      character(len=64) :: tok(100)
      integer :: p, q, k, c

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
         if (n .lt. size(tok)) then
            n = n + 1; tok(n) = s(p:q)
         end if
         p = q + 1
      end do
      if (n .eq. 0) then
         write (error_unit, *) 'ERROR: ', key, ' is empty'
         stop 1
      end if
      allocate (names(n), vals(n))
      do k = 1, n
         c = index(tok(k), ':')
         if (c .le. 1) then
            write (error_unit, *) 'ERROR: ', key, ': entries must look like Element:value, got "', trim(tok(k)), '"'
            stop 1
         end if
         names(k) = tok(k)(1:c - 1); vals(k) = tok(k)(c + 1:)
      end do
   end subroutine split_pairs

   pure function element_case(e) result(r)
      !...'c' / 'C' -> 'C',  'CL' / 'cl' -> 'Cl'
      character(len=*), intent(in) :: e
      character(len=2) :: r
      integer :: k
      r = adjustl(e)
      k = iachar(r(1:1)); if (k .ge. 97 .and. k .le. 122) r(1:1) = achar(k - 32)
      k = iachar(r(2:2)); if (k .ge. 65 .and. k .le. 90) r(2:2) = achar(k + 32)
   end function element_case

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
      if (trim(hamiltonian_model) .eq. 'peptide') then
         if (N_sites .lt. 2) call complain('peptide: need N_sites >= 2')
         if (hel_radius .le. 0.d0) call complain('peptide: hel_radius must be > 0 (alpha-helix: 2.5)')
         if (abs(abs(hand) - 1.d0) .gt. 1.d-12) call complain('peptide: hand must be +1 or -1')
         if (hop_range .lt. 1) call complain('peptide: hop_range must be >= 1')
         if (hop_decay .le. 0.d0) call complain('peptide: hop_decay must be > 0')
         if (trim(soc_model) .ne. 'radial' .and. trim(soc_model) .ne. 'torsion') &
            call complain('peptide: soc_model must be radial or torsion')
      end if
      if (order .lt. 0 .or. order .gt. 2) call complain('order must be 0, 1 or 2')
      if (delta .le. 0.d0) call complain('delta (energy grid step) must be > 0')
      if (delv .le. 0.d0) call complain('delv (voltage step) must be > 0')
      if (Vf .lt. V) call complain('Vf must be >= V')
      if (T .le. 0.d0) call complain('T must be > 0')
      if (nhist .lt. 0) call complain('nhist must be >= 0')
      if (maxscf .lt. 1) call complain('maxscf must be >= 1')
      if (tolf .le. 0.d0) call complain('tolf must be > 0')
      if (tolf_order1 .le. 0.d0) call complain('tolf_order1 must be > 0')
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
