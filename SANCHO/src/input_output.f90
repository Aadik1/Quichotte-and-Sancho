!==============================================================================
! src/input_output.f90  --  module InputOutput   (SANCHO)
!
!   read_input               read input.dat (common/input_parser.f90), check values
!   setup_hubbard            per-site Hubbard U (uniform, or from a file)
!   xyz_input                read atomic positions (EXPERIMENTAL, not used yet)
!   write_hamiltonian_info   H, its eigenvalues/vectors -> info_Hamiltonian.dat
!   print_3matrix            debug dump of an (N,N,n_omega) array
!
! Every input key is listed in read_input, with its default if it has one.
!==============================================================================
module InputOutput
   use, intrinsic :: iso_fortran_env, only: error_unit
   use Constants, only: dp, kb
   use InputParser
   use Parameters
   use Hamiltonian, only: H
   use GreenFunctions, only: U_site, Ev, Eigenvec, coords, elements
   use SCF, only: alpha_mix, nhist, maxscf, tolf, continuation
   implicit none
   private
   public :: read_input, setup_hubbard, xyz_input, write_hamiltonian_info, print_3matrix

   character(len=256) :: hubbard_file = ''   !...optional: per-site U, one value per line

contains

   subroutine read_input()
      logical :: restart_unused
      character(len=8) :: ver_unused

      call input_load([character(len=9) :: 'input.dat'])

      !...voltage sweep, temperature, chemical potential
      call input_get('V', V)
      call input_get('Vf', Vf)
      call input_get('delv', delv)
      call input_get('T', T)
      call input_get('mu', mu)
      !...chain Hamiltonian and interaction
      call input_get('Nsites', Nsites, alias='Natoms')
      call input_get('E_CC', E_CC)
      call input_get('t_hop', t_hop)
      call input_get('U_int', U_int)
      call input_get('hubbard_file', hubbard_file, default='')
      call input_get('order', order)
      !...lead couplings
      call input_get('GammaL', GammaL_coupling, default=1.d0)
      call input_get('GammaR', GammaR_coupling, default=1.d0)
      !...energy grid
      call input_get('dw', dw)
      call input_get('up', up)
      call input_get('delta', delta)
      !...SCF
      call input_get('alpha_mix', alpha_mix, default=0.3d0)
      call input_get('nhist', nhist, default=8)
      call input_get('maxscf', maxscf, default=200)
      call input_get('tolf', tolf, default=1.d-6)
      call input_get('continuation', continuation, default=.true.)
      call input_get('restart', restart_unused, default=.false.)   !...old keys, not used
      call input_get('ver', ver_unused, default='')
      call input_finish()

      !...derived quantities
      beta = 1.d0/(kb*T)
      n_volt_steps = nint((Vf - V)/delv)

      call check_input()

      write (*, *) 'Input read from ', trim(input_file_name)
      write (*, *) 'T:', T, 'V:', V, 'Vf:', Vf, 'delv:', delv, 'mu:', mu
      write (*, *) 'order:', order, 'U_int:', U_int, 'Nsites:', Nsites
      write (*, *) 'GammaL:', GammaL_coupling, 'GammaR:', GammaR_coupling
      write (*, *) 'dw:', dw, 'up:', up, 'delta:', delta
      write (*, *) 'SCF: alpha_mix =', alpha_mix, ' nhist =', nhist, ' maxscf =', maxscf, &
         ' tolf =', tolf, ' continuation = ', continuation
   end subroutine read_input

   subroutine check_input()
      !...stop early, with a clear message, on inputs that cannot work
      logical :: bad

      bad = .false.
      if (Nsites .lt. 1) call complain('Nsites must be >= 1')
      if (order .lt. 0 .or. order .gt. 3) call complain('order must be 0, 1, 2 or 3')
      if (delta .le. 0.d0) call complain('delta (energy grid step) must be > 0')
      if (delv .le. 0.d0) call complain('delv (voltage step) must be > 0')
      if (Vf .lt. V) call complain('Vf must be >= V')
      if (T .le. 0.d0) call complain('T must be > 0')
      if (nhist .lt. 0) call complain('nhist must be >= 0')
      if (maxscf .lt. 1) call complain('maxscf must be >= 1')
      if (tolf .le. 0.d0) call complain('tolf must be > 0')
      if (GammaL_coupling .lt. 0.d0 .or. GammaR_coupling .lt. 0.d0) call complain('GammaL, GammaR must be >= 0')
      if (bad) stop 1

   contains
      subroutine complain(msg)
         character(len=*), intent(in) :: msg
         write (error_unit, *) 'ERROR in input: ', msg
         bad = .true.
      end subroutine complain
   end subroutine check_input

   subroutine setup_hubbard()
      !...U_site(i) = U_int on every site, unless 'hubbard_file' is given in
      !   input.dat, in which case U of site i is read from line i of that file.
      implicit none
      integer :: i, u, ios

      allocate (U_site(Nsites))

      if (len_trim(hubbard_file) .eq. 0) then
         U_site = U_int
      else
         open (newunit=u, file=trim(hubbard_file), status='old', action='read', iostat=ios)
         if (ios .ne. 0) then
            write (error_unit, *) 'ERROR: cannot open hubbard_file ', trim(hubbard_file)
            stop 1
         end if
         do i = 1, Nsites
            read (u, *, iostat=ios) U_site(i)
            if (ios .ne. 0) then
               write (error_unit, *) 'ERROR: ', trim(hubbard_file), ' needs one U per site (', Nsites, ' lines)'
               stop 1
            end if
         end do
         close (u)
         write (*, *) 'Site-dependent U read from ', trim(hubbard_file)
      end if
   end subroutine setup_hubbard

   subroutine xyz_input(xyz_file)
      !...EXPERIMENTAL: read Nsites element symbols and positions from an .xyz
      !   file into 'elements' and 'coords' (for build_hamiltonian_from_xyz).
      implicit none
      character(len=*), intent(in) :: xyz_file
      character(len=256) :: line
      integer :: i, u

      open (newunit=u, file=xyz_file, status='old', action='read')
      read (u, '(A)') line   !...atom count
      read (u, '(A)') line   !...comment
      do i = 1, Nsites
         read (u, *) elements(i), coords(1, i), coords(2, i), coords(3, i)
      end do
      close (u)
   end subroutine xyz_input

   subroutine write_hamiltonian_info()
      !...H, a Hermiticity check, eigenvalues and eigenvectors -> info_Hamiltonian.dat
      implicit none
      integer :: i, j, u
      complex(dp) :: diff

      open (newunit=u, file='info_Hamiltonian.dat', status='replace')

      write (u, '(/a)') '... Hamiltonian:'
      do i = 1, Nsites
         write (u, '(10(a,2f10.5,a))') (' [', H(i, j), '] ', j=1, Nsites)
      end do

      write (u, '(/a)') '... Checking the Hamiltonian is Hermitian:'
      do i = 1, Nsites
         do j = i, Nsites
            diff = H(i, j) - conjg(H(j, i))
            write (u, *) i, j, diff
         end do
      end do

      write (u, '(/a)') '... Eigenvalues:'
      do i = 1, Nsites
         write (u, *) i, '.', Ev(i)
      end do

      write (u, '(/a)') '... Eigenvectors:'
      do j = 1, Nsites
         write (u, '(i3,10(a,2f10.5,a))') j, (' [', Eigenvec(i, j), '] ', i=1, Nsites)
      end do

      close (u)
   end subroutine write_hamiltonian_info

   subroutine print_3matrix(iteration, X, name, Nat, Nw)
      !...debug: dump X(:,:,iw) (first 10x10 block) to <name>_<iteration>.dat
      implicit none
      integer, intent(in) :: iteration, Nat, Nw
      complex(dp), intent(in) :: X(Nat, Nat, Nw)
      character(len=3), intent(in) :: name
      integer :: n, iw, i, j, u
      character(len=12) :: fn

      n = min(Nat, 10)
      write (fn, '(i0)') iteration
      open (newunit=u, file=name//'_'//trim(fn)//'.dat', status='replace')
      do iw = 1, Nw
         write (u, '(/i6)') iw
         do i = 1, n
            write (u, *) i, (X(i, j, iw), j=1, n)
         end do
      end do
      close (u)
      write (*, *) '... written '//name//'_'//trim(fn)//'.dat'
   end subroutine print_3matrix

end module InputOutput
