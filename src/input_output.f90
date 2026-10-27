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
   use Hamiltonian, only: H
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
      !...helix geometry
      call input_get('N_ions', N_ions)
      call input_get('N_turns', N_turns)
      call input_get('hel_length', hel_length)
      call input_get('hel_radius', hel_radius)
      call input_get('hand', hand)
      !...Hamiltonian and interaction
      call input_get('E_CC', E_CC)
      call input_get('t_hop', t_hop)
      call input_get('lamb', lamb)
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

      !...derived quantities
      beta = 1.d0/(kb*T)
      Norb = 2*N_ions*N_turns   !...sites x 2 spins
      n_volt_steps = nint((Vf - V)/delv)

      call check_input()

      write (*, *) 'Input read from ', trim(input_file_name)
      write (*, *) 'T:', T, 'V:', V, 'Vf:', Vf, 'delv:', delv, 'mu:', mu
      write (*, *) 'order:', order, 'U_int:', U_int, 'Norb:', Norb
      write (*, *) 'dw:', dw, 'up:', up, 'delta:', delta
      write (*, *) 'SCF: alpha_mix =', alpha_mix, ' nhist =', nhist, ' maxscf =', maxscf, &
         ' tolf =', tolf, ' continuation = ', continuation
   end subroutine read_input

   subroutine check_input()
      !...stop early, with a clear message, on inputs that cannot work
      logical :: bad

      bad = .false.
      if (N_ions .lt. 2 .or. N_turns .lt. 1) call complain('need N_ions >= 2 and N_turns >= 1')
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
