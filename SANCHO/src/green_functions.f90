!==============================================================================
! src/green_functions.f90  --  module GreenFunctions   (SANCHO)
!
! Storage for the Green's functions on the energy grid and the routines that
! do not involve the interaction:
!   fermi                      Fermi-Dirac distribution
!   noninteracting_gr_ga       G0^r, G0^a on the whole grid
!   noninteracting_gl_gg       G0^<, G0^> on the whole grid
!   complete_input_gf          G^a and G^> of G_in from its G^r and G^<
!   equal_time_lesser          G^<_ii(t=0) = (delta/2pi) sum_w G^<_ii(w)
!   print_input_gf / print_output_gf / print_spectral_debug   debug printing
!
! Naming: G_in  = INPUT Green's functions of the current SCF step
!         G_out = OUTPUT Green's functions built from the self-energy of G_in
!         components %r retarded, %a advanced, %L lesser, %G greater,
!         each (Nsites, Nsites, n_omega).
!==============================================================================
module GreenFunctions
   use Constants, only: dp, hbar, pi, im
   use LinearAlgebra, only: invert_complex, hermitian_conjugate
   use Parameters
   use Hamiltonian, only: H
   implicit none

   type :: GreenFunctionSet
      complex(dp), allocatable, dimension(:, :, :) :: r, a, L, G
   end type GreenFunctionSet
   type(GreenFunctionSet) :: G_in, G_out

   real(dp), allocatable, dimension(:) :: omega    !...energy grid (n_omega)
   real(dp), allocatable, dimension(:) :: U_site   !...Hubbard U per site
   real(dp), allocatable, dimension(:) :: Ev       !...eigenvalues of H
   complex(dp), allocatable, dimension(:, :) :: Eigenvec
   complex(dp), allocatable, dimension(:, :) :: GammaL, GammaR   !...lead coupling matrices
   complex(dp), allocatable, dimension(:) :: G_less_t0          !...G^<_ii(t=0)
   complex(dp), allocatable, dimension(:, :) :: work1, work2, work3
   integer :: INFO

   !...geometry for the EXPERIMENTAL xyz Hamiltonian (see hamiltonian.f90)
   real(dp), allocatable, dimension(:, :) :: coords
   character(len=2), allocatable, dimension(:) :: elements

contains

   !.....fermi-dirac distribution
   real(dp) function fermi(w, V)
      implicit none
      real(dp) ::  arg, w, V

      arg = (w - mu + V/hbar)*beta
      fermi = 1.d0/(exp(hbar*arg) + 1.d0)

   end function fermi

   subroutine noninteracting_gr_ga()
      !............non-interacting G^r and G^a on the whole grid:
      !............G0^r = [hbar*w - H + (i/2hbar)(GammaL+GammaR)]^-1,  G0^a = (G0^r)^dagger
      implicit none
      integer :: j, i
      real(dp) :: w

      do j = 1, n_omega
         work1 = -H + 0.5d0*(im/hbar)*(GammaL + GammaR) !LK <========= must be +
         w = omega(j)
         do i = 1, Nsites
            work1(i, i) = work1(i, i) + hbar*w
         end do

         call invert_complex(Nsites, work1, info)
         call hermitian_conjugate(work1, Nsites, work2)

         G_in%r(:, :, j) = work1; G_in%a(:, :, j) = work2
      end do

   end subroutine noninteracting_gr_ga

   subroutine noninteracting_gl_gg(Volt)
      !............non-interacting G^< and G^> on the whole grid:
      !............G0^< = G0^r [ (i/hbar)(f_L GammaL + f_R GammaR) ] G0^a,  G0^> = G0^< + G0^r - G0^a
      implicit none
      real(dp), intent(in) :: Volt
      real(dp) :: w
      integer :: j

      work1 = (0.d0, 0.d0); work2 = (0.d0, 0.d0); work3 = (0.d0, 0.d0)
      do j = 1, n_omega
         w = omega(j)
         work1 = G_in%r(:, :, j)
         work2 = G_in%a(:, :, j)

         work3 = matmul(matmul(work1, im*(fermi(w, Volt)*GammaL + fermi(w, 0.d0)*GammaR)/hbar), work2)
         G_in%L(:, :, j) = work3
         G_in%G(:, :, j) = G_in%L(:, :, j) + G_in%R(:, :, j) - G_in%A(:, :, j)
      end do

   end subroutine noninteracting_gl_gg

   subroutine equal_time_lesser()
      !.....equal-time lesser GF on the diagonal, G_less_t0(i) = (delta/2pi) sum_w G^<_ii(w)
      implicit none
      integer :: i, k1
      real(dp) :: pp
      complex(dp) :: s

      pp = delta/(2.d0*pi)
      do i = 1, Nsites
         s = (0.d0, 0.d0)
         do k1 = 1, n_omega
            s = s + G_in%L(i, i, k1)
         end do
         G_less_t0(i) = s*pp
      end do

   end subroutine equal_time_lesser

   subroutine print_output_gf()
      implicit none
      integer :: i, j

      write (*, '(/a)') '... GFr:'
      do i = 1, Nsites
         write (*, '(10(a,2f10.5,a))') (' [', G_out%R(i, j, 1), '] ', j=1, Nsites)
      end do

      write (*, '(/a)') '... GFl:'
      do i = 1, Nsites
         write (*, '(10(a,2f10.5,a))') (' [', G_out%L(i, j, 1), '] ', j=1, Nsites)
      end do
   end subroutine print_output_gf

   subroutine print_input_gf()
      implicit none
      integer :: i, j

      write (*, '(/a)') '... GFr:'
      do i = 1, Nsites
         write (*, '(10(a,2f10.5,a))') (' [', G_in%R(i, j, 1), '] ', j=1, Nsites)
      end do

      write (*, '(/a)') '... GFl:'
      do i = 1, Nsites
         write (*, '(10(a,2f10.5,a))') (' [', G_in%L(i, j, 1), '] ', j=1, Nsites)
      end do
   end subroutine print_input_gf

   subroutine print_spectral_debug(iteration)
      !...spectral functions A_jj(w) of the first 10 sites -> sf_<iteration>.dat
      implicit none
      integer :: iteration, n, iw, j
      character :: fn*3

      n = Nsites; if (n .gt. 10) n = 10
      write (fn, '(i0)') iteration
      open (7, file='sf_'//trim(fn)//'.dat', status='unknown')
      do iw = 1, n_omega
         write (7, '(f10.5,x,10(e12.5,x))') omega(iw), (-2.d0*hbar*aimag(G_in%R(j, j, iw)), j=1, n)
      end do
      close (7)
      write (*, *) '... written sf_'//trim(fn)//'.dat'
   end subroutine print_spectral_debug


!.....G_in%a = (G_in%r)^dagger and G_in%G = G_in%L + G_in%r - G_in%a, from the mixed G_in%r, G_in%L
   subroutine complete_input_gf()
      implicit none
      integer :: iw
      do iw = 1, n_omega
         work1 = G_in%r(:, :, iw)
         call hermitian_conjugate(work1, Nsites, work2)
         G_in%a(:, :, iw) = work2
      end do
      G_in%G = G_in%L + G_in%R - G_in%A
   end subroutine complete_input_gf

end module GreenFunctions
