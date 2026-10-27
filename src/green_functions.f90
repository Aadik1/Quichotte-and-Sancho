!==============================================================================
! src/green_functions.f90  --  module GreenFunctions   (QUICHOTTE)
!
! Storage for the Green's functions on the energy grid and the routines that
! do not involve the interaction:
!   fermi                      Fermi-Dirac distribution
!   noninteracting_gr_ga       G0^r, G0^a on the whole grid
!   noninteracting_gl_gg       G0^<, G0^> on the whole grid
!   complete_input_gf          G^a and G^> of G_in from its G^r and G^<
!   equal_time_lesser          G^<(t=0) = (delta/2pi) sum_w G^<(w), on-site blocks
!   print_input_gf / print_output_gf   debug printing
!
! Naming: G_in  = INPUT Green's functions of the current SCF step
!         G_out = OUTPUT Green's functions built from the self-energy of G_in
!         components %r retarded, %a advanced, %L lesser, %G greater,
!         each (Norb, Norb, n_omega).
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
   real(dp), allocatable, dimension(:) :: U_site   !...Hubbard U per spin-orbital (set on the spin-up index)
   real(dp), allocatable, dimension(:) :: Ev       !...eigenvalues of H
   complex(dp), allocatable, dimension(:, :) :: Eigenvec
   complex(dp), allocatable, dimension(:, :) :: GammaL, GammaR   !...lead coupling matrices
   complex(dp), allocatable, dimension(:, :) :: G_less_t0       !...G^<(t=0), on-site 2x2 spin blocks
   complex(dp), allocatable, dimension(:, :) :: work1, work2, work3
   integer :: INFO

contains

   !.....fermi-dirac distribution
   real(dp) function fermi(w, V)
      implicit none
      real(dp) :: arg, w, V

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
         do i = 1, Norb
            work1(i, i) = work1(i, i) + hbar*w
         end do

         call invert_complex(Norb, work1, info)
         call hermitian_conjugate(work1, Norb, work2)

         G_in%r(:, :, j) = work1
         G_in%a(:, :, j) = work2
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

!.....G_in%a = (G_in%r)^dagger and G_in%G = G_in%L + G_in%r - G_in%a, from the mixed G_in%r, G_in%L
   subroutine complete_input_gf()
      implicit none
      integer :: iw
      do iw = 1, n_omega
         work1 = G_in%r(:, :, iw)
         call hermitian_conjugate(work1, Norb, work2)
         G_in%a(:, :, iw) = work2
      end do
      G_in%G = G_in%L + G_in%R - G_in%A
   end subroutine complete_input_gf

   subroutine equal_time_lesser()
      !.....equal-time lesser GF, G_less_t0 = G^<(t=0) = (delta/2pi) sum_w G^<(w),
      !.....only the on-site 2x2 spin blocks (all that Sigma1 and spin_texture need)
      implicit none
      integer :: i, k1, s1, s2
      complex(dp) :: s
      real(dp) :: pp

      pp = delta/(2.d0*pi)
      G_less_t0 = (0.d0, 0.d0)
      do i = 1, Norb, 2
         do s1 = 0, 1
            do s2 = 0, 1

               s = (0.d0, 0.d0)
               do k1 = 1, n_omega
                  s = s + G_in%L(i + s1, i + s2, k1)
               end do
               G_less_t0(i + s1, i + s2) = s*pp

            end do
         end do
      end do

   end subroutine equal_time_lesser

!...Debug print routines (not called by default)

   subroutine print_output_gf()
      implicit none
      integer :: i, j

      write (*, '(/a)') '... GFr:'
      do i = 1, Norb, 2
         write (*, '(10(a,2f10.5,a))') (' [', G_out%R(i, j, 1), '] ', j=1, Norb, 2)
      end do

      write (*, '(/a)') '... GFl:'
      do i = 1, Norb, 2
         write (*, '(10(a,2f10.5,a))') (' [', G_out%L(i, j, 1), '] ', j=1, Norb, 2)
      end do
   end subroutine print_output_gf

   subroutine print_input_gf()
      implicit none
      integer :: i, j

      write (*, '(/a)') '... GFr:'
      do i = 1, Norb, 2
         write (*, '(10(a,2f10.5,a))') (' [', G_in%R(i, j, 1), '] ', j=1, Norb, 2)
      end do

      write (*, '(/a)') '... GFl:'
      do i = 1, Norb, 2
         write (*, '(10(a,2f10.5,a))') (' [', G_in%L(i, j, 1), '] ', j=1, Norb, 2)
      end do
   end subroutine print_input_gf

end module GreenFunctions
