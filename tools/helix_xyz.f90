!==============================================================================
! tools/helix_xyz.f90  --  write the QUICHOTTE helix geometry to helix.xyz
!
! Usage:   helix_xyz.x  [input file]        (default: input.dat)
!
! Reads N_ions, N_turns, hel_length, hel_radius and hand from the same
! input file QUICHOTTE uses, and writes helix.xyz with exactly the
! site positions used by build_helix_hamiltonian (src/hamiltonian.f90). Handy for
! visualising the molecule, or as a starting point for an edited geometry
! to read back with build_hamiltonian_from_xyz.
!
! Build:   gfortran -O2 -o helix_xyz.x helix_xyz.f90   (or: make tools)
!==============================================================================
program generate_helix_xyz
   implicit none
   integer, parameter :: dp = kind(1.d0)
   integer :: N_ions, N_turns, Nsites, i, u, ios, eq
   real(dp) :: hel_radius, hel_length, hand
   real(dp) :: r(3)
   real(dp), parameter :: pi = 4.d0*atan(1.d0)
   character(len=256) :: fname, line, key

   !...defaults (overridden by the input file)
   N_ions = 4; N_turns = 2; hel_length = 5.d0; hel_radius = 1.d0; hand = 1.d0

   fname = 'input.dat'
   if (command_argument_count() .ge. 1) call get_command_argument(1, fname)
   open (newunit=u, file=trim(fname), status='old', action='read', iostat=ios)
   if (ios .eq. 0) then
      do
         read (u, '(A)', iostat=ios) line
         if (ios .ne. 0) exit
         eq = index(line, '=')
         if (eq .eq. 0 .or. adjustl(line(1:1)) .eq. '#') cycle
         key = adjustl(line(:eq - 1))
         select case (trim(key))
         case ('N_ions'); read (line(eq + 1:), *) N_ions
         case ('N_turns'); read (line(eq + 1:), *) N_turns
         case ('hel_length'); read (line(eq + 1:), *) hel_length
         case ('hel_radius'); read (line(eq + 1:), *) hel_radius
         case ('hand'); read (line(eq + 1:), *) hand
         end select
      end do
      close (u)
   else
      write (*, *) 'NOTE: ', trim(fname), ' not found, using built-in defaults'
   end if

   Nsites = N_ions*N_turns

   open (newunit=u, file='helix.xyz', status='replace')
   write (u, '(I0)') Nsites
   write (u, '(A,I0,A,I0,A,F6.2)') 'Helix: N_ions = ', N_ions, ', N_turns = ', N_turns, ', hand = ', hand
   do i = 1, Nsites
      !...same formula as helix_coordinate in src/hamiltonian.f90
      r(1) = hel_radius*cos(hand*(i - 1)*2*pi/N_ions)
      r(2) = hel_radius*sin(hand*(i - 1)*2*pi/N_ions)
      r(3) = hel_length*(i - 1)/(N_turns*(N_ions - 1))
      write (u, '(A2,3F16.8)') 'C ', r(1), r(2), r(3)
   end do
   close (u)

   print *, 'XYZ file generated: helix.xyz with', Nsites, 'sites'
end program generate_helix_xyz
