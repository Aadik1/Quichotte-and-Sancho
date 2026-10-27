!==============================================================================
! tools/stm_junction.f90  --  EXPERIMENTAL: STM tip + surface geometry
!
! Reads stm_junction_input.dat and writes STM_junction.xyz: a 3x3 / 2x2 / 1
! pyramidal tip above a (2*n_surface_x+1) x (2*n_surface_y+1) surface grid.
! Intended for the (unfinished) build_hamiltonian_from_xyz in src/hamiltonian.f90.
! el_bridge / n_bridge are read but not used yet.
!
! Build:   gfortran -O2 -o stm_junction.x stm_junction.f90   (or: make tools)
!==============================================================================
program stm_junction_pyramid
   implicit none
   integer, parameter :: dp = kind(1.d0)
   integer :: n_surface_x, n_surface_y, n_total, i, j
   real(dp) :: dx, dy, dz, x, y, z, hz
   character(len=2) :: el_tip, el_bridge, el_surface
   character(len=200) :: comment
   character(len=256) :: key, line
   integer :: ios, eq_pos

   open (10, file="stm_junction_input.dat", status="old")
   !......read parameters from input file
   do
      read (10, '(A)', iostat=ios) line
      if (ios .ne. 0) exit   ! end of file
      if (trim(line) .eq. '' .or. line(1:1) .eq. '#') cycle  ! skip blank or comment

      eq_pos = index(line, "=")
      if (eq_pos .eq. 0) cycle  ! skip malformed lines

      key = adjustl(trim(line(:eq_pos - 1)))

      select case (trim(key))
      case ("n_surface_x"); read (line(eq_pos + 1:), *) n_surface_x
      case ("n_surface_y"); read (line(eq_pos + 1:), *) n_surface_y
      case ("dx"); read (line(eq_pos + 1:), *) dx
      case ("dy"); read (line(eq_pos + 1:), *) dy
      case ("hz"); read (line(eq_pos + 1:), *) hz !...height of tip
      case ("el_tip"); read (line(eq_pos + 1:), *) el_tip
      case ("el_bridge"); read (line(eq_pos + 1:), *) el_bridge
      case ("el_surface"); read (line(eq_pos + 1:), *) el_surface
      case ("comment"); comment = adjustl(trim(line(eq_pos + 1:)))

      end select
   end do
   close (10)

   dz = 0.7d0*dx !...vertical seperation
   !...number of tip atoms in pyramid = sum over layers
   !...tip pyramid (9 + 4 + 1 atoms) + surface grid
   n_total = 14 + (2*n_surface_x + 1)*(2*n_surface_y + 1)

   open (11, file="STM_junction.xyz", status="unknown")
   write (11, '(I0)') n_total
   write (11, *) trim(comment)

   !...base layer of tip
   do j = 1, 3
      do i = 1, 3
         x = dx*(i - 1)
         y = dy*(j - 1)
         z = hz
         write (11, '(A2,3F10.4)') el_tip, x, y, z
      end do
   end do

   !...second layer of tip
   do j = 1, 3, 2
      do i = 1, 3, 2
         x = dx/2*(i)
         y = dy/2*(j)
         z = hz - dz
         write (11, '(A2,3F10.4)') el_tip, x, y, z
      end do
   end do

   !...tip point
   write (11, '(A2,3F10.4)') el_tip, dx, dy, hz - 2*dz

   !....surface atoms (grid at z=0)
   do j = -n_surface_y, n_surface_y
      do i = -n_surface_x, n_surface_x
         x = dx*(i - 1)
         y = dy*(j - 1)
         z = 0.0d0
         write (11, '(A2,3F10.4)') el_surface, x, y, z
      end do
   end do

   close (11)
   print *, "Wrote STM_junction.xyz with", n_total, "atoms."

end program stm_junction_pyramid
