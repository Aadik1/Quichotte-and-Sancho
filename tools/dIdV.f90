!==============================================================================
! dIdV.f90  --  differential conductance dI/dV from an I-V file
!
! Usage:   dIdV.x  [file]  [column]
!   file    I-V file (default Volt_Current_2.dat); column 1 must be V
!   column  which column holds the current (default 2; e.g. 4 = J_total
!           in QUICHOTTE's Volt_Current_<order>.dat)
!
! Writes dI_dV.dat: V, dI/dV. Central differences inside, one-sided at
! the two ends. Lines starting with '#' and blank lines are skipped.
!
! Build:   gfortran -O2 -o dIdV.x dIdV.f90      (or: make tools)
!==============================================================================
program dIdV
   implicit none
   integer, parameter :: dp = kind(1.d0)
   character(len=512) :: fname, arg, line
   integer :: col, n, i, ios, u
   real(dp), allocatable :: Volt(:), Curr(:), dIdV_(:), row(:)

   fname = 'Volt_Current_2.dat'
   col = 2
   if (command_argument_count() .ge. 1) call get_command_argument(1, fname)
   if (command_argument_count() .ge. 2) then
      call get_command_argument(2, arg)
      read (arg, *) col
   end if

   !...count the data lines
   open (newunit=u, file=trim(fname), status='old', action='read', iostat=ios)
   if (ios .ne. 0) then
      write (*, *) 'ERROR: cannot open ', trim(fname)
      stop 1
   end if
   n = 0
   do
      read (u, '(A)', iostat=ios) line
      if (ios .ne. 0) exit
      if (is_data(line)) n = n + 1
   end do
   if (n .lt. 2) then
      write (*, *) 'ERROR: need at least two data lines in ', trim(fname)
      stop 1
   end if

   !...read V (column 1) and I (column col)
   allocate (Volt(n), Curr(n), dIdV_(n), row(col))
   rewind (u)
   i = 0
   do
      read (u, '(A)', iostat=ios) line
      if (ios .ne. 0) exit
      if (.not. is_data(line)) cycle
      i = i + 1
      read (line, *) row(1:col)
      Volt(i) = row(1); Curr(i) = row(col)
   end do
   close (u)

   !...finite differences
   dIdV_(1) = (Curr(2) - Curr(1))/(Volt(2) - Volt(1))
   do i = 2, n - 1
      dIdV_(i) = (Curr(i + 1) - Curr(i - 1))/(Volt(i + 1) - Volt(i - 1))
   end do
   dIdV_(n) = (Curr(n) - Curr(n - 1))/(Volt(n) - Volt(n - 1))

   open (newunit=u, file='dI_dV.dat', status='replace')
   write (u, '(a,a,a,i0)') '# dI/dV from ', trim(fname), ', current column ', col
   do i = 1, n
      write (u, *) Volt(i), dIdV_(i)
   end do
   close (u)
   write (*, '(a,i0,a)') ' wrote dI_dV.dat (', n, ' points)'

contains

   logical function is_data(l)
      character(len=*), intent(in) :: l
      character(len=len(l)) :: t
      t = adjustl(l)
      is_data = len_trim(t) .gt. 0 .and. t(1:1) .ne. '#'
   end function is_data

end program dIdV
