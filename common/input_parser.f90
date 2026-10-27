!==============================================================================
! common/input_parser.f90  --  module InputParser
!
! Reader for 'key = value' input files, shared by QUICHOTTE and SANCHO.
!
!   call input_load(['input.dat   ', 'inputSOC.dat'])  ! first file that exists
!   call input_get('T', T)                        ! required key
!   call input_get('nhist', nhist, default=8)     ! optional key
!   call input_get('U_int', U, alias='Hubbard')   ! old name still accepted
!   call input_get('N_sites', n, alias='Nsites,Natoms')   ! several old names
!   call input_finish()                           ! report problems, stop if any
!
! File format: one 'key = value' per line, any order. Blank lines and lines
! starting with '#' are skipped; anything after '#' on a line is ignored.
! Keys are case-sensitive.
!
! input_finish reports EVERY problem at once: missing required keys (error),
! keys that were given twice (error), and keys that no input_get asked for
! (warning: usually a typo). Unreadable values stop immediately, with the
! line number.
!
! The files in common/ are IDENTICAL in QUICHOTTE and SANCHO
! (checked by `make check-common`).
!==============================================================================
module InputParser
   use, intrinsic :: iso_fortran_env, only: error_unit
   use Constants, only: dp
   implicit none
   private
   public :: input_load, input_get, input_finish, input_file_name

   integer, parameter :: max_keys = 200, klen = 64, vlen = 256

   character(len=klen) :: keys(max_keys)
   character(len=vlen) :: vals(max_keys)
   integer :: lines(max_keys)
   logical :: used(max_keys)
   integer :: nkeys = 0
   integer :: n_missing = 0
   logical :: duplicate_found = .false.
   character(len=vlen) :: input_file_name = ''

   interface input_get
      module procedure get_real, get_int, get_logical, get_string
   end interface input_get

contains

   subroutine input_load(candidates)
      !...read the first file of 'candidates' that exists into the key table
      character(len=*), intent(in) :: candidates(:)
      character(len=vlen) :: line, key, val
      integer :: k, u, ios, nline, eq, hash, j
      logical :: exists

      input_file_name = ''
      do k = 1, size(candidates)
         inquire (file=trim(candidates(k)), exist=exists)
         if (exists) then
            input_file_name = trim(candidates(k)); exit
         end if
      end do
      if (len_trim(input_file_name) .eq. 0) then
         write (error_unit, '(a)') ' ERROR: no input file found; looked for:'
         do k = 1, size(candidates)
            write (error_unit, '(4x,a)') trim(candidates(k))
         end do
         stop 1
      end if

      nkeys = 0; n_missing = 0; duplicate_found = .false.
      open (newunit=u, file=trim(input_file_name), status='old', action='read')
      nline = 0
      do
         read (u, '(A)', iostat=ios) line
         if (ios .ne. 0) exit
         nline = nline + 1
         line = adjustl(line)
         if (len_trim(line) .eq. 0 .or. line(1:1) .eq. '#') cycle
         hash = index(line, '#')
         if (hash .gt. 0) line = line(:hash - 1)
         eq = index(line, '=')
         if (eq .eq. 0) then
            write (*, '(a,i0,3a)') ' WARNING: line ', nline, ' of ', trim(input_file_name), &
               ' has no "=" (ignored): '//trim(line)
            cycle
         end if
         key = adjustl(line(:eq - 1))
         val = adjustl(line(eq + 1:))

         do j = 1, nkeys
            if (keys(j) .eq. key) then
               write (error_unit, '(3a,i0,a,i0)') ' ERROR: key "', trim(key), '" given twice, lines ', &
                  lines(j), ' and ', nline
               duplicate_found = .true.
            end if
         end do

         if (nkeys .eq. max_keys) then
            write (error_unit, *) 'ERROR: more than', max_keys, 'keys in ', trim(input_file_name)
            stop 1
         end if
         nkeys = nkeys + 1
         keys(nkeys) = key(:klen); vals(nkeys) = val; lines(nkeys) = nline; used(nkeys) = .false.
      end do
      close (u)
   end subroutine input_load

   integer function find_key(key, alias)
      !...index of key (or of one of its aliases) in the table, 0 if absent.
      !...alias may list several old names separated by commas: 'Nsites,Natoms'
      character(len=*), intent(in) :: key
      character(len=*), intent(in), optional :: alias
      integer :: j, s, c

      find_key = 0
      do j = 1, nkeys
         if (keys(j) .eq. key) then
            find_key = j; exit
         end if
      end do
      if (find_key .eq. 0 .and. present(alias)) then
         s = 1
         do while (s .le. len(alias) .and. find_key .eq. 0)
            c = index(alias(s:), ',')
            if (c .eq. 0) then
               c = len(alias) + 1
            else
               c = s + c - 1
            end if
            do j = 1, nkeys
               if (keys(j) .eq. adjustl(alias(s:c - 1))) then
                  find_key = j; exit
               end if
            end do
            s = c + 1
         end do
      end if
      if (find_key .gt. 0) used(find_key) = .true.
   end function find_key

   subroutine missing(key)
      character(len=*), intent(in) :: key
      write (error_unit, '(5a)') ' ERROR: required key "', key, '" missing from ', trim(input_file_name)
      n_missing = n_missing + 1
   end subroutine missing

   subroutine bad_value(j)
      integer, intent(in) :: j
      write (error_unit, '(3a,i0,2a)') ' ERROR: cannot read the value of "', trim(keys(j)), &
         '" on line ', lines(j), ' of ', trim(input_file_name)
      write (error_unit, '(3a)') '        value was: "', trim(vals(j)), '"'
      stop 1
   end subroutine bad_value

   subroutine get_real(key, var, default, alias)
      character(len=*), intent(in) :: key
      real(dp), intent(inout) :: var
      real(dp), intent(in), optional :: default
      character(len=*), intent(in), optional :: alias
      integer :: j, ios

      j = find_key(key, alias)
      if (j .eq. 0) then
         if (present(default)) then
            var = default
         else
            call missing(key)
         end if
         return
      end if
      read (vals(j), *, iostat=ios) var
      if (ios .ne. 0) call bad_value(j)
   end subroutine get_real

   subroutine get_int(key, var, default, alias)
      character(len=*), intent(in) :: key
      integer, intent(inout) :: var
      integer, intent(in), optional :: default
      character(len=*), intent(in), optional :: alias
      integer :: j, ios
      real(dp) :: r

      j = find_key(key, alias)
      if (j .eq. 0) then
         if (present(default)) then
            var = default
         else
            call missing(key)
         end if
         return
      end if
      read (vals(j), *, iostat=ios) r            !...accepts '2', '2.', '2.d0'
      if (ios .ne. 0 .or. abs(r - nint(r)) .gt. 1.d-9) call bad_value(j)
      var = nint(r)
   end subroutine get_int

   subroutine get_logical(key, var, default, alias)
      character(len=*), intent(in) :: key
      logical, intent(inout) :: var
      logical, intent(in), optional :: default
      character(len=*), intent(in), optional :: alias
      integer :: j, ios

      j = find_key(key, alias)
      if (j .eq. 0) then
         if (present(default)) then
            var = default
         else
            call missing(key)
         end if
         return
      end if
      read (vals(j), *, iostat=ios) var          !...T, F, .true., .false.
      if (ios .ne. 0) call bad_value(j)
   end subroutine get_logical

   subroutine get_string(key, var, default, alias)
      character(len=*), intent(in) :: key
      character(len=*), intent(inout) :: var
      character(len=*), intent(in), optional :: default
      character(len=*), intent(in), optional :: alias
      integer :: j

      j = find_key(key, alias)
      if (j .eq. 0) then
         if (present(default)) then
            var = default
         else
            call missing(key)
         end if
         return
      end if
      var = trim(vals(j))                        !...whole value, may contain '/'
   end subroutine get_string

   subroutine input_finish()
      !...warn about keys nobody asked for; stop if anything required was missing
      integer :: j

      do j = 1, nkeys
         if (.not. used(j)) write (*, '(3a,i0,a)') ' WARNING: unknown input key "', trim(keys(j)), &
            '" on line ', lines(j), ' (ignored -- typo?)'
      end do
      if (n_missing .gt. 0 .or. duplicate_found) then
         write (error_unit, '(a)') ' Stopping: fix the input file (see the examples/ folder).'
         stop 1
      end if
   end subroutine input_finish

end module InputParser
