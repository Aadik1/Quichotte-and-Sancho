!==============================================================================
! tests/run_tests.f90  --  regression-test runner (identical in QUICHOTTE and SANCHO)
!
! Usage (normally via `make test` / `make update-references`):
!     run_tests.x  <executable>  <tests directory>  [-u]
!
! <tests directory>/cases.txt lists the cases, one per line:
!     <case name>  <expected exit code>  <reference file> [<reference file> ...]
! ('#' starts a comment). For each case, every plain file in
! cases/<case name>/ (input.dat and any other inputs) is copied into
! work/<case name>/, the executable is run there, and each reference file in
! cases/<case name>/reference/ is compared number by number with the new output:
!
!     |new - ref| <= tol_abs + tol_rel*|ref|      tol_abs = 1e-6, tol_rel = 1e-5
!
! (override with the environment variables TOL_ABS and TOL_REL). Lines that are
! blank or start with '#' are skipped. The SCF tolerance is 1e-6, so changes
! of compiler, BLAS or thread count stay far below this; a real change in the
! physics is far above it.
!
! With -u the reference files are REPLACED by the new output. Do this only
! after a deliberate change to the physics, and commit the new references
! together with that change.
!
! Exit code: 0 if every case passed, 1 otherwise.
!==============================================================================
program run_tests
   implicit none
   integer, parameter :: dp = kind(1.d0), llen = 8192
   character(len=1024) :: exe, tdir, arg, case_name, work, cdir, cmd
   character(len=llen) :: line
   character(len=256) :: refs(50)
   integer :: nref, expected_rc, rc, u, ios, npass, nfail, k, t0, t1, rate
   real(dp) :: tol_abs, tol_rel, worst
   logical :: update, ok, file_ok

   if (command_argument_count() .lt. 2) then
      write (*, '(a)') 'usage: run_tests.x <executable> <tests directory> [-u]'
      stop 1
   end if
   call get_command_argument(1, exe)
   call get_command_argument(2, tdir)
   update = .false.
   if (command_argument_count() .ge. 3) then
      call get_command_argument(3, arg)
      update = (trim(arg) .eq. '-u')
   end if
   tol_abs = env_real('TOL_ABS', 1.d-6)
   tol_rel = env_real('TOL_REL', 1.d-5)

   open (newunit=u, file=trim(tdir)//'/cases.txt', status='old', action='read', iostat=ios)
   if (ios .ne. 0) then
      write (*, '(a)') 'cannot open '//trim(tdir)//'/cases.txt'
      stop 1
   end if

   npass = 0; nfail = 0
   do
      read (u, '(a)', iostat=ios) line
      if (ios .ne. 0) exit
      line = adjustl(line)
      if (len_trim(line) .eq. 0 .or. line(1:1) .eq. '#') cycle
      call parse_case_line(line, case_name, expected_rc, refs, nref)

      cdir = trim(tdir)//'/cases/'//trim(case_name)
      work = trim(tdir)//'/work/'//trim(case_name)
      cmd = 'rm -rf "'//trim(work)//'" && mkdir -p "'//trim(work)//'" && find "'//trim(cdir)// &
            '" -maxdepth 1 -type f -exec cp {} "'//trim(work)//'/" \;'
      call execute_command_line(trim(cmd), exitstat=rc)

      call system_clock(t0, rate)
      cmd = 'cd "'//trim(work)//'" && "'//trim(exe)//'" > stdout.txt 2>&1'
      call execute_command_line(trim(cmd), exitstat=rc)
      call system_clock(t1)

      if (update) then
         do k = 1, nref
            cmd = 'mkdir -p "'//trim(cdir)//'/reference" && cp "'//trim(work)//'/'//trim(refs(k))// &
                  '" "'//trim(cdir)//'/reference/"'
            call execute_command_line(trim(cmd))
         end do
         write (*, '(a,a,a,i0,a,f6.1,a)') 'updated  ', trim(case_name), '  (exit ', rc, ', ', &
            real(t1 - t0)/real(rate), ' s)'
         cycle
      end if

      ok = .true.
      if (rc .ne. expected_rc) then
         write (*, '(a,i0,a,i0,a)') '  exit code ', rc, ', expected ', expected_rc, &
            '  (see '//trim(work)//'/stdout.txt)'
         ok = .false.
      end if
      do k = 1, nref
         call compare_files(trim(work)//'/'//trim(refs(k)), trim(cdir)//'/reference/'//trim(refs(k)), &
                            tol_abs, tol_rel, file_ok, worst)
         if (file_ok) then
            write (*, '(a,a,a,es9.2)') '  ', trim(refs(k)), ':  max |diff| = ', worst
         else
            ok = .false.
         end if
      end do
      if (ok) then
         npass = npass + 1
         write (*, '(a,a,a,f6.1,a)') 'PASS  ', trim(case_name), '  (', real(t1 - t0)/real(rate), ' s)'
      else
         nfail = nfail + 1
         write (*, '(a,a,a)') 'FAIL  ', trim(case_name), '  -> '//trim(work)
      end if
   end do
   close (u)

   if (update) stop
   write (*, '(a)') '----------------------------------------'
   write (*, '(i0,a,i0,a)') npass, ' passed, ', nfail, ' failed'
   if (nfail .gt. 0) error stop 1

contains

   subroutine parse_case_line(l, name, rc_expected, files, nf)
      character(len=*), intent(in) :: l
      character(len=*), intent(out) :: name, files(:)
      integer, intent(out) :: rc_expected, nf
      character(len=256) :: tok(60)
      integer :: nt

      call split(l, tok, nt)
      name = tok(1)
      read (tok(2), *) rc_expected
      nf = nt - 2
      files(1:nf) = tok(3:nt)
   end subroutine parse_case_line

   subroutine split(l, tok, nt)
      !...whitespace-separated tokens of l
      character(len=*), intent(in) :: l
      character(len=*), intent(out) :: tok(:)
      integer, intent(out) :: nt
      integer :: i, s
      logical :: inword

      nt = 0; inword = .false.; s = 1
      do i = 1, len_trim(l) + 1
         if (i .le. len_trim(l) .and. l(i:i) .ne. ' ' .and. l(i:i) .ne. char(9)) then
            if (.not. inword) then
               inword = .true.; s = i
            end if
         else if (inword) then
            inword = .false.
            if (nt .lt. size(tok)) then
               nt = nt + 1; tok(nt) = l(s:i - 1)
            end if
         end if
      end do
   end subroutine split

   logical function data_line(l)
      character(len=*), intent(in) :: l
      character(len=len(l)) :: t
      t = adjustl(l)
      data_line = len_trim(t) .gt. 0 .and. t(1:1) .ne. '#'
   end function data_line

   subroutine next_data_line(unit, l, eof)
      integer, intent(in) :: unit
      character(len=*), intent(out) :: l
      logical, intent(out) :: eof
      integer :: ios
      eof = .false.
      do
         read (unit, '(a)', iostat=ios) l
         if (ios .ne. 0) then
            eof = .true.; return
         end if
         if (data_line(l)) return
      end do
   end subroutine next_data_line

   subroutine compare_files(fnew, fref, ta, tr, ok, worst)
      !...numeric comparison, field by field, of two whitespace-separated files
      character(len=*), intent(in) :: fnew, fref
      real(dp), intent(in) :: ta, tr
      logical, intent(out) :: ok
      real(dp), intent(out) :: worst
      character(len=llen) :: ln, lr
      character(len=64) :: tn(600), tref(600)
      integer :: un, ur, ios, nn, nr, i, nline, nbad
      real(dp) :: a, b, d
      logical :: eofn, eofr

      ok = .true.; worst = 0.d0; nline = 0; nbad = 0
      open (newunit=un, file=fnew, status='old', action='read', iostat=ios)
      if (ios .ne. 0) then
         write (*, '(a)') '  missing output '//fnew
         ok = .false.; return
      end if
      open (newunit=ur, file=fref, status='old', action='read', iostat=ios)
      if (ios .ne. 0) then
         write (*, '(a)') '  missing reference '//fref
         close (un); ok = .false.; return
      end if

      do
         call next_data_line(un, ln, eofn)
         call next_data_line(ur, lr, eofr)
         if (eofn .and. eofr) exit
         nline = nline + 1
         if (eofn .neqv. eofr) then
            write (*, '(a,i0,a)') '  '//fnew//': different number of data lines (after ', nline - 1, ')'
            ok = .false.; exit
         end if
         call split(ln, tn, nn)
         call split(lr, tref, nr)
         if (nn .ne. nr) then
            write (*, '(a,i0,a,i0,a,i0)') '  '//fnew//' line ', nline, ': ', nn, ' columns, expected ', nr
            ok = .false.; exit
         end if
         do i = 1, nn
            read (tn(i), *, iostat=ios) a
            if (ios .ne. 0) cycle                 !...non-numeric field: not compared
            read (tref(i), *, iostat=ios) b
            if (ios .ne. 0) cycle
            d = abs(a - b)
            worst = max(worst, d)
            if (d .gt. ta + tr*abs(b)) then
               nbad = nbad + 1
               if (nbad .le. 3) write (*, '(a,i0,a,i0,a,es16.8,a,es16.8)') '  '//fnew//' line ', nline, &
                  ' col ', i, ': ', a, ' vs ref ', b
               ok = .false.
            end if
         end do
      end do
      close (un); close (ur)
   end subroutine compare_files

   real(dp) function env_real(name, default)
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: default
      character(len=64) :: v
      integer :: st, ios
      call get_environment_variable(name, v, status=st)
      env_real = default
      if (st .eq. 0 .and. len_trim(v) .gt. 0) then
         read (v, *, iostat=ios) env_real
         if (ios .ne. 0) env_real = default
      end if
   end function env_real

end program run_tests
