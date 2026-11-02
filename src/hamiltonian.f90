!==============================================================================
! src/hamiltonian.f90  --  module Hamiltonian   (QUICHOTTE)
!
! Central-region Hamiltonian H (Norb x Norb). The input key 'hamiltonian'
! chooses the model; build_hamiltonian (below) calls the matching routine:
!
!   helix   build_helix_hamiltonian       helical molecule with SOC (default)
!   xyz     build_hamiltonian_from_xyz    same terms, site positions from xyz_file
!   dftb    build_hamiltonian_from_dftb   matrix read from ham_file (soc_dftb block format)
!   custom  build_custom_hamiltonian      YOUR model: edit the routine at the end
!
! Helix model:
!   H = sum_i E_CC n_i
!     + t_hop * sum_<i,i+1>   (spin-diagonal nearest-neighbour hopping)
!     + i*lamb * sum_<i,i+2>  (d_{i,i+1} x d_{i,i+2}) . sigma
! with d the unit bond vectors between sites.
!
! Rules every model must follow (details: manual, "Using your own Hamiltonian"):
!   * basis order 2i-1 = (site i, up), 2i = (site i, down)
!   * Hermitian: fill the upper triangle, then call make_hermitian()
!   * energies in eV; no Hubbard/Hartree terms (U enters via the self-energy)
!   * leads: every orbital of lead_left_atoms / lead_right_atoms (default: site 1, site N_sites)
!==============================================================================
module Hamiltonian
   use Constants, only: dp, pi, im
   use LinearAlgebra, only: cross_product, norm3
   use Parameters
   implicit none

   complex(dp), allocatable, dimension(:, :) :: H   !...central Hamiltonian (Norb x Norb)

contains

   !===========================================================================
   ! Dispatcher: builds the model chosen by the input key 'hamiltonian'.
   ! To add a model of your own: write its routine below, then add ONE case
   ! here, e.g.   case ('ssh');  call build_ssh_hamiltonian()
   ! and select it with  hamiltonian = ssh  (plus N_sites) in input.dat.
   ! H must already be allocated as H(Norb, Norb).
   !===========================================================================
   subroutine build_hamiltonian()
      select case (trim(hamiltonian_model))
      case ('helix')
         call build_helix_hamiltonian()
      case ('xyz')
         call build_hamiltonian_from_xyz(trim(xyz_file))
      case ('dftb')
         call build_hamiltonian_from_dftb(trim(ham_file))
      case ('custom')
         call build_custom_hamiltonian()
      case default
         print *, 'ERROR: unknown hamiltonian = ', trim(hamiltonian_model)
         print *, '       built in: helix, xyz, dftb, custom (see src/hamiltonian.f90)'
         stop 1
      end select
   end subroutine build_hamiltonian

   integer function xyz_site_count(xyz_file_name)
      !...number of atoms (first line) of an .xyz file; used to set N_sites
      character(len=*), intent(in) :: xyz_file_name
      integer :: u, ios

      open (newunit=u, file=trim(xyz_file_name), status='old', action='read', iostat=ios)
      if (ios .ne. 0) then
         print *, 'ERROR: cannot open xyz_file ', trim(xyz_file_name)
         stop 1
      end if
      read (u, *, iostat=ios) xyz_site_count
      close (u)
      if (ios .ne. 0 .or. xyz_site_count .lt. 1) then
         print *, 'ERROR: first line of ', trim(xyz_file_name), ' must be the number of atoms'
         stop 1
      end if
   end function xyz_site_count

   subroutine helix_coordinate(r, i)
      !...Cartesian position of site i = 1 .. N_ions*N_turns on the helix:
      !     x = R cos(hand * 2 pi (i-1)/N_ions),  y = R sin(hand * 2 pi (i-1)/N_ions)
      !     z = (i-1) * dz,   dz = hel_length / (N_turns*(N_ions-1))
      !   i.e. hel_length/N_turns is the rise over the first N_ions-1 bonds of a
      !   turn, so for one turn hel_length is the end-to-end length exactly; for
      !   N_turns turns the end-to-end length is hel_length*(N_sites-1)/(N_turns*(N_ions-1)).
      implicit none
      integer, intent(in) :: i
      real(dp), intent(out) :: r(3)

      r(1) = hel_radius*cos(hand*(i - 1)*2*pi/N_ions)
      r(2) = hel_radius*sin(hand*(i - 1)*2*pi/N_ions)
      r(3) = hel_length*(i - 1)/(N_turns*(N_ions - 1))
   end subroutine helix_coordinate

   subroutine soc_matrix_elements(r, r1, r2, ii, jj)
      !...SOC block between site i (spin-orbitals ii-1, ii) and site i+2
      !   (jj-1, jj), from the positions r, r1, r2 of sites i, i+1, i+2:
      !   H_ij = i*lamb * (d1 x d2) . sigma
      implicit none
      real(dp), intent(in) :: r(3), r1(3), r2(3)
      integer, intent(in) :: ii, jj
      real(dp) :: d1(3), d2(3), v(3)

      d1 = (r - r1)/norm3(r - r1)
      d2 = (r - r2)/norm3(r - r2)
      call cross_product(d1, d2, v)
      H(ii - 1, jj - 1) = im*lamb*v(3)
      H(ii - 1, jj) = im*lamb*(v(1) - im*v(2))
      H(ii, jj - 1) = im*lamb*(v(1) + im*v(2))
      H(ii, jj) = im*lamb*(-v(3))
   end subroutine soc_matrix_elements

   subroutine make_hermitian()
      !...fill the lower triangle of H from the upper triangle
      implicit none
      integer :: i, j

      do i = 1, Norb
         do j = i, Norb
            H(j, i) = conjg(H(i, j))
         end do
      end do
   end subroutine make_hermitian

   subroutine build_helix_hamiltonian()
      !...Hamiltonian of the analytic helix (positions from helix_coordinate)
      implicit none
      integer :: i, ii, jj, Nat
      real(dp) :: r(3), r1(3), r2(3)

      H = (0.0d0, 0.0d0)
      Nat = Norb/2

      do i = 1, Nat
         ii = 2*i

         !...on-site
         H(ii - 1, ii - 1) = E_CC
         H(ii, ii) = E_CC

         !...nearest-neighbour hopping
         if (i .le. Nat - 1) then
            jj = 2*(i + 1)
            H(ii - 1, jj - 1) = t_hop
            H(ii, jj) = t_hop
         end if

         !...next-nearest-neighbour spin-orbit coupling
         if (i .le. Nat - 2) then
            jj = 2*(i + 2)
            call helix_coordinate(r, i)
            call helix_coordinate(r1, i + 1)
            call helix_coordinate(r2, i + 2)
            call soc_matrix_elements(r, r1, r2, ii, jj)
         end if
      end do

      call make_hermitian()
   end subroutine build_helix_hamiltonian

   subroutine build_hamiltonian_from_xyz(xyz_path)
      !...Same model as build_helix_hamiltonian, but site positions are read from an
      !   .xyz file (e.g. one written by tools/helix_xyz.f90).
      !   The file must contain N_sites = Norb/2 atoms (read_input sets N_sites
      !   from it); H must already be allocated. Input: hamiltonian = xyz, xyz_file = ...
      implicit none
      character(len=*), intent(in) :: xyz_path
      integer :: i, ii, jj, Nat, ierr, unit_xyz
      real(dp), allocatable :: coords(:, :)
      character(len=2) :: atom_type
      character(len=200) :: comment_line

      open (newunit=unit_xyz, file=trim(xyz_path), status='old', action='read', iostat=ierr)
      if (ierr /= 0) then
         print *, 'ERROR: cannot open XYZ file: ', trim(xyz_path)
         stop 1
      end if

      read (unit_xyz, *) Nat
      if (2*Nat .ne. Norb) then
         print *, 'ERROR: ', trim(xyz_path), ' has', Nat, 'atoms but the input gives', Norb/2, 'sites'
         stop 1
      end if
      allocate (coords(3, Nat))
      read (unit_xyz, '(A)') comment_line
      do i = 1, Nat
         read (unit_xyz, *) atom_type, coords(1:3, i)
      end do
      close (unit_xyz)

      H = (0.0d0, 0.0d0)
      do i = 1, Nat
         ii = 2*i

         H(ii - 1, ii - 1) = E_CC
         H(ii, ii) = E_CC

         if (i .le. Nat - 1) then
            jj = 2*(i + 1)
            H(ii - 1, jj - 1) = t_hop
            H(ii, jj) = t_hop
         end if

         if (i .le. Nat - 2) then
            jj = 2*(i + 2)
            call soc_matrix_elements(coords(:, i), coords(:, i + 1), coords(:, i + 2), ii, jj)
         end if
      end do

      call make_hermitian()
      deallocate (coords)
   end subroutine build_hamiltonian_from_xyz


   !===========================================================================
   ! DFTB HAMILTONIAN  (input: hamiltonian = dftb, ham_file = ...)
   !
   ! Reads the orthogonal (Lowdin) spin-orbit Hamiltonian written by the
   ! soc_dftb pipeline (dftb_soc_hamiltonian.py, export_fortran_blocks):
   !
   !   # comment lines
   !   BLOCK n | Atom_I: a Atom_J: b | Dist: d A | Dim: (2 n_a)x(2 n_b)
   !     re  im          <- (2 n_a)*(2 n_b) lines, row by row
   !     ...
   !
   ! Rows/columns of a block run over (orbital, spin) of the atom, spin fastest,
   ! orbitals in DFTB+ order (s, p_y, p_z, p_x) -- the same layout as QUICHOTTE
   ! (site = orbital, 2i-1 up, 2i down). Energies in eV. Blocks below the export
   ! cutoff are absent and stay zero. Every atom needs its on-site block (a, a).
   !
   ! N_atoms and the orbitals per atom are read from the block headers by
   ! dftb_blocks_scan (called from read_input); this routine fills H, reports
   ! the Hermiticity error, symmetrises H, and reports time-reversal violation.
   !===========================================================================
   subroutine dftb_blocks_scan(path, n_at, norb)
      !...number of atoms and spatial orbitals per atom, from the BLOCK headers
      character(len=*), intent(in) :: path
      integer, intent(out) :: n_at
      integer, allocatable, intent(out) :: norb(:)
      character(len=512) :: line
      integer :: u, ios, ai, aj, di, dj, pass
      logical :: ok

      n_at = 0
      do pass = 1, 2
         open (newunit=u, file=trim(path), status='old', action='read', iostat=ios)
         if (ios .ne. 0) then
            print *, 'ERROR: cannot open ham_file ', trim(path)
            stop 1
         end if
         if (pass .eq. 2) then
            allocate (norb(n_at)); norb = 0
         end if
         do
            read (u, '(a)', iostat=ios) line
            if (ios .ne. 0) exit
            if (index(adjustl(line), 'BLOCK') .ne. 1) cycle
            call parse_block_header(line, ai, aj, di, dj, ok)
            if (.not. ok) then
               print *, 'ERROR: cannot read block header in ', trim(path), ': ', trim(line)
               stop 1
            end if
            if (mod(di, 2) .ne. 0 .or. mod(dj, 2) .ne. 0) then
               print *, 'ERROR: block dimensions must be even (orbitals x 2 spins): ', trim(line)
               stop 1
            end if
            if (pass .eq. 1) then
               n_at = max(n_at, ai, aj)
            else
               call set_norb(ai, di/2); call set_norb(aj, dj/2)
            end if
         end do
         close (u)
      end do
      if (n_at .lt. 1 .or. any(norb .eq. 0)) then
         print *, 'ERROR: ', trim(path), ' has no BLOCK for some atom (every atom needs its on-site block)'
         stop 1
      end if
   contains
      subroutine set_norb(a, n)
         integer, intent(in) :: a, n
         if (norb(a) .ne. 0 .and. norb(a) .ne. n) then
            print *, 'ERROR: atom', a, ' has inconsistent block dimensions in ', trim(path)
            stop 1
         end if
         norb(a) = n
      end subroutine set_norb
   end subroutine dftb_blocks_scan

   subroutine parse_block_header(line, ai, aj, di, dj, ok)
      !...'BLOCK n | Atom_I: a Atom_J: b | Dist: d A | Dim: 8x8'  ->  a, b, 8, 8
      character(len=*), intent(in) :: line
      integer, intent(out) :: ai, aj, di, dj
      logical, intent(out) :: ok
      integer :: k1, k2, k3, kx, ios

      ok = .false.; ai = 0; aj = 0; di = 0; dj = 0
      k1 = index(line, 'Atom_I:'); k2 = index(line, 'Atom_J:'); k3 = index(line, 'Dim:')
      if (k1 .eq. 0 .or. k2 .eq. 0 .or. k3 .eq. 0) return
      read (line(k1 + 7:), *, iostat=ios) ai; if (ios .ne. 0) return
      read (line(k2 + 7:), *, iostat=ios) aj; if (ios .ne. 0) return
      kx = index(line(k3:), 'x') + k3 - 1
      if (kx .lt. k3) return
      read (line(k3 + 4:kx - 1), *, iostat=ios) di; if (ios .ne. 0) return
      read (line(kx + 1:), *, iostat=ios) dj; if (ios .ne. 0) return
      ok = (ai .ge. 1 .and. aj .ge. 1 .and. di .ge. 1 .and. dj .ge. 1)
   end subroutine parse_block_header

   subroutine build_hamiltonian_from_dftb(path)
      !...fill H from the block file (layout: see the header of this section)
      character(len=*), intent(in) :: path
      character(len=512) :: line
      integer :: u, ios, ai, aj, di, dj, r, c, r0, c0, nblk
      real(dp) :: re, imv, herm
      logical :: ok
      logical, allocatable :: onsite(:)

      allocate (onsite(N_atoms)); onsite = .false.
      H = (0.d0, 0.d0); nblk = 0
      open (newunit=u, file=trim(path), status='old', action='read', iostat=ios)
      if (ios .ne. 0) then
         print *, 'ERROR: cannot open ham_file ', trim(path)
         stop 1
      end if
      do
         read (u, '(a)', iostat=ios) line
         if (ios .ne. 0) exit
         if (index(adjustl(line), 'BLOCK') .ne. 1) cycle
         call parse_block_header(line, ai, aj, di, dj, ok)
         if (.not. ok .or. di .ne. 2*atom_norb(ai) .or. dj .ne. 2*atom_norb(aj)) then
            print *, 'ERROR: bad block header in ', trim(path), ': ', trim(line)
            stop 1
         end if
         r0 = 2*(atom_first(ai) - 1); c0 = 2*(atom_first(aj) - 1)
         do r = 1, di
            do c = 1, dj
               read (u, *, iostat=ios) re, imv
               if (ios .ne. 0) then
                  print *, 'ERROR: block too short in ', trim(path), ': ', trim(line)
                  stop 1
               end if
               H(r0 + r, c0 + c) = cmplx(re, imv, kind=dp)
            end do
         end do
         if (ai .eq. aj) onsite(ai) = .true.
         nblk = nblk + 1
      end do
      close (u)

      if (.not. all(onsite)) then
         print *, 'ERROR: ', trim(path), ' is missing the on-site block of some atom'
         stop 1
      end if
      herm = maxval(abs(H - conjg(transpose(H))))
      H = 0.5d0*(H + conjg(transpose(H)))
      write (*, '(a,i0,a,i0,a,i0,a)') ' DFTB Hamiltonian: ', nblk, ' blocks, ', N_atoms, ' atoms, ', &
         N_sites, ' orbitals'
      write (*, '(a,es10.2,a)') ' max |H - H^dagger| in file = ', herm, ' eV (H symmetrised)'
      call report_time_reversal()
      deallocate (onsite)
   end subroutine build_hamiltonian_from_dftb

   subroutine report_time_reversal()
      !...A non-magnetic H is time-reversal invariant: every 2x2 spin block
      !   h0*1 + h.sigma has h0 real and h purely imaginary. Im(h0) and Re(h)
      !   act like random magnetic fields and can fake spin polarisation, so
      !   report them next to the genuine spin-orbit part Im(h).
      integer :: i, j
      complex(dp) :: b11, b12, b21, b22, h0, hx, hy, hz
      real(dp) :: odd0, oddh, soc

      odd0 = 0.d0; oddh = 0.d0; soc = 0.d0
      do j = 1, N_sites
         do i = 1, N_sites
            b11 = H(2*i - 1, 2*j - 1); b12 = H(2*i - 1, 2*j)
            b21 = H(2*i, 2*j - 1);     b22 = H(2*i, 2*j)
            h0 = 0.5d0*(b11 + b22); hz = 0.5d0*(b11 - b22)
            hx = 0.5d0*(b12 + b21); hy = 0.5d0*im*(b12 - b21)
            odd0 = max(odd0, abs(aimag(h0)))
            oddh = max(oddh, abs(real(hx)), abs(real(hy)), abs(real(hz)))
            soc = max(soc, abs(aimag(hx)), abs(aimag(hy)), abs(aimag(hz)))
         end do
      end do
      write (*, '(a)') ' Time reversal (all spin blocks h0*1 + h.sigma):'
      write (*, '(a,es10.2,a)') '   max |Im h0| =', odd0, ' eV   (should be 0)'
      write (*, '(a,es10.2,a)') '   max |Re h|  =', oddh, ' eV   (should be 0: Zeeman-like)'
      write (*, '(a,es10.2,a)') '   max |Im h|  =', soc, ' eV   (spin-orbit coupling)'
      if (max(odd0, oddh) .gt. 1.d-8) write (*, '(a)') &
         ' WARNING: H breaks time reversal by more than 1e-8 eV; spin polarisation of that size is not physical'
   end subroutine report_time_reversal


   !===========================================================================
   ! YOUR MODEL  (input: hamiltonian = custom, N_sites = ...)
   !
   ! Replace the body with your own Hamiltonian. As shipped it is a plain
   ! spin-degenerate chain (E_CC on site, t_hop between neighbours, no SOC),
   ! so 'hamiltonian = custom' runs out of the box.
   !
   ! Rules: index 2i-1 = (site i, up), 2i = (site i, down), i = 1..N_sites;
   ! fill the UPPER triangle (column >= row) and finish with make_hermitian();
   ! energies in eV; no Hubbard/Hartree terms; leads attach to site 1 and site N_sites.
   ! New parameters: declare them in parameters.f90 and read them in read_input
   ! (input_output.f90) with  call input_get('my_key', my_var[, default=...]).
   !===========================================================================
   subroutine build_custom_hamiltonian()
      integer :: i, ii

      H = (0.d0, 0.d0)
      do i = 1, N_sites
         ii = 2*i                      !...ii-1 = (i,up), ii = (i,down)
         H(ii - 1, ii - 1) = E_CC      !...on-site energy, spin up
         H(ii, ii) = E_CC              !...on-site energy, spin down
         if (i .lt. N_sites) then      !...hop i -> i+1, spin conserving
            H(ii - 1, ii + 1) = t_hop  !...up   -> up
            H(ii, ii + 2) = t_hop      !...down -> down
         end if
      end do
      call make_hermitian()            !...upper triangle -> whole matrix
   end subroutine build_custom_hamiltonian

end module Hamiltonian
