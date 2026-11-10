!==============================================================================
! src/hamiltonian.f90  --  module Hamiltonian   (QUICHOTTE)
!
! Central-region Hamiltonian H (Norb x Norb). The input key 'hamiltonian'
! chooses the model; build_hamiltonian (below) calls the matching routine:
!
!   helix   build_helix_hamiltonian       helical molecule with SOC (default)
!   peptide build_peptide_hamiltonian     alpha-helix, one site per residue, long-range hopping + SOC
!   xyz     build_hamiltonian_from_xyz    same terms, site positions from xyz_file
!   dftb    build_hamiltonian_from_dftb   matrix read from ham_file (soc_dftb block format)
!   dftb_hs build_hamiltonian_from_dftb_hs  H0 and S from DFTB+ (hamsqr1.dat, oversqr.dat),
!                                           atomic spin-orbit coupling, Lowdin orthogonalisation
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
      case ('peptide')
         call build_peptide_hamiltonian()
      case ('xyz')
         call build_hamiltonian_from_xyz(trim(xyz_file))
      case ('dftb')
         call build_hamiltonian_from_dftb(trim(ham_file))
      case ('dftb_hs')
         call build_hamiltonian_from_dftb_hs()
      case ('custom')
         call build_custom_hamiltonian()
      case default
         print *, 'ERROR: unknown hamiltonian = ', trim(hamiltonian_model)
         print *, '       built in: helix, peptide, xyz, dftb, dftb_hs, custom (see src/hamiltonian.f90)'
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

   !===========================================================================
   ! PEPTIDE MODEL  (input: hamiltonian = peptide, N_sites = ...)
   !
   ! alpha-helical peptide, one site per residue (after Guo & Sun, PNAS 111,
   ! 11658 (2014)). Unlike 'helix', the twist need not divide 2*pi (3.6 sites/turn):
   !   r_i   = ( R cos phi_i, R sin phi_i, (i-1)*pep_rise ),  phi_i = hand*(i-1)*pep_twist
   !   H     = sum_i E_CC n_i
   !         + sum_i sum_{j=1..hop_range} t_j  c+_i c_{i+j}  (+ h.c., spin diagonal)
   !         + SOC, soc_model:
   !             radial : sum_i sum_{j=1..hop_range} i*lamb_j c+_i (nu_{i,i+j}.sigma) c_{i+j} + h.c.
   !             torsion: i*lamb sum_<i,i+2> (d_{i,i+1} x d_{i,i+2}).sigma  (exactly as 'helix')
   !   t_j    = t_hop exp(-(l - l_1)/hop_decay),  lamb_j = lamb exp(-(l - l_1)/hop_decay)
   !   l      = |r_{i+j} - r_i|, l_1 the nearest-neighbour distance (t_1 = t_hop, lamb_1 = lamb)
   !   nu     = rhat(phibar) x dhat, phibar = (phi_i + phi_{i+j})/2: radial (confining) field
   !            x bond direction = Guo-Sun's (sin phibar sin th, -cos phibar sin th, cos th).
   !            Their sigma_i + sigma_{i+j} convention means their s_j = lamb_j/2.
   !            phibar is used, not the chord midpoint: for j = 2 (200 deg) the chord
   !            midpoint lies on the far side of the axis and would flip the sign.
   ! alpha-helix: hel_radius = 2.5, pep_twist = 100 (default), pep_rise = 1.5 (default)
   !   -> l_1..l_4 = 4.11, 5.77, 5.15, 6.24 A (note l_3 < l_2).
   ! Writes peptide_geometry.xyz and prints l_j and t_j/t_hop.
   !===========================================================================
   subroutine build_peptide_hamiltonian()
      implicit none
      integer :: i, j, m, ii, jj, Nat, u, jmax
      real(dp) :: dphi, l1, l, fac, lj, phibar
      real(dp) :: d(3), dhat(3), rhat(3), nu(3)
      real(dp), allocatable :: rs(:, :), phi(:)

      Nat = Norb/2
      jmax = min(hop_range, Nat - 1)
      dphi = hand*pep_twist*pi/180.d0
      allocate (rs(3, Nat), phi(Nat))
      do i = 1, Nat
         phi(i) = (i - 1)*dphi
         rs(:, i) = [hel_radius*cos(phi(i)), hel_radius*sin(phi(i)), (i - 1)*pep_rise]
      end do
      l1 = norm3(rs(:, 2) - rs(:, 1))

      H = (0.0d0, 0.0d0)

      !...on-site
      do i = 1, Nat
         ii = 2*i
         H(ii - 1, ii - 1) = E_CC
         H(ii, ii) = E_CC
      end do

      !...torsion SOC first: soc_matrix_elements ASSIGNS its block, the terms below ADD
      if (trim(soc_model) .eq. 'torsion') then
         do i = 1, Nat - 2
            call soc_matrix_elements(rs(:, i), rs(:, i + 1), rs(:, i + 2), 2*i, 2*(i + 2))
         end do
      end if

      !...long-range hopping (and radial SOC) on every bond (i, i+j), j <= hop_range
      do j = 1, jmax
         do i = 1, Nat - j
            m = i + j; ii = 2*i; jj = 2*m
            d = rs(:, m) - rs(:, i); l = norm3(d); dhat = d/l
            fac = exp(-(l - l1)/hop_decay)

            H(ii - 1, jj - 1) = H(ii - 1, jj - 1) + t_hop*fac
            H(ii, jj) = H(ii, jj) + t_hop*fac

            if (trim(soc_model) .eq. 'radial') then
               phibar = 0.5d0*(phi(i) + phi(m))
               rhat = [cos(phibar), sin(phibar), 0.d0]
               call cross_product(rhat, dhat, nu)
               lj = lamb*fac
               H(ii - 1, jj - 1) = H(ii - 1, jj - 1) + im*lj*nu(3)
               H(ii - 1, jj) = H(ii - 1, jj) + im*lj*(nu(1) - im*nu(2))
               H(ii, jj - 1) = H(ii, jj - 1) + im*lj*(nu(1) + im*nu(2))
               H(ii, jj) = H(ii, jj) - im*lj*nu(3)
            end if
         end do
      end do

      call make_hermitian()

      !...geometry (loadable by hamiltonian = xyz) and the hopping table
      open (newunit=u, file='peptide_geometry.xyz', status='replace')
      write (u, '(i0)') Nat
      write (u, '(a)') 'QUICHOTTE peptide sites [Angstrom]'
      do i = 1, Nat
         write (u, '(a,3f14.8)') 'C ', rs(:, i)
      end do
      close (u)
      write (*, '(a,f7.2,a,f6.3,a,f6.3,a,f5.1,a,a)') ' Peptide helix: twist', pep_twist, ' deg, rise', &
         pep_rise, ' A, radius', hel_radius, ' A, hand', hand, ', SOC ', trim(soc_model)
      write (*, '(a,f6.3,a)') '    j    l_j [A]    t_j/t_hop     (hop_decay =', hop_decay, ' A)'
      do j = 1, jmax
         l = norm3(rs(:, 1 + j) - rs(:, 1))
         write (*, '(i5,f11.4,es14.4)') j, l, exp(-(l - l1)/hop_decay)
      end do
      call report_time_reversal()
      deallocate (rs, phi)
   end subroutine build_peptide_hamiltonian

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
   ! DFTB H0 + S  (input: hamiltonian = dftb_hs)
   !
   ! Builds the orthogonal spin-orbit Hamiltonian directly from the DFTB+
   ! matrices of a run with WriteHS = Yes (and Scc = No for a non-interacting H0):
   !
   !   H = (S^-1/2 x 1) [ H0 x 1 + H_SOC ] (S^-1/2 x 1)        (x 1: spin)
   !
   !   H0, S       hamsqr1.dat, oversqr.dat (Hartree -> eV), N x N real
   !   H_SOC       atomic xi L.S on each p shell, L.S = sum_k L_k sigma_k / 2,
   !               (L_k)_ab = -i eps_kab in the real (x,y,z) basis;
   !               soc_form = dual   (default, as DFTB+):  (S Lambda + Lambda S)/2
   !               soc_form = onsite :  Lambda
   !   S^-1/2      symmetric (Lowdin) orthogonalisation
   !
   ! Orbital order per atom as in DFTB+: s, then p_y, p_z, p_x. Elements and
   ! shells from xyz_file and max_angular_momentum; xi per element from soc_xi.
   ! Exact to machine precision and exactly time-reversal symmetric.
   !===========================================================================
   subroutine read_xyz_elements(path, n_at, elem)
      !...element symbols of an .xyz file
      character(len=*), intent(in) :: path
      integer, intent(out) :: n_at
      character(len=2), allocatable, intent(out) :: elem(:)
      integer :: u, ios, a
      character(len=256) :: line

      open (newunit=u, file=trim(path), status='old', action='read', iostat=ios)
      if (ios .ne. 0) then
         print *, 'ERROR: cannot open xyz_file ', trim(path)
         stop 1
      end if
      read (u, *, iostat=ios) n_at
      if (ios .ne. 0 .or. n_at .lt. 1) then
         print *, 'ERROR: first line of ', trim(path), ' must be the number of atoms'
         stop 1
      end if
      read (u, '(a)') line
      allocate (elem(n_at))
      do a = 1, n_at
         read (u, *, iostat=ios) elem(a)
         if (ios .ne. 0) then
            print *, 'ERROR: ', trim(path), ' has fewer than', n_at, ' atom lines'
            stop 1
         end if
      end do
      close (u)
   end subroutine read_xyz_elements

   subroutine read_dftb_square(path, n, A)
      !...a DFTB+ square matrix file (hamsqr1.dat / oversqr.dat): '#' header lines, a line
      !   'REAL NALLORB NKPOINT' values, the k-point index, then n rows of n numbers
      character(len=*), intent(in) :: path
      integer, intent(in) :: n
      real(dp), intent(out) :: A(n, n)
      integer :: u, ios, i, nall, nk
      character(len=8) :: realflag
      character(len=512) :: line

      open (newunit=u, file=trim(path), status='old', action='read', iostat=ios)
      if (ios .ne. 0) then
         print *, 'ERROR: cannot open ', trim(path)
         stop 1
      end if
      !...header: skip comments, read 'T  NALLORB  NKPOINT', skip comments, k-point index, comments
      call next_data_line(u, line)
      read (line, *, iostat=ios) realflag, nall, nk
      if (ios .ne. 0 .or. nall .ne. n) then
         print *, 'ERROR: ', trim(path), ' has', nall, ' orbitals; the geometry and shells give', n
         stop 1
      end if
      if (nk .ne. 1) then
         print *, 'ERROR: ', trim(path), ': only one k-point (a molecule) is supported'
         stop 1
      end if
      call next_data_line(u, line)              !...k-point index
      call next_data_line(u, line)              !...first matrix row (after '# MATRIX')
      backspace (u)
      do i = 1, n
         read (u, *, iostat=ios) A(i, :)
         if (ios .ne. 0) then
            print *, 'ERROR: ', trim(path), ': cannot read row', i
            stop 1
         end if
      end do
      close (u)
   contains
      subroutine next_data_line(uu, l)
         integer, intent(in) :: uu
         character(len=*), intent(out) :: l
         integer :: io
         character(len=len(l)) :: t
         do
            read (uu, '(a)', iostat=io) l
            if (io .ne. 0) then
               print *, 'ERROR: unexpected end of ', trim(path)
               stop 1
            end if
            t = adjustl(l)
            if (len_trim(t) .gt. 0 .and. t(1:1) .ne. '#') exit
         end do
      end subroutine next_data_line
   end subroutine read_dftb_square

   subroutine build_hamiltonian_from_dftb_hs()
      real(dp), parameter :: hartree = 27.211386245988d0   !...eV
      real(dp), allocatable :: S(:, :), H0(:, :), Smh(:, :), ev(:), work(:)
      complex(dp), allocatable :: Lam(:, :), Sx2(:, :), Sm2(:, :), Hao(:, :)
      complex(dp) :: Lk(3, 3, 3), sig(2, 2, 3)
      integer :: N, a, e, i, j, k, p, q, s1, s2, lwork, info, i0
      integer :: pidx(3)
      real(dp) :: xi

      N = N_sites
      allocate (S(N, N), H0(N, N))
      call read_dftb_square(overlap_file, N, S)
      call read_dftb_square(ham_file, N, H0)
      S = 0.5d0*(S + transpose(S)); H0 = 0.5d0*(H0 + transpose(H0))*hartree

      !...S^-1/2 = U diag(1/sqrt(s)) U^T
      allocate (Smh(N, N), ev(N), work(1))
      Smh = S
      call dsyev('V', 'U', N, Smh, N, ev, work, -1, info)
      lwork = int(work(1)); deallocate (work); allocate (work(lwork))
      call dsyev('V', 'U', N, Smh, N, ev, work, lwork, info)
      if (info .ne. 0 .or. minval(ev) .le. 1.d-8) then
         print *, 'ERROR: overlap matrix is not positive definite (min eigenvalue', minval(ev), ')'
         stop 1
      end if
      do i = 1, N
         work(1:N) = Smh(:, i)/sqrt(sqrt(ev(i)))
         Smh(:, i) = work(1:N)
      end do
      Smh = matmul(Smh, transpose(Smh))
      write (*, '(a,es10.2)') ' DFTB H0/S: smallest overlap eigenvalue = ', minval(ev)

      !...angular momentum (L_k)_ab = -i eps_kab, a,b in (x,y,z); Pauli matrices
      Lk = (0.d0, 0.d0)
      Lk(1, 2, 3) = -im; Lk(1, 3, 2) = im     !...L_x
      Lk(2, 3, 1) = -im; Lk(2, 1, 3) = im     !...L_y
      Lk(3, 1, 2) = -im; Lk(3, 2, 1) = im     !...L_z
      sig = (0.d0, 0.d0)
      sig(1, 2, 1) = 1.d0; sig(2, 1, 1) = 1.d0
      sig(1, 2, 2) = -im; sig(2, 1, 2) = im
      sig(1, 1, 3) = 1.d0; sig(2, 2, 3) = -1.d0
      pidx = [2, 3, 1]          !...DFTB+ order p_y, p_z, p_x -> (x,y,z) index 2, 3, 1

      !...Lambda = sum_atoms xi L.S on the p shell (AO basis, spin fastest)
      allocate (Lam(Norb, Norb)); Lam = (0.d0, 0.d0)
      do a = 1, N_atoms
         if (atom_norb(a) .lt. 4) cycle
         e = findloc(elem_name, atom_elem(a), dim=1)
         xi = elem_xi(e)
         if (xi .eq. 0.d0) cycle
         i0 = atom_first(a)              !...s orbital; p_y, p_z, p_x follow
         do p = 1, 3
            do q = 1, 3
               do s1 = 1, 2
                  do s2 = 1, 2
                     i = 2*(i0 + p - 1) + s1; j = 2*(i0 + q - 1) + s2
                     do k = 1, 3
                        Lam(i, j) = Lam(i, j) + xi*Lk(k, pidx(p), pidx(q))*sig(s1, s2, k)/2.d0
                     end do
                  end do
               end do
            end do
         end do
      end do

      !...spin doubling, SOC form, Lowdin transform
      allocate (Sx2(Norb, Norb), Sm2(Norb, Norb), Hao(Norb, Norb))
      Sx2 = (0.d0, 0.d0); Sm2 = (0.d0, 0.d0); Hao = (0.d0, 0.d0)
      do j = 1, N
         do i = 1, N
            do s1 = 1, 2
               Sx2(2*(i - 1) + s1, 2*(j - 1) + s1) = S(i, j)
               Sm2(2*(i - 1) + s1, 2*(j - 1) + s1) = Smh(i, j)
               Hao(2*(i - 1) + s1, 2*(j - 1) + s1) = H0(i, j)
            end do
         end do
      end do
      select case (trim(soc_form))
      case ('dual')
         Hao = Hao + 0.5d0*(matmul(Sx2, Lam) + matmul(Lam, Sx2))
      case ('onsite')
         Hao = Hao + Lam
      case default
         print *, 'ERROR: soc_form must be dual or onsite, not ', trim(soc_form)
         stop 1
      end select
      H = matmul(Sm2, matmul(Hao, Sm2))
      H = 0.5d0*(H + conjg(transpose(H)))

      write (*, '(a,i0,a,i0,a,a)') ' DFTB H0/S Hamiltonian: ', N_atoms, ' atoms, ', N, ' orbitals, SOC form ', &
         trim(soc_form)
      do e = 1, size(elem_name)
         write (*, '(a,a2,a,i0,a,es12.4,a)') '   element ', elem_name(e), ': l_max = ', elem_lmax(e), &
            ', xi(p) = ', elem_xi(e), ' eV'
      end do
      call report_time_reversal()
      deallocate (S, H0, Smh, ev, work, Lam, Sx2, Sm2, Hao)
   end subroutine build_hamiltonian_from_dftb_hs

   subroutine write_hamiltonian_matrix(path)
      !...H in full precision: header, then one line 'i j Re Im' per non-zero element
      character(len=*), intent(in) :: path
      integer :: u, i, j

      open (newunit=u, file=trim(path), status='replace')
      write (u, '(a)') '# QUICHOTTE Hamiltonian [eV]; state index 2*site-1 = up, 2*site = down'
      write (u, '(a,i0)') '# Norb = ', Norb
      write (u, '(a)') '#   i   j   Re H_ij   Im H_ij'
      do j = 1, Norb
         do i = 1, Norb
            if (H(i, j) .ne. (0.d0, 0.d0)) write (u, '(2i6,2es25.16)') i, j, real(H(i, j)), aimag(H(i, j))
         end do
      end do
      close (u)
   end subroutine write_hamiltonian_matrix


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
