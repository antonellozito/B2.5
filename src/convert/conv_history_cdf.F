      subroutine conv_history_cdf (isymm, tria_mode, switch)
*     ..converts the history files of a structured run found in the
*       current folder (b2time.nc, b2batch.nc, b2tallies.nc,
*       b2movies.nc, eirenemovies.nc) to the layout written by the
*       unstructured code; the converted files are named *_us.nc;
*       variables that the run does not have are written as zero, and
*       they are listed together with the variables of the run that
*       the converted files do not have; mapping and geometry are
*       those of b2fgmtry_us, read as the unstructured code reads them
      use b2mod_types
      use b2us_map
      use b2us_geo
      use b2us_io, only : read_b2fgmtry_header, read_b2fgmtry
      use b2mod_switches
      use b2mod_version, only : cfverr
      implicit none
      integer, intent(in) :: isymm, tria_mode
      type (switches), intent(in) :: switch
#ifndef NO_CDF
      type (mapping) :: mpu
      type (geometry) :: gmu
      type (switches) :: sw
      character(10) :: ver
      character(120) :: lbl
      integer :: iun, isy
      logical :: ex
#ifdef LEGACYCOMP
      integer newunit
      external newunit
#endif
      external ipgeti

      call conv_tallies_cdf
      if (tria_mode.eq.0) then
        call conv_eirmovies_cdf (isymm)
      else
        write (*,*) 'eirenemovies.nc not converted: the triangles '//
     .    'were rebuilt (b2us_tria_mode > 0)'
      end if
      inquire (file='b2fgmtry_us', exist=ex)
      if (.not.ex) then
        write (*,*) 'b2time.nc, b2batch.nc, b2movies.nc not '//
     .    'converted: no b2fgmtry_us'
        return
      end if
#ifndef LEGACYCOMP
      open (newunit=iun, file='b2fgmtry_us', status='old',
     .  form='formatted')
#else
      iun = newunit()
      open (unit=iun, file='b2fgmtry_us', status='old',
     .  form='formatted')
#endif
      call cfverr (iun, ver)
      call read_b2fgmtry_header (iun, ver, mpu, lbl, isy)
      call alloc_mapping (mpu)
      call alloc_geometry (gmu, mpu%nCv, mpu%nFc, mpu%nVx, mpu%nFs,
     .  mpu%nFt, mpu%nVmxCv)
      call read_b2fgmtry (iun, mpu, gmu)
      close (iun)
      sw = switch
      sw%use_9pt_stencil = 1
      call ipgeti ('b2mndr_use_9pt_stencil', sw%use_9pt_stencil)
      call init_geometry (mpu, gmu, sw)
      call conv_time_cdf (mpu, gmu, .false.)
      call conv_time_cdf (mpu, gmu, .true.)
      call conv_movies_cdf (mpu)
      call dealloc_geometry (gmu)
      call dealloc_mapping (mpu)
#endif
      return
      end subroutine conv_history_cdf

#ifndef NO_CDF
******************************************************************************

      subroutine conv_time_cdf (mpg, gm, batch)
*     ..b2time.nc (batch false) or b2batch.nc (batch true): the file is
*       defined by b2crtimecdf; traces are copied; target profiles are
*       mapped face by face onto the target face lists and the target
*       fluxes get the outward sign of the face, as written by b2mwti;
*       the midplane profiles are copied when the structured run used
*       the midplane column chosen by b2us (b2us.midplane.parameters);
*       the cell and face lists of the profiles are computed as in
*       read_b2mod_user_namelist and b2mwti
      use b2mod_types
      use b2us_map
      use b2us_geo
      use b2mod_cdf, only : b2crtimecdf
      use b2mod_user_namelist, only : nimp, nomp
      use b2mod_b2cmpa, only : zn, am, is_neutral
      use b2mod_geometry, only : geometryId, GEOMETRY_LIMITER,
     .  GEOMETRY_SN, GEOMETRY_CDN, GEOMETRY_DDN_BOTTOM, GEOMETRY_DDN_TOP
      implicit none
#     include <netcdf.inc>
      type (mapping), intent(in) :: mpg
      type (geometry), intent(in) :: gm
      logical, intent(in) :: batch
*   ..profiles: -1 inner midplane, 0 outer midplane, 1-4 targets
*     (bottom left, top left, top right, bottom right)
      integer, parameter :: mxtgt = 4
      character*256 :: fin, fout
      character*(NF_MAX_NAME) :: name, dname
      character*2 :: tsfx(mxtgt), usfx(mxtgt)
      character*3 :: tdim(mxtgt)
      character*80 :: line
      integer :: nci, nco, iret, nxd, nyd, ns, natm, nmol, nrec, nvars,
     .  ivar, vin, xtype, ndims, natts, irec, n, nin, ninp, nprf, i, j,
     .  k, nk, ntgt, nc, igeo, ix, iy, it, ip, kind, iFc, iCv1, iCv2,
     .  iun, toff, nmx, ismain, ismain0, is, jx_in(-1:0), jx_us(-1:0),
     .  icsep(-1:0),
     .  dimids(NF_MAX_VAR_DIMS), nprof(-1:mxtgt), iy0(-1:mxtgt)
      integer, allocatable :: fc2iy(:), plist(:,:), iyp(:,:)
      real (kind=R8), allocatable :: bi(:), bo(:), sgn(:,:), ds(:),
     .  dsT(:), dsP(:)
      real (kind=R8) :: dval(1), off, fsgn
      logical :: ex, dn, facesgn, found, prof_ok(-1:mxtgt)
#ifdef LEGACYCOMP
      integer newunit
      external newunit
#endif
      external check_cdf_status, ipgeti
      data tsfx / 'l ', 'tl', 'tr', 'r ' /
      data usfx / 'L ', 'TL', 'TR', 'R ' /
      data tdim / 'bl ', 'tl ', 'tr ', 'br ' /

      if (batch) then
        fin = 'b2batch.nc'
        fout = 'b2batch_us.nc'
      else
        fin = 'b2time.nc'
        fout = 'b2time_us.nc'
      end if
      inquire (file=fin, exist=ex)
      if (.not.ex) return
      iret = nf_open (fin, NF_NOWRITE, nci)
      call check_cdf_status (iret)
      if (.not.batch) then
        call cdfh_dimlen (nci, 'nx', nxd)
        call cdfh_dimlen (nci, 'ny', nyd)
        if (nxd.ne.mpg%nx+2 .or. nyd.ne.mpg%ny+2) then
          write (*,*) trim(fin)//' is not on the structured grid of '//
     .      'b2fgmtry: not converted'
          iret = nf_close (nci)
          return
        end if
        call cdfh_dimlen (nci, 'time', nrec)
      else
        call cdfh_dimlen (nci, 'batch', nrec)
      end if
      write (*,*) 'Converting '//trim(fin)
      call cdfh_dimlen (nci, 'ns', ns)
      call cdfh_dimlen (nci, 'natm', natm)
      call cdfh_dimlen (nci, 'nmol', nmol)
      natm = max(natm, 0)
      nmol = max(nmol, 0)

*   ..targets and nc as in b2mwti
      igeo = geometryId (mpg, gm, 1)
      dn = igeo.eq.GEOMETRY_CDN .or. igeo.eq.GEOMETRY_DDN_BOTTOM .or.
     .     igeo.eq.GEOMETRY_DDN_TOP
      ntgt = 0
      if (igeo.eq.GEOMETRY_LIMITER .or. igeo.eq.GEOMETRY_SN .or. dn)
     .  ntgt = mpg%nDiv
      nc = 1
      if (dn) nc = 2
      nprof = 0
      iy0 = 0
      prof_ok = .false.
      nmx = mpg%ny+2
      do it = 1, ntgt
        nmx = max(nmx, mpg%divFcP(it,2))
      end do
      allocate (plist(nmx,-1:mxtgt), iyp(nmx,-1:mxtgt), sgn(nmx,mxtgt))
      plist = 0
      iyp = 0
      sgn = 1.0_R8
      nimp = 0
      nomp = 0

      if (.not.batch) then
*     ..midplane columns: chosen by b2us, used by the structured run
#ifndef LEGACYCOMP
        open (newunit=iun, file='b2us.midplane.parameters')
#else
        iun = newunit()
        open (unit=iun, file='b2us.midplane.parameters')
#endif
        read (iun,'(a)') line
        read (line(index(line,'=')+1:),*) jx_us(0)
        read (iun,'(a)') line
        read (iun,'(a)') line
        read (iun,'(a)') line
        read (line(index(line,'=')+1:),*) jx_us(-1)
        close (iun)
        jx_in = -1
        call cdfh_scalar (nci, 'jxi', jx_in(-1))
        call cdfh_scalar (nci, 'jxa', jx_in(0))
        do ip = -1, 0
*       ..cells of the column ordered by increasing R (cv_intersections)
          n = 0
          do iy = -1, mpg%ny
            if (mpg%imapCv(jx_us(ip),iy).gt.0) then
              n = n+1
              plist(n,ip) = mpg%imapCv(jx_us(ip),iy)
              iyp(n,ip) = iy
            end if
          end do
          if (gm%cvX(plist(1,ip)).gt.gm%cvX(plist(n,ip))) then
            plist(1:n,ip) = plist(n:1:-1,ip)
            iyp(1:n,ip) = iyp(n:1:-1,ip)
          end if
          nprof(ip) = n
*       ..separatrix position in the list (read_b2mod_user_namelist)
          icsep(ip) = 0
          j = 1
          do while (j.le.n .and. icsep(ip).eq.0)
            do i = mpg%cvFcP(plist(j,ip),1),
     .             mpg%cvFcP(plist(j,ip),1)+mpg%cvFcP(plist(j,ip),2)-1
              iFc = mpg%cvFc(i)
              if ((mpg%fcFs(iFc).eq.mpg%iFssep .and.
     .             mpg%iFssep.ne.0) .or.
     .            (mpg%fcFs(iFc).eq.mpg%iFssep2 .and.
     .             mpg%iFssep2.ne.0)) icsep(ip) = j+1
            end do
            j = j+1
          end do
          if (ip.eq.-1) then
            call cdfh_dimlen (nci, 'nyi', nin)
          else
            call cdfh_dimlen (nci, 'nya', nin)
          end if
          call cdfh_first_row (nin, mpg%ny, iy0(ip))
          prof_ok(ip) = jx_in(ip).eq.jx_us(ip) .and. iy0(ip).gt.-2
          if (.not.prof_ok(ip)) write (*,*) trim(fin)//
     .      ': midplane profiles not converted (column ',jx_us(ip),
     .      ' of b2us, column ',jx_in(ip),' of the run)'
        end do
        nimp = nprof(-1)
        nomp = nprof(0)
*     ..target face lists, outward sign and row of each face
        allocate (fc2iy(mpg%nFc))
        fc2iy = -99
        do iy = -1, mpg%ny
          do ix = -1, mpg%nx
            if (mpg%imapFcx(ix,iy).gt.0) fc2iy(mpg%imapFcx(ix,iy)) = iy
          end do
        end do
        do it = 1, ntgt
          ip = it
          if (it.eq.ntgt .and. ntgt.ge.2) ip = 4
          n = mpg%divFcP(it,2)
          nprof(ip) = n
          do j = 1, n
            i = mpg%divFcP(it,1)+j-1
            iFc = mpg%divFc(i)
            plist(j,ip) = iFc
            iyp(j,ip) = fc2iy(iFc)
            sgn(j,ip) = mpg%divFcOr(i)
            if (sgn(j,ip).eq.0.0_R8) then
*           ..as bc_to_struct: +1 if the plasma cell is fcCv(:,1)
              sgn(j,ip) = -1.0_R8
              if (mpg%fcCv(iFc,1).le.mpg%nCi) sgn(j,ip) = 1.0_R8
            end if
          end do
          call cdfh_dimlen (nci, 'ny'//trim(tdim(ip)), nin)
          call cdfh_first_row (nin, mpg%ny, iy0(ip))
          prof_ok(ip) = iy0(ip).gt.-2 .and. minval(iyp(1:n,ip)).ge.0
          if (.not.prof_ok(ip)) write (*,*) trim(fin)//': profiles '//
     .      'of target '//trim(tsfx(ip))//' not converted'
        end do
        deallocate (fc2iy)
      end if

*   ..main plasma and main neutral species as in the driver
      ismain = 1
      call ipgeti ('b2mndr_ismain', ismain)
      ismain0 = 0
      found = .false.
      is = ismain
      do while (is.ge.0 .and. .not.found)
        if (is_neutral(is) .and. nint(zn(is)).eq.nint(zn(ismain)) .and.
     .      nint(am(is)).eq.nint(am(ismain))) then
          ismain0 = is
          found = .true.
        end if
        is = is-1
      end do
      if (.not.found .and. ismain.ne.1) ismain0 = ismain
      call ipgeti ('b2mwti_ismain0', ismain0)
      call b2crtimecdf (fout, mpg, ns, ismain, ismain0, natm, nmol,
     .  ntgt, nc, 0, nco, batch, iret)
      iret = nf_open (fout, NF_WRITE, nco)
      call check_cdf_status (iret)

*   ..fixed lists and abscissae
      if (.not.batch) then
        allocate (bo(nmx))
        do ip = -1, 0
          n = nprof(ip)
          if (n.eq.0) cycle
          bo(1) = icsep(ip)
          if (ip.eq.-1) then
            call cdfh_put (nco, 'icsepimp', bo)
          else
            call cdfh_put (nco, 'icsepomp', bo)
          end if
          bo(1:n) = plist(1:n,ip)
          call cdfh_put (nco, 'cvlist'//merge('i','a',ip.eq.-1), bo)
          bo(1:n) = 0.0_R8
          iret = nf_inq_varid (nci, 'ds'//merge('i','a',ip.eq.-1), vin)
          if (prof_ok(ip) .and. iret.eq.NF_NOERR) then
            call cdfh_recsize (nci, vin, nin)
            allocate (bi(nin))
            call cdfh_rec (nci, vin, 1, bi, nin, .false.)
            do j = 1, n
              bo(j) = bi(iyp(j,ip)-iy0(ip)+1)
            end do
            deallocate (bi)
          end if
          call cdfh_put (nco, 'ds'//merge('i','a',ip.eq.-1), bo)
        end do
        toff = 1
        call ipgeti ('b2mwti_target_offset', toff)
        allocate (ds(nmx), dsT(nmx), dsP(nmx))
        do it = 1, ntgt
          ip = it
          if (it.eq.ntgt .and. ntgt.ge.2) ip = 4
          n = nprof(ip)
          bo(1:n) = plist(1:n,ip)
          call cdfh_put (nco, 'fclist'//trim(tsfx(ip)), bo)
          do j = 1, n
            iCv1 = mpg%fcCv(plist(j,ip),1)
            iCv2 = mpg%fcCv(plist(j,ip),2)
            if ((iCv1.le.mpg%nCi .and. toff.eq.1) .or.
     .          (iCv1.gt.mpg%nCi .and. toff.eq.0)) then
              bo(j) = iCv1
            else
              bo(j) = iCv2
            end if
          end do
          call cdfh_put (nco, 'cvlist'//trim(tsfx(ip)), bo)
          do j = 1, n
            iCv1 = mpg%fcCv(plist(j,ip),1)
            iCv2 = mpg%fcCv(plist(j,ip),2)
            bo(j) = 0.0_R8
            if (iCv1.le.mpg%nCi) bo(j) = iCv1
            if (iCv2.le.mpg%nCi) bo(j) = iCv2
          end do
          call cdfh_put (nco, 'cnlist'//trim(tsfx(ip)), bo)
*       ..as output_ds_fc of b2mwti
          ds(1) = 0.5_R8*gm%fcHt(plist(1,ip))
          do j = 2, n
            ds(j) = ds(j-1)+0.5_R8*(gm%fcHt(plist(j-1,ip))+
     .                              gm%fcHt(plist(j,ip)))
          end do
          i = mpg%ifdiv(it)-1
          if (i.ne.0) then
            off = 0.5_R8*(ds(i)+ds(i+1))
            ds(1:n) = ds(1:n)-off
          end if
          do j = 1, n
            dsT(j) = gm%fcS(plist(j,ip))
            dsP(j) = gm%fcS(plist(j,ip))*abs(gm%fcQalf(plist(j,ip),0))
          end do
          call cdfh_put (nco, 'ds'//trim(tsfx(ip)), ds)
          call cdfh_put (nco, 'ds'//trim(usfx(ip))//'T', dsT)
          call cdfh_put (nco, 'ds'//trim(usfx(ip))//'P', dsP)
        end do
        deallocate (bo, ds, dsT, dsP)
        dval(1) = nrec
        call cdfh_put (nco, 'ntstep', dval)
      end if

*   ..variables per record (time or batch), scalars of b2batch.nc
      iret = nf_inq_nvars (nco, nvars)
      call check_cdf_status (iret)
      do ivar = 1, nvars
        iret = nf_inq_var (nco, ivar, name, xtype, ndims, dimids, natts)
        call check_cdf_status (iret)
        if (ndims.eq.0) then
          iret = nf_inq_varid (nci, trim(name), vin)
          if (batch .and. iret.eq.NF_NOERR) then
            call cdfh_rec (nci, vin, 1, dval, n, .false.)
            call cdfh_rec (nco, ivar, 1, dval, n, .true.)
          end if
          cycle
        end if
        iret = nf_inq_dimname (nco, dimids(ndims), dname)
        call check_cdf_status (iret)
        if (trim(dname).ne.'time' .and. trim(dname).ne.'batch') cycle
*     ..kind: -1 to 4 profile (see above), 98 copy, 99 zero
        kind = 98
        iret = nf_inq_dimname (nco, dimids(1), dname)
        call check_cdf_status (iret)
        select case (trim(dname))
        case ('nyi')
          kind = -1
        case ('nya')
          kind = 0
        case ('nybl')
          kind = 1
        case ('nytl')
          kind = 2
        case ('nytr')
          kind = 3
        case ('nybr')
          kind = 4
        case ('nCv', 'nFc')
          kind = 99
        end select
        if (kind.le.4) then
          if (.not.prof_ok(kind)) kind = 99
        end if
        call cdfh_recsize (nco, ivar, n)
        nin = n
        iret = nf_inq_varid (nci, trim(name), vin)
        if (iret.ne.NF_NOERR) then
          write (*,*) trim(fin)//': '//trim(name)//
     .      ' not in the run: written as zero'
          kind = 99
        else
          call cdfh_recsize (nci, vin, nin)
          if (kind.eq.98 .and. nin.ne.n) then
            write (*,*) trim(fin)//': '//trim(name)//
     .        ' has another size in the run: written as zero'
            kind = 99
          else if (trim(dname).eq.'nCv' .or. trim(dname).eq.'nFc') then
            write (*,*) trim(fin)//': '//trim(name)//
     .        ' is given on the cell or face list: written as zero'
          end if
        end if
*     ..outward sign of the target fluxes and of the target totals
        facesgn = .false.
        fsgn = 1.0_R8
        if (kind.ge.1 .and. kind.le.4) then
          if (name(3:4).eq.'3d') then
            select case (name(1:2))
            case ('fn', 'fe', 'fi', 'fc', 'fl', 'fo', 'ft')
              facesgn = .true.
            end select
          end if
        else if (kind.eq.98 .and. len_trim(name).eq.6) then
          select case (name(1:4))
          case ('fnix', 'feex', 'feix', 'fchx', 'fetx')
            if (name(5:6).eq.'ip' .and. ntgt.ge.1) fsgn = sgn(1,1)
            if (name(5:6).eq.'ap' .and. ntgt.ge.2) fsgn = sgn(1,4)
          end select
        end if
        allocate (bo(n), bi(max(n,nin)))
        do irec = 1, nrec
          select case (kind)
          case (98)
            call cdfh_rec (nci, vin, irec, bo, n, .false.)
            bo(1:n) = fsgn*bo(1:n)
          case (99)
            bo(1:n) = 0.0_R8
          case default
            call cdfh_rec (nci, vin, irec, bi, nin, .false.)
            nprf = nprof(kind)
            nk = n/nprf
            ninp = nin/nk
            do k = 0, nk-1
              do j = 1, nprf
                bo(k*nprf+j) = bi(k*ninp+iyp(j,kind)-iy0(kind)+1)
                if (facesgn) bo(k*nprf+j) = sgn(j,kind)*bo(k*nprf+j)
              end do
            end do
          end select
          call cdfh_rec (nco, ivar, irec, bo, n, .true.)
        end do
        deallocate (bo, bi)
      end do
      deallocate (plist, iyp, sgn)
      call cdfh_dropped (nci, nco, fin)
      iret = nf_close (nco)
      call check_cdf_status (iret)
      iret = nf_close (nci)
      call check_cdf_status (iret)
      return
      end subroutine conv_time_cdf

******************************************************************************

      subroutine cdfh_first_row (nin, ny, iy0)
*     ..first row of a structured profile of nin values: -1 with the
*       guard cells (ny+2 values), 0 without (ny values), -2 otherwise
      implicit none
      integer, intent(in) :: nin, ny
      integer, intent(out) :: iy0

      iy0 = -2
      if (nin.eq.ny+2) iy0 = -1
      if (nin.eq.ny) iy0 = 0
      return
      end subroutine cdfh_first_row

******************************************************************************

      subroutine cdfh_scalar (ncid, name, ival)
*     ..integer value of a scalar variable, unchanged if absent
      use b2mod_types
      implicit none
#     include <netcdf.inc>
      integer, intent(in) :: ncid
      character*(*), intent(in) :: name
      integer, intent(inout) :: ival
      integer :: iret, varid, n
      real (kind=R8) :: val(1)

      iret = nf_inq_varid (ncid, name, varid)
      if (iret.ne.NF_NOERR) return
      call cdfh_rec (ncid, varid, 1, val, n, .false.)
      ival = nint(val(1))
      return
      end subroutine cdfh_scalar

******************************************************************************

      subroutine cdfh_put (ncid, name, val)
*     ..writes a fixed variable whole
      use b2mod_types
      implicit none
#     include <netcdf.inc>
      integer, intent(in) :: ncid
      character*(*), intent(in) :: name
      real (kind=R8), intent(in) :: val(*)
      integer :: iret, varid, n
      real (kind=R8), allocatable :: tmp(:)
      external check_cdf_status

      iret = nf_inq_varid (ncid, name, varid)
      if (iret.ne.NF_NOERR) return
      call cdfh_recsize (ncid, varid, n)
      allocate (tmp(n))
      tmp(1:n) = val(1:n)
      call cdfh_rec (ncid, varid, 1, tmp, n, .true.)
      deallocate (tmp)
      return
      end subroutine cdfh_put

******************************************************************************

      subroutine conv_tallies_cdf
*     ..b2tallies.nc: the x- and y-directed face regions of the
*       structured file become one list of face regions (x regions
*       first, as in b2us.regions.parameters)
      use b2mod_types
      use b2mod_tallies, only : define_cdf_b2mod_tallies
      implicit none
#     include <netcdf.inc>
      character*256 :: fin, fout
      character*(NF_MAX_NAME) :: name
      integer :: nci, nco, iret, nvreg, nxreg, nyreg, nfreg, ns, nrec,
     .  nvars, ivar, vin, vx, vy, xtype, ndims, natts, irec, n, nx, ny,
     .  nk, k, mode, ll, dimids(NF_MAX_VAR_DIMS), nnreg(0:1)
      real (kind=R8), allocatable :: bo(:), bx(:), by(:)
      character*1, allocatable :: txt(:)
      logical :: ex
      external check_cdf_status

      fin = 'b2tallies.nc'
      inquire (file=fin, exist=ex)
      if (.not.ex) return
      iret = nf_open (fin, NF_NOWRITE, nci)
      call check_cdf_status (iret)
      call cdfh_dimlen (nci, 'xreg', nxreg)
      call cdfh_dimlen (nci, 'yreg', nyreg)
      if (nxreg.lt.1 .or. nyreg.lt.1) then
        write (*,*) 'b2tallies.nc is not structured: not converted'
        iret = nf_close (nci)
        return
      end if
      write (*,*) 'Converting b2tallies.nc'
      call cdfh_dimlen (nci, 'vreg', nvreg)
      call cdfh_dimlen (nci, 'ns', ns)
      call cdfh_dimlen (nci, 'time', nrec)
      nnreg(0) = nvreg-1
      nnreg(1) = nxreg+nyreg-2
      nfreg = nnreg(1)+1
      fout = 'b2tallies_us.nc'
      call define_cdf_b2mod_tallies (fout, nco, nnreg, ns)
      iret = nf_open (fout, NF_WRITE, nco)
      call check_cdf_status (iret)
      n = max(nvreg, nfreg)*ns
      allocate (bo(n), bx(nxreg*ns), by(nyreg*ns))
      iret = nf_inq_nvars (nco, nvars)
      call check_cdf_status (iret)
      do ivar = 1, nvars
        iret = nf_inq_var (nco, ivar, name, xtype, ndims, dimids, natts)
        call check_cdf_status (iret)
        if (xtype.eq.NF_CHAR) then
*         ..species names
          iret = nf_inq_varid (nci, trim(name), vin)
          if (iret.eq.NF_NOERR) then
            call cdfh_recsize (nci, vin, ll)
            allocate (txt(ll))
            iret = nf_get_var_text (nci, vin, txt)
            call check_cdf_status (iret)
            iret = nf_put_var_text (nco, ivar, txt)
            call check_cdf_status (iret)
            deallocate (txt)
          end if
          cycle
        end if
        ll = len_trim(name)
        mode = 0
        iret = nf_inq_varid (nci, trim(name), vin)
        if (iret.eq.NF_NOERR) then
          mode = 1
        else if (ll.gt.3) then
          if (name(ll-2:ll).eq.'reg') then
            iret = nf_inq_varid (nci, name(1:ll-3)//'xreg', vx)
            if (iret.eq.NF_NOERR)
     .        iret = nf_inq_varid (nci, name(1:ll-3)//'yreg', vy)
            if (iret.eq.NF_NOERR) mode = 2
          end if
        end if
        if (mode.eq.0) write (*,*) trim(fin)//': '//trim(name)//
     .    ' not in the run: written as zero'
        call cdfh_recsize (nco, ivar, n)
        do irec = 1, nrec
          select case (mode)
          case (1)
            call cdfh_rec (nci, vin, irec, bo, n, .false.)
          case (2)
            call cdfh_rec (nci, vx, irec, bx, nx, .false.)
            call cdfh_rec (nci, vy, irec, by, ny, .false.)
            nk = n/nfreg
            do k = 0, nk-1
              bo(k*nfreg+1) = 0.0_R8
              bo(k*nfreg+2:k*nfreg+nxreg) =
     .          bx(k*nxreg+2:k*nxreg+nxreg)
              bo(k*nfreg+nxreg+1:k*nfreg+nfreg) =
     .          by(k*nyreg+2:k*nyreg+nyreg)
            end do
          case default
            bo(1:n) = 0.0_R8
          end select
          call cdfh_rec (nco, ivar, irec, bo, n, .true.)
        end do
      end do
      deallocate (bo, bx, by)
      call cdfh_dropped (nci, nco, fin)
      iret = nf_close (nco)
      call check_cdf_status (iret)
      iret = nf_close (nci)
      call check_cdf_status (iret)
      return
      end subroutine conv_tallies_cdf

******************************************************************************

      subroutine conv_movies_cdf (mpg)
*     ..b2movies.nc: the converted file has the fields that the
*       unstructured code writes with the cdfmovie_* switches of
*       b2mn.dat (create_cdf_b2movies_fields), so that a run continued
*       from it can append to it; a cell field of the run goes to the
*       cell list, an x-directed flux field to component 0 on the x
*       faces, a y-directed one to component 1 on the y faces (0 on
*       the other faces); a field that the run did not write is zero
      use b2mod_types
      use b2us_map
      use b2mod_movies
      use b2mod_switches, only : switches, set_defaults_switches
      implicit none
#     include <netcdf.inc>
      type (mapping), intent(in) :: mpg
      type (switches) :: sw
      character*256 :: fin, fout
      character*(NF_MAX_NAME) :: dname
      character*(NF_MAX_NAME), allocatable :: vname(:)
      integer :: nci, nco, iret, nxd, nyd, ns, natm, nmol, nrec,
     .  nvars, ivar, ndims, irec, n, nk, nout, vtim, vid, kin,
     .  dimids(NF_MAX_VAR_DIMS)
      integer, allocatable :: nd(:), kk(:), stype(:), idir(:), vin(:)
      real (kind=R8) :: tim(1)
      real (kind=R8), allocatable :: buf(:), out(:)
      logical :: ex, eirene
      external check_cdf_status, ipgeti

      fin = 'b2movies.nc'
      inquire (file=fin, exist=ex)
      if (.not.ex) return
      iret = nf_open (fin, NF_NOWRITE, nci)
      call check_cdf_status (iret)
      call cdfh_dimlen (nci, 'nx', nxd)
      call cdfh_dimlen (nci, 'ny', nyd)
      if (nxd.ne.mpg%nx+2 .or. nyd.ne.mpg%ny+2) then
        write (*,*) 'b2movies.nc is not on the structured grid of '//
     .    'b2fgmtry: not converted'
        iret = nf_close (nci)
        return
      end if
      write (*,*) 'Converting b2movies.nc'
      call cdfh_dimlen (nci, 'ns', ns)
      call cdfh_dimlen (nci, 'natm', natm)
      call cdfh_dimlen (nci, 'nmol', nmol)
      call cdfh_dimlen (nci, 'time', nrec)
      natm = max(natm, 0)
      nmol = max(nmol, 0)
*   ..switches of the fields, as in read_switches
      call set_defaults_switches (sw, ns)
      call ipgeti ('cdfmovie_b25_state_variables',
     .  sw%write_b25_state_variables)
      call ipgeti ('cdfmovie_b25_fluxes', sw%write_b25_fluxes)
      call ipgeti ('cdfmovie_b25_sources', sw%write_b25_sources)
      call ipgeti ('cdfmovie_b25_transport', sw%write_b25_transport)
      call ipgeti ('cdfmovie_eirene_state_variables',
     .  sw%write_eirene_state_variables)
      call ipgeti ('cdfmovie_eirene_fluxes', sw%write_eirene_fluxes)
      call ipgeti ('cdfmovie_eirene_sources', sw%write_eirene_sources)
#ifdef B25_EIRENE
      eirene = .true.
#else
      eirene = .false.
#endif
      fout = 'b2movies_us.nc'
      call create_cdf_start (nco, fout, mpg%nCv, mpg%nFc, ns, natm,
     .  nmol, 0)
      call create_cdf_b2movies_fields (nco, natm, nmol, eirene, sw)
      call create_cdf_end (nco)
*   ..fields of the converted file and the fields of the run that
*     fill them, (nx+2, ny+2[, k], time)
      iret = nf_open (fout, NF_NOWRITE, nco)
      call check_cdf_status (iret)
      iret = nf_inq_nvars (nco, nvars)
      call check_cdf_status (iret)
      allocate (vname(nvars), nd(nvars), kk(nvars), stype(nvars),
     .  idir(nvars), vin(nvars))
      do ivar = 1, nvars
        kk(ivar) = 0
        vin(ivar) = 0
        iret = nf_inq_varname (nco, ivar, vname(ivar))
        call check_cdf_status (iret)
        iret = nf_inq_varndims (nco, ivar, nd(ivar))
        call check_cdf_status (iret)
        if (nd(ivar).lt.2) cycle
        iret = nf_inq_vardimid (nco, ivar, dimids)
        call check_cdf_status (iret)
        kk(ivar) = 1
        if (nd(ivar).eq.3) then
          iret = nf_inq_dimlen (nco, dimids(2), kk(ivar))
          call check_cdf_status (iret)
        end if
        iret = nf_inq_dimname (nco, dimids(1), dname)
        call check_cdf_status (iret)
        stype(ivar) = 1
        idir(ivar) = -1
        if (trim(dname).eq.'nfc') then
          stype(ivar) = 2
          select case (trim(vname(ivar)))
          case ('fnax', 'fhex', 'fhix', 'fchx')
            idir(ivar) = 0
          case ('fnay', 'fhey', 'fhiy', 'fchy')
            idir(ivar) = 1
          end select
        end if
        iret = nf_inq_varid (nci, trim(vname(ivar)), vid)
        if (iret.ne.NF_NOERR) then
          write (*,*) trim(fin)//': '//trim(vname(ivar))//
     .      ' not in the run: written as zero'
          cycle
        end if
        iret = nf_inq_varndims (nci, vid, ndims)
        call check_cdf_status (iret)
        iret = nf_inq_vardimid (nci, vid, dimids)
        call check_cdf_status (iret)
        kin = 1
        if (ndims.eq.4) then
          iret = nf_inq_dimlen (nci, dimids(3), kin)
          call check_cdf_status (iret)
        end if
        if (ndims.ge.3 .and. kin.eq.kk(ivar) .and.
     .      (stype(ivar).eq.1 .or. idir(ivar).ge.0)) then
          vin(ivar) = vid
        else
          write (*,*) trim(fin)//': '//trim(vname(ivar))//
     .      ' has another size in the run: written as zero'
        end if
      end do
      call cdfh_dropped (nci, nco, fin)
      iret = nf_close (nco)
      call check_cdf_status (iret)
*   ..copy the frames
      nk = max(ns, natm, nmol, 1)
      allocate (buf(nxd*nyd*nk), out(max(mpg%nCv, mpg%nFc)*nk))
      iret = nf_inq_varid (nci, 'times', vtim)
      call check_cdf_status (iret)
      do irec = 1, nrec
        call cdfh_rec (nci, vtim, irec, tim, n, .false.)
        call put_cdf_start (nco, fout, tim)
        do ivar = 1, nvars
          if (kk(ivar).eq.0) cycle
          nout = mpg%nCv
          if (stype(ivar).eq.2) nout = mpg%nFc
          if (vin(ivar).eq.0) then
            out(1:nout*kk(ivar)) = 0.0_R8
          else
            call cdfh_rec (nci, vin(ivar), irec, buf, n, .false.)
            if (stype(ivar).eq.1) then
              call cdfh_map_cv (mpg, mpg%nx, mpg%ny, kk(ivar), buf,
     .          out)
            else
              call cdfh_map_fc (mpg, mpg%nx, mpg%ny, kk(ivar),
     .          idir(ivar), buf, out)
            end if
          end if
          if (nd(ivar).eq.2) then
            call put_cdf_field2 (nco, trim(vname(ivar)), out, nout)
          else
            call put_cdf_field3 (nco, trim(vname(ivar)), out, nout,
     .        kk(ivar))
          end if
        end do
        call put_cdf_end (nco)
      end do
      deallocate (vname, nd, kk, stype, idir, vin, buf, out)
      iret = nf_close (nci)
      call check_cdf_status (iret)
      return
      end subroutine conv_movies_cdf

******************************************************************************

      subroutine conv_eirmovies_cdf (isymm)
*     ..eirenemovies.nc: the structured file holds one value per
*       triangle; the unstructured code scores on the plasma cells
*       (cell value = volume-weighted average of its triangles) and on
*       the void triangles (in triangle order), and writes all slots
*       but the last one; triangles and their cells are taken from the
*       triangle files written by b2us; the converted file has the
*       fields that the unstructured code writes
*       (create_cdf_eirenemovies_fields), zero if the run has no such
*       field
      use b2mod_types
      use b2mod_movies
      implicit none
#     include <netcdf.inc>
      integer, intent(in) :: isymm
      character*256 :: fin, fout
      character*(NF_MAX_NAME), allocatable :: vname(:)
      integer :: nci, nco, iret, ntri, natm, nmol, nrec, nvars, ivar,
     .  ndims, irec, n, nk, nco_tri, nel, ncell, nslot, nout, itri, i,
     .  k, iun, idum(9), itfc(3), vtim, vid, kin,
     .  dimids(NF_MAX_VAR_DIMS)
      integer, allocatable :: kk(:), vin(:), elem(:,:), ecv(:),
     .  slot(:)
      real (kind=R8) :: tim(1), xc(3), yc(3)
      real (kind=R8), allocatable :: xco(:), yco(:), w(:), wsum(:),
     .  buf(:,:), out(:,:)
      logical :: ex
#ifdef LEGACYCOMP
      integer newunit
      external newunit
#endif
      external check_cdf_status

      fin = 'eirenemovies.nc'
      inquire (file=fin, exist=ex)
      if (.not.ex) return
      inquire (file='triangle.neighbors', exist=ex)
      if (.not.ex) then
        write (*,*) 'eirenemovies.nc not converted: no triangle files'
        return
      end if
*   ..triangles: nodes, corners, cell of each triangle
#ifndef LEGACYCOMP
      open (newunit=iun, file='triangle.coordinates')
#else
      iun = newunit()
      open (unit=iun, file='triangle.coordinates')
#endif
      read (iun,*) nco_tri
      allocate (xco(nco_tri), yco(nco_tri))
      do i = 1, nco_tri
        read (iun,*) idum(1), xco(i), yco(i)
      end do
      close (iun)
#ifndef LEGACYCOMP
      open (newunit=iun, file='triangle.elements')
#else
      iun = newunit()
      open (unit=iun, file='triangle.elements')
#endif
      read (iun,*) nel
      allocate (elem(3,nel), ecv(nel), slot(nel), w(nel))
      do itri = 1, nel
        read (iun,*) idum(1), elem(1:3,itri)
      end do
      close (iun)
#ifndef LEGACYCOMP
      open (newunit=iun, file='triangle.neighbors')
#else
      iun = newunit()
      open (unit=iun, file='triangle.neighbors')
#endif
      read (iun,*) i
      do itri = 1, nel
        read (iun,*) idum(1), idum(1:9), ecv(itri), itfc(1:3)
      end do
      close (iun)

      iret = nf_open (fin, NF_NOWRITE, nci)
      call check_cdf_status (iret)
      call cdfh_dimlen (nci, 'ntri', ntri)
      if (ntri.ne.nel) then
        write (*,*) 'eirenemovies.nc does not hold one value per '//
     .    'triangle of triangle.elements: not converted'
        iret = nf_close (nci)
        deallocate (xco, yco, elem, ecv, slot, w)
        return
      end if
      write (*,*) 'Converting eirenemovies.nc'
      call cdfh_dimlen (nci, 'natm', natm)
      call cdfh_dimlen (nci, 'nmol', nmol)
      call cdfh_dimlen (nci, 'time', nrec)
      natm = max(natm, 0)
      nmol = max(nmol, 0)

*   ..slots: plasma cells first, then the void triangles in order
      ncell = maxval(ecv)
      nslot = ncell
      do itri = 1, nel
        if (ecv(itri).gt.0) then
          slot(itri) = ecv(itri)
        else
          nslot = nslot+1
          slot(itri) = nslot
        end if
      end do
      nout = nslot-1
*   ..triangle volumes up to a constant factor (see compute_triangle_vol)
      allocate (wsum(nslot))
      wsum = 0.0_R8
      do itri = 1, nel
        xc(1:3) = xco(elem(1:3,itri))
        yc(1:3) = yco(elem(1:3,itri))
        w(itri) = 0.5_R8*abs((xc(2)-xc(1))*(yc(3)-yc(1))-
     .                       (xc(3)-xc(1))*(yc(2)-yc(1)))
        if (isymm.eq.1 .or. isymm.eq.2) then
          w(itri) = w(itri)*(xc(1)+xc(2)+xc(3))/3.0_R8
        else if (isymm.eq.3 .or. isymm.eq.4) then
          w(itri) = w(itri)*(yc(1)+yc(2)+yc(3))/3.0_R8
        end if
        wsum(slot(itri)) = wsum(slot(itri))+w(itri)
      end do

      fout = 'eirenemovies_us.nc'
      call create_cdf_start (nco, fout, 0, 0, 0, natm, nmol, nout)
      call create_cdf_eirenemovies_fields (nco, natm, nmol)
      call create_cdf_end (nco)
*   ..fields of the converted file and the fields of the run that
*     fill them, (ntri, natm or nmol, time)
      iret = nf_open (fout, NF_NOWRITE, nco)
      call check_cdf_status (iret)
      iret = nf_inq_nvars (nco, nvars)
      call check_cdf_status (iret)
      allocate (vname(nvars), kk(nvars), vin(nvars))
      do ivar = 1, nvars
        kk(ivar) = 0
        vin(ivar) = 0
        iret = nf_inq_varname (nco, ivar, vname(ivar))
        call check_cdf_status (iret)
        iret = nf_inq_varndims (nco, ivar, ndims)
        call check_cdf_status (iret)
        if (ndims.ne.3) cycle
        iret = nf_inq_vardimid (nco, ivar, dimids)
        call check_cdf_status (iret)
        iret = nf_inq_dimlen (nco, dimids(2), kk(ivar))
        call check_cdf_status (iret)
        iret = nf_inq_varid (nci, trim(vname(ivar)), vid)
        if (iret.ne.NF_NOERR) then
          write (*,*) trim(fin)//': '//trim(vname(ivar))//
     .      ' not in the run: written as zero'
          cycle
        end if
        iret = nf_inq_varndims (nci, vid, ndims)
        call check_cdf_status (iret)
        iret = nf_inq_vardimid (nci, vid, dimids)
        call check_cdf_status (iret)
        kin = 0
        if (ndims.eq.3) then
          iret = nf_inq_dimlen (nci, dimids(2), kin)
          call check_cdf_status (iret)
        end if
        if (kin.eq.kk(ivar)) then
          vin(ivar) = vid
        else
          write (*,*) trim(fin)//': '//trim(vname(ivar))//
     .      ' has another size in the run: written as zero'
        end if
      end do
      call cdfh_dropped (nci, nco, fin)
      iret = nf_close (nco)
      call check_cdf_status (iret)

      nk = max(natm, nmol, 1)
      allocate (buf(nel,nk), out(nslot,nk))
      iret = nf_inq_varid (nci, 'times', vtim)
      call check_cdf_status (iret)
      do irec = 1, nrec
        call cdfh_rec (nci, vtim, irec, tim, n, .false.)
        call put_cdf_start (nco, fout, tim)
        do ivar = 1, nvars
          if (kk(ivar).eq.0) cycle
          out(:,1:kk(ivar)) = 0.0_R8
          if (vin(ivar).ne.0) then
            call cdfh_rec (nci, vin(ivar), irec, buf, n, .false.)
            do k = 1, kk(ivar)
              do itri = 1, nel
                out(slot(itri),k) = out(slot(itri),k)+
     .            w(itri)*buf(itri,k)
              end do
              do i = 1, nslot
                if (wsum(i).gt.0.0_R8) out(i,k) = out(i,k)/wsum(i)
              end do
            end do
          end if
          call cdfh_put_tri (nco, trim(vname(ivar)), out, nslot, nout,
     .      kk(ivar))
        end do
        call put_cdf_end (nco)
      end do
      deallocate (vname, kk, vin, xco, yco, elem, ecv, slot, w, wsum,
     .  buf, out)
      iret = nf_close (nci)
      call check_cdf_status (iret)
      return
      end subroutine conv_eirmovies_cdf

******************************************************************************

      subroutine cdfh_put_tri (ncid, name, out, nslot, nout, nk)
*     ..writes the first nout slots of each column of out
      use b2mod_types
      use b2mod_movies, only : put_cdf_field_tri
      implicit none
      integer, intent(in) :: ncid, nslot, nout, nk
      character*(*), intent(in) :: name
      real (kind=R8), intent(in) :: out(nslot,nk)
      real (kind=R8) :: outp(nout,nk)

      outp(1:nout,1:nk) = out(1:nout,1:nk)
      call put_cdf_field_tri (ncid, name, outp, nout, nk)
      return
      end subroutine cdfh_put_tri

******************************************************************************

      subroutine cdfh_map_cv (mpg, nx, ny, nk, buf, out)
*     ..structured cell field (guard cells included) to the cell list
      use b2mod_types
      use b2us_map
      implicit none
      type (mapping), intent(in) :: mpg
      integer, intent(in) :: nx, ny, nk
      real (kind=R8), intent(in) :: buf(-1:nx,-1:ny,nk)
      real (kind=R8), intent(out) :: out(mpg%nCv,nk)
      integer :: ix, iy, k

      out = 0.0_R8
      do k = 1, nk
        do iy = -1, ny
          do ix = -1, nx
            if (mpg%imapCv(ix,iy).gt.0)
     .        out(mpg%imapCv(ix,iy),k) = buf(ix,iy,k)
          end do
        end do
      end do
      return
      end subroutine cdfh_map_cv

******************************************************************************

      subroutine cdfh_map_fc (mpg, nx, ny, nk, idir, buf, out)
*     ..structured x- (idir 0) or y-directed (idir 1) face field to the
*       face list, with the face loops of prep_plasma
      use b2mod_types
      use b2us_map
      implicit none
      type (mapping), intent(in) :: mpg
      integer, intent(in) :: nx, ny, nk, idir
      real (kind=R8), intent(in) :: buf(-1:nx,-1:ny,nk)
      real (kind=R8), intent(out) :: out(mpg%nFc,nk)
      integer :: ix, iy, k

      out = 0.0_R8
      do k = 1, nk
        if (idir.eq.0) then
          do iy = 0, ny-1
            do ix = 0, nx
              if (mpg%imapFcx(ix,iy).ne.0)
     .          out(mpg%imapFcx(ix,iy),k) = buf(ix,iy,k)
            end do
          end do
        else
          do ix = 0, nx-1
            do iy = 0, ny
              if (mpg%imapFcy(ix,iy).ne.0)
     .          out(mpg%imapFcy(ix,iy),k) = buf(ix,iy,k)
            end do
          end do
        end if
      end do
      return
      end subroutine cdfh_map_fc

******************************************************************************

      subroutine cdfh_dimlen (ncid, name, n)
*     ..length of a dimension, -1 if the file has no such dimension
      implicit none
#     include <netcdf.inc>
      integer, intent(in) :: ncid
      character*(*), intent(in) :: name
      integer, intent(out) :: n
      integer :: iret, dimid
      external check_cdf_status

      n = -1
      iret = nf_inq_dimid (ncid, name, dimid)
      if (iret.eq.NF_NOERR) then
        iret = nf_inq_dimlen (ncid, dimid, n)
        call check_cdf_status (iret)
      end if
      return
      end subroutine cdfh_dimlen

******************************************************************************

      subroutine cdfh_recsize (ncid, varid, n)
*     ..number of values in one record of a variable
      implicit none
#     include <netcdf.inc>
      integer, intent(in) :: ncid, varid
      integer, intent(out) :: n
      integer :: iret, ndims, unlim, i, len, dimids(NF_MAX_VAR_DIMS)
      external check_cdf_status

      iret = nf_inq_unlimdim (ncid, unlim)
      call check_cdf_status (iret)
      iret = nf_inq_varndims (ncid, varid, ndims)
      call check_cdf_status (iret)
      iret = nf_inq_vardimid (ncid, varid, dimids)
      call check_cdf_status (iret)
      n = 1
      do i = 1, ndims
        if (dimids(i).ne.unlim) then
          iret = nf_inq_dimlen (ncid, dimids(i), len)
          call check_cdf_status (iret)
          n = n*len
        end if
      end do
      return
      end subroutine cdfh_recsize

******************************************************************************

      subroutine cdfh_rec (ncid, varid, irec, buf, n, put)
*     ..reads (put false) or writes (put true) record irec of a
*       variable, or the whole variable if it has no record dimension;
*       n returns the number of values
      use b2mod_types
      implicit none
#     include <netcdf.inc>
      integer, intent(in) :: ncid, varid, irec
      real (kind=R8), intent(inout) :: buf(*)
      integer, intent(out) :: n
      logical, intent(in) :: put
      integer :: iret, ndims, unlim, i, dimids(NF_MAX_VAR_DIMS),
     .  start(NF_MAX_VAR_DIMS), count(NF_MAX_VAR_DIMS)
      external check_cdf_status

      iret = nf_inq_unlimdim (ncid, unlim)
      call check_cdf_status (iret)
      iret = nf_inq_varndims (ncid, varid, ndims)
      call check_cdf_status (iret)
      iret = nf_inq_vardimid (ncid, varid, dimids)
      call check_cdf_status (iret)
      n = 1
      do i = 1, ndims
        start(i) = 1
        if (dimids(i).eq.unlim) then
          start(i) = irec
          count(i) = 1
        else
          iret = nf_inq_dimlen (ncid, dimids(i), count(i))
          call check_cdf_status (iret)
          n = n*count(i)
        end if
      end do
      if (put) then
        iret = nf_put_vara_double (ncid, varid, start, count, buf)
      else
        iret = nf_get_vara_double (ncid, varid, start, count, buf)
      end if
      call check_cdf_status (iret)
      return
      end subroutine cdfh_rec

******************************************************************************

      subroutine cdfh_dropped (nci, nco, fin)
*     ..lists the fields of the run (variables on the record
*       dimension) that the converted file does not have; the x and y
*       region fields of b2tallies.nc are merged, not dropped
      implicit none
#     include <netcdf.inc>
      integer, intent(in) :: nci, nco
      character*(*), intent(in) :: fin
      character*(NF_MAX_NAME) :: name
      integer :: iret, nvars, ivar, ndims, unlim, ll, vid,
     .  dimids(NF_MAX_VAR_DIMS)
      external check_cdf_status

      iret = nf_inq_unlimdim (nci, unlim)
      call check_cdf_status (iret)
      iret = nf_inq_nvars (nci, nvars)
      call check_cdf_status (iret)
      do ivar = 1, nvars
        iret = nf_inq_varname (nci, ivar, name)
        call check_cdf_status (iret)
        iret = nf_inq_varndims (nci, ivar, ndims)
        call check_cdf_status (iret)
        if (ndims.lt.1) cycle
        iret = nf_inq_vardimid (nci, ivar, dimids)
        call check_cdf_status (iret)
        if (dimids(ndims).ne.unlim) cycle
        iret = nf_inq_varid (nco, trim(name), vid)
        if (iret.eq.NF_NOERR) cycle
        ll = len_trim(name)
        if (ll.gt.4) then
          if (name(ll-3:ll).eq.'xreg' .or. name(ll-3:ll).eq.'yreg')
     .      then
            iret = nf_inq_varid (nco, name(1:ll-4)//'reg', vid)
            if (iret.eq.NF_NOERR) cycle
          end if
        end if
        write (*,*) trim(fin)//': '//trim(name)//
     .    ' is not written by the unstructured code: dropped'
      end do
      return
      end subroutine cdfh_dropped
#endif
