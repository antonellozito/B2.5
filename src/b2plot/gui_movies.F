c GUI movie commands guimvinf, guimovie and guimvvec: variables of
c b2movies.nc (file 1, B2.5 cells and faces) and eirenemovies.nc (file
c 2, EIRENE tally cells) for b2plot.py, one frame or the mean of frames

c Open the movie files found by find_file, once per session
      subroutine gui_movie_open
      use b2mod_types
      use b2mod_b2plot_gui
      implicit none
#ifndef NO_CDF
#     include <netcdf.inc>
      character*256 :: fname
      integer :: ifile, iret, unlim
      logical :: ex
      if (gui_mv_opened) return
      gui_mv_opened = .true.
      do ifile = 1, 2
#ifndef B25_EIRENE
        if (ifile.eq.2) cycle
#endif
        if (ifile.eq.1) then
          fname = 'b2movies.nc'
        else
          fname = 'eirenemovies.nc'
        end if
        call find_file(fname,ex)
        if (.not.ex) cycle
        iret = nf_open(fname,NF_NOWRITE,gui_mv_ncid(ifile))
        if (iret.ne.NF_NOERR) then
          gui_mv_ncid(ifile) = -1
          write(*,*) 'gui_movie_open: cannot open ',trim(fname)
          cycle
        end if
        gui_mv_path(ifile) = fname
        iret = nf_inq_unlimdim(gui_mv_ncid(ifile),unlim)
        if (iret.eq.NF_NOERR.and.unlim.gt.0) iret =
     .   nf_inq_dimlen(gui_mv_ncid(ifile),unlim,gui_mv_nt(ifile))
      end do
#endif
      return
      end subroutine gui_movie_open

c Close the movie files and free the running sums
      subroutine gui_movie_close
      use b2mod_types
      use b2mod_b2plot_gui
      implicit none
      integer :: ifile, islot
#ifndef NO_CDF
#     include <netcdf.inc>
      integer :: iret
      do ifile = 1, 2
        if (gui_mv_ncid(ifile).ge.0) iret = nf_close(gui_mv_ncid(ifile))
      end do
#endif
      do ifile = 1, 2
        gui_mv_ncid(ifile) = -1
        gui_mv_nt(ifile) = 0
      end do
      do islot = 1, gui_mv_nslot
        if (allocated(gui_mv_c(islot)%s)) deallocate(gui_mv_c(islot)%s)
        gui_mv_c(islot)%ifile = 0
      end do
      gui_mv_opened = .false.
      return
      end subroutine gui_movie_close

c Location of variable varid of file ifile (1 cells, 2 faces, 3 EIRENE
c tally cells, 0 not usable), size, columns, column dimension, rank
      subroutine gui_movie_var(ifile,varid,vname,loc,n,ncol,dname,nd)
      use b2mod_types
      use b2us_data
      use b2mod_b2plot
      use b2mod_b2plot_gui
#ifdef B25_EIRENE
      use eirmod_ctrig
      use eirmod_cgeom
     , , only: ncltal
#endif
      implicit none
      integer, intent(in) :: ifile, varid
      character*(*), intent(out) :: vname, dname
      integer, intent(out) :: loc, n, ncol, nd
#ifndef NO_CDF
#     include <netcdf.inc>
      integer :: iret, xtype, ndims, natts, unlim, ncid
      integer :: dimids(NF_MAX_VAR_DIMS)
      character*(NF_MAX_NAME) :: d1
#endif
      loc = 0
      n = 0
      ncol = 1
      nd = 0
      vname = ' '
      dname = ' '
#ifndef NO_CDF
      ncid = gui_mv_ncid(ifile)
      if (ncid.lt.0) return
      iret = nf_inq_var(ncid,varid,vname,xtype,ndims,dimids,natts)
      if (iret.ne.NF_NOERR) return
      iret = nf_inq_unlimdim(ncid,unlim)
      if (ndims.lt.2.or.ndims.gt.3) return
      if (dimids(ndims).ne.unlim) return
      nd = ndims
      iret = nf_inq_dimname(ncid,dimids(1),d1)
      iret = nf_inq_dimlen(ncid,dimids(1),n)
      if (ndims.eq.3) then
        iret = nf_inq_dimname(ncid,dimids(2),dname)
        iret = nf_inq_dimlen(ncid,dimids(2),ncol)
      end if
      if (ifile.eq.1.and.ncol.le.size(work,2)) then
        if (d1.eq.'ncv'.and.n.eq.mpg%nCv) loc = 1
        if (d1.eq.'nfc'.and.n.eq.mpg%nFc) loc = 2
      end if
#ifdef B25_EIRENE
      if (ifile.eq.2.and.d1.eq.'ntri'.and.allocated(ncltal)) then
        if (ntrii.gt.0.and.ncol.le.size(eirene_stack,2)) loc = 3
      end if
#endif
#endif
      return
      end subroutine gui_movie_var

c Text attribute of a movie variable (blank if missing)
      subroutine gui_movie_text(ncid,varid,att,text)
      use b2mod_types
      implicit none
      integer, intent(in) :: ncid, varid
      character*(*), intent(in) :: att
      character*(*), intent(out) :: text
#ifndef NO_CDF
#     include <netcdf.inc>
      integer :: iret, alen
#endif
      text = ' '
#ifndef NO_CDF
      iret = nf_inq_attlen(ncid,varid,att,alen)
      if (iret.eq.NF_NOERR.and.alen.le.len(text))
     . iret = nf_get_att_text(ncid,varid,att,text)
#endif
      return
      end subroutine gui_movie_text

c Movie files, frame times and usable variables for the GUI (MOVIEINFO
c per file, then MOVIEVAR per variable)
      subroutine gui_send_movieinfo
      use b2mod_types
      use b2mod_b2plot_gui
      implicit none
      real (kind=R8), allocatable :: times(:)
      character*256 :: vname, dname, long, units
      integer :: ifile, ncid, nt, nvars, varid, nuse, loc, n, ncol, nd
#ifndef NO_CDF
#     include <netcdf.inc>
      integer :: iret
#endif
      external gui_movie_open, gui_movie_var, gui_movie_text

      if (.not.gui_active) return
      call gui_movie_open
      do ifile = 1, 2
        ncid = gui_mv_ncid(ifile)
        if (ncid.lt.0) cycle
        nt = gui_mv_nt(ifile)
        allocate(times(max(nt,1)))
        times = 0.0_R8
        nvars = 0
#ifndef NO_CDF
        iret = nf_inq_varid(ncid,'times',varid)
        if (iret.eq.NF_NOERR.and.nt.gt.0) iret =
     .   nf_get_vara_double(ncid,varid,(/1/),(/nt/),times)
        iret = nf_inq_nvars(ncid,nvars)
#endif
        nuse = 0
        do varid = 1, nvars
          call gui_movie_var(ifile,varid,vname,loc,n,ncol,dname,nd)
          if (loc.gt.0) nuse = nuse+1
        end do
        call gui_begin(79)
        call gui_put_i4(ifile)
        call gui_put_str(gui_mv_path(ifile))
        call gui_put_i4(nt)
        call gui_put_r8v(nt,times)
        call gui_put_i4(nuse)
        call gui_end
        do varid = 1, nvars
          call gui_movie_var(ifile,varid,vname,loc,n,ncol,dname,nd)
          if (loc.eq.0) cycle
          call gui_movie_text(ncid,varid,'long_name',long)
          call gui_movie_text(ncid,varid,'units',units)
          call gui_begin(80)
          call gui_put_i4(ifile)
          call gui_put_str(vname)
          call gui_put_str(long)
          call gui_put_str(units)
          call gui_put_i4(loc)
          call gui_put_i4(ncol)
          call gui_put_str(dname)
          call gui_end
        end do
        deallocate(times)
      end do
      return
      end subroutine gui_send_movieinfo

c Mean over the frames k1 to k2 of variable varid of file ifile, at the
c location of the file; a running sum per variable adds and removes
c only the frames that changed since the last request
      subroutine gui_movie_mean(ifile,varid,vname,k1,k2,n,ncol,nd,mean,
     . ok)
      use b2mod_types
      use b2mod_b2plot_gui
      implicit none
      integer, intent(in) :: ifile, varid, k1, k2, n, ncol, nd
      character*(*), intent(in) :: vname
      real (kind=R8), intent(out) :: mean(n,ncol)
      logical, intent(out) :: ok
      real (kind=R8), allocatable :: frame(:,:)
      integer :: islot, i, k, lo, hi, nmean, nover, nadd, nrem
      integer :: kstart(3), kcount(3)
#ifndef NO_CDF
#     include <netcdf.inc>
      integer :: iret
#endif

      ok = .false.
      mean = 0.0_R8
      islot = 0
      do i = 1, gui_mv_nslot
        if (gui_mv_c(i)%ifile.eq.ifile.and.
     .   gui_mv_c(i)%name.eq.vname) islot = i
      end do
      if (islot.gt.0) then
        if (size(gui_mv_c(islot)%s,1).ne.n.or.
     .   size(gui_mv_c(islot)%s,2).ne.ncol) islot = 0
      end if
      if (islot.eq.0) then
        islot = minloc(gui_mv_c(:)%used,1)
        if (allocated(gui_mv_c(islot)%s)) deallocate(gui_mv_c(islot)%s)
        allocate(gui_mv_c(islot)%s(n,ncol))
        gui_mv_c(islot)%s = 0.0_R8
        gui_mv_c(islot)%name = vname
        gui_mv_c(islot)%ifile = ifile
        gui_mv_c(islot)%lo = 1
        gui_mv_c(islot)%hi = 0
        gui_mv_c(islot)%nupd = 0
      end if
      gui_mv_clock = gui_mv_clock+1
      gui_mv_c(islot)%used = gui_mv_clock
      lo = gui_mv_c(islot)%lo
      hi = gui_mv_c(islot)%hi
      nmean = k2-k1+1
      nover = max(0,min(k2,hi)-max(k1,lo)+1)
      nadd = nmean-nover
      nrem = max(0,hi-lo+1)-nover
      allocate(frame(n,ncol))
      kstart = 1
      kcount(1) = n
      kcount(2) = ncol
      kcount(3) = 1
#ifndef NO_CDF
      if (hi.lt.lo.or.nadd+nrem.ge.nmean.or.
     . gui_mv_c(islot)%nupd.ge.gui_mv_maxupd) then
        gui_mv_c(islot)%s = 0.0_R8
        gui_mv_c(islot)%nupd = 0
        lo = 1
        hi = 0
      else
        gui_mv_c(islot)%nupd = gui_mv_c(islot)%nupd+nadd+nrem
      end if
      do k = min(k1,lo), max(k2,hi)
        if (k.ge.k1.and.k.le.k2.and.(k.lt.lo.or.k.gt.hi)) then
          i = 1
        else if (k.ge.lo.and.k.le.hi.and.(k.lt.k1.or.k.gt.k2)) then
          i = -1
        else
          cycle
        end if
        if (nd.eq.3) then
          kstart(3) = k
        else
          kstart(2) = k
          kcount(2) = 1
        end if
        iret = nf_get_vara_double(gui_mv_ncid(ifile),varid,kstart,
     .   kcount,frame)
        if (iret.ne.NF_NOERR) then
          write(*,*) 'gui_movie_mean: cannot read ',trim(vname),
     .     ' frame ',k
          gui_mv_c(islot)%hi = 0
          gui_mv_c(islot)%lo = 1
          deallocate(frame)
          return
        end if
        gui_mv_c(islot)%s = gui_mv_c(islot)%s+i*frame
      end do
      gui_mv_c(islot)%lo = k1
      gui_mv_c(islot)%hi = k2
      mean = gui_mv_c(islot)%s/nmean
      ok = .true.
#endif
      deallocate(frame)
      return
      end subroutine gui_movie_mean

c Push the mean of frames k1 to k2 of a movie variable on the matrix
c stack (file 1, at the cells) or on the Eirene stack (file 2, at the
c triangles); with yname also its partner on top, both or none
      subroutine gui_movie_push(ifile,xname,yname,k1,k2)
      use b2mod_types
      use b2us_data
      use b2mod_b2plot
      use b2mod_neutrals_namelist
     , , only: textan, textmn, textin
      use b2mod_b2plot_gui
#ifdef B25_EIRENE
      use eirmod_ctrig
      use eirmod_cgeom
     , , only: ncltal
#endif
      implicit none
      integer, intent(in) :: ifile, k1, k2
      character*(*), intent(in) :: xname, yname
      real (kind=R8), allocatable :: mean(:,:,:), cell(:,:)
      character*32, allocatable :: clab(:)
      character*256 :: vname(2), dname(2), long, units, label
      integer :: varid(2), loc(2), n(2), ncol(2), nd(2), nv, iv, is
#ifdef B25_EIRENE
      real (kind=R8), allocatable :: ext(:,:)
      integer :: it, iCv
#endif
      logical :: ok
#ifndef NO_CDF
#     include <netcdf.inc>
      integer :: iret
#endif
      external gui_movie_open, gui_movie_var, gui_movie_mean,
     . gui_movie_text, copy, eirene_copy, intcell

      call gui_movie_open
      nv = 1
      if (yname.ne.' ') nv = 2
      vname(1) = xname
      vname(2) = yname
      if (ifile.lt.1.or.ifile.gt.2) then
        write(*,*) 'guimovie: file must be 1 (b2movies.nc) or 2 ',
     .   '(eirenemovies.nc)'
        return
      end if
      if (gui_mv_ncid(ifile).lt.0) then
        write(*,*) 'guimovie: movie file ',ifile,' not found'
        return
      end if
      if (k1.lt.1.or.k2.lt.k1.or.k2.gt.gui_mv_nt(ifile)) then
        write(*,*) 'guimovie: frames must be 1 <= ',k1,' <= ',k2,
     .   ' <= ',gui_mv_nt(ifile)
        return
      end if
#ifndef NO_CDF
      do iv = 1, nv
        iret = nf_inq_varid(gui_mv_ncid(ifile),vname(iv),varid(iv))
        loc(iv) = 0
        if (iret.eq.NF_NOERR) call gui_movie_var(ifile,varid(iv),
     .   vname(iv),loc(iv),n(iv),ncol(iv),dname(iv),nd(iv))
        if (loc(iv).eq.0) then
          write(*,*) 'guimovie: no usable variable ',trim(vname(iv))
          return
        end if
      end do
      if (nv.eq.2) then
        if (loc(2).ne.loc(1).or.ncol(2).ne.ncol(1)) then
          write(*,*) 'guimvvec: ',trim(vname(1)),' and ',
     .     trim(vname(2)),' differ in location or columns'
          return
        end if
      end if
      allocate(mean(n(1),ncol(1),nv), clab(ncol(1)))
      do iv = 1, nv
        call gui_movie_mean(ifile,varid(iv),vname(iv),k1,k2,n(1),
     .   ncol(1),nd(iv),mean(1,1,iv),ok)
        if (.not.ok) then
          deallocate(mean, clab)
          return
        end if
      end do
      do is = 1, ncol(1)
        clab(is) = ' '
        if (dname(1).eq.'ns'.and.is.le.size(textpl))
     .   clab(is) = textpl(is-1)
        if (dname(1).eq.'natm'.and.is.le.size(textan))
     .   clab(is) = textan(is-1)
        if (dname(1).eq.'nmol'.and.is.le.size(textmn))
     .   clab(is) = textmn(is-1)
        if (dname(1).eq.'nion'.and.is.le.size(textin))
     .   clab(is) = textin(is-1)
      end do
      do iv = 1, nv
        call gui_movie_text(gui_mv_ncid(ifile),varid(iv),'long_name',
     .   long)
        call gui_movie_text(gui_mv_ncid(ifile),varid(iv),'units',units)
        if (long.eq.' ') long = vname(iv)
        if (k1.eq.k2) then
          write(label,'(a,i0,a)') trim(long)//', frame ',k1,
     .     ' ('//trim(units)//')'
        else
          write(label,'(a,i0,a,i0,a)') trim(long)//', frames ',k1,
     .     ' to ',k2,' ('//trim(units)//')'
        end if
        if (loc(1).eq.1) then
          call copy(mean(1,1,iv),1.0_R8,ncol(1),clab)
          graphlabel(iwrk) = label
        else if (loc(1).eq.2) then
          allocate(cell(mpg%nCv,ncol(1)))
          do is = 1, ncol(1)
            if (vname(iv)(len_trim(vname(iv)):len_trim(vname(iv)))
     .       .eq.'y') then
              call intcell(mpg%nFc,mpg%nCv,mpg,mpg%intcellR,
     .         mean(1,is,iv),cell(1,is))
            else
              call intcell(mpg%nFc,mpg%nCv,mpg,mpg%intcellP,
     .         mean(1,is,iv),cell(1,is))
            end if
          end do
          call copy(cell,1.0_R8,ncol(1),clab)
          graphlabel(iwrk) = label
          deallocate(cell)
#ifdef B25_EIRENE
        else
          allocate(ext(ncol(1),ntrii))
          do it = 1, ntrii
            iCv = ncltal(it)
            do is = 1, ncol(1)
              ext(is,it) = 0.0_R8
              if (iCv.ge.1.and.iCv.le.n(1)) ext(is,it) = mean(iCv,is,iv)
            end do
          end do
          call eirene_copy(ext,1.0_R8,ncol(1),clab)
          eirene_label(iewrk) = label
          deallocate(ext)
#endif
        end if
      end do
      deallocate(mean, clab)
#else
      write(*,*) 'guimovie: b2plot is compiled without netCDF'
#endif
      return
      end subroutine gui_movie_push
