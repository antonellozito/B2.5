      module b2mod_b2plot_py
      use b2mod_types
      implicit none

c Message channel of b2plot to its Python graphics: binary messages on
c the FIFO named by B2PLOT_PY_FIFO (see the Python graphics section of
c the b2plot manual)
      integer, parameter :: b2py_i1 = selected_int_kind(2)
      integer, parameter :: b2py_i4 = selected_int_kind(9)
      integer, parameter :: b2py_magic = 1110593607
      logical, save :: b2py_active = .false., b2py_drawn = .false.
c Channel opened by gopks of the NCAR replacement (closed by gclks)
      logical, save :: b2py_by_gks = .false.
c Next trajectories (traj) drawn as an overlay of the GUI (guitraj)
      logical, save :: b2py_traj = .false.
      integer, save :: b2py_unit = -1
      integer(kind=b2py_i1), allocatable, save :: gbuf(:)
      integer, save :: gnbuf = 0, gop = 0
      integer, save :: b2py_tag_kind = 0, b2py_tag_index = 0
      real (kind=R8), save :: b2py_tag_value = 0.0_R8
c Movie files (movie_frames.F): netCDF ids, frame counts and
c running sums of frames and of their squares, a full sum after
c mvf_maxupd updates; file 3 is b2tallies.nc (volume regions)
      integer, parameter :: mvf_nslot = 4, mvf_maxupd = 200
      type mvf_slot
        character*64 :: name = ' '
        integer :: ifile = 0, lo = 1, hi = 0, nupd = 0, used = 0
        real (kind=R8), allocatable :: s(:,:), s2(:,:)
      end type mvf_slot
      logical, save :: mvf_opened = .false.
      integer, save :: mvf_ncid(3) = -1, mvf_nt(3) = 0
      integer, save :: mvf_clock = 0
      character*256, save :: mvf_path(3) = ' '
      type (mvf_slot), save :: mvf_c(mvf_nslot)

      contains

c Open the FIFO and send HELLO and RUNINFO (inactive if not set)
      subroutine b2py_startup(nCv,nCi,nFc,nVx,iscl,nx,ny,ns,natm,nmol,
     . nion,ntrii,nbnd,with_eirene,time,label,nsp,species)
      implicit none
      integer, intent(in) :: nCv,nCi,nFc,nVx,iscl,nx,ny,ns,natm,nmol,
     . nion,ntrii,nbnd,with_eirene,nsp
      real (kind=R8), intent(in) :: time
      character*(*), intent(in) :: label, species(*)
      integer i

      if (.not.b2py_active) call b2py_open('b2plot')
      if (.not.b2py_active) return
      b2py_by_gks = .false.
c RUNINFO
      call b2py_begin(2)
      call b2py_put_i4(nCv)
      call b2py_put_i4(nCi)
      call b2py_put_i4(nFc)
      call b2py_put_i4(nVx)
      call b2py_put_i4(iscl)
      call b2py_put_i4(nx)
      call b2py_put_i4(ny)
      call b2py_put_i4(ns)
      call b2py_put_i4(natm)
      call b2py_put_i4(nmol)
      call b2py_put_i4(nion)
      call b2py_put_i4(ntrii)
      call b2py_put_i4(nbnd)
      call b2py_put_i4(with_eirene)
      call b2py_put_r8(time)
      call b2py_put_str(label)
      call b2py_put_i4(nsp)
      do i = 1, nsp
        call b2py_put_str(species(i))
      end do
      call b2py_end
      call b2py_flush
      return
      end subroutine b2py_startup

c Open the FIFO and send HELLO with the program name (inactive if not
c set)
      subroutine b2py_open(producer)
      implicit none
      character*(*), intent(in) :: producer
      character*512 fifo
      integer ios
#ifndef NO_GETENV
      integer lenval, ierror
#ifndef USE_PXFGETENV
      intrinsic get_environment_variable
#endif
#endif
#ifdef LEGACYCOMP
      integer newunit
      external newunit
#endif

      fifo=' '
#ifndef NO_GETENV
#ifdef USE_PXFGETENV
      call PXFGETENV ('B2PLOT_PY_FIFO', 0, fifo, lenval, ierror)
#else
      call get_environment_variable('B2PLOT_PY_FIFO',
     . status=ierror, length=lenval)
      if (ierror.eq.0) call get_environment_variable('B2PLOT_PY_FIFO',
     . value=fifo)
#endif
#endif
      if (len_trim(fifo).eq.0) return
#ifndef LEGACYCOMP
      open (newunit=b2py_unit, file=trim(fifo), access='stream',
     . form='unformatted', action='write', status='old', iostat=ios)
#else
      b2py_unit = newunit()
      open (unit=b2py_unit, file=trim(fifo), access='stream',
     . form='unformatted', action='write', status='old', iostat=ios)
#endif
      if (ios.ne.0) then
        write(*,'(a)') trim(producer)//
     .    ': cannot open the Python graphics FIFO '//trim(fifo)
        return
      end if
      allocate(gbuf(65536))
      gnbuf = 0
      b2py_active = .true.
c HELLO
      call b2py_begin(1)
      call b2py_put_i4(b2py_magic)
      call b2py_put_str(producer)
      call b2py_end
      call b2py_flush
      return
      end subroutine b2py_open

c Send BYE and close the FIFO
      subroutine b2py_shutdown
      implicit none
      external mvf_close
      call mvf_close
      if (.not.b2py_active) return
      call b2py_begin(4)
      call b2py_end
      call b2py_flush
      if (b2py_active) close(b2py_unit)
      b2py_active = .false.
      b2py_by_gks = .false.
      if (allocated(gbuf)) deallocate(gbuf)
      return
      end subroutine b2py_shutdown

      subroutine b2py_flush
      implicit none
      integer ios
      if (.not.b2py_active) return
      flush(b2py_unit, iostat=ios)
      if (ios.ne.0) call b2py_lost
      return
      end subroutine b2py_flush

c Stop sending after a write error; b2plot goes on as a text tool
      subroutine b2py_lost
      implicit none
      integer ios
      write(*,'(a)') 'b2plot: Python graphics FIFO write failed, '//
     .  'output stopped'
      b2py_active = .false.
      close(b2py_unit, iostat=ios)
      if (allocated(gbuf)) deallocate(gbuf)
      return
      end subroutine b2py_lost

      subroutine b2py_begin(op)
      implicit none
      integer, intent(in) :: op
      gop = op
      gnbuf = 0
      return
      end subroutine b2py_begin

c Make room for n more bytes in the payload buffer
      subroutine b2py_reserve(n)
      implicit none
      integer, intent(in) :: n
      integer(kind=b2py_i1), allocatable :: tmp(:)
      integer nnew
      if (gnbuf+n.le.size(gbuf)) return
      nnew = size(gbuf)
      do while (gnbuf+n.gt.nnew)
        nnew = 2*nnew
      end do
      allocate(tmp(nnew))
      tmp(1:gnbuf) = gbuf(1:gnbuf)
      call move_alloc(tmp, gbuf)
      return
      end subroutine b2py_reserve

      subroutine b2py_put_i4(i)
      implicit none
      integer, intent(in) :: i
      integer(kind=b2py_i1) :: mold
      if (.not.b2py_active) return
      call b2py_reserve(4)
      gbuf(gnbuf+1:gnbuf+4) = transfer(int(i,b2py_i4), mold, 4)
      gnbuf = gnbuf+4
      return
      end subroutine b2py_put_i4

      subroutine b2py_put_r4(x)
      implicit none
      real (kind=R4), intent(in) :: x
      integer(kind=b2py_i1) :: mold
      if (.not.b2py_active) return
      call b2py_reserve(4)
      gbuf(gnbuf+1:gnbuf+4) = transfer(x, mold, 4)
      gnbuf = gnbuf+4
      return
      end subroutine b2py_put_r4

      subroutine b2py_put_r8(x)
      implicit none
      real (kind=R8), intent(in) :: x
      integer(kind=b2py_i1) :: mold
      if (.not.b2py_active) return
      call b2py_reserve(8)
      gbuf(gnbuf+1:gnbuf+8) = transfer(x, mold, 8)
      gnbuf = gnbuf+8
      return
      end subroutine b2py_put_r8

      subroutine b2py_put_i4v(n,a)
      implicit none
      integer, intent(in) :: n, a(*)
      integer(kind=b2py_i1) :: mold
      if (.not.b2py_active.or.n.le.0) return
      call b2py_reserve(4*n)
      gbuf(gnbuf+1:gnbuf+4*n) = transfer(int(a(1:n),b2py_i4), mold, 4*n)
      gnbuf = gnbuf+4*n
      return
      end subroutine b2py_put_i4v

      subroutine b2py_put_r4v(n,a)
      implicit none
      integer, intent(in) :: n
      real (kind=R4), intent(in) :: a(*)
      integer(kind=b2py_i1) :: mold
      if (.not.b2py_active.or.n.le.0) return
      call b2py_reserve(4*n)
      gbuf(gnbuf+1:gnbuf+4*n) = transfer(a(1:n), mold, 4*n)
      gnbuf = gnbuf+4*n
      return
      end subroutine b2py_put_r4v

      subroutine b2py_put_r8v(n,a)
      implicit none
      integer, intent(in) :: n
      real (kind=R8), intent(in) :: a(*)
      integer(kind=b2py_i1) :: mold
      if (.not.b2py_active.or.n.le.0) return
      call b2py_reserve(8*n)
      gbuf(gnbuf+1:gnbuf+8*n) = transfer(a(1:n), mold, 8*n)
      gnbuf = gnbuf+8*n
      return
      end subroutine b2py_put_r8v

c Strings are sent as length and bytes, without trailing blanks
      subroutine b2py_put_str(s)
      implicit none
      character*(*), intent(in) :: s
      integer(kind=b2py_i1) :: mold
      integer n
      if (.not.b2py_active) return
      n = len_trim(s)
      call b2py_put_i4(n)
      if (n.le.0) return
      call b2py_reserve(n)
      gbuf(gnbuf+1:gnbuf+n) = transfer(s(1:n), mold, n)
      gnbuf = gnbuf+n
      return
      end subroutine b2py_put_str

c Write opcode, payload length and payload
      subroutine b2py_end
      implicit none
      integer ios
      if (.not.b2py_active) return
      write(b2py_unit, iostat=ios) int(gop,b2py_i4), int(gnbuf,b2py_i4),
     . gbuf(1:gnbuf)
      if (ios.ne.0) call b2py_lost
c Output primitives mark the current frame as drawn
      if ((gop.ge.40.and.gop.le.49).or.(gop.ge.53.and.gop.le.55)
     . .or.gop.eq.62) b2py_drawn = .true.
      return
      end subroutine b2py_end

c Informative (1) or warning (2) text for the Python graphics
      subroutine b2py_message(level,text)
      implicit none
      integer, intent(in) :: level
      character*(*), intent(in) :: text
      if (.not.b2py_active) return
      call b2py_begin(6)
      call b2py_put_i4(level)
      call b2py_put_str(text)
      call b2py_end
      return
      end subroutine b2py_message

c Pending cell tag for the next FILL (1 B2 cell, 2 EIRENE triangle,
c 3 wall segment, 4 target element)
      subroutine b2py_tag(kind,index,value)
      implicit none
      integer, intent(in) :: kind, index
      real (kind=R8), intent(in) :: value
      b2py_tag_kind = kind
      b2py_tag_index = index
      b2py_tag_value = value
      return
      end subroutine b2py_tag

c Handshake at the end of a command batch
      subroutine b2py_sync(id,depths,toplabel,topnz)
      implicit none
      integer, intent(in) :: id, depths(7), topnz
      character*(*), intent(in) :: toplabel
      if (.not.b2py_active) return
      call b2py_begin(3)
      call b2py_put_i4(id)
      call b2py_put_i4v(7,depths)
      call b2py_put_str(toplabel)
      call b2py_put_i4(topnz)
      call b2py_end
      call b2py_flush
      return
      end subroutine b2py_sync

c Start and end of an overlay drawn by a b2plot routine (role 1 device
c structure, 2 separatrices, 3 grid boundary); the lines come between
      subroutine gui_overlay_begin(role,name)
      implicit none
      integer, intent(in) :: role
      character*(*), intent(in) :: name
      if (.not.b2py_active) return
      call b2py_begin(71)
      call b2py_put_i4(role)
      call b2py_put_str(name)
      call b2py_end
      return
      end subroutine gui_overlay_begin

      subroutine gui_overlay_end
      implicit none
      if (.not.b2py_active) return
      call b2py_begin(72)
      call b2py_end
      return
      end subroutine gui_overlay_end

c Page and panel description, sent by init_page
      subroutine b2py_page_meta(title,linlog,physical,aspect,colorbar)
      use b2mod_graphics
      implicit none
      character*(*), intent(in) :: title
      integer, intent(in) :: linlog
      logical, intent(in) :: physical, aspect, colorbar
      if (.not.b2py_active) return
      call b2py_begin(11)
      call b2py_put_str(title)
      call b2py_put_str(xlabel)
      call b2py_put_str(ylabel)
      call b2py_put_str(extralabel)
      call b2py_put_str(globalheader)
      call b2py_put_i4(linlog)
      call b2py_put_i4(merge(1,0,physical))
      call b2py_put_i4(merge(1,0,aspect))
      call b2py_put_i4(merge(1,0,colorbar))
      call b2py_put_i4(ix_page)
      call b2py_put_i4(iy_page)
      call b2py_put_i4(nx_page)
      call b2py_put_i4(ny_page)
      call b2py_put_r4(xwidth)
      call b2py_put_r4(yheight)
      call b2py_put_r4(xlo)
      call b2py_put_r4(xhi)
      call b2py_put_r4(ylo)
      call b2py_put_r4(yhi)
      call b2py_end
      return
      end subroutine b2py_page_meta

c Print the b2plot catalogue and send CATALOG to the Python graphics
      subroutine b2py_help_catalog
      use b2mod_b2plot_help
      implicit none
      integer, parameter :: ncat = 9
      character*24 catname(ncat)
      character*64 desc
      integer icat, i, n
      data catname / 'Matrix quantities', 'Matrix functions',
     . 'Eirene quantities', 'Eirene functions', 'Wall quantities',
     . 'Wall functions', 'Target functions', 'Parameters',
     . 'Constants' /

      if (b2py_active) then
        call b2py_begin(5)
        call b2py_put_i4(ncat)
      end if
      do icat = 1, ncat
        n = b2py_help_size(icat)
        write(*,'(/,a,a,i4,a)') trim(catname(icat)),' (',n,')'
        if (b2py_active) then
          call b2py_put_str(catname(icat))
          call b2py_put_i4(n)
        end if
        do i = 1, b2py_help_last(icat)
          if (len_trim(b2py_help_name(icat,i)).eq.0) cycle
          desc = b2py_help_desc(icat,i)
          write(*,'(2x,a8,i3,2x,a)') b2py_help_name(icat,i),
     .     b2py_help_dim(icat,i), trim(desc)
          if (b2py_active) then
            call b2py_put_str(b2py_help_name(icat,i))
            call b2py_put_str(desc)
            call b2py_put_i4(b2py_help_dim(icat,i))
          end if
        end do
      end do
      if (b2py_active) then
        call b2py_end
        call b2py_flush
      end if
      return
      end subroutine b2py_help_catalog

c Number of array entries of catalogue category icat
      integer function b2py_help_last(icat)
      use b2mod_b2plot_help
      implicit none
      integer, intent(in) :: icat
      select case (icat)
      case (1)
        b2py_help_last = nlabels
      case (2)
        b2py_help_last = nfns
      case (3)
        b2py_help_last = nelabels
      case (4)
        b2py_help_last = nefns
      case (5)
        b2py_help_last = nwlabels
      case (6)
        b2py_help_last = nwfns
      case (7)
        b2py_help_last = ntrgfns
      case (8)
        b2py_help_last = nsparameters
      case default
        b2py_help_last = nconstants
      end select
      return
      end function b2py_help_last

c Number of non-blank names of catalogue category icat
      integer function b2py_help_size(icat)
      implicit none
      integer, intent(in) :: icat
      integer i
      b2py_help_size = 0
      do i = 1, b2py_help_last(icat)
        if (len_trim(b2py_help_name(icat,i)).gt.0)
     .   b2py_help_size = b2py_help_size+1
      end do
      return
      end function b2py_help_size

      character*8 function b2py_help_name(icat,i)
      use b2mod_b2plot_help
      implicit none
      integer, intent(in) :: icat, i
#ifndef WG_TODO
c Commands compiled out of b2plot.F on these grids: not listed
      character*8, parameter :: absent(19) = (/ 'bvis    ',
     . 'niso    ', 'lix     ', 'liy     ', 'rix     ', 'riy     ',
     . 'bix     ', 'biy     ', 'tix     ', 'tiy     ', 'crx     ',
     . 'cry     ', 'saw     ', 'oson    ', 'uedgegrd', 'trgnorth',
     . 'trgsouth', 'trgwest ', 'trgeast ' /)
#endif
      select case (icat)
      case (1)
        b2py_help_name = labels(i)
      case (2)
        b2py_help_name = fns(i)
      case (3)
        b2py_help_name = elabels(i)
      case (4)
        b2py_help_name = efns(i)
      case (5)
        b2py_help_name = wlabels(i)
      case (6)
        b2py_help_name = wfns(i)
      case (7)
        b2py_help_name = trgfns(i)
      case (8)
        b2py_help_name = sparameters(i)
      case default
        b2py_help_name = constant_name(i)
      end select
#ifndef WG_TODO
      if (any(absent.eq.b2py_help_name)) b2py_help_name = ' '
#endif
      return
      end function b2py_help_name

      character*64 function b2py_help_desc(icat,i)
      use b2mod_b2plot_help
      implicit none
      integer, intent(in) :: icat, i
      select case (icat)
      case (1)
        b2py_help_desc = descriptions(i)
      case (3)
        b2py_help_desc = edescriptions(i)
      case (5)
        b2py_help_desc = wdescriptions(i)
      case (9)
        write(b2py_help_desc,'(1pe14.7)') constant_value(i)
      case default
        b2py_help_desc = ' '
      end select
      return
      end function b2py_help_desc

      integer function b2py_help_dim(icat,i)
      use b2mod_b2plot_help
      implicit none
      integer, intent(in) :: icat, i
      select case (icat)
      case (1)
        b2py_help_dim = dimensions(i)
      case (3)
        b2py_help_dim = edimensions(i)
      case (5)
        b2py_help_dim = wdimensions(i)
      case default
        b2py_help_dim = 0
      end select
      return
      end function b2py_help_dim

      end module b2mod_b2plot_py
