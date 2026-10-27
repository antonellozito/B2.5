      module b2mod_b2plot_gui
      use b2mod_types
      implicit none

c GUI message channel of b2plot: binary messages on the FIFO named by
c B2PLOT_GUI_FIFO (see the b2plot.py documentation in the manual)
      integer, parameter :: gui_i1 = selected_int_kind(2)
      integer, parameter :: gui_i4 = selected_int_kind(9)
      integer, parameter :: gui_protocol_version = 1
      integer, parameter :: gui_magic = 1110593607
      logical, save :: gui_active = .false., gui_drawn = .false.
      integer, save :: gui_unit = -1
      integer(kind=gui_i1), allocatable, save :: gbuf(:)
      integer, save :: gnbuf = 0, gop = 0
      integer, save :: gui_tag_kind = 0, gui_tag_index = 0
      real (kind=R8), save :: gui_tag_value = 0.0_R8

      contains

c Open the FIFO and send HELLO and RUNINFO (inactive if not set)
      subroutine gui_startup(nCv,nCi,nFc,nVx,iscl,nx,ny,ns,natm,nmol,
     . nion,ntrii,nbnd,with_eirene,time,label,nsp,species)
      implicit none
      integer, intent(in) :: nCv,nCi,nFc,nVx,iscl,nx,ny,ns,natm,nmol,
     . nion,ntrii,nbnd,with_eirene,nsp
      real (kind=R8), intent(in) :: time
      character*(*), intent(in) :: label, species(*)
      character*512 fifo
      integer ios, i
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
      call PXFGETENV ('B2PLOT_GUI_FIFO', 0, fifo, lenval, ierror)
#else
      call get_environment_variable('B2PLOT_GUI_FIFO',
     . status=ierror, length=lenval)
      if (ierror.eq.0) call get_environment_variable('B2PLOT_GUI_FIFO',
     . value=fifo)
#endif
#endif
      if (len_trim(fifo).eq.0) return
#ifndef LEGACYCOMP
      open (newunit=gui_unit, file=trim(fifo), access='stream',
     . form='unformatted', action='write', status='old', iostat=ios)
#else
      gui_unit = newunit()
      open (unit=gui_unit, file=trim(fifo), access='stream',
     . form='unformatted', action='write', status='old', iostat=ios)
#endif
      if (ios.ne.0) then
        write(*,'(a)') 'b2plot: cannot open the GUI FIFO '//trim(fifo)
        return
      end if
      allocate(gbuf(65536))
      gnbuf = 0
      gui_active = .true.
c HELLO
      call gui_begin(1)
      call gui_put_i4(gui_magic)
      call gui_put_i4(gui_protocol_version)
      call gui_put_str('b2plot')
      call gui_end
c RUNINFO
      call gui_begin(2)
      call gui_put_i4(nCv)
      call gui_put_i4(nCi)
      call gui_put_i4(nFc)
      call gui_put_i4(nVx)
      call gui_put_i4(iscl)
      call gui_put_i4(nx)
      call gui_put_i4(ny)
      call gui_put_i4(ns)
      call gui_put_i4(natm)
      call gui_put_i4(nmol)
      call gui_put_i4(nion)
      call gui_put_i4(ntrii)
      call gui_put_i4(nbnd)
      call gui_put_i4(with_eirene)
      call gui_put_r8(time)
      call gui_put_str(label)
      call gui_put_i4(nsp)
      do i = 1, nsp
        call gui_put_str(species(i))
      end do
      call gui_end
      call gui_flush
      return
      end subroutine gui_startup

c Send BYE and close the FIFO
      subroutine gui_shutdown
      implicit none
      if (.not.gui_active) return
      call gui_begin(4)
      call gui_end
      call gui_flush
      if (gui_active) close(gui_unit)
      gui_active = .false.
      if (allocated(gbuf)) deallocate(gbuf)
      return
      end subroutine gui_shutdown

      subroutine gui_flush
      implicit none
      integer ios
      if (.not.gui_active) return
      flush(gui_unit, iostat=ios)
      if (ios.ne.0) call gui_lost
      return
      end subroutine gui_flush

c Stop sending after a write error; b2plot goes on as a text tool
      subroutine gui_lost
      implicit none
      integer ios
      write(*,'(a)') 'b2plot: GUI FIFO write failed, GUI output stopped'
      gui_active = .false.
      close(gui_unit, iostat=ios)
      if (allocated(gbuf)) deallocate(gbuf)
      return
      end subroutine gui_lost

      subroutine gui_begin(op)
      implicit none
      integer, intent(in) :: op
      gop = op
      gnbuf = 0
      return
      end subroutine gui_begin

c Make room for n more bytes in the payload buffer
      subroutine gui_reserve(n)
      implicit none
      integer, intent(in) :: n
      integer(kind=gui_i1), allocatable :: tmp(:)
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
      end subroutine gui_reserve

      subroutine gui_put_i4(i)
      implicit none
      integer, intent(in) :: i
      integer(kind=gui_i1) :: mold
      if (.not.gui_active) return
      call gui_reserve(4)
      gbuf(gnbuf+1:gnbuf+4) = transfer(int(i,gui_i4), mold, 4)
      gnbuf = gnbuf+4
      return
      end subroutine gui_put_i4

      subroutine gui_put_r4(x)
      implicit none
      real (kind=R4), intent(in) :: x
      integer(kind=gui_i1) :: mold
      if (.not.gui_active) return
      call gui_reserve(4)
      gbuf(gnbuf+1:gnbuf+4) = transfer(x, mold, 4)
      gnbuf = gnbuf+4
      return
      end subroutine gui_put_r4

      subroutine gui_put_r8(x)
      implicit none
      real (kind=R8), intent(in) :: x
      integer(kind=gui_i1) :: mold
      if (.not.gui_active) return
      call gui_reserve(8)
      gbuf(gnbuf+1:gnbuf+8) = transfer(x, mold, 8)
      gnbuf = gnbuf+8
      return
      end subroutine gui_put_r8

      subroutine gui_put_i4v(n,a)
      implicit none
      integer, intent(in) :: n, a(*)
      integer(kind=gui_i1) :: mold
      if (.not.gui_active.or.n.le.0) return
      call gui_reserve(4*n)
      gbuf(gnbuf+1:gnbuf+4*n) = transfer(int(a(1:n),gui_i4), mold, 4*n)
      gnbuf = gnbuf+4*n
      return
      end subroutine gui_put_i4v

      subroutine gui_put_r4v(n,a)
      implicit none
      integer, intent(in) :: n
      real (kind=R4), intent(in) :: a(*)
      integer(kind=gui_i1) :: mold
      if (.not.gui_active.or.n.le.0) return
      call gui_reserve(4*n)
      gbuf(gnbuf+1:gnbuf+4*n) = transfer(a(1:n), mold, 4*n)
      gnbuf = gnbuf+4*n
      return
      end subroutine gui_put_r4v

      subroutine gui_put_r8v(n,a)
      implicit none
      integer, intent(in) :: n
      real (kind=R8), intent(in) :: a(*)
      integer(kind=gui_i1) :: mold
      if (.not.gui_active.or.n.le.0) return
      call gui_reserve(8*n)
      gbuf(gnbuf+1:gnbuf+8*n) = transfer(a(1:n), mold, 8*n)
      gnbuf = gnbuf+8*n
      return
      end subroutine gui_put_r8v

c Strings are sent as length and bytes, without trailing blanks
      subroutine gui_put_str(s)
      implicit none
      character*(*), intent(in) :: s
      integer(kind=gui_i1) :: mold
      integer n
      if (.not.gui_active) return
      n = len_trim(s)
      call gui_put_i4(n)
      if (n.le.0) return
      call gui_reserve(n)
      gbuf(gnbuf+1:gnbuf+n) = transfer(s(1:n), mold, n)
      gnbuf = gnbuf+n
      return
      end subroutine gui_put_str

c Write opcode, payload length and payload
      subroutine gui_end
      implicit none
      integer ios
      if (.not.gui_active) return
      write(gui_unit, iostat=ios) int(gop,gui_i4), int(gnbuf,gui_i4),
     . gbuf(1:gnbuf)
      if (ios.ne.0) call gui_lost
c Output primitives mark the current frame as drawn
      if ((gop.ge.40.and.gop.le.49).or.(gop.ge.53.and.gop.le.55)
     . .or.gop.eq.62) gui_drawn = .true.
      return
      end subroutine gui_end

c Informative (1) or warning (2) text for the GUI
      subroutine gui_message(level,text)
      implicit none
      integer, intent(in) :: level
      character*(*), intent(in) :: text
      if (.not.gui_active) return
      call gui_begin(6)
      call gui_put_i4(level)
      call gui_put_str(text)
      call gui_end
      return
      end subroutine gui_message

c Pending cell tag for the next FILL (1 B2 cell, 2 EIRENE triangle,
c 3 wall segment, 4 target element)
      subroutine gui_tag(kind,index,value)
      implicit none
      integer, intent(in) :: kind, index
      real (kind=R8), intent(in) :: value
      gui_tag_kind = kind
      gui_tag_index = index
      gui_tag_value = value
      return
      end subroutine gui_tag

c Handshake at the end of a command batch
      subroutine gui_sync(id,depths,toplabel,topnz)
      implicit none
      integer, intent(in) :: id, depths(7), topnz
      character*(*), intent(in) :: toplabel
      if (.not.gui_active) return
      call gui_begin(3)
      call gui_put_i4(id)
      call gui_put_i4v(7,depths)
      call gui_put_str(toplabel)
      call gui_put_i4(topnz)
      call gui_end
      call gui_flush
      return
      end subroutine gui_sync

c Page and panel description, sent by init_page
      subroutine gui_page_meta(title,linlog,physical,aspect,colorbar)
      use b2mod_graphics
      implicit none
      character*(*), intent(in) :: title
      integer, intent(in) :: linlog
      logical, intent(in) :: physical, aspect, colorbar
      if (.not.gui_active) return
      call gui_begin(11)
      call gui_put_str(title)
      call gui_put_str(xlabel)
      call gui_put_str(ylabel)
      call gui_put_str(extralabel)
      call gui_put_str(globalheader)
      call gui_put_i4(linlog)
      call gui_put_i4(merge(1,0,physical))
      call gui_put_i4(merge(1,0,aspect))
      call gui_put_i4(merge(1,0,colorbar))
      call gui_put_i4(ix_page)
      call gui_put_i4(iy_page)
      call gui_put_i4(nx_page)
      call gui_put_i4(ny_page)
      call gui_put_r4(xwidth)
      call gui_put_r4(yheight)
      call gui_put_r4(xlo)
      call gui_put_r4(xhi)
      call gui_put_r4(ylo)
      call gui_put_r4(yhi)
      call gui_end
      return
      end subroutine gui_page_meta

c Print the b2plot catalogue and, for the GUI, send CATALOG
      subroutine gui_help_catalog
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

      if (gui_active) then
        call gui_begin(5)
        call gui_put_i4(ncat)
      end if
      do icat = 1, ncat
        n = gui_help_size(icat)
        write(*,'(/,a,a,i4,a)') trim(catname(icat)),' (',n,')'
        if (gui_active) then
          call gui_put_str(catname(icat))
          call gui_put_i4(n)
        end if
        do i = 1, gui_help_last(icat)
          if (len_trim(gui_help_name(icat,i)).eq.0) cycle
          desc = gui_help_desc(icat,i)
          write(*,'(2x,a8,i3,2x,a)') gui_help_name(icat,i),
     .     gui_help_dim(icat,i), trim(desc)
          if (gui_active) then
            call gui_put_str(gui_help_name(icat,i))
            call gui_put_str(desc)
            call gui_put_i4(gui_help_dim(icat,i))
          end if
        end do
      end do
      if (gui_active) then
        call gui_end
        call gui_flush
      end if
      return
      end subroutine gui_help_catalog

c Number of array entries of catalogue category icat
      integer function gui_help_last(icat)
      use b2mod_b2plot_help
      implicit none
      integer, intent(in) :: icat
      select case (icat)
      case (1)
        gui_help_last = nlabels
      case (2)
        gui_help_last = nfns
      case (3)
        gui_help_last = nelabels
      case (4)
        gui_help_last = nefns
      case (5)
        gui_help_last = nwlabels
      case (6)
        gui_help_last = nwfns
      case (7)
        gui_help_last = ntrgfns
      case (8)
        gui_help_last = nsparameters
      case default
        gui_help_last = nconstants
      end select
      return
      end function gui_help_last

c Number of non-blank names of catalogue category icat
      integer function gui_help_size(icat)
      implicit none
      integer, intent(in) :: icat
      integer i
      gui_help_size = 0
      do i = 1, gui_help_last(icat)
        if (len_trim(gui_help_name(icat,i)).gt.0)
     .   gui_help_size = gui_help_size+1
      end do
      return
      end function gui_help_size

      character*8 function gui_help_name(icat,i)
      use b2mod_b2plot_help
      implicit none
      integer, intent(in) :: icat, i
      select case (icat)
      case (1)
        gui_help_name = labels(i)
      case (2)
        gui_help_name = fns(i)
      case (3)
        gui_help_name = elabels(i)
      case (4)
        gui_help_name = efns(i)
      case (5)
        gui_help_name = wlabels(i)
      case (6)
        gui_help_name = wfns(i)
      case (7)
        gui_help_name = trgfns(i)
      case (8)
        gui_help_name = sparameters(i)
      case default
        gui_help_name = constant_name(i)
      end select
      return
      end function gui_help_name

      character*64 function gui_help_desc(icat,i)
      use b2mod_b2plot_help
      implicit none
      integer, intent(in) :: icat, i
      select case (icat)
      case (1)
        gui_help_desc = descriptions(i)
      case (3)
        gui_help_desc = edescriptions(i)
      case (5)
        gui_help_desc = wdescriptions(i)
      case (9)
        write(gui_help_desc,'(1pe14.7)') constant_value(i)
      case default
        gui_help_desc = ' '
      end select
      return
      end function gui_help_desc

      integer function gui_help_dim(icat,i)
      use b2mod_b2plot_help
      implicit none
      integer, intent(in) :: icat, i
      select case (icat)
      case (1)
        gui_help_dim = dimensions(i)
      case (3)
        gui_help_dim = edimensions(i)
      case (5)
        gui_help_dim = wdimensions(i)
      case default
        gui_help_dim = 0
      end select
      return
      end function gui_help_dim

      end module b2mod_b2plot_gui
