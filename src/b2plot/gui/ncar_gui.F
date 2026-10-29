c NCAR Graphics replacement for b2pl_py.exe: each NCAR and GKS entry
c point used by b2plot sends a message on the GUI FIFO instead of
c drawing (protocol in the b2plot.py documentation of the manual)

      module ncar_gui_state
      use b2mod_types
      implicit none
c Parameter tables: 1 Conpackt, 2 Conpack
      integer, parameter :: gui_ct = 1, gui_cp = 2, gui_nlev = 256
      integer, save :: tpai(2) = 1, tncl(2) = 0, tcls(2) = 0,
     . tllp(2) = 0, tllo(2) = 0, tmap(2) = 0, tset(2) = 1,
     . tpic(2) = 0, tgic(2) = 3, tgil(2) = 3, tgis(2) = 4,
     . tnvs(2) = 1, tilp(2) = 0
      integer, save :: tclu(gui_nlev,2) = 0, taia(gui_nlev,2) = 0,
     . taib(gui_nlev,2) = 0
      real (kind=R4), save :: tclv(gui_nlev,2) = 0.0_R4,
     . tcll(gui_nlev,2) = 0.0_R4, tt2d(2) = 0.0_R4,
     . tspv(2) = 0.0_R4, txc1(2) = 0.0_R4, txcm(2) = 0.0_R4,
     . tyc1(2) = 0.0_R4, tycn(2) = 0.0_R4, tilx(2) = 0.0_R4,
     . tily(2) = 0.0_R4
      character*256, save :: tilt(2) = ' ', tcft(2) = ' '
c Contour source of the current plot (gui_ct or gui_cp)
      integer, save :: current_source = 0
c Colour capture while probing the area fill callback
      logical, save :: probing = .false.
      integer, save :: probe_color = -1
c Open workstations (GFLASH stores excluded)
      integer, save :: nwkopen = 0, iwkopen(16) = 0
c Current polyline colour and stored LABMOD arguments
      integer, save :: cur_plci = 1
      character*32, save :: lm_fmtx = ' ', lm_fmty = ' '
      integer, save :: lm_numx = 0, lm_numy = 0, lm_iszx = 0,
     . lm_iszy = 0, lm_ixdc = 0, lm_iydc = 0, lm_ixor = 0
c Unknown parameter names already reported
      integer, save :: nunknown = 0
      character*3, save :: unknown(64) = ' '
c Sink for NCAR arguments the GUI does not need
      integer, save :: idum = 0
      real (kind=R4), save :: rdum = 0.0_R4

      contains

c Lowest open workstation, 0 if none is open
      integer function gui_lowest_wk()
      implicit none
      if (nwkopen.eq.0) then
        gui_lowest_wk = 0
      else
        gui_lowest_wk = minval(iwkopen(1:nwkopen))
      end if
      return
      end function gui_lowest_wk

c First three characters of a parameter name, in upper case
      character*3 function gui_pkey(pnam)
      implicit none
      character*(*), intent(in) :: pnam
      integer i, ic
      gui_pkey = pnam
      do i = 1, 3
        ic = iachar(gui_pkey(i:i))
        if (ic.ge.iachar('a').and.ic.le.iachar('z'))
     .   gui_pkey(i:i) = achar(ic-32)
      end do
      return
      end function gui_pkey

c Report a parameter name that the tables do not know, once
      subroutine gui_unknown(key)
      use b2mod_b2plot_gui
      implicit none
      character*3, intent(in) :: key
      integer i
      do i = 1, nunknown
        if (unknown(i).eq.key) return
      end do
      if (nunknown.lt.size(unknown)) then
        nunknown = nunknown+1
        unknown(nunknown) = key
      end if
      write(*,'(a)') 'ncar_gui: parameter '//key//' not emulated'
      call gui_message(2, 'ncar_gui: parameter '//key//' not emulated')
      return
      end subroutine gui_unknown

c Store a Conpack or Conpackt parameter (kind 1 int, 2 real, 3 char)
      subroutine gui_pset(it,pnam,kind,ival,rval,cval)
      implicit none
      integer, intent(in) :: it, kind, ival
      real (kind=R4), intent(in) :: rval
      character*(*), intent(in) :: pnam, cval
      character*3 key
      integer iv, ip
      real (kind=R4) :: rv
      logical lev
      key = gui_pkey(pnam)
      if (kind.eq.2) then
        iv = nint(rval)
        rv = rval
      else
        iv = ival
        rv = real(ival,R4)
      end if
      ip = tpai(it)
      lev = ip.ge.1.and.ip.le.gui_nlev
      select case (key)
      case ('PAI')
        tpai(it) = iv
      case ('CLV')
        if (lev) tclv(ip,it) = rv
      case ('CLU')
        if (lev) tclu(ip,it) = iv
      case ('AIA')
        if (lev) taia(ip,it) = iv
      case ('AIB')
        if (lev) taib(ip,it) = iv
      case ('CLL')
        if (lev) tcll(ip,it) = rv
      case ('NCL')
        tncl(it) = max(0,min(iv,gui_nlev))
      case ('CLS')
        tcls(it) = iv
      case ('LLP')
        tllp(it) = iv
      case ('LLO')
        tllo(it) = iv
      case ('MAP')
        tmap(it) = iv
      case ('SET')
        tset(it) = iv
      case ('PIC')
        tpic(it) = iv
      case ('GIC')
        tgic(it) = iv
      case ('GIL')
        tgil(it) = iv
      case ('GIS')
        tgis(it) = iv
      case ('NVS')
        tnvs(it) = iv
      case ('ILP')
        tilp(it) = iv
      case ('T2D')
        tt2d(it) = rv
      case ('SPV')
        tspv(it) = rv
      case ('XC1')
        txc1(it) = rv
      case ('XCM')
        txcm(it) = rv
      case ('YC1')
        tyc1(it) = rv
      case ('YCN')
        tycn(it) = rv
      case ('ILX')
        tilx(it) = rv
      case ('ILY')
        tily(it) = rv
      case ('ILT')
        if (kind.eq.3) tilt(it) = cval
      case ('CFT')
        if (kind.eq.3) tcft(it) = cval
      case ('NEL','NET','NEU','NLS','NLZ','NOF','NSD')
c numeric label format, not needed by the GUI
      case default
        call gui_unknown(key)
      end select
      return
      end subroutine gui_pset

c Return a real-valued parameter of table it
      real (kind=R4) function gui_pget(it,pnam)
      implicit none
      integer, intent(in) :: it
      character*(*), intent(in) :: pnam
      integer ip
      ip = max(1,min(tpai(it),gui_nlev))
      select case (gui_pkey(pnam))
      case ('CLV')
        gui_pget = tclv(ip,it)
      case ('CLL')
        gui_pget = tcll(ip,it)
      case ('T2D')
        gui_pget = tt2d(it)
      case ('SPV')
        gui_pget = tspv(it)
      case ('XC1')
        gui_pget = txc1(it)
      case ('XCM')
        gui_pget = txcm(it)
      case ('YC1')
        gui_pget = tyc1(it)
      case ('YCN')
        gui_pget = tycn(it)
      case ('ILX')
        gui_pget = tilx(it)
      case ('ILY')
        gui_pget = tily(it)
      case ('NCL')
        gui_pget = real(tncl(it),R4)
      case default
        gui_pget = 0.0_R4
        call gui_unknown(gui_pkey(pnam))
      end select
      return
      end function gui_pget

c Send CONTOUR_LINES or CONTOUR_INFO for table it
      subroutine gui_contour_lines(it)
      use b2mod_b2plot_gui
      implicit none
      integer, intent(in) :: it
      integer n
      if (.not.gui_active) return
      n = tncl(it)
      call gui_begin(54)
      call gui_put_i4(it)
      call gui_put_i4(n)
      call gui_put_r4v(n,tclv(1,it))
      call gui_put_i4v(n,tclu(1,it))
      call gui_put_r4v(n,tcll(1,it))
      call gui_put_i4(cur_plci)
      call gui_put_i4(tllp(it))
      call gui_put_i4(tllo(it))
      call gui_put_r4(txc1(it))
      call gui_put_r4(txcm(it))
      call gui_put_r4(tyc1(it))
      call gui_put_r4(tycn(it))
      call gui_put_r4(tspv(it))
      call gui_put_i4(tmap(it))
      call gui_end
      return
      end subroutine gui_contour_lines

      subroutine gui_contour_info(it)
      use b2mod_b2plot_gui
      implicit none
      integer, intent(in) :: it
      if (.not.gui_active) return
      call gui_begin(55)
      call gui_put_i4(it)
      call gui_put_str(tilt(it))
      call gui_put_r4(tilx(it))
      call gui_put_r4(tily(it))
      call gui_put_i4(tilp(it))
      call gui_put_i4(tllp(it))
      call gui_put_i4(tllo(it))
      call gui_end
      return
      end subroutine gui_contour_info

c Send a polyline or polygon of n points
      subroutine gui_points(n,x,y)
      use b2mod_b2plot_gui
      implicit none
      integer, intent(in) :: n
      real (kind=R4), intent(in) :: x(*), y(*)
      call gui_put_i4(max(n,0))
      call gui_put_r4v(n,x)
      call gui_put_r4v(n,y)
      return
      end subroutine gui_points

      end module ncar_gui_state

c ===================== GKS control and workstations ===================
      subroutine gopks(errfil,bufa)
      use ncar_gui_state
      implicit none
      integer errfil, bufa
      if (.false.) idum = errfil+bufa
      nwkopen = 0
      return
      end subroutine gopks

      subroutine gopwk(wkid,conid,wtype)
      use ncar_gui_state
      implicit none
      integer wkid, conid, wtype
      integer i
      if (.false.) idum = conid
      if (wtype.eq.3) return
      do i = 1, nwkopen
        if (iwkopen(i).eq.wkid) return
      end do
      if (nwkopen.lt.size(iwkopen)) then
        nwkopen = nwkopen+1
        iwkopen(nwkopen) = wkid
      end if
      return
      end subroutine gopwk

      subroutine gclwk(wkid)
      use ncar_gui_state
      implicit none
      integer wkid
      integer i, j
      j = 0
      do i = 1, nwkopen
        if (iwkopen(i).ne.wkid) then
          j = j+1
          iwkopen(j) = iwkopen(i)
        end if
      end do
      nwkopen = j
      return
      end subroutine gclwk

      subroutine gacwk(wkid)
      use ncar_gui_state
      implicit none
      integer wkid
      if (.false.) idum = wkid
      return
      end subroutine gacwk

      subroutine gdawk(wkid)
      use ncar_gui_state
      implicit none
      integer wkid
      if (.false.) idum = wkid
      return
      end subroutine gdawk

      subroutine gclks
      implicit none
      return
      end subroutine gclks

      subroutine opngks
      implicit none
      return
      end subroutine opngks

      subroutine clsgks
      implicit none
      return
      end subroutine clsgks

      subroutine ngseti(pnam,ival)
      use ncar_gui_state
      implicit none
      character*(*) pnam
      integer ival
      if (.false.) idum = len(pnam)+ival
      return
      end subroutine ngseti

      subroutine ngsetc(pnam,cval)
      use ncar_gui_state
      implicit none
      character*(*) pnam, cval
      if (.false.) idum = len(pnam)+len(cval)
      return
      end subroutine ngsetc

c ======================== Frames and mapping ==========================
c FRAME only after something was drawn, as NCAR clears only used frames
      subroutine frame
      use b2mod_b2plot_gui
      implicit none
      if (.not.gui_active.or..not.gui_drawn) return
      call gui_begin(10)
      call gui_end
      call gui_flush
      gui_drawn = .false.
      return
      end subroutine frame

c One FRAME per page: end_page calls this once per workstation
      subroutine ngpict(wkid,action)
      use b2mod_b2plot_gui
      use ncar_gui_state
      implicit none
      integer wkid, action
      if (.false.) idum = action
      if (.not.gui_active.or..not.gui_drawn) return
      if (wkid.ne.gui_lowest_wk()) return
      call gui_begin(10)
      call gui_end
      call gui_flush
      gui_drawn = .false.
      return
      end subroutine ngpict

      subroutine set(vl,vr,vb,vt,wl,wr,wb,wt,lf)
      use b2mod_b2plot_gui
      implicit none
      real (kind=R4) :: vl, vr, vb, vt, wl, wr, wb, wt
      integer lf
      if (.not.gui_active) return
      call gui_begin(12)
      call gui_put_r4(vl)
      call gui_put_r4(vr)
      call gui_put_r4(vb)
      call gui_put_r4(vt)
      call gui_put_r4(wl)
      call gui_put_r4(wr)
      call gui_put_r4(wb)
      call gui_put_r4(wt)
      call gui_put_i4(lf)
      call gui_end
      return
      end subroutine set

c ===================== Colours and attributes =========================
      subroutine gscr(wkid,ci,cr,cg,cb)
      use b2mod_b2plot_gui
      use ncar_gui_state
      implicit none
      integer wkid, ci
      real (kind=R4) :: cr, cg, cb
      integer iw
      if (.not.gui_active) return
      iw = gui_lowest_wk()
      if (iw.ne.0.and.wkid.ne.iw) return
      call gui_begin(20)
      call gui_put_i4(ci)
      call gui_put_r4(cr)
      call gui_put_r4(cg)
      call gui_put_r4(cb)
      call gui_end
      return
      end subroutine gscr

c Integer attribute message
      subroutine gui_attr_i(op,ival)
      use b2mod_b2plot_gui
      implicit none
      integer op, ival
      if (.not.gui_active) return
      call gui_begin(op)
      call gui_put_i4(ival)
      call gui_end
      return
      end subroutine gui_attr_i

c Real attribute message
      subroutine gui_attr_r(op,rval)
      use b2mod_b2plot_gui
      implicit none
      integer op
      real (kind=R4) :: rval
      if (.not.gui_active) return
      call gui_begin(op)
      call gui_put_r4(rval)
      call gui_end
      return
      end subroutine gui_attr_r

c Character attribute message (1 CHH, 2 CHSP, 3 CHXP, 4 CHUP)
      subroutine gui_chattr(which,a,b)
      use b2mod_b2plot_gui
      implicit none
      integer which
      real (kind=R4) :: a, b
      if (.not.gui_active) return
      call gui_begin(34)
      call gui_put_i4(which)
      call gui_put_r4(a)
      call gui_put_r4(b)
      call gui_end
      return
      end subroutine gui_chattr

      subroutine gsplci(coli)
      use ncar_gui_state
      implicit none
      integer coli
      cur_plci = coli
      call gui_attr_i(21,coli)
      return
      end subroutine gsplci

      subroutine gslwsc(lwidth)
      use b2mod_types
      implicit none
      real (kind=R4) :: lwidth
      call gui_attr_r(22,lwidth)
      return
      end subroutine gslwsc

      subroutine gsln(ltype)
      implicit none
      integer ltype
      call gui_attr_i(23,ltype)
      return
      end subroutine gsln

c While probing the fill callback only the colour is captured
      subroutine gsfaci(coli)
      use ncar_gui_state
      implicit none
      integer coli
      if (probing) then
        probe_color = coli
        return
      end if
      call gui_attr_i(25,coli)
      return
      end subroutine gsfaci

      subroutine gsfais(ints)
      implicit none
      integer ints
      call gui_attr_i(26,ints)
      return
      end subroutine gsfais

      subroutine gspmci(coli)
      implicit none
      integer coli
      call gui_attr_i(27,coli)
      return
      end subroutine gspmci

      subroutine gsmk(mtype)
      implicit none
      integer mtype
      call gui_attr_i(28,mtype)
      return
      end subroutine gsmk

      subroutine gsmksc(mszsf)
      use b2mod_types
      implicit none
      real (kind=R4) :: mszsf
      call gui_attr_r(29,mszsf)
      return
      end subroutine gsmksc

      subroutine gstxfp(font,prec)
      use b2mod_b2plot_gui
      implicit none
      integer font, prec
      if (.not.gui_active) return
      call gui_begin(30)
      call gui_put_i4(font)
      call gui_put_i4(prec)
      call gui_end
      return
      end subroutine gstxfp

      subroutine gschh(chh)
      use b2mod_types
      implicit none
      real (kind=R4) :: chh
      call gui_chattr(1,chh,0.0_R4)
      return
      end subroutine gschh

      subroutine gschsp(chsp)
      use b2mod_types
      implicit none
      real (kind=R4) :: chsp
      call gui_chattr(2,chsp,0.0_R4)
      return
      end subroutine gschsp

      subroutine gschxp(chxp)
      use b2mod_types
      implicit none
      real (kind=R4) :: chxp
      call gui_chattr(3,chxp,0.0_R4)
      return
      end subroutine gschxp

      subroutine gschup(chux,chuy)
      use b2mod_types
      implicit none
      real (kind=R4) :: chux, chuy
      call gui_chattr(4,chux,chuy)
      return
      end subroutine gschup

c ======================= Lines, markers, fills ========================
      subroutine dashdc(ipat,jcrt,jsize)
      use b2mod_b2plot_gui
      implicit none
      character*(*) ipat
      integer jcrt, jsize
      if (.not.gui_active) return
      call gui_begin(24)
      call gui_put_str(ipat)
      call gui_put_i4(jcrt)
      call gui_put_i4(jsize)
      call gui_end
      return
      end subroutine dashdc

c Argument order (xa,ya,xb,yb) as used by b2plot
      subroutine lined(xa,ya,xb,yb)
      use b2mod_b2plot_gui
      implicit none
      real (kind=R4) :: xa, ya, xb, yb
      if (.not.gui_active) return
      call gui_begin(40)
      call gui_put_r4(xa)
      call gui_put_r4(ya)
      call gui_put_r4(xb)
      call gui_put_r4(yb)
      call gui_end
      return
      end subroutine lined

      subroutine curved(x,y,n)
      use b2mod_b2plot_gui
      use ncar_gui_state
      implicit none
      integer n
      real (kind=R4) :: x(*), y(*)
      if (.not.gui_active) return
      call gui_begin(41)
      call gui_points(n,x,y)
      call gui_end
      return
      end subroutine curved

      subroutine gpl(n,px,py)
      use b2mod_b2plot_gui
      use ncar_gui_state
      implicit none
      integer n
      real (kind=R4) :: px(*), py(*)
      if (.not.gui_active) return
      call gui_begin(41)
      call gui_points(n,px,py)
      call gui_end
      return
      end subroutine gpl

      subroutine lastd
      use b2mod_b2plot_gui
      implicit none
      if (.not.gui_active) return
      call gui_begin(42)
      call gui_end
      return
      end subroutine lastd

      subroutine reset
      implicit none
      return
      end subroutine reset

      subroutine gpm(n,px,py)
      use b2mod_b2plot_gui
      use ncar_gui_state
      implicit none
      integer n
      real (kind=R4) :: px(*), py(*)
      if (.not.gui_active) return
      call gui_begin(44)
      call gui_points(n,px,py)
      call gui_end
      return
      end subroutine gpm

c FILL message with colour ci (-1 current) and the pending tag
      subroutine gui_fill(n,x,y,ci)
      use b2mod_b2plot_gui
      use ncar_gui_state
      implicit none
      integer n, ci
      real (kind=R4) :: x(*), y(*)
      if (gui_active) then
        call gui_begin(43)
        call gui_points(n,x,y)
        call gui_put_i4(ci)
        call gui_put_i4(gui_tag_kind)
        call gui_put_i4(gui_tag_index)
        call gui_put_r8(gui_tag_value)
        call gui_end
      end if
      call gui_tag(0,0,0.0_R8)
      return
      end subroutine gui_fill

      subroutine gfa(n,px,py)
      use ncar_gui_state
      implicit none
      integer n
      real (kind=R4) :: px(*), py(*)
      if (probing) return
      call gui_fill(n,px,py,-1)
      return
      end subroutine gfa

      subroutine sfsgfa(xra,yra,nra,dst,nst,ind,nnd,ici)
      use ncar_gui_state
      implicit none
      integer nra, nst, ind(*), nnd, ici
      real (kind=R4) :: xra(*), yra(*), dst(*)
      if (.false.) then
        idum = nst+nnd+ind(1)
        rdum = dst(1)
      end if
      call gui_fill(nra,xra,yra,ici)
      return
      end subroutine sfsgfa

      subroutine sfseti(cnp,ivp)
      use b2mod_b2plot_gui
      implicit none
      character*(*) cnp
      integer ivp
      if (.not.gui_active) return
      call gui_begin(33)
      call gui_put_str(cnp)
      call gui_put_i4(ivp)
      call gui_end
      return
      end subroutine sfseti

      subroutine gca(xp1,yp1,xp2,yp2,nxclr,nyclr,ixclr,iyclr,nxdo,nydo,
     . clria)
      use b2mod_b2plot_gui
      implicit none
      real (kind=R4) :: xp1, yp1, xp2, yp2
      integer nxclr, nyclr, ixclr, iyclr, nxdo, nydo, clria(*)
      if (.not.gui_active) return
      call gui_begin(49)
      call gui_put_r4(xp1)
      call gui_put_r4(yp1)
      call gui_put_r4(xp2)
      call gui_put_r4(yp2)
      call gui_put_i4(nxclr)
      call gui_put_i4(nyclr)
      call gui_put_i4(ixclr)
      call gui_put_i4(iyclr)
      call gui_put_i4(nxdo)
      call gui_put_i4(nydo)
      call gui_put_i4v(nxclr*nyclr,clria)
      call gui_end
      return
      end subroutine gca

c ============================== Text ==================================
      subroutine plchhq(xpos,ypos,chrs,size,angd,cntr)
      use b2mod_b2plot_gui
      implicit none
      real (kind=R4) :: xpos, ypos, size, angd, cntr
      character*(*) chrs
      if (.not.gui_active) return
      call gui_begin(45)
      call gui_put_r4(xpos)
      call gui_put_r4(ypos)
      call gui_put_r4(size)
      call gui_put_r4(angd)
      call gui_put_r4(cntr)
      call gui_put_str(chrs)
      call gui_end
      return
      end subroutine plchhq

      subroutine gtx(px,py,chars)
      use b2mod_b2plot_gui
      implicit none
      real (kind=R4) :: px, py
      character*(*) chars
      if (.not.gui_active) return
      call gui_begin(46)
      call gui_put_r4(px)
      call gui_put_r4(py)
      call gui_put_str(chars)
      call gui_end
      return
      end subroutine gtx

c Plotchar (1 int, 3 char) and Gridall (1 int, 2 real) parameters
      subroutine gui_named_param(op,pnam,kind,ival,rval,cval)
      use b2mod_b2plot_gui
      implicit none
      integer op, kind, ival
      real (kind=R4) :: rval
      character*(*) pnam, cval
      if (.not.gui_active) return
      call gui_begin(op)
      call gui_put_str(pnam)
      call gui_put_i4(kind)
      call gui_put_i4(ival)
      if (op.eq.31) then
        call gui_put_str(cval)
      else
        call gui_put_r4(rval)
      end if
      call gui_end
      return
      end subroutine gui_named_param

      subroutine pcseti(pnam,ival)
      use b2mod_types
      implicit none
      character*(*) pnam
      integer ival
      call gui_named_param(31,pnam,1,ival,0.0_R4,' ')
      return
      end subroutine pcseti

      subroutine pcsetc(pnam,cval)
      use b2mod_types
      implicit none
      character*(*) pnam, cval
      call gui_named_param(31,pnam,3,0,0.0_R4,cval)
      return
      end subroutine pcsetc

      subroutine gaseti(pnam,ival)
      use b2mod_types
      implicit none
      character*(*) pnam
      integer ival
      call gui_named_param(32,pnam,1,ival,0.0_R4,' ')
      return
      end subroutine gaseti

      subroutine gasetr(pnam,rval)
      use b2mod_types
      implicit none
      character*(*) pnam
      real (kind=R4) :: rval
      call gui_named_param(32,pnam,2,0,rval,' ')
      return
      end subroutine gasetr

c ======================= Axes and colour bars =========================
      subroutine labmod(fmtx,fmty,numx,numy,iszx,iszy,ixdc,iydc,ixor)
      use ncar_gui_state
      implicit none
      character*(*) fmtx, fmty
      integer numx, numy, iszx, iszy, ixdc, iydc, ixor
      integer i, j
c Formats are read up to their closing parenthesis, as in NCAR
      i = index(fmtx,')')
      if (i.eq.0) i = len(fmtx)
      j = index(fmty,')')
      if (j.eq.0) j = len(fmty)
      lm_fmtx = fmtx(1:i)
      lm_fmty = fmty(1:j)
      lm_numx = numx
      lm_numy = numy
      lm_iszx = iszx
      lm_iszy = iszy
      lm_ixdc = ixdc
      lm_iydc = iydc
      lm_ixor = ixor
      return
      end subroutine labmod

      subroutine gridal(mjrx,mnrx,mjry,mnry,ixlb,iylb,igph,xint,yint)
      use b2mod_b2plot_gui
      use ncar_gui_state
      implicit none
      integer mjrx, mnrx, mjry, mnry, ixlb, iylb, igph
      real (kind=R4) :: xint, yint
      if (.not.gui_active) return
      call gui_begin(47)
      call gui_put_i4(mjrx)
      call gui_put_i4(mnrx)
      call gui_put_i4(mjry)
      call gui_put_i4(mnry)
      call gui_put_i4(ixlb)
      call gui_put_i4(iylb)
      call gui_put_i4(igph)
      call gui_put_r4(xint)
      call gui_put_r4(yint)
      call gui_put_str(lm_fmtx)
      call gui_put_str(lm_fmty)
      call gui_put_i4(lm_numx)
      call gui_put_i4(lm_numy)
      call gui_put_i4(lm_iszx)
      call gui_put_i4(lm_iszy)
      call gui_put_i4(lm_ixdc)
      call gui_put_i4(lm_iydc)
      call gui_put_i4(lm_ixor)
      call gui_end
      return
      end subroutine gridal

      subroutine lblbar(ihov,xleb,xreb,ybeb,yteb,nbox,wsfb,hsfb,lfin,
     . iftp,llbs,nlbs,lbab)
      use b2mod_b2plot_gui
      implicit none
      integer ihov, nbox, lfin(*), iftp, nlbs, lbab
      real (kind=R4) :: xleb, xreb, ybeb, yteb, wsfb, hsfb
      character*(*) llbs(*)
      integer i
      if (.not.gui_active) return
      call gui_begin(48)
      call gui_put_i4(ihov)
      call gui_put_r4(xleb)
      call gui_put_r4(xreb)
      call gui_put_r4(ybeb)
      call gui_put_r4(yteb)
      call gui_put_i4(nbox)
      call gui_put_r4(wsfb)
      call gui_put_r4(hsfb)
      call gui_put_i4v(nbox,lfin)
      call gui_put_i4(iftp)
      call gui_put_i4(max(nlbs,0))
      do i = 1, nlbs
        call gui_put_str(llbs(i))
      end do
      call gui_put_i4(lbab)
      call gui_end
      return
      end subroutine lblbar

c ===================== Conpack and Conpackt ===========================
      subroutine cpseti(pnam,ival)
      use ncar_gui_state
      implicit none
      character*(*) pnam
      integer ival
      call gui_pset(gui_cp,pnam,1,ival,0.0_R4,' ')
      return
      end subroutine cpseti

      subroutine cpsetr(pnam,rval)
      use ncar_gui_state
      implicit none
      character*(*) pnam
      real (kind=R4) :: rval
      call gui_pset(gui_cp,pnam,2,0,rval,' ')
      return
      end subroutine cpsetr

      subroutine cpsetc(pnam,cval)
      use ncar_gui_state
      implicit none
      character*(*) pnam, cval
      call gui_pset(gui_cp,pnam,3,0,0.0_R4,cval)
      return
      end subroutine cpsetc

      subroutine cpgetr(pnam,rval)
      use ncar_gui_state
      implicit none
      character*(*) pnam
      real (kind=R4) :: rval
      rval = gui_pget(gui_cp,pnam)
      return
      end subroutine cpgetr

      subroutine ctseti(pnam,ival)
      use ncar_gui_state
      implicit none
      character*(*) pnam
      integer ival
      call gui_pset(gui_ct,pnam,1,ival,0.0_R4,' ')
      return
      end subroutine ctseti

      subroutine ctsetr(pnam,rval)
      use ncar_gui_state
      implicit none
      character*(*) pnam
      real (kind=R4) :: rval
      call gui_pset(gui_ct,pnam,2,0,rval,' ')
      return
      end subroutine ctsetr

      subroutine ctsetc(pnam,cval)
      use ncar_gui_state
      implicit none
      character*(*) pnam, cval
      call gui_pset(gui_ct,pnam,3,0,0.0_R4,cval)
      return
      end subroutine ctsetc

      subroutine ctmesh(rpnt,npnt,lopn,iedg,nedg,loen,itri,ntri,lotn,
     . rwrk,lrwk,iwrk,liwk)
      use b2mod_b2plot_gui
      use ncar_gui_state
      implicit none
      integer npnt, lopn, iedg(*), nedg, loen, itri(*), ntri, lotn,
     . lrwk, iwrk(*), liwk
      real (kind=R4) :: rpnt(*), rwrk(*)
      if (.false.) then
        idum = lrwk+liwk+iwrk(1)
        rdum = rwrk(1)
      end if
      current_source = gui_ct
      if (.not.gui_active) return
      call gui_begin(50)
      call gui_put_i4(npnt)
      call gui_put_i4(lopn)
      call gui_put_i4(nedg)
      call gui_put_i4(loen)
      call gui_put_i4(ntri)
      call gui_put_i4(lotn)
      call gui_put_r4v(npnt,rpnt)
      call gui_put_i4v(nedg,iedg)
      call gui_put_i4v(ntri,itri)
      call gui_end
      return
      end subroutine ctmesh

      subroutine cprect(zdat,kzdt,mzdt,nzdt,rwrk,lrwk,iwrk,liwk)
      use b2mod_b2plot_gui
      use ncar_gui_state
      implicit none
      integer kzdt, mzdt, nzdt, lrwk, iwrk(*), liwk
      real (kind=R4) :: zdat(*), rwrk(*)
      if (.false.) then
        idum = lrwk+liwk+iwrk(1)
        rdum = rwrk(1)
      end if
      current_source = gui_cp
      if (.not.gui_active) return
      call gui_begin(51)
      call gui_put_i4(kzdt)
      call gui_put_i4(mzdt)
      call gui_put_i4(nzdt)
      call gui_put_r4v(kzdt*nzdt,zdat)
      call gui_end
      return
      end subroutine cprect

c Contour lines and labels are drawn by Python from the tables
      subroutine cpclam(zdat,rwrk,iwrk,iama)
      use ncar_gui_state
      implicit none
      real (kind=R4) :: zdat(*), rwrk(*)
      integer iwrk(*), iama(*)
      if (.false.) then
        idum = iwrk(1)+iama(1)
        rdum = zdat(1)+rwrk(1)
      end if
      return
      end subroutine cpclam

      subroutine ctclam(rpnt,iedg,itri,rwrk,iwrk,iama)
      use ncar_gui_state
      implicit none
      real (kind=R4) :: rpnt(*), rwrk(*)
      integer iedg(*), itri(*), iwrk(*), iama(*)
      if (.false.) then
        idum = iedg(1)+itri(1)+iwrk(1)+iama(1)
        rdum = rpnt(1)+rwrk(1)
      end if
      return
      end subroutine ctclam

      subroutine cplbam(zdat,rwrk,iwrk,iama)
      use ncar_gui_state
      implicit none
      real (kind=R4) :: zdat(*), rwrk(*)
      integer iwrk(*), iama(*)
      if (.false.) then
        idum = iwrk(1)+iama(1)
        rdum = zdat(1)+rwrk(1)
      end if
      return
      end subroutine cplbam

      subroutine ctlbam(rpnt,iedg,itri,rwrk,iwrk,iama)
      use ncar_gui_state
      implicit none
      real (kind=R4) :: rpnt(*), rwrk(*)
      integer iedg(*), itri(*), iwrk(*), iama(*)
      if (.false.) then
        idum = iedg(1)+itri(1)+iwrk(1)+iama(1)
        rdum = rpnt(1)+rwrk(1)
      end if
      return
      end subroutine ctlbam

      subroutine cpcldm(zdat,rwrk,iwrk,iama,rtpl)
      use ncar_gui_state
      implicit none
      real (kind=R4) :: zdat(*), rwrk(*)
      integer iwrk(*), iama(*)
      external rtpl
      if (.false.) then
        idum = iwrk(1)+iama(1)
        rdum = zdat(1)+rwrk(1)
        call rtpl
      end if
      call gui_contour_lines(gui_cp)
      return
      end subroutine cpcldm

      subroutine ctcldm(rpnt,iedg,itri,rwrk,iwrk,iama,rtpl)
      use ncar_gui_state
      implicit none
      real (kind=R4) :: rpnt(*), rwrk(*)
      integer iedg(*), itri(*), iwrk(*), iama(*)
      external rtpl
      if (.false.) then
        idum = iedg(1)+itri(1)+iwrk(1)+iama(1)
        rdum = rpnt(1)+rwrk(1)
        call rtpl
      end if
      call gui_contour_lines(gui_ct)
      return
      end subroutine ctcldm

      subroutine cplbdr(zdat,rwrk,iwrk)
      use ncar_gui_state
      implicit none
      real (kind=R4) :: zdat(*), rwrk(*)
      integer iwrk(*)
      if (.false.) then
        idum = iwrk(1)
        rdum = zdat(1)+rwrk(1)
      end if
      call gui_contour_info(gui_cp)
      return
      end subroutine cplbdr

      subroutine ctlbdr(rpnt,iedg,itri,rwrk,iwrk)
      use ncar_gui_state
      implicit none
      real (kind=R4) :: rpnt(*), rwrk(*)
      integer iedg(*), itri(*), iwrk(*)
      if (.false.) then
        idum = iedg(1)+itri(1)+iwrk(1)
        rdum = rpnt(1)+rwrk(1)
      end if
      call gui_contour_info(gui_ct)
      return
      end subroutine ctlbdr

      subroutine cpdrpl(xcs,ycs,ncs,iai,iag,nai)
      use ncar_gui_state
      implicit none
      real (kind=R4) :: xcs(*), ycs(*)
      integer ncs, iai(*), iag(*), nai
      if (.false.) then
        idum = ncs+iai(1)+iag(1)+nai
        rdum = xcs(1)+ycs(1)
      end if
      return
      end subroutine cpdrpl

      subroutine cpcnrc(zdat,kzdt,mzdt,nzdt,flow,fhgh,finc,kset,nhgh,
     . ndsh)
      use b2mod_b2plot_gui
      use ncar_gui_state
      implicit none
      real (kind=R4) :: zdat(*), flow, fhgh, finc
      integer kzdt, mzdt, nzdt, kset, nhgh, ndsh
      if (.false.) then
        idum = kzdt+mzdt+nzdt+kset+nhgh+ndsh
        rdum = zdat(1)+flow+fhgh+finc
      end if
      call gui_message(2, 'ncar_gui: cpcnrc not supported')
      return
      end subroutine cpcnrc

c ================================ Areas ===============================
      subroutine arinam(iama,lama)
      use b2mod_b2plot_gui
      use ncar_gui_state
      implicit none
      integer iama(*), lama
      if (.false.) idum = iama(1)+lama
      if (.not.gui_active) return
      call gui_begin(56)
      call gui_end
      return
      end subroutine arinam

      subroutine aredam(iama,xca,yca,lca,igi,idl,idr)
      use b2mod_b2plot_gui
      use ncar_gui_state
      implicit none
      integer iama(*), lca, igi, idl, idr
      real (kind=R4) :: xca(*), yca(*)
      if (.false.) idum = iama(1)
      if (.not.gui_active) return
      call gui_begin(52)
      call gui_put_i4(igi)
      call gui_put_i4(idl)
      call gui_put_i4(idr)
      call gui_points(lca,xca,yca)
      call gui_end
      return
      end subroutine aredam

      subroutine arpram(iama,if1,if2,if3)
      use ncar_gui_state
      implicit none
      integer iama(*), if1, if2, if3
      if (.false.) idum = iama(1)+if1+if2+if3
      return
      end subroutine arpram

c Contour fill: the colour of each area identifier is obtained by
c calling the user fill routine apr in probing mode
      subroutine arscam(iama,xcs,ycs,mcs,iai,iag,mai,apr)
      use b2mod_b2plot_gui
      use ncar_gui_state
      implicit none
      integer iama(*), mcs, iai(*), iag(*), mai
      real (kind=R4) :: xcs(*), ycs(*)
      external apr
      integer it, k, n, nid, iarea(2), igrp(2), icol(gui_nlev+1)
      real (kind=R4) :: xp(4), yp(4)
      if (.false.) then
        idum = iama(1)+mcs+iai(1)+iag(1)+mai
        rdum = xcs(1)+ycs(1)
      end if
      if (.not.gui_active) return
      it = current_source
      if (it.ne.gui_ct) it = gui_cp
      n = tncl(it)
      nid = n+1
      xp = 0.0_R4
      yp = 0.0_R4
      igrp(1) = 1
      igrp(2) = 2
      do k = 1, nid
        iarea(1) = k
        iarea(2) = 1
        probe_color = -1
        probing = .true.
        call apr(xp,yp,4,iarea,igrp,2)
        probing = .false.
        icol(k) = probe_color
      end do
      call gui_begin(53)
      call gui_put_i4(it)
      call gui_put_i4(n)
      call gui_put_r4v(n,tclv(1,it))
      call gui_put_i4v(n,taia(1,it))
      call gui_put_i4v(n,taib(1,it))
      call gui_put_i4(nid)
      call gui_put_i4v(nid,icol)
      call gui_put_r4(txc1(it))
      call gui_put_r4(txcm(it))
      call gui_put_r4(tyc1(it))
      call gui_put_r4(tycn(it))
      call gui_put_r4(tspv(it))
      call gui_put_i4(tmap(it))
      call gui_end
      return
      end subroutine arscam

      subroutine ardbpa(iama,igip,label)
      use ncar_gui_state
      implicit none
      integer iama(*), igip
      character*(*) label
      if (.false.) idum = iama(1)+igip+len(label)
      return
      end subroutine ardbpa

c ========================== Error handling ============================
      subroutine entsr(irold,irnew)
      use ncar_gui_state
      implicit none
      integer irold, irnew
      if (.false.) idum = irnew
      irold = 0
      return
      end subroutine entsr

      subroutine retsr(irold)
      use ncar_gui_state
      implicit none
      integer irold
      if (.false.) idum = irold
      return
      end subroutine retsr

      subroutine errof
      implicit none
      return
      end subroutine errof

      subroutine eprin
      implicit none
      return
      end subroutine eprin

      integer function nerro(nerrf)
      implicit none
      integer nerrf
      nerrf = 0
      nerro = 0
      return
      end function nerro

c =================== GFLASH segments (trajectories) ===================
      subroutine gflas1(ib)
      implicit none
      integer ib
      call gui_attr_i(60,ib)
      return
      end subroutine gflas1

      subroutine gflas2
      use b2mod_b2plot_gui
      implicit none
      if (.not.gui_active) return
      call gui_begin(61)
      call gui_end
      return
      end subroutine gflas2

      subroutine gflas3(ib)
      implicit none
      integer ib
      call gui_attr_i(62,ib)
      return
      end subroutine gflas3

c ======================= Colour conversion ============================
c HLS to RGB with the NCAR convention: hue 0 is blue, L and S in percent
      subroutine hlsrgb(h,l,s,r,g,b)
      use b2mod_types
      implicit none
      real (kind=R4) :: h, l, s, r, g, b
      real (kind=R4) :: xl, xs, m1, m2, gui_hls_value
      xl = l/100.0_R4
      xs = s/100.0_R4
      if (xs.le.0.0_R4) then
        r = xl
        g = xl
        b = xl
        return
      end if
      if (xl.le.0.5_R4) then
        m2 = xl*(1.0_R4+xs)
      else
        m2 = xl+xs-xl*xs
      end if
      m1 = 2.0_R4*xl-m2
      r = gui_hls_value(m1,m2,h)
      g = gui_hls_value(m1,m2,h-120.0_R4)
      b = gui_hls_value(m1,m2,h-240.0_R4)
      return
      end subroutine hlsrgb

      real (kind=R4) function gui_hls_value(m1,m2,hue)
      use b2mod_types
      implicit none
      real (kind=R4) :: m1, m2, hue
      real (kind=R4) :: hh
      hh = modulo(hue,360.0_R4)
      if (hh.lt.60.0_R4) then
        gui_hls_value = m1+(m2-m1)*hh/60.0_R4
      else if (hh.lt.180.0_R4) then
        gui_hls_value = m2
      else if (hh.lt.240.0_R4) then
        gui_hls_value = m1+(m2-m1)*(240.0_R4-hh)/60.0_R4
      else
        gui_hls_value = m1
      end if
      return
      end function gui_hls_value
