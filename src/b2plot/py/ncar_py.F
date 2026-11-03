c NCAR Graphics replacement for b2pl_py.exe: each NCAR and GKS entry
c point used by b2plot sends a message on the b2plot.py FIFO instead of
c drawing (protocol in the b2plot.py documentation of the manual)

      module ncar_py_state
      use b2mod_types
      implicit none
c Parameter tables: 1 Conpackt, 2 Conpack
      integer, parameter :: ncpy_ct = 1, ncpy_cp = 2, ncpy_nlev = 256
      integer, save :: tpai(2) = 1, tncl(2) = 0, tcls(2) = 0,
     . tllp(2) = 0, tllo(2) = 0, tmap(2) = 0, tset(2) = 1,
     . tpic(2) = 0, tgic(2) = 3, tgil(2) = 3, tgis(2) = 4,
     . tnvs(2) = 1, tilp(2) = 0
      integer, save :: tclu(ncpy_nlev,2) = 0, taia(ncpy_nlev,2) = 0,
     . taib(ncpy_nlev,2) = 0
      real (kind=R4), save :: tclv(ncpy_nlev,2) = 0.0_R4,
     . tcll(ncpy_nlev,2) = 0.0_R4, tt2d(2) = 0.0_R4,
     . tspv(2) = 0.0_R4, txc1(2) = 0.0_R4, txcm(2) = 0.0_R4,
     . tyc1(2) = 0.0_R4, tycn(2) = 0.0_R4, tilx(2) = 0.0_R4,
     . tily(2) = 0.0_R4
      character*256, save :: tilt(2) = ' ', tcft(2) = ' '
c Contour source of the current plot (ncpy_ct or ncpy_cp)
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
c Sink for NCAR arguments b2plot.py does not need
      integer, save :: idum = 0
      real (kind=R4), save :: rdum = 0.0_R4

      contains

c Lowest open workstation, 0 if none is open
      integer function ncpy_lowest_wk()
      implicit none
      if (nwkopen.eq.0) then
        ncpy_lowest_wk = 0
      else
        ncpy_lowest_wk = minval(iwkopen(1:nwkopen))
      end if
      return
      end function ncpy_lowest_wk

c First three characters of a parameter name, in upper case
      character*3 function ncpy_pkey(pnam)
      implicit none
      character*(*), intent(in) :: pnam
      integer i, ic
      ncpy_pkey = pnam
      do i = 1, 3
        ic = iachar(ncpy_pkey(i:i))
        if (ic.ge.iachar('a').and.ic.le.iachar('z'))
     .   ncpy_pkey(i:i) = achar(ic-32)
      end do
      return
      end function ncpy_pkey

c Report a parameter name that the tables do not know, once
      subroutine ncpy_unknown(key)
      use b2mod_b2plot_py
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
      write(*,'(a)') 'ncar_py: parameter '//key//' not emulated'
      call b2py_message(2, 'ncar_py: parameter '//key//' not emulated')
      return
      end subroutine ncpy_unknown

c Store a Conpack or Conpackt parameter (kind 1 int, 2 real, 3 char)
      subroutine ncpy_pset(it,pnam,kind,ival,rval,cval)
      implicit none
      integer, intent(in) :: it, kind, ival
      real (kind=R4), intent(in) :: rval
      character*(*), intent(in) :: pnam, cval
      character*3 key
      integer iv, ip
      real (kind=R4) :: rv
      logical lev
      key = ncpy_pkey(pnam)
      if (kind.eq.2) then
        iv = nint(rval)
        rv = rval
      else
        iv = ival
        rv = real(ival,R4)
      end if
      ip = tpai(it)
      lev = ip.ge.1.and.ip.le.ncpy_nlev
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
        tncl(it) = max(0,min(iv,ncpy_nlev))
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
c numeric label format, not needed by b2plot.py
      case default
        call ncpy_unknown(key)
      end select
      return
      end subroutine ncpy_pset

c Return a real-valued parameter of table it
      real (kind=R4) function ncpy_pget(it,pnam)
      implicit none
      integer, intent(in) :: it
      character*(*), intent(in) :: pnam
      integer ip
      ip = max(1,min(tpai(it),ncpy_nlev))
      select case (ncpy_pkey(pnam))
      case ('CLV')
        ncpy_pget = tclv(ip,it)
      case ('CLL')
        ncpy_pget = tcll(ip,it)
      case ('T2D')
        ncpy_pget = tt2d(it)
      case ('SPV')
        ncpy_pget = tspv(it)
      case ('XC1')
        ncpy_pget = txc1(it)
      case ('XCM')
        ncpy_pget = txcm(it)
      case ('YC1')
        ncpy_pget = tyc1(it)
      case ('YCN')
        ncpy_pget = tycn(it)
      case ('ILX')
        ncpy_pget = tilx(it)
      case ('ILY')
        ncpy_pget = tily(it)
      case ('NCL')
        ncpy_pget = real(tncl(it),R4)
      case default
        ncpy_pget = 0.0_R4
        call ncpy_unknown(ncpy_pkey(pnam))
      end select
      return
      end function ncpy_pget

c Send CONTOUR_LINES or CONTOUR_INFO for table it
      subroutine ncpy_contour_lines(it)
      use b2mod_b2plot_py
      implicit none
      integer, intent(in) :: it
      integer n
      if (.not.b2py_active) return
      n = tncl(it)
      call b2py_begin(54)
      call b2py_put_i4(it)
      call b2py_put_i4(n)
      call b2py_put_r4v(n,tclv(1,it))
      call b2py_put_i4v(n,tclu(1,it))
      call b2py_put_r4v(n,tcll(1,it))
      call b2py_put_i4(cur_plci)
      call b2py_put_i4(tllp(it))
      call b2py_put_i4(tllo(it))
      call b2py_put_r4(txc1(it))
      call b2py_put_r4(txcm(it))
      call b2py_put_r4(tyc1(it))
      call b2py_put_r4(tycn(it))
      call b2py_put_r4(tspv(it))
      call b2py_put_i4(tmap(it))
      call b2py_end
      return
      end subroutine ncpy_contour_lines

      subroutine ncpy_contour_info(it)
      use b2mod_b2plot_py
      implicit none
      integer, intent(in) :: it
      if (.not.b2py_active) return
      call b2py_begin(55)
      call b2py_put_i4(it)
      call b2py_put_str(tilt(it))
      call b2py_put_r4(tilx(it))
      call b2py_put_r4(tily(it))
      call b2py_put_i4(tilp(it))
      call b2py_put_i4(tllp(it))
      call b2py_put_i4(tllo(it))
      call b2py_end
      return
      end subroutine ncpy_contour_info

c Send a polyline or polygon of n points
      subroutine ncpy_points(n,x,y)
      use b2mod_b2plot_py
      implicit none
      integer, intent(in) :: n
      real (kind=R4), intent(in) :: x(*), y(*)
      call b2py_put_i4(max(n,0))
      call b2py_put_r4v(n,x)
      call b2py_put_r4v(n,y)
      return
      end subroutine ncpy_points

      end module ncar_py_state

c ===================== GKS control and workstations ===================
      subroutine gopks(errfil,bufa)
      use ncar_py_state
      implicit none
      integer errfil, bufa
      if (.false.) idum = errfil+bufa
      nwkopen = 0
      return
      end subroutine gopks

      subroutine gopwk(wkid,conid,wtype)
      use ncar_py_state
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
      use ncar_py_state
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
      use ncar_py_state
      implicit none
      integer wkid
      if (.false.) idum = wkid
      return
      end subroutine gacwk

      subroutine gdawk(wkid)
      use ncar_py_state
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
      use ncar_py_state
      implicit none
      character*(*) pnam
      integer ival
      if (.false.) idum = len(pnam)+ival
      return
      end subroutine ngseti

      subroutine ngsetc(pnam,cval)
      use ncar_py_state
      implicit none
      character*(*) pnam, cval
      if (.false.) idum = len(pnam)+len(cval)
      return
      end subroutine ngsetc

c ======================== Frames and mapping ==========================
c FRAME only after something was drawn, as NCAR clears only used frames
      subroutine frame
      use b2mod_b2plot_py
      implicit none
      if (.not.b2py_active.or..not.b2py_drawn) return
      call b2py_begin(10)
      call b2py_end
      call b2py_flush
      b2py_drawn = .false.
      return
      end subroutine frame

c One FRAME per page: end_page calls this once per workstation
      subroutine ngpict(wkid,action)
      use b2mod_b2plot_py
      use ncar_py_state
      implicit none
      integer wkid, action
      if (.false.) idum = action
      if (.not.b2py_active.or..not.b2py_drawn) return
      if (wkid.ne.ncpy_lowest_wk()) return
      call b2py_begin(10)
      call b2py_end
      call b2py_flush
      b2py_drawn = .false.
      return
      end subroutine ngpict

      subroutine set(vl,vr,vb,vt,wl,wr,wb,wt,lf)
      use b2mod_b2plot_py
      implicit none
      real (kind=R4) :: vl, vr, vb, vt, wl, wr, wb, wt
      integer lf
      if (.not.b2py_active) return
      call b2py_begin(12)
      call b2py_put_r4(vl)
      call b2py_put_r4(vr)
      call b2py_put_r4(vb)
      call b2py_put_r4(vt)
      call b2py_put_r4(wl)
      call b2py_put_r4(wr)
      call b2py_put_r4(wb)
      call b2py_put_r4(wt)
      call b2py_put_i4(lf)
      call b2py_end
      return
      end subroutine set

c ===================== Colours and attributes =========================
      subroutine gscr(wkid,ci,cr,cg,cb)
      use b2mod_b2plot_py
      use ncar_py_state
      implicit none
      integer wkid, ci
      real (kind=R4) :: cr, cg, cb
      integer iw
      if (.not.b2py_active) return
      iw = ncpy_lowest_wk()
      if (iw.ne.0.and.wkid.ne.iw) return
      call b2py_begin(20)
      call b2py_put_i4(ci)
      call b2py_put_r4(cr)
      call b2py_put_r4(cg)
      call b2py_put_r4(cb)
      call b2py_end
      return
      end subroutine gscr

c Integer attribute message
      subroutine ncpy_attr_i(op,ival)
      use b2mod_b2plot_py
      implicit none
      integer op, ival
      if (.not.b2py_active) return
      call b2py_begin(op)
      call b2py_put_i4(ival)
      call b2py_end
      return
      end subroutine ncpy_attr_i

c Real attribute message
      subroutine ncpy_attr_r(op,rval)
      use b2mod_b2plot_py
      implicit none
      integer op
      real (kind=R4) :: rval
      if (.not.b2py_active) return
      call b2py_begin(op)
      call b2py_put_r4(rval)
      call b2py_end
      return
      end subroutine ncpy_attr_r

c Character attribute message (1 CHH, 2 CHSP, 3 CHXP, 4 CHUP)
      subroutine ncpy_chattr(which,a,b)
      use b2mod_b2plot_py
      implicit none
      integer which
      real (kind=R4) :: a, b
      if (.not.b2py_active) return
      call b2py_begin(34)
      call b2py_put_i4(which)
      call b2py_put_r4(a)
      call b2py_put_r4(b)
      call b2py_end
      return
      end subroutine ncpy_chattr

      subroutine gsplci(coli)
      use ncar_py_state
      implicit none
      integer coli
      cur_plci = coli
      call ncpy_attr_i(21,coli)
      return
      end subroutine gsplci

      subroutine gslwsc(lwidth)
      use b2mod_types
      implicit none
      real (kind=R4) :: lwidth
      call ncpy_attr_r(22,lwidth)
      return
      end subroutine gslwsc

      subroutine gsln(ltype)
      implicit none
      integer ltype
      call ncpy_attr_i(23,ltype)
      return
      end subroutine gsln

c While probing the fill callback only the colour is captured
      subroutine gsfaci(coli)
      use ncar_py_state
      implicit none
      integer coli
      if (probing) then
        probe_color = coli
        return
      end if
      call ncpy_attr_i(25,coli)
      return
      end subroutine gsfaci

      subroutine gsfais(ints)
      implicit none
      integer ints
      call ncpy_attr_i(26,ints)
      return
      end subroutine gsfais

      subroutine gspmci(coli)
      implicit none
      integer coli
      call ncpy_attr_i(27,coli)
      return
      end subroutine gspmci

      subroutine gsmk(mtype)
      implicit none
      integer mtype
      call ncpy_attr_i(28,mtype)
      return
      end subroutine gsmk

      subroutine gsmksc(mszsf)
      use b2mod_types
      implicit none
      real (kind=R4) :: mszsf
      call ncpy_attr_r(29,mszsf)
      return
      end subroutine gsmksc

      subroutine gstxfp(font,prec)
      use b2mod_b2plot_py
      implicit none
      integer font, prec
      if (.not.b2py_active) return
      call b2py_begin(30)
      call b2py_put_i4(font)
      call b2py_put_i4(prec)
      call b2py_end
      return
      end subroutine gstxfp

      subroutine gschh(chh)
      use b2mod_types
      implicit none
      real (kind=R4) :: chh
      call ncpy_chattr(1,chh,0.0_R4)
      return
      end subroutine gschh

      subroutine gschsp(chsp)
      use b2mod_types
      implicit none
      real (kind=R4) :: chsp
      call ncpy_chattr(2,chsp,0.0_R4)
      return
      end subroutine gschsp

      subroutine gschxp(chxp)
      use b2mod_types
      implicit none
      real (kind=R4) :: chxp
      call ncpy_chattr(3,chxp,0.0_R4)
      return
      end subroutine gschxp

      subroutine gschup(chux,chuy)
      use b2mod_types
      implicit none
      real (kind=R4) :: chux, chuy
      call ncpy_chattr(4,chux,chuy)
      return
      end subroutine gschup

c ======================= Lines, markers, fills ========================
      subroutine dashdc(ipat,jcrt,jsize)
      use b2mod_b2plot_py
      implicit none
      character*(*) ipat
      integer jcrt, jsize
      if (.not.b2py_active) return
      call b2py_begin(24)
      call b2py_put_str(ipat)
      call b2py_put_i4(jcrt)
      call b2py_put_i4(jsize)
      call b2py_end
      return
      end subroutine dashdc

c Argument order (xa,ya,xb,yb) as used by b2plot
      subroutine lined(xa,ya,xb,yb)
      use b2mod_b2plot_py
      implicit none
      real (kind=R4) :: xa, ya, xb, yb
      if (.not.b2py_active) return
      call b2py_begin(40)
      call b2py_put_r4(xa)
      call b2py_put_r4(ya)
      call b2py_put_r4(xb)
      call b2py_put_r4(yb)
      call b2py_end
      return
      end subroutine lined

      subroutine curved(x,y,n)
      use b2mod_b2plot_py
      use ncar_py_state
      implicit none
      integer n
      real (kind=R4) :: x(*), y(*)
      if (.not.b2py_active) return
      call b2py_begin(41)
      call ncpy_points(n,x,y)
      call b2py_end
      return
      end subroutine curved

      subroutine gpl(n,px,py)
      use b2mod_b2plot_py
      use ncar_py_state
      implicit none
      integer n
      real (kind=R4) :: px(*), py(*)
      if (.not.b2py_active) return
      call b2py_begin(41)
      call ncpy_points(n,px,py)
      call b2py_end
      return
      end subroutine gpl

      subroutine lastd
      use b2mod_b2plot_py
      implicit none
      if (.not.b2py_active) return
      call b2py_begin(42)
      call b2py_end
      return
      end subroutine lastd

      subroutine reset
      implicit none
      return
      end subroutine reset

      subroutine gpm(n,px,py)
      use b2mod_b2plot_py
      use ncar_py_state
      implicit none
      integer n
      real (kind=R4) :: px(*), py(*)
      if (.not.b2py_active) return
      call b2py_begin(44)
      call ncpy_points(n,px,py)
      call b2py_end
      return
      end subroutine gpm

c FILL message with colour ci (-1 current) and the pending tag
      subroutine ncpy_fill(n,x,y,ci)
      use b2mod_b2plot_py
      use ncar_py_state
      implicit none
      integer n, ci
      real (kind=R4) :: x(*), y(*)
      if (b2py_active) then
        call b2py_begin(43)
        call ncpy_points(n,x,y)
        call b2py_put_i4(ci)
        call b2py_put_i4(b2py_tag_kind)
        call b2py_put_i4(b2py_tag_index)
        call b2py_put_r8(b2py_tag_value)
        call b2py_end
      end if
      call b2py_tag(0,0,0.0_R8)
      return
      end subroutine ncpy_fill

      subroutine gfa(n,px,py)
      use ncar_py_state
      implicit none
      integer n
      real (kind=R4) :: px(*), py(*)
      if (probing) return
      call ncpy_fill(n,px,py,-1)
      return
      end subroutine gfa

      subroutine sfsgfa(xra,yra,nra,dst,nst,ind,nnd,ici)
      use ncar_py_state
      implicit none
      integer nra, nst, ind(*), nnd, ici
      real (kind=R4) :: xra(*), yra(*), dst(*)
      if (.false.) then
        idum = nst+nnd+ind(1)
        rdum = dst(1)
      end if
      call ncpy_fill(nra,xra,yra,ici)
      return
      end subroutine sfsgfa

      subroutine sfseti(cnp,ivp)
      use b2mod_b2plot_py
      implicit none
      character*(*) cnp
      integer ivp
      if (.not.b2py_active) return
      call b2py_begin(33)
      call b2py_put_str(cnp)
      call b2py_put_i4(ivp)
      call b2py_end
      return
      end subroutine sfseti

      subroutine gca(xp1,yp1,xp2,yp2,nxclr,nyclr,ixclr,iyclr,nxdo,nydo,
     . clria)
      use b2mod_b2plot_py
      implicit none
      real (kind=R4) :: xp1, yp1, xp2, yp2
      integer nxclr, nyclr, ixclr, iyclr, nxdo, nydo, clria(*)
      if (.not.b2py_active) return
      call b2py_begin(49)
      call b2py_put_r4(xp1)
      call b2py_put_r4(yp1)
      call b2py_put_r4(xp2)
      call b2py_put_r4(yp2)
      call b2py_put_i4(nxclr)
      call b2py_put_i4(nyclr)
      call b2py_put_i4(ixclr)
      call b2py_put_i4(iyclr)
      call b2py_put_i4(nxdo)
      call b2py_put_i4(nydo)
      call b2py_put_i4v(nxclr*nyclr,clria)
      call b2py_end
      return
      end subroutine gca

c ============================== Text ==================================
      subroutine plchhq(xpos,ypos,chrs,size,angd,cntr)
      use b2mod_b2plot_py
      implicit none
      real (kind=R4) :: xpos, ypos, size, angd, cntr
      character*(*) chrs
      if (.not.b2py_active) return
      call b2py_begin(45)
      call b2py_put_r4(xpos)
      call b2py_put_r4(ypos)
      call b2py_put_r4(size)
      call b2py_put_r4(angd)
      call b2py_put_r4(cntr)
      call b2py_put_str(chrs)
      call b2py_end
      return
      end subroutine plchhq

      subroutine gtx(px,py,chars)
      use b2mod_b2plot_py
      implicit none
      real (kind=R4) :: px, py
      character*(*) chars
      if (.not.b2py_active) return
      call b2py_begin(46)
      call b2py_put_r4(px)
      call b2py_put_r4(py)
      call b2py_put_str(chars)
      call b2py_end
      return
      end subroutine gtx

c Plotchar (1 int, 3 char) and Gridall (1 int, 2 real) parameters
      subroutine ncpy_named_param(op,pnam,kind,ival,rval,cval)
      use b2mod_b2plot_py
      implicit none
      integer op, kind, ival
      real (kind=R4) :: rval
      character*(*) pnam, cval
      if (.not.b2py_active) return
      call b2py_begin(op)
      call b2py_put_str(pnam)
      call b2py_put_i4(kind)
      call b2py_put_i4(ival)
      if (op.eq.31) then
        call b2py_put_str(cval)
      else
        call b2py_put_r4(rval)
      end if
      call b2py_end
      return
      end subroutine ncpy_named_param

      subroutine pcseti(pnam,ival)
      use b2mod_types
      implicit none
      character*(*) pnam
      integer ival
      call ncpy_named_param(31,pnam,1,ival,0.0_R4,' ')
      return
      end subroutine pcseti

      subroutine pcsetc(pnam,cval)
      use b2mod_types
      implicit none
      character*(*) pnam, cval
      call ncpy_named_param(31,pnam,3,0,0.0_R4,cval)
      return
      end subroutine pcsetc

      subroutine gaseti(pnam,ival)
      use b2mod_types
      implicit none
      character*(*) pnam
      integer ival
      call ncpy_named_param(32,pnam,1,ival,0.0_R4,' ')
      return
      end subroutine gaseti

      subroutine gasetr(pnam,rval)
      use b2mod_types
      implicit none
      character*(*) pnam
      real (kind=R4) :: rval
      call ncpy_named_param(32,pnam,2,0,rval,' ')
      return
      end subroutine gasetr

c ======================= Axes and colour bars =========================
      subroutine labmod(fmtx,fmty,numx,numy,iszx,iszy,ixdc,iydc,ixor)
      use ncar_py_state
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
      use b2mod_b2plot_py
      use ncar_py_state
      implicit none
      integer mjrx, mnrx, mjry, mnry, ixlb, iylb, igph
      real (kind=R4) :: xint, yint
      if (.not.b2py_active) return
      call b2py_begin(47)
      call b2py_put_i4(mjrx)
      call b2py_put_i4(mnrx)
      call b2py_put_i4(mjry)
      call b2py_put_i4(mnry)
      call b2py_put_i4(ixlb)
      call b2py_put_i4(iylb)
      call b2py_put_i4(igph)
      call b2py_put_r4(xint)
      call b2py_put_r4(yint)
      call b2py_put_str(lm_fmtx)
      call b2py_put_str(lm_fmty)
      call b2py_put_i4(lm_numx)
      call b2py_put_i4(lm_numy)
      call b2py_put_i4(lm_iszx)
      call b2py_put_i4(lm_iszy)
      call b2py_put_i4(lm_ixdc)
      call b2py_put_i4(lm_iydc)
      call b2py_put_i4(lm_ixor)
      call b2py_end
      return
      end subroutine gridal

      subroutine lblbar(ihov,xleb,xreb,ybeb,yteb,nbox,wsfb,hsfb,lfin,
     . iftp,llbs,nlbs,lbab)
      use b2mod_b2plot_py
      implicit none
      integer ihov, nbox, lfin(*), iftp, nlbs, lbab
      real (kind=R4) :: xleb, xreb, ybeb, yteb, wsfb, hsfb
      character*(*) llbs(*)
      integer i
      if (.not.b2py_active) return
      call b2py_begin(48)
      call b2py_put_i4(ihov)
      call b2py_put_r4(xleb)
      call b2py_put_r4(xreb)
      call b2py_put_r4(ybeb)
      call b2py_put_r4(yteb)
      call b2py_put_i4(nbox)
      call b2py_put_r4(wsfb)
      call b2py_put_r4(hsfb)
      call b2py_put_i4v(nbox,lfin)
      call b2py_put_i4(iftp)
      call b2py_put_i4(max(nlbs,0))
      do i = 1, nlbs
        call b2py_put_str(llbs(i))
      end do
      call b2py_put_i4(lbab)
      call b2py_end
      return
      end subroutine lblbar

c ===================== Conpack and Conpackt ===========================
      subroutine cpseti(pnam,ival)
      use ncar_py_state
      implicit none
      character*(*) pnam
      integer ival
      call ncpy_pset(ncpy_cp,pnam,1,ival,0.0_R4,' ')
      return
      end subroutine cpseti

      subroutine cpsetr(pnam,rval)
      use ncar_py_state
      implicit none
      character*(*) pnam
      real (kind=R4) :: rval
      call ncpy_pset(ncpy_cp,pnam,2,0,rval,' ')
      return
      end subroutine cpsetr

      subroutine cpsetc(pnam,cval)
      use ncar_py_state
      implicit none
      character*(*) pnam, cval
      call ncpy_pset(ncpy_cp,pnam,3,0,0.0_R4,cval)
      return
      end subroutine cpsetc

      subroutine cpgetr(pnam,rval)
      use ncar_py_state
      implicit none
      character*(*) pnam
      real (kind=R4) :: rval
      rval = ncpy_pget(ncpy_cp,pnam)
      return
      end subroutine cpgetr

      subroutine ctseti(pnam,ival)
      use ncar_py_state
      implicit none
      character*(*) pnam
      integer ival
      call ncpy_pset(ncpy_ct,pnam,1,ival,0.0_R4,' ')
      return
      end subroutine ctseti

      subroutine ctsetr(pnam,rval)
      use ncar_py_state
      implicit none
      character*(*) pnam
      real (kind=R4) :: rval
      call ncpy_pset(ncpy_ct,pnam,2,0,rval,' ')
      return
      end subroutine ctsetr

      subroutine ctsetc(pnam,cval)
      use ncar_py_state
      implicit none
      character*(*) pnam, cval
      call ncpy_pset(ncpy_ct,pnam,3,0,0.0_R4,cval)
      return
      end subroutine ctsetc

      subroutine ctmesh(rpnt,npnt,lopn,iedg,nedg,loen,itri,ntri,lotn,
     . rwrk,lrwk,iwrk,liwk)
      use b2mod_b2plot_py
      use ncar_py_state
      implicit none
      integer npnt, lopn, iedg(*), nedg, loen, itri(*), ntri, lotn,
     . lrwk, iwrk(*), liwk
      real (kind=R4) :: rpnt(*), rwrk(*)
      if (.false.) then
        idum = lrwk+liwk+iwrk(1)
        rdum = rwrk(1)
      end if
      current_source = ncpy_ct
      if (.not.b2py_active) return
      call b2py_begin(50)
      call b2py_put_i4(npnt)
      call b2py_put_i4(lopn)
      call b2py_put_i4(nedg)
      call b2py_put_i4(loen)
      call b2py_put_i4(ntri)
      call b2py_put_i4(lotn)
      call b2py_put_r4v(npnt,rpnt)
      call b2py_put_i4v(nedg,iedg)
      call b2py_put_i4v(ntri,itri)
      call b2py_end
      return
      end subroutine ctmesh

      subroutine cprect(zdat,kzdt,mzdt,nzdt,rwrk,lrwk,iwrk,liwk)
      use b2mod_b2plot_py
      use ncar_py_state
      implicit none
      integer kzdt, mzdt, nzdt, lrwk, iwrk(*), liwk
      real (kind=R4) :: zdat(*), rwrk(*)
      if (.false.) then
        idum = lrwk+liwk+iwrk(1)
        rdum = rwrk(1)
      end if
      current_source = ncpy_cp
      if (.not.b2py_active) return
      call b2py_begin(51)
      call b2py_put_i4(kzdt)
      call b2py_put_i4(mzdt)
      call b2py_put_i4(nzdt)
      call b2py_put_r4v(kzdt*nzdt,zdat)
      call b2py_end
      return
      end subroutine cprect

c Contour lines and labels are drawn by Python from the tables
      subroutine cpclam(zdat,rwrk,iwrk,iama)
      use ncar_py_state
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
      use ncar_py_state
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
      use ncar_py_state
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
      use ncar_py_state
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
      use ncar_py_state
      implicit none
      real (kind=R4) :: zdat(*), rwrk(*)
      integer iwrk(*), iama(*)
      external rtpl
      if (.false.) then
        idum = iwrk(1)+iama(1)
        rdum = zdat(1)+rwrk(1)
        call rtpl
      end if
      call ncpy_contour_lines(ncpy_cp)
      return
      end subroutine cpcldm

      subroutine ctcldm(rpnt,iedg,itri,rwrk,iwrk,iama,rtpl)
      use ncar_py_state
      implicit none
      real (kind=R4) :: rpnt(*), rwrk(*)
      integer iedg(*), itri(*), iwrk(*), iama(*)
      external rtpl
      if (.false.) then
        idum = iedg(1)+itri(1)+iwrk(1)+iama(1)
        rdum = rpnt(1)+rwrk(1)
        call rtpl
      end if
      call ncpy_contour_lines(ncpy_ct)
      return
      end subroutine ctcldm

      subroutine cplbdr(zdat,rwrk,iwrk)
      use ncar_py_state
      implicit none
      real (kind=R4) :: zdat(*), rwrk(*)
      integer iwrk(*)
      if (.false.) then
        idum = iwrk(1)
        rdum = zdat(1)+rwrk(1)
      end if
      call ncpy_contour_info(ncpy_cp)
      return
      end subroutine cplbdr

      subroutine ctlbdr(rpnt,iedg,itri,rwrk,iwrk)
      use ncar_py_state
      implicit none
      real (kind=R4) :: rpnt(*), rwrk(*)
      integer iedg(*), itri(*), iwrk(*)
      if (.false.) then
        idum = iedg(1)+itri(1)+iwrk(1)
        rdum = rpnt(1)+rwrk(1)
      end if
      call ncpy_contour_info(ncpy_ct)
      return
      end subroutine ctlbdr

      subroutine cpdrpl(xcs,ycs,ncs,iai,iag,nai)
      use ncar_py_state
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
      use b2mod_b2plot_py
      use ncar_py_state
      implicit none
      real (kind=R4) :: zdat(*), flow, fhgh, finc
      integer kzdt, mzdt, nzdt, kset, nhgh, ndsh
      if (.false.) then
        idum = kzdt+mzdt+nzdt+kset+nhgh+ndsh
        rdum = zdat(1)+flow+fhgh+finc
      end if
      call b2py_message(2, 'ncar_py: cpcnrc not supported')
      return
      end subroutine cpcnrc

c ================================ Areas ===============================
      subroutine arinam(iama,lama)
      use b2mod_b2plot_py
      use ncar_py_state
      implicit none
      integer iama(*), lama
      if (.false.) idum = iama(1)+lama
      if (.not.b2py_active) return
      call b2py_begin(56)
      call b2py_end
      return
      end subroutine arinam

      subroutine aredam(iama,xca,yca,lca,igi,idl,idr)
      use b2mod_b2plot_py
      use ncar_py_state
      implicit none
      integer iama(*), lca, igi, idl, idr
      real (kind=R4) :: xca(*), yca(*)
      if (.false.) idum = iama(1)
      if (.not.b2py_active) return
      call b2py_begin(52)
      call b2py_put_i4(igi)
      call b2py_put_i4(idl)
      call b2py_put_i4(idr)
      call ncpy_points(lca,xca,yca)
      call b2py_end
      return
      end subroutine aredam

      subroutine arpram(iama,if1,if2,if3)
      use ncar_py_state
      implicit none
      integer iama(*), if1, if2, if3
      if (.false.) idum = iama(1)+if1+if2+if3
      return
      end subroutine arpram

c Contour fill: the colour of each area identifier is obtained by
c calling the user fill routine apr in probing mode
      subroutine arscam(iama,xcs,ycs,mcs,iai,iag,mai,apr)
      use b2mod_b2plot_py
      use ncar_py_state
      implicit none
      integer iama(*), mcs, iai(*), iag(*), mai
      real (kind=R4) :: xcs(*), ycs(*)
      external apr
      integer it, k, n, nid, iarea(2), igrp(2), icol(ncpy_nlev+1)
      real (kind=R4) :: xp(4), yp(4)
      if (.false.) then
        idum = iama(1)+mcs+iai(1)+iag(1)+mai
        rdum = xcs(1)+ycs(1)
      end if
      if (.not.b2py_active) return
      it = current_source
      if (it.ne.ncpy_ct) it = ncpy_cp
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
      call b2py_begin(53)
      call b2py_put_i4(it)
      call b2py_put_i4(n)
      call b2py_put_r4v(n,tclv(1,it))
      call b2py_put_i4v(n,taia(1,it))
      call b2py_put_i4v(n,taib(1,it))
      call b2py_put_i4(nid)
      call b2py_put_i4v(nid,icol)
      call b2py_put_r4(txc1(it))
      call b2py_put_r4(txcm(it))
      call b2py_put_r4(tyc1(it))
      call b2py_put_r4(tycn(it))
      call b2py_put_r4(tspv(it))
      call b2py_put_i4(tmap(it))
      call b2py_end
      return
      end subroutine arscam

      subroutine ardbpa(iama,igip,label)
      use ncar_py_state
      implicit none
      integer iama(*), igip
      character*(*) label
      if (.false.) idum = iama(1)+igip+len(label)
      return
      end subroutine ardbpa

c ========================== Error handling ============================
      subroutine entsr(irold,irnew)
      use ncar_py_state
      implicit none
      integer irold, irnew
      if (.false.) idum = irnew
      irold = 0
      return
      end subroutine entsr

      subroutine retsr(irold)
      use ncar_py_state
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
      call ncpy_attr_i(60,ib)
      return
      end subroutine gflas1

      subroutine gflas2
      use b2mod_b2plot_py
      implicit none
      if (.not.b2py_active) return
      call b2py_begin(61)
      call b2py_end
      return
      end subroutine gflas2

      subroutine gflas3(ib)
      implicit none
      integer ib
      call ncpy_attr_i(62,ib)
      return
      end subroutine gflas3

c ======================= Colour conversion ============================
c HLS to RGB with the NCAR convention: hue 0 is blue, L and S in percent
      subroutine hlsrgb(h,l,s,r,g,b)
      use b2mod_types
      implicit none
      real (kind=R4) :: h, l, s, r, g, b
      real (kind=R4) :: xl, xs, m1, m2, ncpy_hls_value
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
      r = ncpy_hls_value(m1,m2,h)
      g = ncpy_hls_value(m1,m2,h-120.0_R4)
      b = ncpy_hls_value(m1,m2,h-240.0_R4)
      return
      end subroutine hlsrgb

      real (kind=R4) function ncpy_hls_value(m1,m2,hue)
      use b2mod_types
      implicit none
      real (kind=R4) :: m1, m2, hue
      real (kind=R4) :: hh
      hh = modulo(hue,360.0_R4)
      if (hh.lt.60.0_R4) then
        ncpy_hls_value = m1+(m2-m1)*hh/60.0_R4
      else if (hh.lt.180.0_R4) then
        ncpy_hls_value = m2
      else if (hh.lt.240.0_R4) then
        ncpy_hls_value = m1+(m2-m1)*(240.0_R4-hh)/60.0_R4
      else
        ncpy_hls_value = m1
      end if
      return
      end function ncpy_hls_value
