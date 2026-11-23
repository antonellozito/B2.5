      subroutine rlcomp
      use b2mod_types
      use b2mod_subsys
      use b2mod_rdload
      implicit none
c
c======================================================================
c*** Computation of the radiation loads.
c*** Returns the total radiation power coming to each of the specified
c*** wall segments, with shadowing, in toroidal geometry (rad3d): each
c*** source is a ring around the symmetry axis, each segment a surface
c*** of revolution.
c***
c*** Input:
c***  xysrc   R and Z coordinates of the sources
c***  ssrcm   their power (W)
c***  xp,yp   coordinates of the boundary polygon P
c***  nbndry  number of points in shadowing structure S
c***  lbnd    indices of the S segments (assumed closed) in P
c***  ntrg    number of segments representing the radiation target T
c***  ltrg    their indices in P (must belong to S!)
c***
c*** Output
c***  radld   radiation power coming to the segments of T
c***  ttsrc   the total power radiated from the sources
c***  tttrg   the total power coming to T
c======================================================================
      integer i, j, k
      integer, allocatable :: kt(:)
      real (kind=R8), allocatable :: es(:,:)
      intrinsic mod
      external rad3d, xerrab
c======================================================================
      call subini ('rlcomp')
      allocate(es(4,max(1,nbndry)), kt(max(1,ntrg)))
      do k=1,nbndry
        i=lbnd(k)
        j=lbnd(mod(k,nbndry)+1)
        es(1,k)=xp(i)
        es(2,k)=yp(i)
        es(3,k)=xp(j)
        es(4,k)=yp(j)
      end do
      do j=1,ntrg
        kt(j)=0
        do k=1,nbndry
          if(lbnd(k) .eq. ltrg(j)) then
            kt(j)=k
            exit
          end if
        end do
        if(kt(j) .eq. 0) then
          write(*,*) 'Error in rlcomp: target not in structure',
     ,      j,ltrg(j)
          call xerrab ('Target not in structure')
        end if
      end do
      call rad3d(nsrc,xysrc,ssrcm,nbndry,es,ntrg,kt,radld)
      ttsrc=sum(ssrcm(1:nsrc))
      tttrg=sum(radld(1:ntrg))
      deallocate(es, kt)
      call subend ()
      return
      end subroutine rlcomp
c
c======================================================================
c*** Radiation power reaching wall segments in toroidal geometry.
c*** Source i, of power ps(i), is a ring of radius xys(1,i) at height
c*** xys(2,i) radiating isotropically. Segment k of the wall contour,
c*** from (es(1,k),es(2,k)) to (es(3,k),es(4,k)) in (R,Z), is a surface
c*** of revolution. A point of a ring is seen from a point of a target
c*** segment when the straight line between them, whose trace in the
c*** (R,Z) plane is a hyperbola, crosses no other segment. Target j is
c*** segment kt(j); the power reaching it is
c***   pw(j) = sum_i ps(i)/(4 pi^2) int dA int_0^pi |cos|/d^2 dphi
c*** with phi the toroidal angle between the two points, d their
c*** distance and cos the angle of the line to the segment normal.
c*** The toroidal angle is integrated with nphi points, mapped by a
c*** tangent so that they are dense where the ring is closest; the
c*** segment with up to nsubm points, more when sources are close.
c*** Sources further than 4h+L from the segment centre (L its length)
c*** are merged within the squares of side h of a grid, at their
c*** centre weighted with the absolute power.
c======================================================================
      subroutine rad3d(ns,xys,ps,ne,es,nt,kt,pw)
      use b2mod_types
      use b2mod_constants, only : pi
      use b2mod_subsys
      implicit none
      integer, intent(in) :: ns, ne, nt, kt(nt)
      real (kind=R8), intent(in) :: xys(2,ns), ps(ns), es(4,ne)
      real (kind=R8), intent(out) :: pw(nt)
      integer nphi, nsubm
      real (kind=R8) :: h, teps, tend
      parameter (nphi=24, nsubm=24, h=0.03_R8)
      parameter (teps=1.0e-9_R8, tend=1.0_R8-1.0e-7_R8)
      integer i, j, k, kk, ib, ic, nc, nbr, nbz, nsl, nsub, isub, iph,
     .  nk
      integer, allocatable :: head(:), nxt(:), lc(:), lk(:)
      real (kind=R8), allocatable :: cp(:), ca(:), cr(:), cz(:),
     .  sr(:), sz(:), sp(:), erl(:), erh(:), ezl(:), ezh(:)
      real (kind=R8) :: rmn, rmx, zmn, zmx, ra, za, rb, zb, el, xm, zm,
     .  dmin, xn, zn, f, rw, zw, rs, zs, dz, wd, tmax, tt, phi, dph,
     .  x0, y0, vx, vy, vz, d2, cs, sm, e, aa, bb, cc, rmin, rh, zl, zh,
     .  ts
      logical seen
      intrinsic abs, atan, cos, sin, tan, sqrt, min, max,
     .  minval, maxval, int, ceiling, real
c======================================================================
      call subini ('rad3d')
      pw(1:nt)=0.0_R8
      if(ns.le.0 .or. nt.le.0) then
        call subend ()
        return
      end if
c
c*** Squares of the grid merging the far sources: their power, absolute
c*** power, centre, and the list of their sources
c
      rmn=minval(xys(1,1:ns))
      rmx=maxval(xys(1,1:ns))
      zmn=minval(xys(2,1:ns))
      zmx=maxval(xys(2,1:ns))
      nbr=int((rmx-rmn)/h)+1
      nbz=int((zmx-zmn)/h)+1
      allocate(head(nbr*nbz), nxt(ns), lc(nbr*nbz), cp(nbr*nbz),
     .  ca(nbr*nbz), cr(nbr*nbz), cz(nbr*nbz))
      head=0
      cp=0.0_R8
      ca=0.0_R8
      cr=0.0_R8
      cz=0.0_R8
      do i=1,ns
        nxt(i)=0
        if(ps(i).eq.0.0_R8) cycle
        ib=1+min(nbr-1,int((xys(1,i)-rmn)/h))+
     .    nbr*min(nbz-1,int((xys(2,i)-zmn)/h))
        nxt(i)=head(ib)
        head(ib)=i
        cp(ib)=cp(ib)+ps(i)
        ca(ib)=ca(ib)+abs(ps(i))
        cr(ib)=cr(ib)+abs(ps(i))*xys(1,i)
        cz(ib)=cz(ib)+abs(ps(i))*xys(2,i)
      end do
      nc=0
      do ib=1,nbr*nbz
        if(ca(ib).le.0.0_R8) cycle
        nc=nc+1
        lc(nc)=ib
        cr(ib)=cr(ib)/ca(ib)
        cz(ib)=cz(ib)/ca(ib)
      end do
c
c*** Bounding boxes of the segments
c
      allocate(sr(ns+nc), sz(ns+nc), sp(ns+nc), lk(ne),
     .  erl(ne), erh(ne), ezl(ne), ezh(ne))
      do k=1,ne
        erl(k)=min(es(1,k),es(3,k))
        erh(k)=max(es(1,k),es(3,k))
        ezl(k)=min(es(2,k),es(4,k))
        ezh(k)=max(es(2,k),es(4,k))
      end do
c
c*** Loop over the targets
c
      do j=1,nt
        k=kt(j)
        ra=es(1,k)
        za=es(2,k)
        rb=es(3,k)
        zb=es(4,k)
        el=sqrt((rb-ra)**2+(zb-za)**2)
        if(el.le.0.0_R8) cycle
        xm=0.5_R8*(ra+rb)
        zm=0.5_R8*(za+zb)
        xn=(zb-za)/el
        zn=-(rb-ra)/el
c Sources: the near ones one by one, the far ones merged
        nsl=0
        do ic=1,nc
          ib=lc(ic)
          if((cr(ib)-xm)**2+(cz(ib)-zm)**2.lt.(4.0_R8*h+el)**2)
     .      then
            i=head(ib)
            do while (i.gt.0)
              nsl=nsl+1
              sr(nsl)=xys(1,i)
              sz(nsl)=xys(2,i)
              sp(nsl)=ps(i)
              i=nxt(i)
            end do
          else
            nsl=nsl+1
            sr(nsl)=cr(ib)
            sz(nsl)=cz(ib)
            sp(nsl)=cp(ib)
          end if
        end do
        if(nsl.eq.0) cycle
        dmin=sqrt(minval((sr(1:nsl)-xm)**2+(sz(1:nsl)-zm)**2))
        nsub=min(nsubm,max(1,ceiling(3.0_R8*el/max(dmin,1.0e-6_R8))))
c Points of the segment
        do isub=1,nsub
          f=(real(isub,R8)-0.5_R8)/real(nsub,R8)
          rw=ra+f*(rb-ra)
          zw=za+f*(zb-za)
          if(rw.le.0.0_R8) cycle
          e=0.0_R8
          do i=1,nsl
            rs=sr(i)
            zs=sz(i)
            if(rs.le.0.0_R8 .or. sp(i).eq.0.0_R8) cycle
            dz=zw-zs
            wd=min(pi,max(1.0e-5_R8,sqrt(((rs-rw)**2+dz*dz)/
     .        (rs*rw))))
            tmax=atan(pi/wd)
c Segments that may cross the line: within its range of Z, below its
c largest R
            zl=min(zs,zw)
            zh=max(zs,zw)
            rh=max(rs,rw)
            nk=0
            do kk=1,ne
              if(kk.eq.k) cycle
              if(ezh(kk).lt.zl .or. ezl(kk).gt.zh .or.
     .           erl(kk).gt.rh) cycle
              nk=nk+1
              lk(nk)=kk
            end do
            sm=0.0_R8
            do iph=1,nphi
              tt=(real(iph,R8)-0.5_R8)/real(nphi,R8)*tmax
              phi=wd*tan(tt)
              dph=wd/cos(tt)**2*tmax/real(nphi,R8)
              x0=rs*cos(phi)
              y0=rs*sin(phi)
              vx=x0-rw
              vy=y0
              vz=zs-zw
              d2=vx*vx+vy*vy+vz*vz
              cs=abs(vx*xn+vz*zn)/sqrt(d2)
              if(cs.le.0.0_R8) cycle
c Line from the ring point (t=0) to the wall point (t=1): R(t)**2 =
c aa*t**2+bb*t+cc, Z(t) = zs+t*dz; smallest R on it
              aa=vx*vx+vy*vy
              bb=-2.0_R8*(x0*vx+y0*vy)
              cc=rs*rs
              rmin=min(rs,rw)
              if(aa.gt.0.0_R8) then
                ts=-bb/(2.0_R8*aa)
                if(ts.gt.0.0_R8 .and. ts.lt.1.0_R8)
     .            rmin=sqrt(max(0.0_R8,cc-0.25_R8*bb*bb/aa))
              end if
              seen=.true.
              do ic=1,nk
                kk=lk(ic)
                if(erh(kk).lt.rmin) cycle
                if(crosses(es(1,kk),es(2,kk),es(3,kk),es(4,kk))) then
                  seen=.false.
                  exit
                end if
              end do
              if(seen) sm=sm+cs/d2*dph
            end do
            e=e+sp(i)*sm
          end do
          pw(j)=pw(j)+e*rw*el/(2.0_R8*pi*real(nsub,R8))
        end do
      end do
      deallocate(head, nxt, lc, cp, ca, cr, cz, sr, sz, sp, lk,
     .  erl, erh, ezl, ezh)
      call subend ()
      return

      contains

c True when the line crosses the segment from (r1,z1) to (r2,z2): a
c point of the segment, R = al+be*t at the height of the line, with
c R(t)**2 = (al+be*t)**2 for 0 < t < 1
      logical function crosses(r1,z1,r2,z2)
      real (kind=R8), intent(in) :: r1, z1, r2, z2
      real (kind=R8) :: dze, u0, u1, al, be, a2, b2, c2, disc, sq, t, r
      integer m
      crosses=.false.
      dze=z2-z1
      if(abs(dze).gt.1.0e-12_R8) then
        u0=(zs-z1)/dze
        u1=dz/dze
        al=r1+(r2-r1)*u0
        be=(r2-r1)*u1
        a2=aa-be*be
        b2=bb-2.0_R8*al*be
        c2=cc-al*al
        if(abs(a2).le.1.0e-14_R8*(abs(b2)+1.0e-30_R8)) then
          if(b2.ne.0.0_R8) crosses=hit(-c2/b2,u0,u1,al,be)
        else
          disc=b2*b2-4.0_R8*a2*c2
          if(disc.ge.0.0_R8) then
            sq=sqrt(disc)
            do m=-1,1,2
              t=(-b2+real(m,R8)*sq)/(2.0_R8*a2)
              if(hit(t,u0,u1,al,be)) then
                crosses=.true.
                return
              end if
            end do
          end if
        end if
      else if(dz.ne.0.0_R8) then
        t=(z1-zs)/dz
        if(t.gt.teps .and. t.lt.tend) then
          r=sqrt(max(0.0_R8,(aa*t+bb)*t+cc))
          crosses=r.ge.min(r1,r2) .and. r.le.max(r1,r2)
        end if
      end if
      return
      end function crosses

c True when t is inside the line and the point of the segment at t
c is on it, at positive R
      logical function hit(t,u0,u1,al,be)
      real (kind=R8), intent(in) :: t, u0, u1, al, be
      real (kind=R8) :: u
      u=u0+u1*t
      hit=t.gt.teps .and. t.lt.tend .and. u.ge.0.0_R8 .and.
     .  u.le.1.0_R8 .and. al+be*t.ge.0.0_R8
      return
      end function hit

      end subroutine rad3d
