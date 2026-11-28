c Toroidal tiles of the B2.5 wall elements: the effective tile of an
c element (TILE_SECTORS tiles around the torus, each tilted toroidally
c by TILE_TILT degrees, separated by toroidal gaps of TILE_GAP m, each
c raised by TILE_MISALIGN m above its upstream neighbour), the split of
c the loads of the element over the tile and the heat conduction in the
c tile (toroidal position x depth).
c
c The cut: the plane of the toroidal direction and the normal of the
c element (the tiles are uniform along the element: no poloidal gaps).
c zeta runs along the toroidal motion of the particles reaching the
c element, y along the normal of the element into the plasma. Each tile
c is a rectangle, PLATE_THICKNESS deep, of toroidal extent
c L = per - TILE_GAP (per: the toroidal length of the element over
c TILE_SECTORS), tilted by TILE_TILT (positive: its top rising along
c zeta, facing the incoming particles); corners E (bottom, leading),
c A (top, leading), B (top, trailing), C (bottom, trailing), tile j is
c tile 0 moved by j (per, TILE_MISALIGN). The surface is the polyline
c E A B C E' ... : edge 1 the leading side (E A), 2 the top (A B), 3 the
c trailing side (B C), 4 the floor of the gap (C E'), which belongs to
c the next tile. Points on the tile: along the top from A, down the
c leading side from A, down the trailing side from B.
c
c Loads, per unit length of the element and per tile (W/m): the
c field-aligned heat flux of the element (W/m^2 of the element)
c deposited where the field lines (electrons, tile_fieldline) or the ion
c orbits (tile_gyro) first cross the surface, from one period of the
c array mapped onto tile 0 (exact power: q0 (per + TILE_MISALIGN /
c tan theta) per tile, theta the angle of the field lines to the
c element in the cut), the other loads spread evenly on the top
c (tile_loads). Conduction: tile2d_steady, tile2d_transient.

c..................................................................
c Shape of the tiles of element ibnd in the cut (torl: toroidal length
c of the element): corners vz, vy (E, A, B, C) of tile 0, its period
c (tp, tm), the length of its top tlen and its depth.
      subroutine tile_shape(ibnd,torl,vz,vy,tp,tm,tlen,depth)
      use b2mod_types
      use b2mod_constants
      use b2mod_wall
      implicit none
      integer, intent(in) :: ibnd
      real(kind=R8), intent(in) :: torl
      real(kind=R8), intent(out) :: vz(4), vy(4), tp, tm, tlen, depth

      real(kind=R8) :: dl, cd, sd, h
      external xertst

      call xertst(tile_sectors(ibnd).gt.0,
     . 'tile_shape: TILE_SECTORS must be positive')
      tp = torl/real(tile_sectors(ibnd),R8)
      tm = tile_misalign(ibnd)
      dl = tp-tile_gap(ibnd)
      depth = plate_thickness(ibnd)
      call xertst(dl.gt.0.0_R8 .and. tile_gap(ibnd).ge.0.0_R8,
     . 'tile_shape: TILE_GAP outside 0 .. toroidal length/sectors')
      call xertst(depth.gt.0.0_R8,
     . 'tile_shape: PLATE_THICKNESS must be positive')
      call xertst(abs(tile_tilt(ibnd)).lt.45.0_R8,
     . 'tile_shape: TILE_TILT must be below 45 degrees')
      cd = cos(tile_tilt(ibnd)*pi/180.0_R8)
      sd = sin(tile_tilt(ibnd)*pi/180.0_R8)
      tlen = dl/cd
      h = dl*sd/cd-tm
      call xertst(tile_gap(ibnd)*cd-h*sd.ge.0.0_R8,
     . 'tile_shape: the tiles overlap (TILE_GAP too small)')
      vz(2) = 0.0_R8
      vy(2) = 0.0_R8
      vz(3) = dl
      vy(3) = dl*sd/cd
      vz(1) = vz(2)+depth*sd
      vy(1) = vy(2)-depth*cd
      vz(4) = vz(3)+depth*sd
      vy(4) = vy(3)-depth*cd

      return
      end subroutine tile_shape

c..................................................................
c Direction of the particles reaching an element along the magnetic
c field, in the cut. In: the normal (rn, zn) of the element into the
c plasma, the field (br, bz, bphi), in the right-handed (R, phi, Z).
c Out: sina0 = |b.n| (zero: no incidence along the field); the cut
c components of the unit direction d = -sign(b.n) b into the element:
c dz = |d_phi| (zeta along d_phi), dn = d.n = -sina0; the cut
c components of the unit vectors u1, u2 normal to b with (u1, u2, b)
c right-handed: an ion (charge > 0) is at rho (sin psi u1 + cos psi u2)
c from its guiding centre, psi = psi0 + Omega t.
      subroutine tile_direction(rn,zn,br,bz,bphi,sina0,dz,dn,
     .  u1z,u1n,u2z,u2n)
      use b2mod_types
      implicit none
      real(kind=R8), intent(in) :: rn, zn, br, bz, bphi
      real(kind=R8), intent(out) :: sina0, dz, dn, u1z, u1n, u2z, u2n

      real(kind=R8) :: b(3), d(3), en(3), ez(3), w(3), u1(3), u2(3),
     .  bm, bn, ww, sg
      external xertst

c     Cartesian at the element: x = R, y = phi, z = Z
      bm = sqrt(br**2+bz**2+bphi**2)
      call xertst(bm.gt.0.0_R8,'tile_direction: no magnetic field')
      b(1) = br/bm
      b(2) = bphi/bm
      b(3) = bz/bm
      en(1) = rn/sqrt(rn**2+zn**2)
      en(2) = 0.0_R8
      en(3) = zn/sqrt(rn**2+zn**2)
      bn = b(1)*en(1)+b(3)*en(3)
      sina0 = abs(bn)
      sg = -sign(1.0_R8,bn)
      d = sg*b
      ez = 0.0_R8
      ez(2) = sign(1.0_R8,d(2))
      dz = d(2)*ez(2)
      dn = d(1)*en(1)+d(3)*en(3)
      w = ez-(ez(2)*b(2))*b
      ww = sqrt(w(1)**2+w(2)**2+w(3)**2)
      if (ww.lt.1.0e-6_R8) then
        w = en-bn*b
        ww = sqrt(w(1)**2+w(2)**2+w(3)**2)
      endif
      u1 = w/ww
      u2(1) = b(2)*u1(3)-b(3)*u1(2)
      u2(2) = b(3)*u1(1)-b(1)*u1(3)
      u2(3) = b(1)*u1(2)-b(2)*u1(1)
      u1z = u1(2)*ez(2)
      u1n = u1(1)*en(1)+u1(3)*en(3)
      u2z = u2(2)*ez(2)
      u2n = u2(1)*en(1)+u2(3)*en(3)

      return
      end subroutine tile_direction

c..................................................................
c Number of nodes along the top of a tile (tile_grid).
      function tile_nx()
      use b2mod_wall
      implicit none
      integer tile_nx

      tile_nx = 2*min(50,max(10,ndepth_nml-1))+1

      return
      end function tile_nx

c..................................................................
c Nodes xg(nx) along the top of the tiles of element ibnd, 0 .. tlen
c (nx = tile_nx()): symmetric; from each end the first interval of the
c depth points growing by 1.15 up to a uniform spacing in the middle
c (geometric throughout if that cannot fill the tile, uniform if the
c first interval is not finer).
      subroutine tile_grid(ibnd,tlen,nx,xg)
      use b2mod_types
      use b2mod_wall
      implicit none
      integer, intent(in) :: ibnd, nx
      real(kind=R8), intent(in) :: tlen
      real(kind=R8), intent(out) :: xg(nx)

      real(kind=R8) :: rg
      parameter (rg=1.15_R8)
      integer i, it, nh, ng
      real(kind=R8) :: e0, half, r, rlo, rhi, s, dlo, dhi, dm, sp

      nh = (nx-1)/2
      half = 0.5_R8*tlen
      e0 = target_depth(ibnd,2)-target_depth(ibnd,1)
      r = rg
      if (e0*real(nh,R8).ge.half .or. e0.le.0.0_R8) then
c       uniform
        e0 = half/real(nh,R8)
        ng = 0
        dm = e0
      else if (e0*(rg**nh-1.0_R8)/(rg-1.0_R8).le.half) then
c       geometric throughout, its ratio to fill the half
        rlo = rg
        rhi = 2.0_R8*rg
        do while (e0*(rhi**nh-1.0_R8)/(rhi-1.0_R8).lt.half)
          rhi = 2.0_R8*rhi
        enddo
        do it = 1, 200
          r = 0.5_R8*(rlo+rhi)
          s = e0*(r**nh-1.0_R8)/(r-1.0_R8)
          if (s.lt.half) then
            rlo = r
          else
            rhi = r
          endif
        enddo
        r = 0.5_R8*(rlo+rhi)
        ng = nh
        dm = 0.0_R8
      else
c       geometric up to the spacing dm, then uniform: dm by bisection
        dlo = e0
        dhi = half
        do it = 1, 200
          dm = 0.5_R8*(dlo+dhi)
          ng = min(nh,ceiling(log(dm/e0)/log(rg)))
          s = e0*(rg**ng-1.0_R8)/(rg-1.0_R8)+real(nh-ng,R8)*dm
          if (s.lt.half) then
            dlo = dm
          else
            dhi = dm
          endif
        enddo
        dm = 0.5_R8*(dlo+dhi)
        ng = min(nh,ceiling(log(dm/e0)/log(rg)))
      endif
      xg(1) = 0.0_R8
      do i = 1, nh-1
        if (i.le.ng) then
          sp = e0*r**(i-1)
        else
          sp = dm
        endif
        xg(i+1) = xg(i)+sp
      enddo
      xg(nh+1) = half
      do i = 1, nh
        xg(nx+1-i) = tlen-xg(i)
      enddo

      return
      end subroutine tile_grid

c..................................................................
c Adds the power pw (W/m), spread evenly over sa .. sb (a point if they
c are equal), to the nodes of grid g(n): node k holds the interval
c between the midpoints to its neighbours (g(1) .. g(n) in all).
      subroutine tile_deposit(sa,sb,pw,n,g,q)
      use b2mod_types
      implicit none
      integer, intent(in) :: n
      real(kind=R8), intent(in) :: sa, sb, pw, g(n)
      real(kind=R8), intent(inout) :: q(n)

      integer k
      real(kind=R8) :: a, b, lo, hi, ov

      a = max(g(1),min(g(n),min(sa,sb)))
      b = max(g(1),min(g(n),max(sa,sb)))
      do k = 1, n
        if (k.eq.1) then
          lo = g(1)
        else
          lo = 0.5_R8*(g(k-1)+g(k))
        endif
        if (k.eq.n) then
          hi = g(n)
        else
          hi = 0.5_R8*(g(k)+g(k+1))
        endif
        if (b.gt.a) then
          ov = max(0.0_R8,min(b,hi)-max(a,lo))
          q(k) = q(k)+pw*ov/(b-a)
        else if (a.ge.lo .and. (a.lt.hi .or. k.eq.n)) then
          q(k) = q(k)+pw
          return
        endif
      enddo

      return
      end subroutine tile_deposit

c..................................................................
c Puts a hit on edge e of a tile, at the fraction f from its first
c corner, onto tile 0: power pw (W/m) to pw(e) and to the nodes (top:
c along the top from A; sides: down from A or B; the floor: no nodes).
      subroutine tile_hit(e,f,fb,pwr,tlen,depth,nx,xg,nz,zg,qtop,
     .  qlead,qtrail,pw)
      use b2mod_types
      implicit none
      integer, intent(in) :: e, nx, nz
      real(kind=R8), intent(in) :: f, fb, pwr, tlen, depth, xg(nx),
     .  zg(nz)
      real(kind=R8), intent(inout) :: qtop(nx), qlead(nz), qtrail(nz),
     .  pw(4)
      external tile_deposit

      pw(e) = pw(e)+pwr
      if (e.eq.1) then
        call tile_deposit((1.0_R8-f)*depth,(1.0_R8-fb)*depth,pwr,
     .    nz,zg,qlead)
      else if (e.eq.2) then
        call tile_deposit(f*tlen,fb*tlen,pwr,nx,xg,qtop)
      else if (e.eq.3) then
        call tile_deposit(f*depth,fb*depth,pwr,nz,zg,qtrail)
      endif

      return
      end subroutine tile_hit

c..................................................................
c First crossing of the surface by the field line c (across the cut
c direction, unit vectors dh along and nh across it): edge e and tile j
c of the crossing, fraction f along the edge; e = 0 if none (cv, lv:
c c and position along dh of the corners 1..5, corner 5 = E of tile 1;
c dc, dl: the same for the period).
      subroutine tile_first(c,cv,lv,dc,dl,e,j,f)
      use b2mod_types
      implicit none
      real(kind=R8), intent(in) :: c, cv(5), lv(5), dc, dl
      integer, intent(out) :: e, j
      real(kind=R8), intent(out) :: f

      integer ie, jj, jlo, jhi
      real(kind=R8) :: cmin, cmax, fr, lam, best

      e = 0
      j = 0
      f = 0.0_R8
      best = huge(1.0_R8)
      do ie = 1, 4
        cmin = min(cv(ie),cv(ie+1))
        cmax = max(cv(ie),cv(ie+1))
        if (cmax-cmin.le.1.0e-14_R8*dc) cycle
        jlo = ceiling((c-cmax)/dc)
        jhi = floor((c-cmin)/dc)
        do jj = jlo, jhi
          fr = (c-cv(ie)-real(jj,R8)*dc)/(cv(ie+1)-cv(ie))
          if (fr.lt.0.0_R8 .or. fr.gt.1.0_R8) cycle
          lam = lv(ie)+real(jj,R8)*dl+fr*(lv(ie+1)-lv(ie))
          if (lam.lt.best) then
            best = lam
            e = ie
            j = jj
            f = fr
          endif
        enddo
      enddo

      return
      end subroutine tile_first

c..................................................................
c Exact split of a field-aligned load q0 (W/m^2 of the element) over
c the tile: the field lines of one period of the array, c in
c [c(B) - dc, c(B)), each to its first crossing of the surface, mapped
c onto tile 0. Between the c of the corners the same edge is crossed
c and the crossing moves linearly with c: the power (q0 |d_cut| / sina0
c per unit c) is put exactly on uniform pieces. Adds to qtop(nx),
c qlead(nz), qtrail(nz) (W/m into each node) and pw(4) (W/m: leading
c side, top, trailing side, floor of the gap); info 1: the tiles shadow
c each other entirely (nothing deposited).
      subroutine tile_fieldline(vz,vy,tp,tm,dz,dn,sina0,q0,tlen,depth,
     .  nx,xg,nz,zg,qtop,qlead,qtrail,pw,info)
      use b2mod_types
      implicit none
      integer, intent(in) :: nx, nz
      real(kind=R8), intent(in) :: vz(4), vy(4), tp, tm, dz, dn, sina0,
     .  q0, tlen, depth, xg(nx), zg(nz)
      real(kind=R8), intent(inout) :: qtop(nx), qlead(nz), qtrail(nz),
     .  pw(4)
      integer, intent(out) :: info

      integer i, k, nb, e, j
      real(kind=R8) :: dd, dhz, dhn, nhz, nhn, cv(5), lv(5), dc, dl,
     .  c0, cb(6), t, ca, cc, fa, fb, f, pwr, wz(5), wy(5)
      external tile_first, tile_hit, xertst

      info = 0
      call xertst(sina0.gt.0.0_R8,
     . 'tile_fieldline: no incidence along the field')
      dd = sqrt(dz**2+dn**2)
      dhz = dz/dd
      dhn = dn/dd
      nhz = -dhn
      nhn = dhz
      wz(1:4) = vz
      wy(1:4) = vy
      wz(5) = vz(1)+tp
      wy(5) = vy(1)+tm
      do i = 1, 5
        cv(i) = nhz*wz(i)+nhn*wy(i)
        lv(i) = dhz*wz(i)+dhn*wy(i)
      enddo
      dc = nhz*tp+nhn*tm
      dl = dhz*tp+dhn*tm
      if (dc.le.0.0_R8) then
        info = 1
        return
      endif
      if (q0.eq.0.0_R8) return
c     sub-bands: the c of the corners within the period
      c0 = cv(3)-dc
      nb = 2
      cb(1) = c0
      cb(2) = c0+dc
      do i = 1, 4
        t = c0+modulo(cv(i)-c0,dc)
        nb = nb+1
        cb(nb) = t
      enddo
      do i = 2, nb
        t = cb(i)
        k = i-1
        do while (k.ge.1)
          if (cb(k).le.t) exit
          cb(k+1) = cb(k)
          k = k-1
        enddo
        cb(k+1) = t
      enddo
      do i = 1, nb-1
        ca = cb(i)
        cc = cb(i+1)
        if (cc-ca.le.1.0e-13_R8*dc) cycle
        call tile_first(0.5_R8*(ca+cc),cv,lv,dc,dl,e,j,f)
        call xertst(e.gt.0,'tile_fieldline: a field line misses')
        fa = (ca-cv(e)-real(j,R8)*dc)/(cv(e+1)-cv(e))
        fb = (cc-cv(e)-real(j,R8)*dc)/(cv(e+1)-cv(e))
        fa = max(0.0_R8,min(1.0_R8,fa))
        fb = max(0.0_R8,min(1.0_R8,fb))
        pwr = q0*dd/sina0*(cc-ca)
        call tile_hit(e,fa,fb,pwr,tlen,depth,nx,xg,nz,zg,qtop,
     .    qlead,qtrail,pw)
      enddo

      return
      end subroutine tile_fieldline

c..................................................................
c Radical inverse of i in base nb (Halton sequence).
      function tile_halton(i,nb)
      use b2mod_types
      implicit none
      integer, intent(in) :: i, nb
      real(kind=R8) :: tile_halton

      integer k
      real(kind=R8) :: f

      tile_halton = 0.0_R8
      f = 1.0_R8/real(nb,R8)
      k = i
      do while (k.gt.0)
        tile_halton = tile_halton+f*real(mod(k,nb),R8)
        k = k/nb
        f = f/real(nb,R8)
      enddo

      return
      end function tile_halton

c..................................................................
c Ion part of a field-aligned load q0 (W/m^2 of the element) over the
c tile, with the orbits about the field lines (no electric field): ions
c of mass amass (proton masses) and charge zion, a Maxwellian at ti
c (eV) weighted by its energy flux along the field (energy E/ti of
c density x^2 exp(-x); pitch angle beta to b of density sin(2 beta)),
c gyrophase uniform, guiding centres on the field lines of one period
c of the array (as tile_fieldline), each orbit (helix, field bmag, in
c the cut) to its first crossing of the surface, mapped onto tile 0.
c nsamp orbits of equal power (Halton points in guiding centre,
c gyrophase, energy and pitch): the total is exact (q0 |d_cut| / sina0
c times the period across the field lines), its spread sampled. Adds
c to qtop, qlead, qtrail, pw as tile_fieldline; nfall: orbits put on
c their field line (no crossing found).
      subroutine tile_gyro(vz,vy,tp,tm,dz,dn,u1z,u1n,u2z,u2n,sina0,
     .  q0,ti,amass,zion,bmag,nsamp,tlen,depth,nx,xg,nz,zg,qtop,qlead,
     .  qtrail,pw,nfall)
      use b2mod_types
      use b2mod_constants
      implicit none
      integer, intent(in) :: nsamp, nx, nz
      real(kind=R8), intent(in) :: vz(4), vy(4), tp, tm, dz, dn, u1z,
     .  u1n, u2z, u2n, sina0, q0, ti, amass, zion, bmag, tlen, depth,
     .  xg(nx), zg(nz)
      real(kind=R8), intent(inout) :: qtop(nx), qlead(nz), qtrail(nz),
     .  pw(4)
      integer, intent(out) :: nfall

      integer nfine
      parameter (nfine=8)
      integer is, i, k, it, e, j, jlo, jhi, jj, ie, ehit, nstep
      real(kind=R8) :: dd, dhz, dhn, nhz, nhn, cv(5), lv(5), dc, dl,
     .  wz(5), wy(5), eta(5), etamax, vt, u, x, xlo, xhi, fx, v, beta,
     .  rho, cotb, az, an, ucut, gz, gn, psi, psi0, cs, lam, etad,
     .  zmin, zmax, reach, kmin, ke, kl, kx, kd, sl, sx, rate, tr, a1z,
     .  a1n, a2z, a2n, tz, tn, ez, en, ln, pz1, pn1, d1, d2,
     .  tbest, fr, frhit, psa, psb, psm, sm, pwr, pmz, pmn, hz, hn, sa,
     .  ca, cb, cr, ph, psn, pl, pr
      real(kind=R8) :: tile_halton
      logical found
      external tile_halton, tile_first, tile_hit, xertst

      nfall = 0
      if (q0.eq.0.0_R8) return
      call xertst(sina0.gt.0.0_R8 .and. ti.gt.0.0_R8 .and.
     . amass.gt.0.0_R8 .and. zion.gt.0.0_R8 .and. bmag.gt.0.0_R8,
     . 'tile_gyro: faulty input')
      dd = sqrt(dz**2+dn**2)
      dhz = dz/dd
      dhn = dn/dd
      nhz = -dhn
      nhn = dhz
      wz(1:4) = vz
      wy(1:4) = vy
      wz(5) = vz(1)+tp
      wy(5) = vy(1)+tm
      do i = 1, 5
        cv(i) = nhz*wz(i)+nhn*wy(i)
        lv(i) = dhz*wz(i)+dhn*wy(i)
        eta(i) = wy(i)-tm/tp*wz(i)
      enddo
      dc = nhz*tp+nhn*tm
      dl = dhz*tp+dhn*tm
      if (dc.le.0.0_R8) return
      etamax = maxval(eta)
      etad = dhn-tm/tp*dhz
      zmin = minval(wz)
      zmax = maxval(wz)
      ucut = sqrt(u1z**2+u1n**2+u2z**2+u2n**2)
      vt = sqrt(2.0_R8*ti*ev/(amass*mp))
      pwr = q0*dd/sina0*dc/real(nsamp,R8)
      do is = 1, nsamp
c       energy: x^2 exp(-x), its distribution 1-exp(-x)(1+x+x^2/2) = u
        u = tile_halton(is,5)
        xlo = 0.0_R8
        xhi = 60.0_R8
        do it = 1, 100
          x = 0.5_R8*(xlo+xhi)
          fx = 1.0_R8-exp(-x)*(1.0_R8+x+0.5_R8*x**2)
          if (fx.lt.u) then
            xlo = x
          else
            xhi = x
          endif
        enddo
        x = 0.5_R8*(xlo+xhi)
        v = vt*sqrt(x)
        beta = asin(sqrt(tile_halton(is,7)))
        rho = amass*mp*v*sin(beta)/(zion*qe*bmag)
        cotb = cos(beta)/max(sin(beta),1.0e-300_R8)
c       guiding centre motion per radian of gyrophase
        az = rho*cotb*dz
        an = rho*cotb*dn
        a1z = rho*u1z
        a1n = rho*u1n
        a2z = rho*u2z
        a2n = rho*u2n
c       start above the surface (eta: height above the staircase)
        cs = cv(3)-dc+dc*tile_halton(is,2)
        psi0 = 2.0_R8*pi*tile_halton(is,3)
        lam = (etamax+rho*ucut*(1.0_R8+abs(tm/tp))+1.0e-6_R8*tp
     .    -(cs*nhn-tm/tp*cs*nhz))/etad
        gz = cs*nhz+lam*dhz
        gn = cs*nhn+lam*dhn
        psi = psi0
        found = .false.
        nstep = 0
        do while (.not.found .and. nstep.lt.2000000)
c         tiles within reach of the next gyration
          reach = rho*ucut+2.0_R8*pi*sqrt(az**2+an**2)+tp
          jlo = ceiling((gz-reach-zmax)/tp)
          jhi = floor((gz+reach-zmin)/tp)
c         whole gyrations while every edge stays clear of the orbit
          kmin = huge(1.0_R8)
          do jj = jlo, jhi
            do ie = 1, 4
              ln = sqrt((wz(ie+1)-wz(ie))**2+(wy(ie+1)-wy(ie))**2)
              tz = (wz(ie+1)-wz(ie))/ln
              tn = (wy(ie+1)-wy(ie))/ln
              ez = -tn
              en = tz
              hz = wz(ie)+real(jj,R8)*tp
              hn = wy(ie)+real(jj,R8)*tm
              sl = ez*(gz-hz)+en*(gn-hn)
     .          -sqrt((ez*a1z+en*a1n)**2+(ez*a2z+en*a2n)**2)
              kl = 0.0_R8
              if (sl.gt.0.0_R8) then
                rate = 2.0_R8*pi*(ez*az+en*an)
                if (rate.lt.0.0_R8) then
                  kl = sl/(-rate)
                else
                  kl = huge(1.0_R8)
                endif
              endif
              pmz = hz+0.5_R8*ln*tz
              pmn = hn+0.5_R8*ln*tn
              tr = tz*(gz-pmz)+tn*(gn-pmn)
              sx = abs(tr)-0.5_R8*ln
     .          -sqrt((tz*a1z+tn*a1n)**2+(tz*a2z+tn*a2n)**2)
              kx = 0.0_R8
              if (sx.gt.0.0_R8) then
                rate = 2.0_R8*pi*(tz*az+tn*an)*sign(1.0_R8,tr)
                if (rate.lt.0.0_R8) then
                  kx = sx/(-rate)
                else
                  kx = huge(1.0_R8)
                endif
              endif
c             distance of the guiding centre to the edge
              tr = max(-0.5_R8*ln,min(0.5_R8*ln,tr))
              sx = sqrt((gz-pmz-tr*tz)**2+(gn-pmn-tr*tn)**2)-rho*ucut
              kd = 0.0_R8
              if (sx.gt.0.0_R8) kd = sx/(2.0_R8*pi*sqrt(az**2+an**2)
     .          +1.0e-300_R8)
              ke = max(kl,kx,kd)
              kmin = min(kmin,ke)
            enddo
          enddo
          if (kmin.ge.1.0_R8) then
            if (az.gt.0.0_R8) kmin = min(kmin,tp/(2.0_R8*pi*az))
            kmin = min(kmin,1.0e9_R8)
            if (kmin.ge.1.0_R8) then
              k = int(kmin)
              gz = gz+2.0_R8*pi*real(k,R8)*az
              gn = gn+2.0_R8*pi*real(k,R8)*an
              nstep = nstep+1
              cycle
            endif
          endif
c         one gyration in steps of 1/nfine: the first crossing of each
c         edge on the orbit, from the analytic minimum of the distance
c         (drift plus sinusoid) in the step
          do i = 1, nfine
            psa = psi+2.0_R8*pi*real(i-1,R8)/real(nfine,R8)
            psb = psi+2.0_R8*pi*real(i,R8)/real(nfine,R8)
            tbest = huge(1.0_R8)
            ehit = 0
            do jj = jlo, jhi
              do ie = 1, 4
                ln = sqrt((wz(ie+1)-wz(ie))**2+(wy(ie+1)-wy(ie))**2)
                tz = (wz(ie+1)-wz(ie))/ln
                tn = (wy(ie+1)-wy(ie))/ln
                ez = -tn
                en = tz
                hz = wz(ie)+real(jj,R8)*tp
                hn = wy(ie)+real(jj,R8)*tm
                d1 = ez*(gz+az*(psa-psi)+a1z*sin(psa)+a2z*cos(psa)-hz)
     .            +en*(gn+an*(psa-psi)+a1n*sin(psa)+a2n*cos(psa)-hn)
                if (d1.lt.0.0_R8) cycle
                d2 = ez*(gz+az*(psb-psi)+a1z*sin(psb)+a2z*cos(psb)-hz)
     .            +en*(gn+an*(psb-psi)+a1n*sin(psb)+a2n*cos(psb)-hn)
                sa = ez*az+en*an
                ca = ez*a1z+en*a1n
                cb = ez*a2z+en*a2n
                cr = sqrt(ca**2+cb**2)
                psn = psb+1.0_R8
                if (cr.gt.abs(sa)) then
                  ph = atan2(cb,ca)+acos(-sa/cr)
                  psm = -ph+2.0_R8*pi*real(ceiling((psa+ph)/
     .              (2.0_R8*pi)),R8)
                  if (psm.gt.psa .and. psm.lt.psb) then
                    sm = ez*(gz+az*(psm-psi)+a1z*sin(psm)
     .                +a2z*cos(psm)-hz)+en*(gn+an*(psm-psi)
     .                +a1n*sin(psm)+a2n*cos(psm)-hn)
                    if (sm.lt.0.0_R8) psn = psm
                  endif
                endif
                if (psn.gt.psb .and. d2.lt.0.0_R8) psn = psb
                if (psn.gt.psb) cycle
c               the root in psa .. psn, by bisection
                pl = psa
                pr = psn
                do it = 1, 60
                  psm = 0.5_R8*(pl+pr)
                  pz1 = gz+az*(psm-psi)+a1z*sin(psm)+a2z*cos(psm)
                  pn1 = gn+an*(psm-psi)+a1n*sin(psm)+a2n*cos(psm)
                  sm = ez*(pz1-hz)+en*(pn1-hn)
                  if (sm.ge.0.0_R8) then
                    pl = psm
                  else
                    pr = psm
                  endif
                enddo
                fr = (tz*(pz1-hz)+tn*(pn1-hn))/ln
                if (fr.ge.0.0_R8 .and. fr.le.1.0_R8 .and.
     .              psm.lt.tbest) then
                  tbest = psm
                  ehit = ie
                  frhit = fr
                endif
              enddo
            enddo
            nstep = nstep+1
            if (ehit.gt.0) then
              call tile_hit(ehit,frhit,frhit,pwr,tlen,depth,nx,xg,
     .          nz,zg,qtop,qlead,qtrail,pw)
              found = .true.
              exit
            endif
          enddo
          if (.not.found) then
            gz = gz+2.0_R8*pi*az
            gn = gn+2.0_R8*pi*an
          endif
        enddo
        if (.not.found) then
c         no crossing found: the orbit on its field line
          nfall = nfall+1
          call tile_first(cs,cv,lv,dc,dl,e,j,frhit)
          call xertst(e.gt.0,'tile_gyro: a field line misses')
          call tile_hit(e,frhit,frhit,pwr,tlen,depth,nx,xg,nz,zg,
     .      qtop,qlead,qtrail,pw)
        endif
      enddo

      return
      end subroutine tile_gyro

c..................................................................
c Loads of the tile of element ibnd (W/m per tile, per unit length of
c the element): qe and qi (W/m^2 of the element) the electron and ion
c parts of the heat flux along the field, qo the rest (radiation,
c neutrals, transport across the field; spread evenly on the top).
c Electrons on the field lines (tile_fieldline); ions with their orbits
c (tile_gyro: ti, amass, zion, bmag, nsamp orbits) if nsamp > 0, else
c on the field lines. Without incidence along the field (sina0 = 0)
c everything is spread on the top. Out: qtop(nx) at the nodes xg along
c the top, qlead and qtrail at the depth points on the sides, pw(5):
c leading side, top, trailing side, floor of the gap (lost to the
c support), and the power of the element over one tile, (qe+qi+qo) per.
c info 1: the tile is entirely in the shadow of its neighbour (only qo
c on it).
      subroutine tile_loads(ibnd,torl,rn,zn,br,bz,bphi,qe,qi,qo,ti,
     .  amass,zion,nsamp,nx,xg,qtop,qlead,qtrail,pw,nfall,info)
      use b2mod_types
      use b2mod_wall
      implicit none
      integer, intent(in) :: ibnd, nsamp, nx
      real(kind=R8), intent(in) :: torl, rn, zn, br, bz, bphi, qe, qi,
     .  qo, ti, amass, zion
      real(kind=R8), intent(out) :: xg(nx), qtop(nx),
     .  qlead(ndepth_nml), qtrail(ndepth_nml), pw(5)
      integer, intent(out) :: nfall, info

      integer i, nz
      real(kind=R8) :: vz(4), vy(4), tp, tm, tlen, depth, sina0, dz,
     .  dn, u1z, u1n, u2z, u2n, zg(ndepth_nml), qs, w, bmag
      external tile_shape, tile_direction, tile_grid, tile_fieldline,
     .  tile_gyro

      nz = ndepth_nml
      zg(1:nz) = target_depth(ibnd,1:nz)
      call tile_shape(ibnd,torl,vz,vy,tp,tm,tlen,depth)
      call tile_grid(ibnd,tlen,nx,xg)
      call tile_direction(rn,zn,br,bz,bphi,sina0,dz,dn,u1z,u1n,u2z,
     .  u2n)
      qtop = 0.0_R8
      qlead = 0.0_R8
      qtrail = 0.0_R8
      pw = 0.0_R8
      nfall = 0
      info = 0
      if (sina0.gt.1.0e-10_R8) then
        qs = qo
        call tile_fieldline(vz,vy,tp,tm,dz,dn,sina0,qe,tlen,depth,
     .    nx,xg,nz,zg,qtop,qlead,qtrail,pw,info)
        if (nsamp.gt.0) then
          bmag = sqrt(br**2+bz**2+bphi**2)
          call tile_gyro(vz,vy,tp,tm,dz,dn,u1z,u1n,u2z,u2n,sina0,qi,
     .      ti,amass,zion,bmag,nsamp,tlen,depth,nx,xg,nz,zg,qtop,
     .      qlead,qtrail,pw,nfall)
        else
          call tile_fieldline(vz,vy,tp,tm,dz,dn,sina0,qi,tlen,depth,
     .      nx,xg,nz,zg,qtop,qlead,qtrail,pw,info)
        endif
      else
        qs = qo+qe+qi
      endif
c     the rest evenly on the top
      do i = 1, nx
        if (i.eq.1) then
          w = 0.5_R8*(xg(2)-xg(1))
        else if (i.eq.nx) then
          w = 0.5_R8*(xg(nx)-xg(nx-1))
        else
          w = 0.5_R8*(xg(i+1)-xg(i-1))
        endif
        qtop(i) = qtop(i)+qs*tp*w/tlen
      enddo
      pw(2) = pw(2)+qs*tp
      pw(5) = (qe+qi+qo)*tp

      return
      end subroutine tile_loads

c..................................................................
c Net heat inflow f(nx,nz) (W/m) into the nodes of the tile of element
c ibnd at temperatures t (nodes xg along the top, the depth points of
c the element below): conduction between neighbours (along the top
c Kappa_perp, in depth Kappa: the harmonic mean of the two nodes over
c their distance, through the control length of the nodes across), the
c loads qtop, qlead, qtrail, radiation and evaporation from the top
c (qrad, qevap) and the contact at the back (heatdiff_contact), as
c heatdiff1D; the sides insulated but for their loads.
      subroutine tile2d_fluxes(ibnd,ns,nx,nz,xg,t,qtop,qlead,qtrail,f)
      use b2mod_types
      use b2mod_wall
      implicit none
      integer, intent(in) :: ibnd, ns, nx, nz
      real(kind=R8), intent(in) :: xg(nx), t(nx,nz), qtop(nx),
     .  qlead(nz), qtrail(nz)
      real(kind=R8), intent(out) :: f(nx,nz)

      integer i, j
      real(kind=R8) :: wx(nx), wz(nz), zg(nz), k(nx,nz), kp(nx,nz),
     .  fl, tc
      real(kind=R8) :: Kappa, Kappa_perp, qrad, qevap, heatdiff_contact
      external Kappa, Kappa_perp, qrad, qevap, heatdiff_contact

      zg(1:nz) = target_depth(ibnd,1:nz)
      tc = backplate_temp(ibnd)
      do i = 1, nx
        if (i.eq.1) then
          wx(i) = 0.5_R8*(xg(2)-xg(1))
        else if (i.eq.nx) then
          wx(i) = 0.5_R8*(xg(nx)-xg(nx-1))
        else
          wx(i) = 0.5_R8*(xg(i+1)-xg(i-1))
        endif
      enddo
      do j = 1, nz
        if (j.eq.1) then
          wz(j) = 0.5_R8*(zg(2)-zg(1))
        else if (j.eq.nz) then
          wz(j) = 0.5_R8*(zg(nz)-zg(nz-1))
        else
          wz(j) = 0.5_R8*(zg(j+1)-zg(j-1))
        endif
        do i = 1, nx
          k(i,j) = Kappa(ibnd,j,t(i,j))
          kp(i,j) = Kappa_perp(ibnd,j,t(i,j))
        enddo
      enddo
      f = 0.0_R8
      do j = 1, nz
        do i = 1, nx-1
          fl = 2.0_R8*kp(i,j)*kp(i+1,j)/((kp(i,j)+kp(i+1,j))*
     .      (xg(i+1)-xg(i)))*(t(i,j)-t(i+1,j))*wz(j)
          f(i,j) = f(i,j)-fl
          f(i+1,j) = f(i+1,j)+fl
        enddo
      enddo
      do j = 1, nz-1
        do i = 1, nx
          fl = 2.0_R8*k(i,j)*k(i,j+1)/((k(i,j)+k(i,j+1))*
     .      (zg(j+1)-zg(j)))*(t(i,j)-t(i,j+1))*wx(i)
          f(i,j) = f(i,j)-fl
          f(i,j+1) = f(i,j+1)+fl
        enddo
      enddo
      do i = 1, nx
        f(i,1) = f(i,1)+qtop(i)-wx(i)*(qrad(ibnd,t(i,1),.true.)
     .    +qevap(ibnd,ns,t(i,1)))
        f(i,nz) = f(i,nz)-wx(i)*heatdiff_contact(ibnd,t(i,nz))*
     .    (t(i,nz)-tc)
      enddo
      do j = 1, nz
        f(1,j) = f(1,j)+qlead(j)
        f(nx,j) = f(nx,j)+qtrail(j)
      enddo

      return
      end subroutine tile2d_fluxes

c..................................................................
c Jacobian -df/dt of the tile at t (tile2d_fluxes) by differences, five
c colours of nodes, in LAPACK band storage ab(3 kb+1, nx nz) for the
c band LU; kb = min(nx,nz), the shorter direction inner.
      subroutine tile2d_jac(ibnd,ns,nx,nz,xg,t,qtop,qlead,qtrail,kb,ab)
      use b2mod_types
      implicit none
      integer, intent(in) :: ibnd, ns, nx, nz, kb
      real(kind=R8), intent(in) :: xg(nx), t(nx,nz), qtop(nx),
     .  qlead(nz), qtrail(nz)
      real(kind=R8), intent(out) :: ab(3*kb+1,nx*nz)

      integer i, j, m, ii, jj, di, dj, row, col
      real(kind=R8) :: f0(nx,nz), f1(nx,nz), tp(nx,nz), del(nx,nz)
      logical zin
      external tile2d_fluxes

      zin = nz.le.nx
      ab = 0.0_R8
      call tile2d_fluxes(ibnd,ns,nx,nz,xg,t,qtop,qlead,qtrail,f0)
      do m = 0, 4
        tp = t
        do j = 1, nz
          do i = 1, nx
            if (mod(i+2*j,5).eq.m) then
              del(i,j) = 1.0e-6_R8*max(abs(t(i,j)),1.0_R8)
              tp(i,j) = t(i,j)+del(i,j)
            endif
          enddo
        enddo
        call tile2d_fluxes(ibnd,ns,nx,nz,xg,tp,qtop,qlead,qtrail,f1)
        do j = 1, nz
          do i = 1, nx
            if (mod(i+2*j,5).ne.m) cycle
            col = merge(j+(i-1)*nz, i+(j-1)*nx, zin)
            do dj = -1, 1
              do di = -1, 1
                if (abs(di)+abs(dj).gt.1) cycle
                ii = i+di
                jj = j+dj
                if (ii.lt.1 .or. ii.gt.nx .or. jj.lt.1 .or. jj.gt.nz)
     .            cycle
                row = merge(jj+(ii-1)*nz, ii+(jj-1)*nx, zin)
                ab(2*kb+1+row-col,col) = -(f1(ii,jj)-f0(ii,jj))/del(i,j)
              enddo
            enddo
          enddo
        enddo
      enddo

      return
      end subroutine tile2d_jac

c..................................................................
c Adds the storage rv Cp(t) to the diagonal of the band Jacobian ab
c (tile2d_jac) and factorises it (sgbtrf).
      subroutine tile2d_factor(ibnd,nx,nz,t,rv,kb,ab,ipiv,ok)
      use b2mod_types
      implicit none
      integer, intent(in) :: ibnd, nx, nz, kb
      real(kind=R8), intent(in) :: t(nx,nz), rv(nx,nz)
      real(kind=R8), intent(inout) :: ab(3*kb+1,nx*nz)
      integer, intent(out) :: ipiv(nx*nz)
      logical, intent(out) :: ok

      integer i, j, col, info
      logical zin
      real(kind=R8) :: Cp
      external Cp, sgbtrf

      zin = nz.le.nx
      do j = 1, nz
        do i = 1, nx
          col = merge(j+(i-1)*nz, i+(j-1)*nx, zin)
          ab(2*kb+1,col) = ab(2*kb+1,col)+rv(i,j)*Cp(ibnd,j,t(i,j))
        enddo
      enddo
      call sgbtrf(nx*nz,nx*nz,kb,kb,ab,3*kb+1,ipiv,info)
      ok = info.eq.0

      return
      end subroutine tile2d_factor

c..................................................................
c Newton correction dt = -A^-1 r with the factors of tile2d_factor.
      subroutine tile2d_solve(nx,nz,kb,ab,ipiv,r,dt,ok)
      use b2mod_types
      implicit none
      integer, intent(in) :: nx, nz, kb, ipiv(nx*nz)
      real(kind=R8), intent(in) :: ab(3*kb+1,nx*nz), r(nx,nz)
      real(kind=R8), intent(out) :: dt(nx,nz)
      logical, intent(out) :: ok

      integer i, j, info
      real(kind=R8) :: rhs(nx*nz)
      logical zin
      external sgbtrs

      zin = nz.le.nx
      do j = 1, nz
        do i = 1, nx
          rhs(merge(j+(i-1)*nz, i+(j-1)*nx, zin)) = -r(i,j)
        enddo
      enddo
      call sgbtrs('N',nx*nz,kb,kb,1,ab,3*kb+1,ipiv,rhs,nx*nz,info)
      ok = info.eq.0
      do j = 1, nz
        do i = 1, nx
          dt(i,j) = rhs(merge(j+(i-1)*nz, i+(j-1)*nx, zin))
        enddo
      enddo

      return
      end subroutine tile2d_solve

c..................................................................
c Newton step of the tile: dt = -J^-1 r for the residual
c r = rv (H(t) - H(t0)) - f(t), J = rv Cp - df/dt (tile2d_jac,
c tile2d_factor, tile2d_solve); rv(nx,nz): rho V / h of the nodes,
c zero when steady (r = -f).
      subroutine tile2d_newton(ibnd,ns,nx,nz,xg,t,rv,qtop,qlead,
     .  qtrail,r,dt,ok)
      use b2mod_types
      implicit none
      integer, intent(in) :: ibnd, ns, nx, nz
      real(kind=R8), intent(in) :: xg(nx), t(nx,nz),
     .  rv(nx,nz), qtop(nx), qlead(nz), qtrail(nz), r(nx,nz)
      real(kind=R8), intent(out) :: dt(nx,nz)
      logical, intent(out) :: ok

      integer kb
      integer, allocatable :: ipiv(:)
      real(kind=R8), allocatable :: ab(:,:)
      external tile2d_jac, tile2d_factor, tile2d_solve

      kb = min(nx,nz)
      allocate(ab(3*kb+1,nx*nz), ipiv(nx*nz))
      call tile2d_jac(ibnd,ns,nx,nz,xg,t,qtop,qlead,qtrail,kb,ab)
      call tile2d_factor(ibnd,nx,nz,t,rv,kb,ab,ipiv,ok)
      if (ok) call tile2d_solve(nx,nz,kb,ab,ipiv,r,dt,ok)
      deallocate(ab, ipiv)

      return
      end subroutine tile2d_newton

c..................................................................
c Control volumes rho V (kg/m, per unit length of the element) of the
c nodes of the tile.
      subroutine tile2d_mass(ibnd,nx,nz,xg,rm)
      use b2mod_types
      use b2mod_wall
      implicit none
      integer, intent(in) :: ibnd, nx, nz
      real(kind=R8), intent(in) :: xg(nx)
      real(kind=R8), intent(out) :: rm(nx,nz)

      integer i, j
      real(kind=R8) :: wx, wz
      real(kind=R8) :: mass_density
      external mass_density

      do j = 1, nz
        if (j.eq.1) then
          wz = 0.5_R8*(target_depth(ibnd,2)-target_depth(ibnd,1))
        else if (j.eq.nz) then
          wz = 0.5_R8*(target_depth(ibnd,nz)-target_depth(ibnd,nz-1))
        else
          wz = 0.5_R8*(target_depth(ibnd,j+1)-target_depth(ibnd,j-1))
        endif
        do i = 1, nx
          if (i.eq.1) then
            wx = 0.5_R8*(xg(2)-xg(1))
          else if (i.eq.nx) then
            wx = 0.5_R8*(xg(nx)-xg(nx-1))
          else
            wx = 0.5_R8*(xg(i+1)-xg(i-1))
          endif
          rm(i,j) = mass_density(ibnd,j)*wx*wz
        enddo
      enddo

      return
      end subroutine tile2d_mass

c..................................................................
c One implicit step h of the tile from t0 (backward Euler, enthalpy
c form: rho V (H(t) - H(t0))/h = f(t) at every node, as heatdiff_be):
c Newton iterations with the Jacobian of the start of the step,
c renewed every ten iterations; ok false if they do not converge.
      subroutine tile2d_be(ibnd,ns,nx,nz,xg,qtop,qlead,qtrail,t0,h,t,
     .  ok)
      use b2mod_types
      implicit none
      integer, intent(in) :: ibnd, ns, nx, nz
      real(kind=R8), intent(in) :: xg(nx), qtop(nx), qlead(nz),
     .  qtrail(nz), t0(nx,nz), h
      real(kind=R8), intent(out) :: t(nx,nz)
      logical, intent(out) :: ok

      integer i, j, it, kb
      integer, allocatable :: ipiv(:)
      real(kind=R8), allocatable :: ab(:,:)
      real(kind=R8) :: rv(nx,nz), r(nx,nz), f(nx,nz), dt(nx,nz), step
      logical oknt
      real(kind=R8) :: heatdiff_enthalpy
      external heatdiff_enthalpy, tile2d_fluxes, tile2d_jac,
     .  tile2d_factor, tile2d_solve, tile2d_mass

      kb = min(nx,nz)
      allocate(ab(3*kb+1,nx*nz), ipiv(nx*nz))
      call tile2d_mass(ibnd,nx,nz,xg,rv)
      rv = rv/h
      t = t0
      ok = .false.
      do it = 1, 60
        if (mod(it-1,10).eq.0) then
          call tile2d_jac(ibnd,ns,nx,nz,xg,t,qtop,qlead,qtrail,kb,ab)
          call tile2d_factor(ibnd,nx,nz,t,rv,kb,ab,ipiv,oknt)
          if (.not.oknt) exit
        endif
        call tile2d_fluxes(ibnd,ns,nx,nz,xg,t,qtop,qlead,qtrail,f)
        do j = 1, nz
          do i = 1, nx
            r(i,j) = rv(i,j)*heatdiff_enthalpy(ibnd,j,t0(i,j),t(i,j))
     .        -f(i,j)
          enddo
        enddo
        call tile2d_solve(nx,nz,kb,ab,ipiv,r,dt,oknt)
        if (.not.oknt) exit
        step = maxval(abs(dt))
        if (.not.(step.lt.1.0e30_R8)) exit
        if (step.gt.500.0_R8) dt = dt*500.0_R8/step
        t = t+dt
        if (minval(t).le.0.0_R8) exit
        if (step.le.1.0e-10_R8*max(maxval(abs(t)),1.0_R8)) then
          ok = .true.
          exit
        endif
      enddo
      deallocate(ab, ipiv)

      return
      end subroutine tile2d_be

c..................................................................
c Temperatures t(nx,nz) of the tile of element ibnd after dt (s) under
c the loads qtop, qlead, qtrail (tile_loads), from t: implicit steps
c as heatdiff_implicit (each step whole and in two halves, the next
c step from their difference to 1e-2 K, the enthalpies extrapolated:
c second order, the energy exact), capped at the vaporisation
c temperature.
      subroutine tile2d_transient(ibnd,ns,nx,xg,qtop,qlead,qtrail,dt,
     .  t,output)
      use b2mod_types
      use b2mod_wall
      implicit none
      integer, intent(in) :: ibnd, ns, nx, output
      real(kind=R8), intent(in) :: xg(nx), qtop(nx), qlead(ndepth_nml),
     .  qtrail(ndepth_nml), dt
      real(kind=R8), intent(inout) :: t(nx,ndepth_nml)

      real(kind=R8) :: tol
      parameter (tol=1.0e-2_R8)
      integer i, j, it, nz, nstep, nrej
      real(kind=R8) :: uf(nx,ndepth_nml), uh(nx,ndepth_nml),
     .  u2(nx,ndepth_nml), tt, h, err, hexp, alpha, dh, ue, du, dxm
      logical ok1, ok2, ok3
      real(kind=R8) :: Kappa, Kappa_perp, Cp, mass_density,
     .  heatdiff_enthalpy
      external Kappa, Kappa_perp, Cp, mass_density, heatdiff_enthalpy,
     .  tile2d_be, xertst

      nz = ndepth_nml
      if (dt.le.0.0_R8) return
      hexp = dt
      do j = 1, nz-1
        do i = 1, nx-1
          alpha = max(Kappa(ibnd,j,t(i,j)),Kappa_perp(ibnd,j,t(i,j)))/
     .      (mass_density(ibnd,j)*Cp(ibnd,j,t(i,j)))
          dxm = min(xg(i+1)-xg(i),
     .      target_depth(ibnd,j+1)-target_depth(ibnd,j))
          hexp = min(hexp,0.25_R8*dxm**2/max(alpha,1.0e-30_R8))
        enddo
      enddo
      h = min(dt,10.0_R8*hexp)
      tt = 0.0_R8
      nstep = 0
      nrej = 0
      do while (dt-tt.gt.1.0e-12_R8*dt)
        h = min(h,dt-tt)
        ok3 = .false.
        call tile2d_be(ibnd,ns,nx,nz,xg,qtop,qlead,qtrail,t,h,uf,ok1)
        call tile2d_be(ibnd,ns,nx,nz,xg,qtop,qlead,qtrail,t,
     .    0.5_R8*h,uh,ok2)
        if (ok2) call tile2d_be(ibnd,ns,nx,nz,xg,qtop,qlead,qtrail,uh,
     .    0.5_R8*h,u2,ok3)
        if (.not.(ok1.and.ok2.and.ok3)) then
          h = 0.25_R8*h
          nrej = nrej+1
          call xertst(h.gt.1.0e-14_R8*dt .and. nrej.lt.100000,
     .     'tile2d_transient: no convergence')
          cycle
        endif
        err = maxval(abs(u2-uf))
        if (err.le.tol .or. h.le.1.0e-14_R8*dt) then
          do j = 1, nz
            do i = 1, nx
              dh = 2.0_R8*heatdiff_enthalpy(ibnd,j,t(i,j),u2(i,j))
     .          -heatdiff_enthalpy(ibnd,j,t(i,j),uf(i,j))
              ue = 2.0_R8*u2(i,j)-uf(i,j)
              do it = 1, 30
                du = (dh-heatdiff_enthalpy(ibnd,j,t(i,j),ue))/
     .            Cp(ibnd,j,ue)
                ue = ue+du
                if (abs(du).le.1.0e-12_R8*abs(ue)) exit
              enddo
              if (ue.gt.0.0_R8 .and. abs(du).le.1.0e-9_R8*abs(ue))
     .          then
                t(i,j) = ue
              else
                t(i,j) = u2(i,j)
              endif
            enddo
          enddo
          tt = tt+h
          nstep = nstep+1
          h = h*min(4.0_R8,max(0.5_R8,
     .      0.9_R8*sqrt(tol/max(err,1.0e-300_R8))))
        else
          h = h*max(0.1_R8,0.9_R8*sqrt(tol/err))
          nrej = nrej+1
        endif
      enddo
      if (maxval(t).gt.matsurf_vapor_temp(ibnd)) write(*,'(a,i5,a)')
     . ' tile2d_transient: element ',ibnd,
     . ' above the vaporisation temperature (capped)'
      t = min(t,matsurf_vapor_temp(ibnd))
      if (output.ge.10) write(*,*) 'tile2d_transient: element ',ibnd,
     . ', steps ',nstep,', rejected ',nrej

      return
      end subroutine tile2d_transient

c..................................................................
c Steady temperatures t(nx,nz) of the tile of element ibnd under the
c loads qtop, qlead, qtrail (tile_loads): Newton iterations from t if
c it holds temperatures (t(1,1) > 0), else from the steady depth
c profile of heatdiff1D_steady for the mean load of the top; steps of
c at most 500 K, backtracking. Capped at the vaporisation temperature
c (a message says so).
      subroutine tile2d_steady(ibnd,ns,nx,xg,qtop,qlead,qtrail,t,
     .  output)
      use b2mod_types
      use b2mod_wall
      implicit none
      integer, intent(in) :: ibnd, ns, nx, output
      real(kind=R8), intent(in) :: xg(nx), qtop(nx), qlead(ndepth_nml),
     .  qtrail(ndepth_nml)
      real(kind=R8), intent(inout) :: t(nx,ndepth_nml)

      integer i, nz, it, ib
      real(kind=R8) :: rv(nx,ndepth_nml), r(nx,ndepth_nml),
     .  f(nx,ndepth_nml), dt(nx,ndepth_nml), tn(nx,ndepth_nml),
     .  t1(ndepth_nml), q, rnorm, rnew, lam, step, tol
      logical ok
      external heatdiff1D_steady, tile2d_fluxes, tile2d_newton

      nz = ndepth_nml
      if (t(1,1).le.0.0_R8) then
        q = (sum(qtop)+sum(qlead)+sum(qtrail))/xg(nx)
        t1 = 0.0_R8
        call heatdiff1D_steady(ibnd,ns,q,t1,0)
        do i = 1, nx
          t(i,1:nz) = t1(1:nz)
        enddo
      endif
      rv = 0.0_R8
      tol = 1.0e-10_R8*max(sum(abs(qtop))+sum(abs(qlead))
     .  +sum(abs(qtrail)),1.0e-30_R8)
      call tile2d_fluxes(ibnd,ns,nx,nz,xg,t,qtop,qlead,qtrail,f)
      r = -f
      rnorm = maxval(abs(r))
      ok = rnorm.le.tol
      it = 0
      do while (.not.ok .and. it.lt.200)
        it = it+1
        call tile2d_newton(ibnd,ns,nx,nz,xg,t,rv,qtop,qlead,qtrail,r,
     .    dt,ok)
        if (.not.ok) exit
        ok = .false.
        step = maxval(abs(dt))
        if (step.le.1.0e-11_R8*maxval(abs(t))) then
          t = t+dt
          call tile2d_fluxes(ibnd,ns,nx,nz,xg,t,qtop,qlead,qtrail,f)
          r = -f
          rnorm = maxval(abs(r))
          ok = .true.
          exit
        endif
        lam = min(1.0_R8,500.0_R8/max(step,1.0e-300_R8))
        do ib = 1, 40
          tn = max(t+lam*dt,1.0_R8)
          call tile2d_fluxes(ibnd,ns,nx,nz,xg,tn,qtop,qlead,qtrail,f)
          rnew = maxval(abs(f))
          if (rnew.lt.rnorm .or. rnew.le.tol) exit
          lam = 0.5_R8*lam
        enddo
        t = tn
        r = -f
        rnorm = rnew
        ok = rnorm.le.tol .or. lam*step.le.1.0e-12_R8*maxval(t)
      enddo
      if (.not.ok) write(*,'(a,i5,a,1p,e12.4,a)')
     . ' tile2d_steady: element ',ibnd,
     . ' not converged (residual ',rnorm,' W/m)'
      if (maxval(t).gt.matsurf_vapor_temp(ibnd)) write(*,'(a,i5,a)')
     . ' tile2d_steady: element ',ibnd,
     . ' above the vaporisation temperature (capped)'
      t = min(t,matsurf_vapor_temp(ibnd))
      if (output.ge.10) write(*,*) 'tile2d_steady: element ',ibnd,
     . ', iterations ',it

      return
      end subroutine tile2d_steady
