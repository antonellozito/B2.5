c======================================================================
c*** Wall interaction of each incident species on the wall elements,
c*** estimated for runs without the B2.5 wall model (b2stbr_sput_dst
c*** negative) from the plasma of the run, with the data and functions
c*** of B2.5 (b2mod_sputter): physical sputtering yield and energy
c*** yield, reflection of particles and of energy, chemical sputtering
c*** and radiation enhanced sublimation of carbon. The incident energy
c*** is that of the wall label ei, the incidence angle that of the TRIM
c*** column b2stbr_sput_phys_col (default 3, 30 degrees), every
c*** species counts with weight one, without prompt redeposition.
c*** The material and the surface temperature of each element are those
c*** of its EIRENE surface (eirene.input.json, fort.35): the non-default
c*** standard surface of its face if it is a wall (ILIIN > 0), else the
c*** first wall met along the outward normal of the face (the additional
c*** surfaces and the elements already assigned).
c*** Runs with the wall model keep the values they saved. The mean
c*** implantation depth of each incident species (wall_depths) is
c*** computed for every run, on the materials of the elements.
c======================================================================
      subroutine wall_estimates(mpg,geo,ns,rza)
      use b2mod_types
      use b2mod_sputter, only : sput_dst
      use b2us_map
      use b2us_geo
      implicit none
      type (mapping), intent(in) :: mpg
      type (geometry), intent(in) :: geo
      integer, intent(in) :: ns
      real (kind=R8), intent(in) :: rza(mpg%nCv,0:ns-1)
      logical, save :: done = .false.
      external wall_yields_eirene, wall_depths
      if (done) return
      done = .true.
      if (sput_dst.lt.0) call wall_yields_eirene(mpg,geo,ns,rza)
      call wall_depths(ns)
      return
      end subroutine wall_estimates

c======================================================================
c*** Materials, temperatures and yields of the elements from EIRENE
c======================================================================
      subroutine wall_yields_eirene(mpg,geo,ns,rza)
      use b2mod_types
      use b2mod_constants, only : ev, KBolt
      use b2mod_subsys
      use b2mod_b2cmpa, only : zn, am, is_neutral
      use b2mod_wall
      use b2mod_sputter, only : read_sputter_data,
     .  sputter_t, sputter2_t, reflect, reflect2, chmsput, ressput_t,
     .  target_index
      use b2mod_b2plot, only : incident_energy, wall_flux, wall_temp
      use b2us_map
      use b2us_geo
      implicit none
      type (mapping), intent(in) :: mpg
      type (geometry), intent(in) :: geo
      integer, intent(in) :: ns
      real (kind=R8), intent(in) :: rza(mpg%nCv,0:ns-1)
      integer, parameter :: lu = 97, msurf = 4096
      integer nstd, nadd, ib, iFc, is, k, itg, iz, nown, nray, nnone,
     .  col, smodel, rmodel, ios, ntri, it, idum(14)
      integer std_iliin(msurf), add_iliin(msurf), std_z(msurf),
     .  add_z(msurf)
      real (kind=R8) :: std_a(msurf), add_a(msurf), std_t(msurf),
     .  add_t(msurf), add_p(4,msurf)
      integer, allocatable :: fsurf(:), ez(:)
      real (kind=R8), allocatable :: ea(:), et(:)
      real (kind=R8) :: e, zf, fl, y, erel, h1, h2
      logical ok
      character*256 line, fname
      logical open_file
      external open_file, ipgeti, xertst, find_file
c======================================================================
      call subini ('wall_yields_eirene')
      write(*,*) 'Wall interaction estimated from the plasma of the',
     .  ' run (no B2.5 wall model), materials from EIRENE'
c*** EIRENE surfaces
      call eirene_surfaces(msurf,nstd,std_iliin,std_z,std_a,std_t,
     .  nadd,add_iliin,add_z,add_a,add_t,add_p,ok)
      if (.not.ok) then
        write(*,*) 'Wall interaction not estimated: no EIRENE input',
     .    ' in JSON form (eirene.input.json)'
        call subend ()
        return
      end if
c*** Non-default standard surface of each B2.5 face (fort.35)
      allocate(fsurf(mpg%nFc), ez(nbnd), ea(nbnd), et(nbnd))
      fsurf = 0
      fname = 'fort.35'
      call find_file(fname,ok)
      open(lu,file=trim(fname),status='old',action='read',iostat=ios)
      if (ios.eq.0) then
        read(lu,*,iostat=ios) ntri
        do it = 1, ntri
          if (ios.ne.0) exit
          read(lu,*,iostat=ios) idum
          if (ios.ne.0) exit
          do k = 1, 3
            iFc = abs(idum(11+k))
            if (iFc.gt.0 .and. iFc.le.mpg%nFc .and. idum(3*k+1).lt.0)
     .        fsurf(iFc) = -idum(3*k+1)
          end do
        end do
        close(lu)
      else
        write(*,*) 'Wall interaction: no fort.35, materials only',
     .    ' along the outward normals'
      end if
c*** Material and temperature of each element
      ez = 0
      ea = 0.0_R8
      et = 0.0_R8
      nown = 0
      do ib = 1, nbnd
        k = fsurf(imapfc(ib))
        if (k.lt.1 .or. k.gt.nstd) cycle
        if (std_iliin(k).le.0 .or. std_z(k).le.0) cycle
        ez(ib) = std_z(k)
        ea(ib) = std_a(k)
        et(ib) = std_t(k)
        nown = nown+1
      end do
      call wall_along_normals(nbnd,r1,z1,r2,z2,ez,ea,et,nadd,
     .  add_iliin,add_z,add_a,add_t,add_p,nray,nnone)
      write(*,'(a,i6,a,i6,a,i6,a)') ' Wall interaction: ',nown,
     .  ' elements on their EIRENE wall,',nray,
     .  ' facing one,',nnone,' without a material'
c*** Material constituents of the elements, data of B2.5 (the arrays
c*** of the constituents are allocated by the wall model otherwise)
      if (.not.allocated(matsurf_nzconstituents)) then
        allocate(matsurf_nzconstituents(nwall,1),
     .    matsurf_nrelconstituents(nwall,1),
     .    matsurf_relconstituents(nwall,1))
        matsurf_nzconstituents = 0
        matsurf_nrelconstituents = 0
        matsurf_relconstituents = 0.0_R8
      end if
      do ib = 1, nbnd
        if (ez(ib).gt.0) then
          matsurf_nconstituents(ib) = 1
          matsurf_nzconstituents(ib,1) = ez(ib)
          matsurf_nrelconstituents(ib,1) = 1
          matsurf_relconstituents(ib,1) = 1.0_R8
        else
          matsurf_nconstituents(ib) = 0
        end if
        layer_nconstituents(ib) = matsurf_nconstituents(ib)
        layer_nzconstituents(ib,1) = ez(ib)
        layer_nrelconstituents(ib,1) = 1
        layer_relconstituents(ib,1) = 1.0_R8
        wall_temp(ib) = et(ib)
        target_temp(ib,1) = et(ib)
      end do
      erel = 0.0_R8
      if (any(ez(1:nbnd).eq.6)) then
        if (open_file(lu,'Surface_properties/C',
     .      'data.local|modules/B2.5/Database')) then
          do k = 1, 4
            read(lu,'(a)',iostat=ios) line
          end do
          read(lu,*,iostat=ios) erel
          close(lu)
        end if
      end if
      where (ez(1:nbnd).eq.6) matsurf_chemsput_erel(1:nbnd) = erel
      col = 3
      smodel = 1
      rmodel = 1
      call ipgeti ('b2stbr_sput_phys_col', col)
      call ipgeti ('b2stbr_sput_phys_model', smodel)
      call ipgeti ('b2stbr_refl_model', rmodel)
      call read_sputter_data(ns,col,smodel,rmodel)
c*** Yields of each incident species
      physical_sputtering(1:nbnd,0:ns-1) = 0.0_R8
      physical_sputtering_energy(1:nbnd,0:ns-1) = 0.0_R8
      backscattering(1:nbnd,0:ns-1) = 0.0_R8
      backscattering_energy(1:nbnd,0:ns-1) = 0.0_R8
      chemical_sputtering(1:nbnd,0:ns-1) = 0.0_R8
      RES_sputtering(1:nbnd,0:ns-1) = 0.0_R8
      do ib = 1, nbnd
        iz = ez(ib)
        if (iz.le.0) cycle
        itg = target_index(iz)
        if (itg.le.0) cycle
        iFc = imapfc(ib)
        h1 = geo%fcHc(iFc,1)
        h2 = geo%fcHc(iFc,2)
        do is = 0, ns-1
          e = incident_energy(ib,is)
          if (e.le.0.0_R8) cycle
          zf = (rza(mpg%fcCv(iFc,1),is)*h2+rza(mpg%fcCv(iFc,2),is)*h1)/
     .      (h1+h2)
          physical_sputtering(ib,is) =
     .      sputter_t(is,zf,itg,real(iz,R8),ea(ib),e*ev)
          physical_sputtering_energy(ib,is) = sputter2_t(is,itg,e*ev)
          if (.not.is_neutral(is)) then
            backscattering(ib,is) = reflect(is,e*ev,ib)
            backscattering_energy(ib,is) = reflect2(is,e*ev,ib)
          end if
          if (iz.eq.6) then
            fl = max(0.0_R8,wall_flux(ib,is))
            if (nint(zn(is)).eq.1) then
              call chmsput(ib,e,zf,zn(is),6.0_R8,am(is),ea(ib),et(ib),
     .          fl,y)
              chemical_sputtering(ib,is) = y
            end if
            call ressput_t(e,zf,zn(is),6.0_R8,am(is),ea(ib),et(ib),fl,
     .        itg,y)
            RES_sputtering(ib,is) = y
          end if
        end do
      end do
      deallocate(fsurf, ez, ea, et)
      call subend ()
      return
      end subroutine wall_yields_eirene

c======================================================================
c*** Elements without a material take that of the first wall met along
c*** the outward normal of their face (the plasma is on the right of
c*** each element): the additional surfaces of EIRENE that are walls
c*** (ILIIN > 0, coordinates in cm) and the elements already assigned.
c======================================================================
      subroutine wall_along_normals(nb,r1,z1,r2,z2,ez,ea,et,nadd,
     .  add_iliin,add_z,add_a,add_t,add_p,nray,nnone)
      use b2mod_types
      implicit none
      integer, intent(in) :: nb, nadd, add_iliin(*), add_z(*)
      real (kind=R8), intent(in) :: r1(nb), z1(nb), r2(nb), z2(nb),
     .  add_a(*), add_t(*), add_p(4,*)
      integer, intent(inout) :: ez(nb)
      real (kind=R8), intent(inout) :: ea(nb), et(nb)
      integer, intent(out) :: nray, nnone
      integer ib, k, kbest, kind0
      integer, allocatable :: own(:)
      real (kind=R8) :: xm, ym, ux, uy, el, tbest, t, ax, ay, bx, by
      logical hit
      allocate(own(nb))
      own = ez
      nray = 0
      nnone = 0
      do ib = 1, nb
        if (own(ib).gt.0) cycle
        el = sqrt((r2(ib)-r1(ib))**2+(z2(ib)-z1(ib))**2)
        if (el.le.0.0_R8) then
          nnone = nnone+1
          cycle
        end if
        xm = 0.5_R8*(r1(ib)+r2(ib))
        ym = 0.5_R8*(z1(ib)+z2(ib))
        ux = -(z2(ib)-z1(ib))/el
        uy = (r2(ib)-r1(ib))/el
        tbest = huge(1.0_R8)
        kbest = 0
        kind0 = 0
        do k = 1, nadd
          if (add_iliin(k).le.0 .or. add_z(k).le.0) cycle
          ax = 0.01_R8*add_p(1,k)
          ay = 0.01_R8*add_p(2,k)
          bx = 0.01_R8*add_p(3,k)
          by = 0.01_R8*add_p(4,k)
          call ray_hit(xm,ym,ux,uy,ax,ay,bx,by,t,hit)
          if (hit .and. t.lt.tbest) then
            tbest = t
            kbest = k
            kind0 = 1
          end if
        end do
        do k = 1, nb
          if (own(k).le.0 .or. k.eq.ib) cycle
          call ray_hit(xm,ym,ux,uy,r1(k),z1(k),r2(k),z2(k),t,hit)
          if (hit .and. t.lt.tbest) then
            tbest = t
            kbest = k
            kind0 = 2
          end if
        end do
        if (kind0.eq.1) then
          ez(ib) = add_z(kbest)
          ea(ib) = add_a(kbest)
          et(ib) = add_t(kbest)
          nray = nray+1
        else if (kind0.eq.2) then
          ez(ib) = own(kbest)
          ea(ib) = ea(kbest)
          et(ib) = et(kbest)
          nray = nray+1
        else
          nnone = nnone+1
        end if
      end do
      deallocate(own)
      return
      end subroutine wall_along_normals

c Distance t > 0 along the ray (x,y)+t*(ux,uy) to the segment from
c (ax,ay) to (bx,by), hit false if it misses it
      subroutine ray_hit(x,y,ux,uy,ax,ay,bx,by,t,hit)
      use b2mod_types
      implicit none
      real (kind=R8), intent(in) :: x, y, ux, uy, ax, ay, bx, by
      real (kind=R8), intent(out) :: t
      logical, intent(out) :: hit
      real (kind=R8) :: dx, dy, den, s
      dx = bx-ax
      dy = by-ay
      den = ux*dy-uy*dx
      hit = .false.
      t = 0.0_R8
      if (abs(den).le.1.0e-14_R8) return
      t = ((ax-x)*dy-(ay-y)*dx)/den
      s = ((ax-x)*uy-(ay-y)*ux)/den
      hit = t.gt.1.0e-6_R8 .and. s.ge.0.0_R8 .and. s.le.1.0_R8
      return
      end subroutine ray_hit

c======================================================================
c*** Surfaces of the EIRENE input in JSON form: for the non-default
c*** standard surfaces and the additional surfaces (also from the file
c*** INCLUDE names), their type ILIIN, the atomic number and mass of
c*** their material (ZNML = 100*mass+number, directly or through their
c*** SURFMOD) and their temperature (K, from EWALL in eV). The
c*** additional surfaces also give their ends P1, P2 (cm).
c======================================================================
      subroutine eirene_surfaces(msurf,nstd,std_iliin,std_z,std_a,
     .  std_t,nadd,add_iliin,add_z,add_a,add_t,add_p,ok)
      use b2mod_types
      use b2mod_constants, only : ev, KBolt
      implicit none
      integer, intent(in) :: msurf
      integer, intent(out) :: nstd, nadd, std_iliin(msurf),
     .  std_z(msurf), add_iliin(msurf), add_z(msurf)
      real (kind=R8), intent(out) :: std_a(msurf), std_t(msurf),
     .  add_a(msurf), add_t(msurf), add_p(4,msurf)
      logical, intent(out) :: ok
      character(len=:), allocatable :: buf, inc
      integer, parameter :: mmod = 256
      character*80 mname(mmod)
      real (kind=R8) :: mznml(mmod), mewall(mmod)
      integer nmod, a, b, va, vb, i, n, ia(msurf), ib(msurf)
      real (kind=R8) :: x(3), p4(4)
      character*256 fname
      logical found
      external find_file

      nstd = 0
      nadd = 0
c the files are looked for as by find_file (here, one level up, ...)
      fname = 'eirene.input.json'
      call find_file(fname,found)
      call json_read(trim(fname),buf,ok)
      if (.not.ok) return
      a = index(buf,'{')
      b = len(buf)
c surface models
      nmod = 0
      call json_value(buf,a,b,'REFLECTION_MODELS',va,vb)
      if (va.gt.0) then
        call json_value(buf,va,vb,'SURFMODS',i,n)
        if (i.gt.0) then
          call json_items(buf,i,n,mmod,ia,ib,nmod)
          do i = 1, nmod
            call json_string(buf,ia(i),ib(i),'NAME',mname(i))
            call json_number(buf,ia(i),ib(i),'ZNML',mznml(i),found)
            if (.not.found) mznml(i) = 0.0_R8
            call json_number(buf,ia(i),ib(i),'EWALL',mewall(i),found)
            if (.not.found) mewall(i) = 0.0_R8
          end do
        end if
      end if
c non-default standard surfaces
      call json_value(buf,a,b,'NONDEF_STD_SURFACES',va,vb)
      if (va.gt.0) then
        call json_value(buf,va,vb,'SURFACES',i,n)
        if (i.gt.0) then
          call json_items(buf,i,n,msurf,ia,ib,nstd)
          do i = 1, nstd
            call surface(buf,ia(i),ib(i),std_iliin(i),std_z(i),
     .        std_a(i),std_t(i),p4)
          end do
        end if
      end if
c additional surfaces, in this file or the one it includes
      call json_value(buf,a,b,'ADDITIONAL_SURFACES',va,vb)
      if (va.gt.0) then
        call json_string(buf,va,vb,'INCLUDE',fname)
        if (fname.ne.' ') then
          call find_file(fname,found)
          call json_read(trim(fname),inc,found)
          if (found) then
            call json_value(inc,index(inc,'{'),len(inc),
     .        'ADDITIONAL_SURFACES',va,vb)
            if (va.gt.0) then
              call json_value(inc,va,vb,'SURFACES',i,n)
              if (i.gt.0) then
                call json_items(inc,i,n,msurf,ia,ib,nadd)
                do i = 1, nadd
                  call surface(inc,ia(i),ib(i),add_iliin(i),add_z(i),
     .              add_a(i),add_t(i),add_p(1,i))
                end do
              end if
            end if
          end if
        else
          call json_value(buf,va,vb,'SURFACES',i,n)
          if (i.gt.0) then
            call json_items(buf,i,n,msurf,ia,ib,nadd)
            do i = 1, nadd
              call surface(buf,ia(i),ib(i),add_iliin(i),add_z(i),
     .          add_a(i),add_t(i),add_p(1,i))
            end do
          end if
        end if
      end if
      return

      contains

c Type, material, temperature and ends of the surface object buf(i:j)
      subroutine surface(txt,i,j,iliin,iz,am2,tk,p)
      character(len=*), intent(in) :: txt
      integer, intent(in) :: i, j
      integer, intent(out) :: iliin, iz
      real (kind=R8), intent(out) :: am2, tk, p(4)
      real (kind=R8) :: v, znml, ewall
      character*80 smod
      character*200 num
      integer k, kk, l, m
      logical f
      call json_number(txt,i,j,'ILIIN',v,f)
      iliin = 0
      if (f) iliin = nint(v)
      znml = 0.0_R8
      ewall = 0.0_R8
      call json_string(txt,i,j,'SURFMOD',smod)
      if (smod.ne.' ') then
        do k = 1, nmod
          if (trim(adjustl(mname(k))).eq.trim(adjustl(smod))) then
            znml = mznml(k)
            ewall = mewall(k)
            exit
          end if
        end do
      else
        call json_number(txt,i,j,'ZNML',znml,f)
        if (.not.f) znml = 0.0_R8
        call json_number(txt,i,j,'EWALL',ewall,f)
        if (.not.f) ewall = 0.0_R8
      end if
      iz = mod(nint(znml),100)
      am2 = real(nint(znml)/100,R8)
      tk = abs(ewall)*ev/KBolt
      p = 0.0_R8
      do kk = 1, 2
        if (kk.eq.1) then
          call json_value(txt,i,j,'P1',l,m)
        else
          call json_value(txt,i,j,'P2',l,m)
        end if
        if (l.gt.0 .and. m-l.gt.1) then
          num = ' '
          num = txt(l+1:min(m-1,l+len(num)))
          do k = 1, len_trim(num)
            if (iachar(num(k:k)).lt.32) num(k:k) = ' '
          end do
          read(num,*,iostat=k) x(1:2)
          if (k.eq.0) p(2*kk-1:2*kk) = x(1:2)
        end if
      end do
      return
      end subroutine surface

c Whole text file into a string, ok false if it cannot be read
      subroutine json_read(fname,buf,ok)
      character(len=*), intent(in) :: fname
      character(len=:), allocatable, intent(out) :: buf
      logical, intent(out) :: ok
      integer, parameter :: lu = 96
      integer n, ios
      ok = .false.
      inquire(file=fname,exist=ok,size=n)
      if (.not.ok .or. n.le.0) then
        ok = .false.
        return
      end if
      allocate(character(len=n) :: buf)
      open(lu,file=fname,access='stream',form='unformatted',
     .  status='old',action='read',iostat=ios)
      if (ios.eq.0) read(lu,iostat=ios) buf
      close(lu)
      ok = ios.eq.0
      return
      end subroutine json_read

      end subroutine eirene_surfaces


c Extent va:vb of the value of key among the members of the object or
c the elements of the array txt(a:b) (a at its opening bracket), not
c in the objects inside it; va zero if absent
      subroutine json_value(txt,a,b,key,va,vb)
      implicit none
      character(len=*), intent(in) :: txt, key
      integer, intent(in) :: a, b
      integer, intent(out) :: va, vb
      integer i, d, k
      va = 0
      vb = 0
      d = 0
      i = a+1
      do while (i.lt.b)
        select case (txt(i:i))
        case ('"')
          k = index(txt(i+1:b),'"')
          if (k.eq.0) return
          if (d.eq.0 .and. txt(i+1:i+k-1).eq.key) then
            i = i+k+1
            do while (i.lt.b .and. (txt(i:i).eq.' ' .or.
     .        txt(i:i).eq.char(10) .or. txt(i:i).eq.char(13) .or.
     .        txt(i:i).eq.char(9)))
              i = i+1
            end do
            if (txt(i:i).ne.':') cycle
            i = i+1
            do while (i.lt.b .and. (txt(i:i).eq.' ' .or.
     .        txt(i:i).eq.char(10) .or. txt(i:i).eq.char(13) .or.
     .        txt(i:i).eq.char(9)))
              i = i+1
            end do
            va = i
            call json_end(txt,i,b,vb)
            return
          end if
          i = i+k+1
        case ('{','[')
          d = d+1
          i = i+1
        case ('}',']')
          d = d-1
          i = i+1
        case default
          i = i+1
        end select
      end do
      return
      end subroutine json_value

c End vb of the value starting at txt(va) (string, bracket or number)
      subroutine json_end(txt,va,b,vb)
      implicit none
      character(len=*), intent(in) :: txt
      integer, intent(in) :: va, b
      integer, intent(out) :: vb
      integer i, d, k
      character c
      c = txt(va:va)
      if (c.eq.'"') then
        k = index(txt(va+1:b),'"')
        vb = va+k
      else if (c.eq.'{' .or. c.eq.'[') then
        d = 0
        i = va
        do while (i.le.b)
          select case (txt(i:i))
          case ('"')
            k = index(txt(i+1:b),'"')
            if (k.eq.0) exit
            i = i+k
          case ('{','[')
            d = d+1
          case ('}',']')
            d = d-1
            if (d.eq.0) exit
          end select
          i = i+1
        end do
        vb = min(i,b)
      else
        i = va
        do while (i.lt.b .and. index(',}] '//char(10)//char(13),
     .    txt(i:i)).eq.0)
          i = i+1
        end do
        vb = i-1
      end if
      return
      end subroutine json_end

c Extents of the elements of the array txt(a:b)
      subroutine json_items(txt,a,b,mx,ia,ib,n)
      implicit none
      character(len=*), intent(in) :: txt
      integer, intent(in) :: a, b, mx
      integer, intent(out) :: ia(mx), ib(mx), n
      integer i, e
      n = 0
      i = a+1
      do while (i.lt.b .and. n.lt.mx)
        if (index(', '//char(10)//char(13)//char(9),txt(i:i)).gt.0)
     .    then
          i = i+1
          cycle
        end if
        call json_end(txt,i,b,e)
        n = n+1
        ia(n) = i
        ib(n) = e
        i = e+1
      end do
      return
      end subroutine json_items

c Number value of key in the object txt(a:b), found false if absent
      subroutine json_number(txt,a,b,key,v,found)
      use b2mod_types
      implicit none
      character(len=*), intent(in) :: txt, key
      integer, intent(in) :: a, b
      real (kind=R8), intent(out) :: v
      logical, intent(out) :: found
      integer va, vb, ios
      v = 0.0_R8
      call json_value(txt,a,b,key,va,vb)
      found = va.gt.0
      if (.not.found) return
      read(txt(va:vb),*,iostat=ios) v
      found = ios.eq.0
      return
      end subroutine json_number

c String value of key in the object txt(a:b), blank if absent
      subroutine json_string(txt,a,b,key,s)
      implicit none
      character(len=*), intent(in) :: txt, key
      integer, intent(in) :: a, b
      character(len=*), intent(out) :: s
      integer va, vb
      s = ' '
      call json_value(txt,a,b,key,va,vb)
      if (va.gt.0 .and. txt(va:va).eq.'"' .and. vb.gt.va+1)
     .  s = txt(va+1:vb-1)
      return
      end subroutine json_string

c======================================================================
c*** Mean implantation depth (m) of each incident species in the wall
c*** elements: first moment of the depth distribution of the TRIM range
c*** tables (trim.data/range.data, <projectile><target>.1m, in Angstrom)
c*** at the incident energy of ei and the angle column
c*** b2stbr_sput_phys_col, weighted with the constituents of the
c*** element; as the sputtered energy yield, zero below the energies of
c*** the table, its last value above them, zero without a table.
c======================================================================
      subroutine wall_depths(ns)
      use b2mod_types
      use b2mod_subsys
      use b2mod_wall, only : nbnd, layer_nconstituents,
     .  layer_nzconstituents, layer_relconstituents
      use b2mod_elements, only : elements, is_codes
      use b2mod_b2plot, only : incident_energy, implant_depth
      implicit none
      integer, intent(in) :: ns
      integer, parameter :: mtab = 128, mrow = 64, lu = 95
      character*2, save :: tproj(mtab)
      integer, save :: tz(mtab), tn(mtab), ntab = 0
      real (kind=R8), save :: te(mrow,mtab), td(mrow,mtab)
      integer ib, ic, is, iz, k, col, i
      real (kind=R8) :: e, d, f
      external ipgeti
      call subini ('wall_depths')
      col = 3
      call ipgeti ('b2stbr_sput_phys_col', col)
      implant_depth(1:nbnd,0:ns-1) = 0.0_R8
      do ib = 1, nbnd
        do ic = 1, layer_nconstituents(ib)
          iz = layer_nzconstituents(ib,ic)
          if (iz.le.0 .or. iz.gt.size(elements)) cycle
          do is = 0, ns-1
            e = incident_energy(ib,is)
            if (e.le.0.0_R8) cycle
            k = 0
            do i = 1, ntab
              if (tproj(i).eq.is_codes(is) .and. tz(i).eq.iz) k = i
            end do
            if (k.eq.0) then
              if (ntab.ge.mtab) cycle
              ntab = ntab+1
              k = ntab
              tproj(k) = is_codes(is)
              tz(k) = iz
              call range_table(is_codes(is),elements(iz),col,mrow,
     .          te(1,k),td(1,k),tn(k))
            end if
            if (tn(k).le.0) cycle
            if (e.lt.te(1,k)) then
              d = 0.0_R8
            else if (e.ge.te(tn(k),k)) then
              d = td(tn(k),k)
            else
              i = 1
              do while (te(i+1,k).le.e)
                i = i+1
              end do
              f = (e-te(i,k))/(te(i+1,k)-te(i,k))
              d = td(i,k)+f*(td(i+1,k)-td(i,k))
            end if
            implant_depth(ib,is) = implant_depth(ib,is)+
     .        layer_relconstituents(ib,ic)*d*1.0e-10_R8
          end do
        end do
      end do
      call subend ()
      return

      contains

c Energies (eV) and depths (Angstrom) of the column col of the range
c table of projectile proj on target targ, n rows (0 without table)
      subroutine range_table(proj,targ,col,mx,ev,dp,n)
      character*2, intent(in) :: proj, targ
      integer, intent(in) :: col, mx
      real (kind=R8), intent(out) :: ev(mx), dp(mx)
      integer, intent(out) :: n
      character*256 fname, line
      real (kind=R8) :: v(10)
      integer ios
      logical open_file
      external open_file, lower_case
      n = 0
      fname = trim(proj)//trim(targ)//'.1m'
      call lower_case(fname)
      if (.not.open_file(lu,fname,'data.local/trim.data/range.data|'//
     .    'modules/B2.5/Database/trim.data/range.data')) return
      do
        read(lu,'(a)',iostat=ios) line
        if (ios.ne.0) exit
        if (line(1:1).eq.'c' .or. line(1:1).eq.'C') cycle
        read(line,*,iostat=ios) v
        if (ios.ne.0 .or. n.ge.mx) cycle
        n = n+1
        ev(n) = v(1)
        dp(n) = v(1+min(max(col,1),9))
      end do
      close(lu)
      return
      end subroutine range_table

      end subroutine wall_depths
