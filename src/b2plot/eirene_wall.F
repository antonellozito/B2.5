c======================================================================
c*** EIRENE wall: the wall stack elements extended beyond the B2.5 wall
c*** elements. Elements 1 to nbnd are the B2.5 wall elements (the same
c*** faces, so that EIRENE and B2.5 wall quantities match element by
c*** element), followed by the faces of the non-default standard
c*** surfaces of EIRENE that are not wall elements (the core boundary),
c*** then by the additional surfaces of EIRENE (one element for the
c*** additional surfaces with the same ends). A face of a standard
c*** surface takes the values of its resolved segment of fort.44
c*** (fnmti_b2); the quantities EIRENE does not resolve, and the
c*** additional surfaces, take the total of their EIRENE surfaces over
c*** the area. The extra elements get their ends in r1, z1, r2, z2 after
c*** nbnd, their area in wall_area, imapx and imapy -2.
c======================================================================
      subroutine eirene_wall_init
      use b2mod_types
      use b2mod_subsys
      use b2mod_wall, only : nbnd, nwall, imapfc, imapx, imapy,
     .  r1, z1, r2, z2
      use b2mod_b2plot, only : b2_fnmti, wall_area, p1, p2, nlimi,
     .  nstsi, n_ewall, ew_res, ew_srf, ew_add
      use b2us_data, only : mpg, geo
#ifdef B25_EIRENE
      use eirmod_parmmod, only : nlim
      use eirmod_wneutrals, only : eirdiag_nds_ind, sarea_res, wlarea
      use eirmod_eirbra, only : fnmti_b2
#endif
      implicit none
      logical, save :: done = .false.
      integer n, n0, nshare, e, ib, is, i, iFc, ires, ix
      logical, allocatable :: onwall(:)
      external xertst
      if (done) return
      done = .true.
      n_ewall = 0
#ifdef B25_EIRENE
      call subini ('eirene_wall_init')
      if (.not.allocated(wlarea) .or. nlimi+nstsi.le.0) then
        call subend ()
        return
      end if
      if (allocated(ew_res)) deallocate(ew_res, ew_srf)
      allocate(ew_res(nwall), ew_srf(nwall), onwall(mpg%nFc))
      ew_res = 0
      ew_srf = 0
      onwall = .false.
c the B2.5 wall elements, with their segment and standard surface
      do ib = 1, nbnd
        iFc = imapfc(ib)
        onwall(iFc) = .true.
        ires = 0
        if (allocated(sarea_res)) ires = b2_fnmti(iFc)
        if (ires.le.0) cycle
        do is = 1, nstsi
          if (eirdiag_nds_ind(is).le.0) cycle
          if (ires.gt.eirdiag_nds_ind(is) .and.
     .        ires.lt.eirdiag_nds_ind(is+1)) then
            ew_res(ib) = ires
            ew_srf(ib) = nlim+is
          end if
        end do
      end do
      n = nbnd
c the other faces of the standard surfaces
      if (allocated(sarea_res) .and. allocated(fnmti_b2)) then
        do is = 1, nstsi
          if (eirdiag_nds_ind(is).le.0) cycle
          do i = eirdiag_nds_ind(is)+1, eirdiag_nds_ind(is+1)-1
            iFc = fnmti_b2(i)
            if (iFc.le.0 .or. iFc.gt.mpg%nFc) cycle
            if (onwall(iFc)) cycle
            n = n+1
            call xertst (n.le.nwall,
     .       'eirene_wall_init: more EIRENE wall elements than nwall')
            ew_res(n) = i
            ew_srf(n) = nlim+is
            r1(n) = geo%vxX(mpg%fcVx(iFc,1))
            z1(n) = geo%vxY(mpg%fcVx(iFc,1))
            r2(n) = geo%vxX(mpg%fcVx(iFc,2))
            z2(n) = geo%vxY(mpg%fcVx(iFc,2))
            wall_area(n) = sarea_res(i)
            imapx(n) = -2
            imapy(n) = -2
          end do
        end do
      end if
c the additional surfaces; one with the same ends as an earlier one
c shares its element, whose area is the larger of their areas
      if (allocated(ew_add)) deallocate(ew_add)
      allocate(ew_add(max(1,nlimi)))
      ew_add = 0
      n0 = n
      nshare = 0
      do ix = 1, nlimi
        do e = n0+1, n
          if (same(e,p1(1,ix),p1(2,ix),p2(1,ix),p2(2,ix)) .or.
     .        same(e,p2(1,ix),p2(2,ix),p1(1,ix),p1(2,ix))) then
            ew_add(ix) = e
            wall_area(e) = max(wall_area(e),wlarea(ix))
            nshare = nshare+1
            exit
          end if
        end do
        if (ew_add(ix).gt.0) cycle
        n = n+1
        call xertst (n.le.nwall,
     .   'eirene_wall_init: more EIRENE wall elements than nwall')
        ew_add(ix) = n
        ew_srf(n) = ix
        r1(n) = p1(1,ix)
        z1(n) = p1(2,ix)
        r2(n) = p2(1,ix)
        z2(n) = p2(2,ix)
        wall_area(n) = wlarea(ix)
        imapx(n) = -2
        imapy(n) = -2
      end do
      n_ewall = n
      deallocate(onwall)
      write(*,'(a,i6,a,i6,a,i6,a,i6,a)') ' EIRENE wall:',n_ewall,
     . ' elements (',nbnd,' B2.5 wall elements,',nlimi,
     . ' additional surfaces,',nshare,' sharing an element)'
      call subend ()
#endif
      return
      contains
      logical function same(e,ra,za,rb,zb)
      integer, intent(in) :: e
      real (kind=R8), intent(in) :: ra, za, rb, zb
      real (kind=R8), parameter :: tol = 1.0e-6_R8
      same = abs(r1(e)-ra).lt.tol .and. abs(z1(e)-za).lt.tol .and.
     .       abs(r2(e)-rb).lt.tol .and. abs(z2(e)-zb).lt.tol
      end function same
      end subroutine eirene_wall_init

c======================================================================
c*** Values of the EIRENE wall quantity key on the EIRENE wall elements,
c*** as densities per element area: w(1:n_ewall,0:nz-1), igrp the kind
c*** of its columns (1 atoms, 2 molecules, 3 one total, 4 atoms,
c*** molecules and test ions); nsp the number of B2.5 species. Flux
c*** densities in m^-2 s^-1, energy flux densities in W m^-2. EIRENE
c*** runs every stratum with its full source strength, B2.5 keeps the
c*** fraction of its recycling coefficient (recyceir of the stratum, or
c*** with eirene_use_recyceir 0 recyc of the first ion of the element of
c*** the column, as b2mod_eirene scales the sources): every value of a
c*** stratum is multiplied by it. Totals of EIRENE surfaces are the sums
c*** over the strata (as the resolved segments); the incident energy on
c*** an additional surface is the sum over the strata of the flux times
c*** its mean energy. With the fields of each stratum (fort.44 from
c*** 20261126, wlstra) every value is the sum over the strata of its
c*** values by stratum. Without them a resolved segment holds the sum
c*** over the strata, split among them by the shares of the strata in
c*** the total of its surface (by their incident particles where the
c*** surface has no values by stratum or where they change sign), and a
c*** quantity without values by stratum is scaled only if all the strata
c*** have the same coefficient. ewldk and ewldp are the net kinetic and
c*** the potential parts of ewldt. Additional surfaces sharing an
c*** element add their totals. nz is 0 for an unknown key.
c======================================================================
      subroutine eirene_wall_field(key,w,nz,igrp,nsp)
      use b2mod_types
      use b2mod_constants, only : ev
      use b2mod_wall, only : nwall
      use b2mod_b2plot, only : natmi, nmoli, nioni, nlimi, n_ewall,
     .  ew_res, ew_srf, ew_add, wall_area, nstrai_eir, ew_stratum
#ifdef B25_EIRENE
      use b2mod_b2cmpa, only : is_neutral, lnext
      use b2mod_neutrals_namelist, only : recyc, recyceir, b2espcr,
     .  eb2atcr, lmolscl, lionscl, species_start
      use eirmod_parmmod, only : nlim
      use eirmod_wneutrals, only : wldna, wldnm, wldra, wldrm, ewlda,
     .  ewldm, wldnek, wldnep, wldpa, wldpm, wldpeb, wldspta, wldsptm,
     .  wlpump, wlarea, wldna_res, wldnm_res, ewlda_res, ewldm_res,
     .  ewldea_res, ewldem_res, ewldt_res, ewldrp_res, ewldmr_res,
     .  wldspta_res, wldsptm_res, wlpump_res, wldra_res, wldrm_res,
     .  wldpa_res, wldpm_res, wlstra, wlstra_res, wlstra_field,
     .  nwlstra_res
#endif
      implicit none
      character*(*), intent(in) :: key
      real (kind=R8), intent(out) :: w(nwall,0:*)
      integer, intent(out) :: nz, igrp
      integer, intent(in) :: nsp
      integer e, i, ix, is, ns, kk, k
      logical res, split, same_c, strat
      real (kind=R8) :: a
      real (kind=R8), allocatable :: c(:,:), sh(:)
      integer, save :: use_recyceir = -1
      external ipgeti
      kk = ew_stratum                   ! wstrat: one stratum, once
      ew_stratum = 0
      nz = 0
      igrp = 3
      select case (key)
      case ('wldna','wldra','ewlda','ewldea','wldpa','wldspta')
        nz = natmi
        igrp = 1
      case ('wldnm','wldrm','ewldm','ewldem','wldpm','wldsptm')
        nz = nmoli
        igrp = 2
      case ('ewldt','ewldk','ewldp','ewldrp')
        nz = 1
        igrp = 3
      case ('wlpump')
        nz = natmi+nmoli+nioni
        igrp = 4
      case default
        return
      end select
      if (nz.le.0) return
      w(1:nwall,0:nz-1) = 0.0_R8
#ifdef B25_EIRENE
      ns = max(1,nstrai_eir)
c the fields of each stratum (fort.44 from 20261126)
      strat = allocated(wlstra)
c by stratum: not the quantities without values by stratum
      if (kk.lt.0 .or. kk.gt.ns) return
      if (kk.gt.0 .and. .not.strat .and. (key.eq.'ewldea' .or.
     .    key.eq.'ewldem' .or. key.eq.'wlpump')) return
      if (use_recyceir.lt.0) then
        use_recyceir = 1
        call ipgeti ('eirene_use_recyceir', use_recyceir)
      end if
c the recycling coefficient of each column and stratum
      allocate(c(0:nz-1,ns), sh(ns))
      do k = 1, ns
        do is = 0, nz-1
          c(is,k) = coefficient(column_atom(is,k),k)
        end do
      end do
c quantities EIRENE resolves on the faces of the standard surfaces
      select case (key)
      case ('wldna','wldnm','ewlda','ewldm','ewldea','ewldem','ewldt',
     .      'ewldk','ewldp','ewldrp','wldspta','wldsptm','wlpump')
        res = allocated(ewldt_res)
      case ('wldra','wldrm','wldpa','wldpm')
        res = allocated(wldra_res)
      case default
        res = .false.
      end select
c with the fields of each stratum, the resolved segments of each one
      if (res .and. strat) then
        do k = 1, ns
          if (kk.gt.0 .and. k.ne.kk) cycle
          call wlstra_res(k,.false.)
          do e = 1, n_ewall
            i = ew_res(e)
            if (i.le.0) cycle
            do is = 0, nz-1
              w(e,is) = w(e,is)+face_value(i,is)*c(is,k)
            end do
          end do
        end do
        call wlstra_res(0,.false.)
      end if
c faces of the standard surfaces: their resolved segment, split among
c the strata, or the total of their surface over its area
      do e = 1, n_ewall
        i = ew_res(e)
        ix = ew_srf(e)
        if (res .and. i.gt.0) then
          if (strat) cycle
          do is = 0, nz-1
            split = shares(ix,is)
            same_c = all(c(is,1:ns).eq.c(is,1))
            if (kk.gt.0) then
              if (split) then
                w(e,is) = face_value(i,is)*c(is,kk)*sh(kk)
              else
                w(e,is) = 0.0_R8
              end if
            else if (split) then
              w(e,is) = face_value(i,is)*sum(c(is,1:ns)*sh(1:ns))
            else if (same_c) then
              w(e,is) = face_value(i,is)*c(is,1)
            else
              w(e,is) = face_value(i,is)
              call not_scaled
            end if
          end do
        else if (ix.gt.nlim) then
          a = wlarea(ix)
          if (a.le.0.0_R8) cycle
          do is = 0, nz-1
            w(e,is) = surface_total(ix,is,kk)/a
          end do
        end if
      end do
c additional surfaces: their totals over the area of their element
      if (allocated(ew_add)) then
        do ix = 1, min(nlimi,size(ew_add))
          e = ew_add(ix)
          if (e.le.0) cycle
          a = wall_area(e)
          if (a.le.0.0_R8) cycle
          do is = 0, nz-1
            w(e,is) = w(e,is)+surface_total(ix,is,kk)/a
          end do
        end do
      end if
      deallocate(c, sh)
#endif
      return
#ifdef B25_EIRENE
      contains
c the EIRENE atom whose recycling coefficient column is of stratum k
c takes (0: none)
      integer function column_atom(is,k)
      integer, intent(in) :: is, k
      integer m
      column_atom = 0
      select case (igrp)
      case (1)
        column_atom = is+1
      case (2)
        m = lmolscl(is+1)
        if (m.ge.1 .and. m.le.natmi) column_atom = m
      case (3)                  ! the element the stratum launches
        m = species_start(k)
        if (m.ge.0 .and. m.lt.nsp) column_atom = b2espcr(m)
      case (4)
        if (is.lt.natmi) then
          column_atom = is+1
        else if (is.lt.natmi+nmoli) then
          m = lmolscl(is-natmi+1)
          if (m.ge.1 .and. m.le.natmi) column_atom = m
        else
          m = lionscl(is-natmi-nmoli+1)
          if (m.ge.1 .and. m.le.natmi) column_atom = m
        end if
      end select
      if (column_atom.gt.natmi) column_atom = 0
      end function column_atom
c recycling coefficient of the EIRENE atom ia in stratum k: recyceir, or
c recyc of the first ion of its element (the next species for a
c duplicate atomic species), as eirene_mc of b2mod_eirene; 1 if none
      real (kind=R8) function coefficient(ia,k)
      integer, intent(in) :: ia, k
      integer j, iref
      if (use_recyceir.gt.0) then
        coefficient = recyceir(k)
        return
      end if
      coefficient = 1.0_R8
      if (ia.lt.1) return
      iref = -1
      do j = 0, nsp-1
        if (is_neutral(j)) cycle
        if (b2espcr(j).eq.ia) then
          iref = j
          exit
        end if
      end do
      if (iref.lt.0) then
        j = eb2atcr(ia)
        if (j.ge.0 .and. j+1.lt.nsp) then
          if (lnext(j,j+1)) iref = j+1
        end if
      end if
      if (iref.ge.0) coefficient = recyc(iref,k)
      end function coefficient
c a quantity without values by stratum, whose strata have different
c recycling coefficients: left as EIRENE gives it
      subroutine not_scaled
      logical, save :: told = .false.
      if (told) return
      told = .true.
      write(*,*) 'eirene_wall: ',trim(key),' not scaled by the ',
     .  'recycling coefficients (no values by stratum, the strata ',
     .  'have different ones)'
      end subroutine not_scaled
c value of column is on the resolved segment i (all the strata)
      real (kind=R8) function face_value(i,is)
      integer, intent(in) :: i, is
      face_value = 0.0_R8
      select case (key)
      case ('wldna')
        face_value = wldna_res(is+1,i)
      case ('wldnm')
        face_value = wldnm_res(is+1,i)
      case ('ewlda')
        face_value = ewlda_res(is+1,i)
      case ('ewldm')
        face_value = ewldm_res(is+1,i)
      case ('ewldea')
        face_value = ewldea_res(is+1,i)
      case ('ewldem')
        face_value = ewldem_res(is+1,i)
      case ('ewldt')
        face_value = ewldt_res(i)
      case ('ewldk')
        face_value = ewldt_res(i)-sum(ewldmr_res(1:nmoli,i))
      case ('ewldp')
        face_value = sum(ewldmr_res(1:nmoli,i))
      case ('ewldrp')
        face_value = ewldrp_res(i)
      case ('wldspta')
        face_value = wldspta_res(is+1,i)/ev
      case ('wldsptm')
        face_value = wldsptm_res(is+1,i)/ev
      case ('wlpump')
        face_value = wlpump_res(is+1,i)/ev
      case ('wldra')
        face_value = wldra_res(is+1,i)
      case ('wldrm')
        face_value = wldrm_res(is+1,i)
      case ('wldpa')
        face_value = wldpa_res(is+1,i)
      case ('wldpm')
        face_value = wldpm_res(is+1,i)
      end select
      end function face_value
c value of column is on the EIRENE surface ix in stratum k, as EIRENE
c gives it (zero for the quantities without values by stratum, before
c fort.44 has the fields of each stratum); with a
c single stratum the totals of fort.44 (ngread copies only the first
c ones to the stratum)
      real (kind=R8) function stratum_value(ix,is,k)
      integer, intent(in) :: ix, is, k
      integer kp
      kp = k
      if (ns.eq.1) kp = 0
      stratum_value = 0.0_R8
      select case (key)
      case ('wldna')
        stratum_value = wldna(ix,is+1,k)
      case ('wldnm')
        stratum_value = wldnm(ix,is+1,k)
      case ('wldra')
        stratum_value = wldra(ix,is+1,k)
      case ('wldrm')
        stratum_value = wldrm(ix,is+1,k)
      case ('ewlda')
        stratum_value = wldna(ix,is+1,k)*ewlda(ix,is+1,k)*ev
      case ('ewldm')
        stratum_value = wldnm(ix,is+1,k)*ewldm(ix,is+1,k)*ev
      case ('ewldt')
        stratum_value = wldnek(ix,k)+wldnep(ix,k)
      case ('ewldk')
        stratum_value = wldnek(ix,k)
      case ('ewldp')
        stratum_value = wldnep(ix,k)
      case ('wldpa')
        stratum_value = wldpa(ix,is+1,kp)
      case ('wldpm')
        stratum_value = wldpm(ix,is+1,kp)
      case ('ewldrp')
        stratum_value = wldpeb(ix,kp)
      case ('wldspta')
        stratum_value = wldspta(ix,is+1,kp)
      case ('wldsptm')
        stratum_value = wldsptm(ix,is+1,kp)
      case ('ewldea')
        if (strat) stratum_value = on_surface(1,is+1,ix,k)
      case ('ewldem')
        if (strat) stratum_value = on_surface(2,is+1,ix,k)
      case ('wlpump')
        if (strat) then
          if (is.lt.natmi) then
            stratum_value = on_surface(4,is+1,ix,k)
          else if (is.lt.natmi+nmoli) then
            stratum_value = on_surface(5,is-natmi+1,ix,k)
          else
            stratum_value = on_surface(6,is-natmi-nmoli+1,ix,k)
          end if
        end if
      end select
      end function stratum_value
c value of species j of the m-th field on the surfaces of wlstra (1
c ewldea, 2 ewldem, 3 ewldmr, 4 to 6 the pumped atoms, molecules, test
c ions) on the EIRENE surface ix in stratum k
      real (kind=R8) function on_surface(m,j,ix,k)
      integer, intent(in) :: m, j, ix, k
      integer nsp, nel, off
      call wlstra_field(nwlstra_res+m,nsp,nel,off)
      on_surface = wlstra(off+j+(ix-1)*nsp,k)
      end function on_surface
c total of column is on the EIRENE surface ix, of all the strata (ks
c zero) or of stratum ks, each scaled by its recycling coefficient (zero
c for the energy of the backscattered particles before fort.44 gives it
c on the surfaces, from 20261126)
      real (kind=R8) function surface_total(ix,is,ks)
      integer, intent(in) :: ix, is, ks
      integer k
      surface_total = 0.0_R8
      if (key.eq.'wlpump' .and. .not.strat) then
        if (all(c(is,1:ns).eq.c(is,1))) then
          surface_total = wlpump(is+1,ix)*c(is,1)
        else
          surface_total = wlpump(is+1,ix)
          call not_scaled
        end if
      else if (ks.gt.0) then
        surface_total = stratum_value(ix,is,ks)*c(is,ks)
      else
        do k = 1, ns
          surface_total = surface_total+stratum_value(ix,is,k)*c(is,k)
        end do
      end if
      end function surface_total
c incident particles of column is on the EIRENE surface ix in stratum k:
c the atoms or molecules of the column, or all of them for a total
      real (kind=R8) function particles(ix,is,k)
      integer, intent(in) :: ix, is, k
      particles = 0.0_R8
      select case (igrp)
      case (1)
        particles = wldna(ix,is+1,k)
      case (2)
        particles = wldnm(ix,is+1,k)
      case (3)
        particles = sum(wldna(ix,1:natmi,k))+sum(wldnm(ix,1:nmoli,k))
      end select
      end function particles
c shares sh of the strata in column is on the EIRENE surface ix: their
c values (all of one sign), else their incident particles; false if
c neither gives any
      logical function shares(ix,is)
      integer, intent(in) :: ix, is
      integer k
      real (kind=R8) :: t
      shares = .false.
      sh = 0.0_R8
      if (key.eq.'ewldea' .or. key.eq.'ewldem' .or. key.eq.'wlpump')
     .  return
      do k = 1, ns
        sh(k) = stratum_value(ix,is,k)
      end do
      if (any(sh.gt.0.0_R8) .and. any(sh.lt.0.0_R8)) sh = 0.0_R8
      t = sum(sh)
      if (t.eq.0.0_R8) then
        do k = 1, ns
          sh(k) = particles(ix,is,k)
        end do
        t = sum(sh)
      end if
      if (t.eq.0.0_R8) return
      sh = sh/t
      shares = .true.
      end function shares
#endif
      end subroutine eirene_wall_field

c======================================================================
c*** Push w(1:n_ewall,0:nz-1) on the wall stack (as wall_copy for the
c*** B2.5 wall elements)
c======================================================================
      subroutine eirene_wall_copy(w,nz,texttmp)
      use b2mod_types
      use b2mod_wall, only : nwall
      use b2mod_b2plot, only : wall_stack, nwdim, iwwrk, wklng,
     .  textwwk, n_ewall
      use b2mod_b2plot_values, only : values
      implicit none
      integer, intent(in) :: nz
      real (kind=R8), intent(in) :: w(nwall,0:nz-1)
      character*(*), intent(in) :: texttmp(0:nz-1)
      external dropwallstackbottom
      if (iwwrk.eq.wklng) call dropwallstackbottom()
      iwwrk = iwwrk+1
      wall_stack(1:n_ewall,0:nz-1,iwwrk) = w(1:n_ewall,0:nz-1)
      values(14) = minval(wall_stack(1:n_ewall,0:nz-1,iwwrk))
      values(15) = maxval(wall_stack(1:n_ewall,0:nz-1,iwwrk))
      textwwk(0:nz-1,iwwrk) = texttmp(0:nz-1)
      nwdim(1,iwwrk) = n_ewall
      nwdim(2,iwwrk) = nz
      return
      end subroutine eirene_wall_copy
