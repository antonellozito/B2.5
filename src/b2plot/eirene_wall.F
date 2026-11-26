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
c*** molecules and test ions). Flux densities in m^-2 s^-1, energy flux
c*** densities in W m^-2. Totals of EIRENE surfaces are the sums over
c*** the strata of fort.44 (as the resolved segments); the incident
c*** energy on an additional surface is the sum over the strata of the
c*** flux times its mean energy. Additional surfaces sharing an element
c*** add their totals. nz is 0 for an unknown key.
c======================================================================
      subroutine eirene_wall_field(key,w,nz,igrp)
      use b2mod_types
      use b2mod_constants, only : ev
      use b2mod_wall, only : nwall
      use b2mod_b2plot, only : natmi, nmoli, nioni, nlimi, n_ewall,
     .  ew_res, ew_srf, ew_add, wall_area, nstrai_eir, ew_stratum
#ifdef B25_EIRENE
      use eirmod_parmmod, only : nlim
      use eirmod_wneutrals, only : wldna, wldnm, wldra, wldrm, ewlda,
     .  ewldm, wldnek, wldnep, wldpa, wldpm, wldpeb, wldspta, wldsptm,
     .  wlpump, wlarea, wldna_res, wldnm_res, ewlda_res, ewldm_res,
     .  ewldea_res, ewldem_res, ewldt_res, ewldrp_res, wldspta_res,
     .  wldsptm_res, wlpump_res
#endif
      implicit none
      character*(*), intent(in) :: key
      real (kind=R8), intent(out) :: w(nwall,0:*)
      integer, intent(out) :: nz, igrp
      integer e, i, ix, is, ns, kk, k
      logical res
      real (kind=R8) :: a, all
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
      case ('ewldt','ewldrp')
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
c by stratum: not the quantities EIRENE only resolves, nor the pumping
      if (kk.ne.0 .and. (kk.lt.0 .or. kk.gt.ns .or.
     .    key.eq.'ewldea' .or. key.eq.'ewldem' .or. key.eq.'wlpump'))
     .  return
c quantities EIRENE resolves on the faces of the standard surfaces
      select case (key)
      case ('wldna','wldnm','ewlda','ewldm','ewldea','ewldem','ewldt',
     .      'ewldrp','wldspta','wldsptm','wlpump')
        res = allocated(ewldt_res)
      case default
        res = .false.
      end select
c faces of the standard surfaces: their resolved segment, or the total
c of their surface over its area; one stratum (kk): the resolved value
c times the share of the stratum in the total of the surface (EIRENE does
c not resolve the strata on the faces; for an energy flux whose surface
c has no energies by stratum, the share of its incident particles), or
c its total over the area
      do e = 1, n_ewall
        i = ew_res(e)
        ix = ew_srf(e)
        if (res .and. i.gt.0) then
          do is = 0, nz-1
            w(e,is) = face_value(i,is)
            if (kk.gt.0) then
              all = 0.0_R8
              do k = 1, ns
                all = all+surface_total(ix,is,k)
              end do
              if (all.ne.0.0_R8) then
                w(e,is) = w(e,is)*surface_total(ix,is,kk)/all
              else              ! no energies by stratum: by particles
                do k = 1, ns
                  all = all+particles(ix,is,k)
                end do
                if (all.ne.0.0_R8) then
                  w(e,is) = w(e,is)*particles(ix,is,kk)/all
                else
                  w(e,is) = 0.0_R8
                end if
              end if
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
#endif
      return
#ifdef B25_EIRENE
      contains
c value of column is on the resolved segment i
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
      case ('ewldrp')
        face_value = ewldrp_res(i)
      case ('wldspta')
        face_value = wldspta_res(is+1,i)/ev
      case ('wldsptm')
        face_value = wldsptm_res(is+1,i)/ev
      case ('wlpump')
        face_value = wlpump_res(is+1,i)/ev
      end select
      end function face_value
c total of column is on the EIRENE surface ix (zero for the energy of
c the backscattered particles, which EIRENE only resolves), of all the
c strata (ks zero) or of stratum ks
      real (kind=R8) function surface_total(ix,is,ks)
      integer, intent(in) :: ix, is, ks
      integer k, k1, k2, k0
      surface_total = 0.0_R8
      k1 = 1                            ! the strata summed
      k2 = ns
      k0 = 0                            ! the totals EIRENE writes
      if (ks.gt.0) then
        k1 = ks
        k2 = ks
        k0 = ks
      end if
      select case (key)
      case ('wldna')
        surface_total = sum(wldna(ix,is+1,k1:k2))
      case ('wldnm')
        surface_total = sum(wldnm(ix,is+1,k1:k2))
      case ('wldra')
        surface_total = sum(wldra(ix,is+1,k1:k2))
      case ('wldrm')
        surface_total = sum(wldrm(ix,is+1,k1:k2))
      case ('ewlda')
        do k = k1, k2
          surface_total = surface_total+
     .      wldna(ix,is+1,k)*ewlda(ix,is+1,k)*ev
        end do
      case ('ewldm')
        do k = k1, k2
          surface_total = surface_total+
     .      wldnm(ix,is+1,k)*ewldm(ix,is+1,k)*ev
        end do
      case ('ewldt')
        surface_total = sum(wldnek(ix,k1:k2))+sum(wldnep(ix,k1:k2))
      case ('wldpa')
        surface_total = wldpa(ix,is+1,k0)
      case ('wldpm')
        surface_total = wldpm(ix,is+1,k0)
      case ('ewldrp')
        surface_total = wldpeb(ix,k0)
      case ('wldspta')
        surface_total = wldspta(ix,is+1,k0)
      case ('wldsptm')
        surface_total = wldsptm(ix,is+1,k0)
      case ('wlpump')
        surface_total = wlpump(is+1,ix)
      end select
      end function surface_total
c incident particles of an energy flux on the EIRENE surface ix, of
c stratum ks (zero for the other quantities)
      real (kind=R8) function particles(ix,is,ks)
      integer, intent(in) :: ix, is, ks
      particles = 0.0_R8
      select case (key)
      case ('ewlda')
        particles = wldna(ix,is+1,ks)
      case ('ewldm')
        particles = wldnm(ix,is+1,ks)
      end select
      end function particles
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
