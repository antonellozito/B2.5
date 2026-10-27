c GUI data commands guigeom, guifield, guiefld and guivect: raw grid,
c overlays and fields for the b2plot.py GUI (inactive without it)

c Grid, X-points, overlays and EIRENE mesh (GEOMETRY, OVERLAY, TRIMESH)
      subroutine gui_send_geometry
      use b2mod_types
      use b2us_data
      use b2mod_b2plot
      use b2mod_graphics
      use b2mod_b2plot_gui
#ifdef B25_EIRENE
      use eirmod_ctrig
#endif
      implicit none
      integer, allocatable :: cells(:), cnv(:), vl(:), creg(:),
     . cix(:), ciy(:)
      real (kind=R4), allocatable :: buf(:)
      real (kind=R4) :: gwin(4)
      real (kind=R8) :: psisep
      integer :: iCv, k, ncell, nlist, ityp
      logical :: ylog, drawn
      external gui_cell_polygon, pltvessel, pltsptrix, pltoutline

      if (.not.gui_active) return
      ncell = count(mpg%cvReg(1:mpg%nCv).ne.0)
      allocate(cells(ncell), cnv(ncell), creg(ncell), cix(ncell),
     . ciy(ncell))
      allocate(vl(sum(mpg%cvFcP(1:mpg%nCv,2))))
      allocate(buf(max(mpg%nVx,ncell)))
      k = 0
      nlist = 0
      do iCv = 1, mpg%nCv
        if (mpg%cvReg(iCv).eq.0) cycle
        k = k+1
        cells(k) = iCv
        call gui_cell_polygon(iCv,cnv(k),vl(nlist+1))
        nlist = nlist+cnv(k)
        creg(k) = mpg%cvReg(iCv)
        cix(k) = -1
        ciy(k) = -1
        if (mpg%isClassicalGrid.eq.1.and.allocated(cvIx)) then
          cix(k) = cvIx(iCv)
          ciy(k) = cvIy(iCv)
        end if
      end do
      psisep = 0.0_R8
      if (mpg%iFssep.gt.0) psisep = geo%fsPsi(mpg%iFssep)

      call gui_begin(70)
      call gui_put_i4(mpg%nCv)
      call gui_put_i4(mpg%nVx)
      call gui_put_i4(mpg%isClassicalGrid)
      buf(1:mpg%nVx) = real(geo%vxX(1:mpg%nVx),R4)
      call gui_put_r4v(mpg%nVx,buf)
      buf(1:mpg%nVx) = real(geo%vxY(1:mpg%nVx),R4)
      call gui_put_r4v(mpg%nVx,buf)
c Normalised poloidal flux as in the psi profiles of b2plot (fy.F)
      call gui_put_i4(merge(1,0,psisep.ne.0.0_R8))
      buf(1:mpg%nVx) = 0.0_R4
      if (psisep.ne.0.0_R8) buf(1:mpg%nVx) =
     . real(geo%vxFpsi(1:mpg%nVx)/psisep,R4)
      call gui_put_r4v(mpg%nVx,buf)
      call gui_put_i4(ncell)
      call gui_put_i4v(ncell,cells)
      call gui_put_i4v(ncell,cnv)
      call gui_put_i4(nlist)
      call gui_put_i4v(nlist,vl)
      buf(1:ncell) = real(geo%cvX(cells),R4)
      call gui_put_r4v(ncell,buf)
      buf(1:ncell) = real(geo%cvY(cells),R4)
      call gui_put_r4v(ncell,buf)
      call gui_put_i4v(ncell,creg)
      call gui_put_i4v(ncell,cix)
      call gui_put_i4v(ncell,ciy)
      call gui_put_i4(mpg%nXpt)
      do k = 1, mpg%nXpt
        call gui_put_r4(real(geo%vxX(mpg%Xpt(k)),R4))
        call gui_put_r4(real(geo%vxY(mpg%Xpt(k)),R4))
      end do
      call gui_end

c Overlays drawn by the b2plot routines, without the window clipping
      gwin = (/ xlo, xhi, ylo, yhi /)
      ylog = logy
      ityp = typeofplot
      drawn = gui_drawn
      xlo = -1.0e30_R4
      xhi = 1.0e30_R4
      ylo = -1.0e30_R4
      yhi = 1.0e30_R4
      logy = .false.
      typeofplot = 1
      call gui_overlay_begin(1,'Device structure')
      call pltvessel
      call gui_overlay_end
      call gui_overlay_begin(2,'Separatrices')
      call pltsptrix(.true.)
      call gui_overlay_end
      call gui_overlay_begin(3,'Grid boundary')
      call pltoutline(.true.)
      call gui_overlay_end
      xlo = gwin(1)
      xhi = gwin(2)
      ylo = gwin(3)
      yhi = gwin(4)
      logy = ylog
      typeofplot = ityp
      gui_drawn = drawn

#ifdef B25_EIRENE
      if (ntrii.gt.0.and.allocated(xtrian)) then
        deallocate(buf)
        allocate(buf(nrknot))
        call gui_begin(73)
        call gui_put_i4(nrknot)
        buf(1:nrknot) = real(xtrian(1:nrknot),R4)
        call gui_put_r4v(nrknot,buf)
        buf(1:nrknot) = real(ytrian(1:nrknot),R4)
        call gui_put_r4v(nrknot,buf)
        call gui_put_i4(ntrii)
        call gui_put_i4v(3*ntrii,necke)
        call gui_end
      end if
#endif
      call gui_flush
      deallocate(cells, cnv, creg, cix, ciy, vl, buf)
      return
      end subroutine gui_send_geometry

c Vertices of cell iCv in boundary order, chained along its faces
      subroutine gui_cell_polygon(iCv,nv,vl)
      use b2us_data
      implicit none
      integer, intent(in) :: iCv
      integer, intent(out) :: nv, vl(*)
      integer :: nf, i, j, iFc, inext
      logical :: used(mpg%cvFcP(iCv,2)), found

      nf = mpg%cvFcP(iCv,2)
      iFc = mpg%cvFc(mpg%cvFcP(iCv,1))
      vl(1) = mpg%fcVx(iFc,1)
      vl(2) = mpg%fcVx(iFc,2)
      nv = 2
      used = .false.
      used(1) = .true.
      do j = 2, nf
        found = .false.
        do i = 2, nf
          if (used(i)) cycle
          iFc = mpg%cvFc(mpg%cvFcP(iCv,1)+i-1)
          if (mpg%fcVx(iFc,1).eq.vl(nv)) then
            inext = mpg%fcVx(iFc,2)
          else if (mpg%fcVx(iFc,2).eq.vl(nv)) then
            inext = mpg%fcVx(iFc,1)
          else
            cycle
          end if
          used(i) = .true.
          found = .true.
          if (inext.ne.vl(1)) then
            nv = nv+1
            vl(nv) = inext
          end if
          exit
        end do
        if (.not.found) exit
      end do
      return
      end subroutine gui_cell_polygon

c Top matrix stack field with its values at the vertices (FIELD)
      subroutine gui_send_field
      use b2mod_types
      use b2us_data
      use b2mod_b2plot
      use b2mod_b2plot_gui
      implicit none
      real (kind=R8), allocatable :: fvx(:)
      real (kind=R4), allocatable :: buf(:)
      integer :: iz, nz
      external intvertex

      if (.not.gui_active.or.iwrk.lt.1) return
      nz = ndim(2,iwrk)
      allocate(fvx(mpg%nVx), buf(max(mpg%nCv,mpg%nVx)))
      call gui_begin(74)
      call gui_put_str(graphlabel(iwrk))
      call gui_put_i4(1)
      call gui_put_i4(mpg%nCv)
      call gui_put_i4(nz)
      do iz = 0, nz-1
        call gui_put_str(textwk(iz,iwrk))
      end do
      do iz = 0, nz-1
        buf(1:mpg%nCv) = real(min(1.0e30_R8,max(-1.0e30_R8,
     .   work(1:mpg%nCv,iz,iwrk))),R4)
        call gui_put_r4v(mpg%nCv,buf)
      end do
      call gui_put_i4(mpg%nVx)
      do iz = 0, nz-1
        call intvertex(mpg%nCv,mpg%nVx,mpg,geo%vxVol,work(1,iz,iwrk),
     .   fvx)
        buf(1:mpg%nVx) = real(min(1.0e30_R8,max(-1.0e30_R8,fvx)),R4)
        call gui_put_r4v(mpg%nVx,buf)
      end do
      call gui_end
      call gui_flush
      deallocate(fvx, buf)
      return
      end subroutine gui_send_field

c Top EIRENE stack field on the triangles (FIELD)
      subroutine gui_send_efield
      use b2mod_types
      use b2mod_b2plot
      use b2mod_b2plot_gui
      implicit none
#ifdef B25_EIRENE
      real (kind=R4), allocatable :: buf(:)
      integer :: iz, n, nz

      if (.not.gui_active.or.iewrk.lt.1) return
      n = nedim(1,iewrk)
      nz = nedim(2,iewrk)
      allocate(buf(n))
      call gui_begin(74)
      call gui_put_str(eirene_label(iewrk))
      call gui_put_i4(2)
      call gui_put_i4(n)
      call gui_put_i4(nz)
      do iz = 0, nz-1
        call gui_put_str(textewk(iz,iewrk))
      end do
      do iz = 0, nz-1
        buf(1:n) = real(min(1.0e30_R8,max(-1.0e30_R8,
     .   eirene_stack(1:n,iz,iewrk))),R4)
        call gui_put_r4v(n,buf)
      end do
      call gui_put_i4(0)
      call gui_end
      call gui_flush
      deallocate(buf)
#endif
      return
      end subroutine gui_send_efield

c Top two matrix stack fields (poloidal below, radial on top) as an
c (R,Z) vector at the cell centres and at the vertices (VECTOR)
      subroutine gui_send_vector
      use b2mod_types
      use b2us_data
      use b2mod_b2plot
      use b2mod_b2plot_gui
      implicit none
      real (kind=R8), allocatable :: vr(:,:), vz(:,:), fvx(:)
      real (kind=R4), allocatable :: buf(:)
      integer :: iCv, iz, nz, ip, ir
      external arrow_rz, intvertex

      if (.not.gui_active.or.iwrk.lt.2) return
      ip = iwrk-1
      ir = iwrk
      nz = min(ndim(2,ip),ndim(2,ir))
      allocate(vr(mpg%nCv,0:nz-1), vz(mpg%nCv,0:nz-1), fvx(mpg%nVx))
      allocate(buf(max(mpg%nCv,mpg%nVx)))
      vr = 0.0_R8
      vz = 0.0_R8
      do iz = 0, nz-1
        do iCv = 1, mpg%nCi
          call arrow_rz(iCv,work(iCv,iz,ip),work(iCv,iz,ir),
     .     vr(iCv,iz),vz(iCv,iz))
        end do
      end do
      call gui_begin(75)
      call gui_put_str(graphlabel(ip))
      call gui_put_str(graphlabel(ir))
      call gui_put_i4(mpg%nCv)
      call gui_put_i4(nz)
      do iz = 0, nz-1
        call gui_put_str(textwk(iz,ip))
      end do
      do iz = 0, nz-1
        buf(1:mpg%nCv) = real(min(1.0e30_R8,max(-1.0e30_R8,
     .   vr(:,iz))),R4)
        call gui_put_r4v(mpg%nCv,buf)
      end do
      do iz = 0, nz-1
        buf(1:mpg%nCv) = real(min(1.0e30_R8,max(-1.0e30_R8,
     .   vz(:,iz))),R4)
        call gui_put_r4v(mpg%nCv,buf)
      end do
      call gui_put_i4(mpg%nVx)
      do iz = 0, nz-1
        call intvertex(mpg%nCv,mpg%nVx,mpg,geo%vxVol,vr(1,iz),fvx)
        buf(1:mpg%nVx) = real(min(1.0e30_R8,max(-1.0e30_R8,fvx)),R4)
        call gui_put_r4v(mpg%nVx,buf)
      end do
      do iz = 0, nz-1
        call intvertex(mpg%nCv,mpg%nVx,mpg,geo%vxVol,vz(1,iz),fvx)
        buf(1:mpg%nVx) = real(min(1.0e30_R8,max(-1.0e30_R8,fvx)),R4)
        call gui_put_r4v(mpg%nVx,buf)
      end do
      call gui_end
      call gui_flush
      deallocate(vr, vz, fvx, buf)
      return
      end subroutine gui_send_vector
