c GUI data commands guigeom, guifield, guiefld, guivect and guievect:
c raw grid, overlays and fields for the b2plot.py GUI (inactive without)

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
      external gui_cell_polygon, pltvessel, pltsptrix, pltoutline,
     . gui_send_structure, gui_send_psimap

      if (.not.gui_active) return
      ncell = count(mpg%cvReg(1:mpg%nCv).ne.0)
      allocate(cells(ncell), cnv(ncell), creg(ncell), cix(ncell),
     . ciy(ncell))
      allocate(vl(sum(max(mpg%cvFcP(1:mpg%nCv,2),2))+1))
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
      call gui_send_structure
      xlo = gwin(1)
      xhi = gwin(2)
      ylo = gwin(3)
      yhi = gwin(4)
      logy = ylog
      typeofplot = ityp
      gui_drawn = drawn

      call gui_send_psimap
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
      nv = 0
      if (nf.lt.1) return
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
      call gui_put_i4(1)
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

c Top two EIRENE stack fields (R component below, Z on top) at the
c triangle centroids, as plot_trig_arrows draws them (VECTOR)
      subroutine gui_send_evector
      use b2mod_types
      use b2mod_b2plot
      use b2mod_b2plot_gui
      implicit none
#ifdef B25_EIRENE
      real (kind=R4), allocatable :: buf(:)
      integer :: iz, n, nz, ip, ir

      if (.not.gui_active.or.iewrk.lt.2) return
      ip = iewrk-1
      ir = iewrk
      n = nedim(1,ip)
      nz = min(nedim(2,ip),nedim(2,ir))
      allocate(buf(n))
      call gui_begin(75)
      call gui_put_str(eirene_label(ip))
      call gui_put_str(eirene_label(ir))
      call gui_put_i4(2)
      call gui_put_i4(n)
      call gui_put_i4(nz)
      do iz = 0, nz-1
        call gui_put_str(textewk(iz,ip))
      end do
      do iz = 0, nz-1
        buf(1:n) = real(min(1.0e30_R8,max(-1.0e30_R8,
     .   eirene_stack(1:n,iz,ip))),R4)
        call gui_put_r4v(n,buf)
      end do
      do iz = 0, nz-1
        buf(1:n) = real(min(1.0e30_R8,max(-1.0e30_R8,
     .   eirene_stack(1:n,iz,ir))),R4)
        call gui_put_r4v(n,buf)
      end do
      call gui_put_i4(0)
      call gui_end
      call gui_flush
      deallocate(buf)
#endif
      return
      end subroutine gui_send_evector

c Open a text input file found by find_file (ok false if absent)
      subroutine gui_open_text(name,iu,ok)
      implicit none
      character*(*), intent(in) :: name
      integer, intent(out) :: iu
      logical, intent(out) :: ok
      character*256 filename
      integer ios
#ifdef LEGACYCOMP
      integer newunit
      external newunit
#endif
      external find_file

      filename = name
      call find_file(filename,ok)
      if (.not.ok) return
#ifndef LEGACYCOMP
      open(newunit=iu,file=trim(filename),status='old',action='read',
     . iostat=ios)
#else
      iu = newunit()
      open(unit=iu,file=trim(filename),status='old',action='read',
     . iostat=ios)
#endif
      ok = ios.eq.0
      return
      end subroutine gui_open_text

c Skip lines of unit iu up to and including the first containing key
      subroutine gui_skip_to(iu,key,ios)
      implicit none
      integer, intent(in) :: iu
      character*(*), intent(in) :: key
      integer, intent(out) :: ios
      character*256 line

      ios = 0
      do while (ios.eq.0)
        read(iu,'(a)',iostat=ios) line
        if (ios.eq.0.and.index(line,key).gt.0) return
      end do
      return
      end subroutine gui_skip_to

c Device structure from structure.dat (one polyline per structure) as
c the overlay of role 4; nothing if the file is absent
      subroutine gui_send_structure
      use b2mod_types
      use b2mod_b2plot_gui
      implicit none
      real (kind=R4), allocatable :: x(:), y(:)
      integer :: iu, ios, nstr, nel, n, i, k
      logical :: ok
      external gui_open_text, gui_skip_to

      call gui_open_text('structure.dat',iu,ok)
      if (.not.ok) return
      read(iu,*,iostat=ios) nstr
      if (ios.eq.0) call gui_skip_to(iu,'$structures',ios)
      if (ios.ne.0) then
        call gui_message(2,'structure.dat: no $structures section')
        close(iu)
        return
      end if
      call gui_overlay_begin(4,'Device structure')
      do k = 1, nstr
        call gui_skip_to(iu,'Structure',ios)
        if (ios.eq.0) read(iu,*,iostat=ios) nel
        if (ios.ne.0) exit
        n = abs(nel)
        allocate(x(n), y(n))
        read(iu,*,iostat=ios) (x(i),y(i),i=1,n)
        if (ios.eq.0.and.n.ge.2) then
          call gui_begin(41)
          call gui_put_i4(n)
          call gui_put_r4v(n,x)
          call gui_put_r4v(n,y)
          call gui_end
        end if
        deallocate(x, y)
        if (ios.ne.0) exit
      end do
      call gui_overlay_end
      if (ios.ne.0) call gui_message(2,'structure.dat: read error')
      close(iu)
      return
      end subroutine gui_send_structure

c Normalised poloidal flux on the box of rzpsi.dat (PSIMAP):
c (psi - psi_O)/(psi_S - psi_O), with psi_O at the O-point (extremum of
c psi in the box of the separatrix, refined by a Newton step) and psi_S
c at the primary X-point (or a separatrix vertex without X-point)
      subroutine gui_send_psimap
      use b2mod_types
      use b2us_data
      use b2mod_b2plot_gui
      implicit none
      real (kind=R8), allocatable :: r(:), z(:), psi(:,:)
      real (kind=R4), allocatable :: buf(:)
      real (kind=R8) :: xmin, xmax, ymin, ymax, best, gr, gz, hrr, hzz,
     . hrz, det, du, dv, rs, zs, po, ps, ro, zo, tr, tz
      integer :: iu, ios, nr, nz, i, j, k, ib, jb, iFs, iFc, iVx, n
      logical :: ok
      external gui_open_text, gui_skip_to

      call gui_open_text('rzpsi.dat',iu,ok)
      if (.not.ok) return
      read(iu,*,iostat=ios) nr, nz
      if (ios.ne.0.or.nr.lt.3.or.nz.lt.3) then
        call gui_message(2,'rzpsi.dat: bad dimensions')
        close(iu)
        return
      end if
      allocate(r(nr), z(nz), psi(nr,nz))
      call gui_skip_to(iu,'nr=',ios)
      if (ios.eq.0) read(iu,*,iostat=ios) (r(i),i=1,nr)
      if (ios.eq.0) call gui_skip_to(iu,'nz=',ios)
      if (ios.eq.0) read(iu,*,iostat=ios) (z(j),j=1,nz)
      if (ios.eq.0) call gui_skip_to(iu,'$psi',ios)
      if (ios.eq.0) read(iu,*,iostat=ios) ((psi(i,j),i=1,nr),j=1,nz)
      close(iu)
      if (ios.ne.0) then
        call gui_message(2,'rzpsi.dat: read error')
        deallocate(r, z, psi)
        return
      end if
c Separatrix surface and its point used for psi_S
      if (mpg%nXpt.gt.0) then
        k = 1
        if (allocated(mpg%isPrimaryXpt)) then
          if (size(mpg%isPrimaryXpt).ge.mpg%nXpt) then
            do i = mpg%nXpt, 1, -1
              if (mpg%isPrimaryXpt(i).ne.0) k = i
            end do
          end if
        end if
        iVx = mpg%Xpt(k)
        iFs = mpg%vxFs(iVx)
      else
        iFs = mpg%iFssep
        iVx = 0
      end if
      if (iFs.le.0) then
        deallocate(r, z, psi)
        return
      end if
      xmin = huge(1.0_R8)
      xmax = -huge(1.0_R8)
      ymin = huge(1.0_R8)
      ymax = -huge(1.0_R8)
      do n = 1, mpg%fsFcP(iFs,2)
        iFc = mpg%fsFc(mpg%fsFcP(iFs,1)+n-1)
        do k = 1, 2
          xmin = min(xmin,geo%vxX(mpg%fcVx(iFc,k)))
          xmax = max(xmax,geo%vxX(mpg%fcVx(iFc,k)))
          ymin = min(ymin,geo%vxY(mpg%fcVx(iFc,k)))
          ymax = max(ymax,geo%vxY(mpg%fcVx(iFc,k)))
        end do
        if (iVx.eq.0) iVx = mpg%fcVx(iFc,1)
      end do
c O-point: positive Hessian determinant and smallest gradient
      best = huge(1.0_R8)
      ib = 0
      jb = 0
      do j = 2, nz-1
        do i = 2, nr-1
          if (r(i).lt.xmin.or.r(i).gt.xmax.or.
     .        z(j).lt.ymin.or.z(j).gt.ymax) cycle
          call gui_psi_derivs(nr,nz,r,z,psi,i,j,gr,gz,hrr,hzz,hrz)
          det = hrr*hzz-hrz**2
          if (det.gt.0.0_R8.and.gr**2+gz**2.lt.best) then
            best = gr**2+gz**2
            ib = i
            jb = j
          end if
        end do
      end do
      if (ib.eq.0) then
        call gui_message(2,'rzpsi.dat: no O-point found')
        deallocate(r, z, psi)
        return
      end if
      call gui_psi_derivs(nr,nz,r,z,psi,ib,jb,gr,gz,hrr,hzz,hrz)
      det = hrr*hzz-hrz**2
      du = -(hzz*gr-hrz*gz)/det
      dv = -(hrr*gz-hrz*gr)/det
      du = max(-(r(ib+1)-r(ib-1))/2.0_R8,
     . min((r(ib+1)-r(ib-1))/2.0_R8,du))
      dv = max(-(z(jb+1)-z(jb-1))/2.0_R8,
     . min((z(jb+1)-z(jb-1))/2.0_R8,dv))
      ro = r(ib)+du
      zo = z(jb)+dv
      po = psi(ib,jb)+gr*du+gz*dv
     . +0.5_R8*(hrr*du**2+2.0_R8*hrz*du*dv+hzz*dv**2)
c psi_S by bilinear interpolation at the separatrix point
      rs = geo%vxX(iVx)
      zs = geo%vxY(iVx)
      i = 1
      do while (i.lt.nr-1.and.r(i+1).lt.rs)
        i = i+1
      end do
      j = 1
      do while (j.lt.nz-1.and.z(j+1).lt.zs)
        j = j+1
      end do
      tr = (rs-r(i))/(r(i+1)-r(i))
      tz = (zs-z(j))/(z(j+1)-z(j))
      ps = (1.0_R8-tr)*(1.0_R8-tz)*psi(i,j)+tr*(1.0_R8-tz)*psi(i+1,j)
     . +(1.0_R8-tr)*tz*psi(i,j+1)+tr*tz*psi(i+1,j+1)
      if (ps.eq.po) then
        deallocate(r, z, psi)
        return
      end if
      allocate(buf(max(nr,nz,nr*nz)))
      call gui_begin(77)
      call gui_put_i4(nr)
      call gui_put_i4(nz)
      buf(1:nr) = real(r,R4)
      call gui_put_r4v(nr,buf)
      buf(1:nz) = real(z,R4)
      call gui_put_r4v(nz,buf)
      buf(1:nr*nz) = real(reshape((psi-po)/(ps-po),(/nr*nz/)),R4)
      call gui_put_r4v(nr*nz,buf)
      call gui_put_r4(real(ro,R4))
      call gui_put_r4(real(zo,R4))
      call gui_put_r4(real(rs,R4))
      call gui_put_r4(real(zs,R4))
      call gui_end
      deallocate(r, z, psi, buf)
      return
      end subroutine gui_send_psimap

c Central differences of psi at node (i,j) of the rzpsi.dat grid
      subroutine gui_psi_derivs(nr,nz,r,z,psi,i,j,gr,gz,hrr,hzz,hrz)
      use b2mod_types
      implicit none
      integer, intent(in) :: nr, nz, i, j
      real (kind=R8), intent(in) :: r(nr), z(nz), psi(nr,nz)
      real (kind=R8), intent(out) :: gr, gz, hrr, hzz, hrz
      real (kind=R8) :: hr, hz

      hr = 0.5_R8*(r(i+1)-r(i-1))
      hz = 0.5_R8*(z(j+1)-z(j-1))
      gr = (psi(i+1,j)-psi(i-1,j))/(2.0_R8*hr)
      gz = (psi(i,j+1)-psi(i,j-1))/(2.0_R8*hz)
      hrr = (psi(i+1,j)-2.0_R8*psi(i,j)+psi(i-1,j))/hr**2
      hzz = (psi(i,j+1)-2.0_R8*psi(i,j)+psi(i,j-1))/hz**2
      hrz = (psi(i+1,j+1)-psi(i+1,j-1)-psi(i-1,j+1)+psi(i-1,j-1))
     . /(4.0_R8*hr*hz)
      return
      end subroutine gui_psi_derivs

c Contour levels set by clvs, clvd, clvn or clvp (kind 1) or by their
c Eirene versions (kind 2), for the GUI (LEVELS)
      subroutine gui_send_levels(kind)
      use b2mod_types
      use b2mod_graphics
      use b2mod_b2plot_gui
      implicit none
      integer, intent(in) :: kind

      if (.not.gui_active) return
      call gui_begin(78)
      call gui_put_i4(kind)
      if (kind.eq.1) then
        call gui_put_i4(nofcls)
        call gui_put_r4v(nofcls,clval)
      else
        call gui_put_i4(nofecls)
        call gui_put_r4v(nofecls,eclval)
      end if
      call gui_end
      return
      end subroutine gui_send_levels
