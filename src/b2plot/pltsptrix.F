c
      subroutine pltsptrix(physical)
      use b2mod_types
      use b2mod_b2plot
      use b2us_data
      implicit none
      logical physical
      integer i,j,iCv1,iCv2,iFc,iFs,iVx1,iVx2,iXpt
      logical, allocatable :: onsep(:)

#ifndef PLPLOT
      call lastd
      if(scol .gt. 1) call gsplci(scol)
      call gslwsc_d(sepwidth)
#endif
c     the flux surfaces of the separatrices: every one with a face at an
c     X-point (an X-point is on several, e.g. on the two halves of a
c     double null grid with flux surfaces of their own), the second
c     separatrix, without X-point the separatrix
      allocate(onsep(mpg%nFs))
      onsep = .false.
      do iXpt = 1, mpg%nXpt
        iVx1 = mpg%Xpt(iXpt)
        if (mpg%vxFs(iVx1).gt.0) onsep(mpg%vxFs(iVx1)) = .true.
        do j = 1, mpg%vxFcP(iVx1,2)
          iFc = mpg%vxFc(mpg%vxFcP(iVx1,1)+j-1)
          if (mpg%fcFs(iFc).gt.0) onsep(mpg%fcFs(iFc)) = .true.
        end do
      end do
      if (mpg%nXpt.ge.2.and.mpg%iFssep2.gt.0)
     &  onsep(mpg%iFssep2) = .true.
      if (mpg%nXpt.eq.0.and.mpg%iFssep.gt.0) onsep(mpg%iFssep) = .true.
      if(physical) then
        do iFs = 1, mpg%nFs
          if (.not.onsep(iFs)) cycle
          do i = 1, mpg%fsFcP(iFs,2)
            iFc = mpg%fsFc(mpg%fsFcP(iFs,1)+i-1)
            iVx1 = mpg%fcVx(iFc,1)
            iVx2 = mpg%fcVx(iFc,2)
            call lined_d(geo%vxX(iVx1), geo%vxY(iVx1),
     &                   geo%vxX(iVx2), geo%vxY(iVx2))
          end do
        end do
      else if(mpg%isClassicalGrid.eq.1) then
        do iFs = 1, mpg%nFs
          if (.not.onsep(iFs)) cycle
          do i = 1, mpg%fsFcP(iFs,2)
            iFc = mpg%fsFc(mpg%fsFcP(iFs,1)+i-1)
            iVx1 = mpg%fcVx(iFc,1)
            iVx2 = mpg%fcVx(iFc,2)
            call lined_d(real(vxIx(iVx1),R8)-0.5_R8,
     &                   real(vxIy(iVx1),R8)-0.5_R8,
     &                   real(vxIx(iVx2),R8)-0.5_R8,
     &                   real(vxIy(iVx2),R8)-0.5_R8)
          end do
        end do
        do iFc = 1, mpg%nFc
          iCv1 = mpg%fcCv(iFc,1)
          iCv2 = mpg%fcCv(iFc,2)
          if (fcxIx(iFc).ne.-2.and.fcxIy(iFc).ne.-2) then
            if (cvIx(iCv1).eq.-2.or.cvIy(iCv1).eq.-2.or.
     &          cvIx(iCv2).eq.-2.or.cvIy(iCv2).eq.-2.or.
     &          cvIy(iCv1).ne.cvIy(iCv2).or.(.not.
     &         (cvIx(iCv1).eq.cvIx(iCv2)+1.or.
     &          cvIx(iCv1).eq.cvIx(iCv2)-1))) then
              call lined_d(real(fcxIx(iFc),R8)-0.5_R8,
     &                     real(fcxIy(iFc),R8)-0.5_R8,
     &                     real(fcxIx(iFc),R8)-0.5_R8,
     &                     real(fcxIy(iFc),R8)+0.5_R8)
              if (fcxIx(iFc).eq.cvIx(iCv1)) then
                call lined_d(real(cvIx(iCv2),R8)-0.5_R8,
     &                       real(cvIy(iCv2),R8)-0.5_R8,
     &                       real(cvIx(iCv2),R8)-0.5_R8,
     &                       real(cvIy(iCv2),R8)+0.5_R8)
              end if
              if (fcxIx(iFc).eq.cvIx(iCv2)) then
                call lined_d(real(cvIx(iCv1),R8)+0.5_R8,
     &                       real(cvIy(iCv1),R8)-0.5_R8,
     &                       real(cvIx(iCv1),R8)+0.5_R8,
     &                       real(cvIy(iCv1),R8)+0.5_R8)
              end if
            end if
          end if
          if (fcyIx(iFc).ne.-2.and.fcyIy(iFc).ne.-2) then
            if (cvIx(iCv1).eq.-2.or.cvIy(iCv1).eq.-2.or.
     &          cvIx(iCv2).eq.-2.or.cvIy(iCv2).eq.-2.or.
     &          cvIx(iCv1).ne.cvIx(iCv2).or.(.not.
     &         (cvIy(iCv1).eq.cvIy(iCv2)+1.or.
     &          cvIy(iCv1).eq.cvIy(iCv2)-1))) then
              call lined_d(real(fcyIx(iFc),R8)-0.5_R8,
     &                     real(fcyIy(iFc),R8)-0.5_R8,
     &                     real(fcyIx(iFc),R8)+0.5_R8,
     &                     real(fcyIy(iFc),R8)-0.5_R8)
              if (fcyIy(iFc).eq.cvIy(iCv1)) then
                call lined_d(real(cvIx(iCv2),R8)-0.5_R8,
     &                       real(cvIy(iCv2),R8)-0.5_R8,
     &                       real(cvIx(iCv2),R8)+0.5_R8,
     &                       real(cvIy(iCv2),R8)-0.5_R8)
              end if
              if (fcyIy(iFc).eq.cvIy(iCv2)) then
                call lined_d(real(cvIx(iCv1),R8)-0.5_R8,
     &                       real(cvIy(iCv1),R8)+0.5_R8,
     &                       real(cvIx(iCv1),R8)+0.5_R8,
     &                       real(cvIy(iCv1),R8)+0.5_R8)
              end if
            end if
          end if
        end do
      end if
      deallocate(onsep)
#ifndef PLPLOT
      call gslwsc_d(lwidth)
      call lastd
      if(scol .gt. 1) call gsplci(1)
#endif
      return
      end subroutine pltsptrix
