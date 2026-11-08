c
      subroutine map_chords(filename)

      use b2mod_types
      use b2us_data
      use b2mod_b2plot
      use b2mod_b2plot_values
      implicit none

      character*(*) filename

      real (kind=R8) ::
     *   mach(mpg%nCv)
      character*80 dataname,ccol
      logical open_file,found
      real (kind=R8) ::
     *   rcol,r,ro,
     1   x,x1,x2,xo,y,y1,y2,yo,z,z1,z2,zo
      integer npoints,i,iCv

      if(.not.open_file(90,filename,
     1 'data.local/chords|data/chords/')) return
      read(90,*,err=30,end=30) dataname,npoints
      npoints=abs(npoints)

   10 read(90,*,err=30,end=40) x1,y1,z1,x2,y2,z2,rcol
      write(ccol,'(1p,g12.5)') rcol
      x=x1
      y=y1
      z=z1
      r=sqrt(x**2+y**2)
      mach=0.0_R8
      call init_page(
     1 trim(dataname)//' '//trim(ccol),
     1 values(5),values(6),values(7),values(8),1,.true.,.true.,.true.)
      do i=1,npoints
        xo=x
        yo=y
        zo=z
        ro=r
        x=x1*real(npoints-i,R8)/real(npoints,R8)
     .   +x2*real(i,R8)/real(npoints,R8)
        y=y1*real(npoints-i,R8)/real(npoints,R8)
     .   +y2*real(i,R8)/real(npoints,R8)
        z=z1*real(npoints-i,R8)/real(npoints,R8)
     .   +z2*real(i,R8)/real(npoints,R8)
        r=sqrt(x**2+y**2)
        call chord_find(mpg,geo,
     1     (x+xo)/2.0_R8,(y+yo)/2.0_R8,(z+zo)/2.0_R8,
     2     found,iCv)
        if(found) then
          mach(iCv)=mach(iCv)+1.0_R8/real(npoints)
        endif
        call lined_d(ro,zo,r,z)
      enddo
      call plotpatch(mpg,geo,0,mach,
     1 0.0_R8,1.0_R8,.true.,linlog,1)
      call data_copy(0)
      call additional_elements(mpg%nCv,.true.)
      call end_page
      goto 10
   30 write(*,*) 'Error in reading ',trim(filename)
      return
   40 continue
      close(90)
      return
      end subroutine map_chords

c Contribution of chord ich of set iset (the *.chords files read as for
c b2chords.nc) to each cell (tria false, w(nCv)) or Eirene triangle
c (tria true, w(ntrii)), sampled as chmp and echmp draw it
      subroutine chord_weights(iset,ich,tria,n,w)
      use b2mod_types
      use b2us_data
      use b2mod_chords, only : read_chord_files
#ifdef B25_EIRENE
      use eirmod_ctrig
#endif
      implicit none
      integer iset,ich,n
      logical tria
      real (kind=R8) :: w(n)
      character*80, allocatable :: setname(:), setlabel(:)
      integer, allocatable :: setnp(:), cset(:), cnum(:)
      real (kind=R8), allocatable :: cxyz(:,:,:), cval(:)
      real (kind=R8) :: p(3), q(3), m(3)
      integer :: nset, nch, k, i, np, ic
      logical :: found

      w = 0.0_R8
      call read_chord_files(nset, setname, setlabel, setnp, nch,
     . cset, cnum, cxyz, cval)
      do k = 1, nch
        if (cset(k).ne.iset .or. cnum(k).ne.ich) cycle
        np = max(abs(setnp(iset)),1)
        p = cxyz(1:3,1,k)
        do i = 1, np
          q = cxyz(1:3,1,k)+real(i,R8)/real(np,R8)*
     .     (cxyz(1:3,2,k)-cxyz(1:3,1,k))
          m = 0.5_R8*(p+q)
          found = .false.
          ic = 0
          if (tria) then
#ifdef B25_EIRENE
            call chord_tria_find(ntrii,nrknot,necke,xtrian,ytrian,
     .       m(1),m(2),m(3),found,ic)
#endif
          else
            call chord_find(mpg,geo,m(1),m(2),m(3),found,ic)
          end if
          if (found.and.ic.ge.1.and.ic.le.n) w(ic) = w(ic)+
     .     1.0_R8/real(np,R8)
          p = q
        end do
      end do
      if (allocated(setname)) deallocate(setname, setlabel, setnp)
      if (allocated(cset)) deallocate(cset, cnum, cxyz, cval)
      return
      end subroutine chord_weights
