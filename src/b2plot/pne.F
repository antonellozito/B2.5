c
      subroutine pne()
      use b2mod_types
      use b2us_data
      use b2mod_b2plot
      use b2mod_b2plot_values
      implicit none
      integer is
      character*256 stmp

      do is=0,ndim(2,iwrk)-1
        work(:,is,iwrk)=work(:,is,iwrk)/state%dv%ne(:)
      enddo
      values(9)=minval(work(:,0:ndim(2,iwrk)-1,iwrk))
      values(10)=maxval(work(:,0:ndim(2,iwrk)-1,iwrk))
      stmp=trim(graphlabel(iwrk))//'/(m^{-3}) (FRACTIONAL)'
      call customise_label(stmp,graphlabel(iwrk))
      return
      end subroutine pne
c
      subroutine pwne()
      use b2mod_types
      use b2mod_wall
      use b2us_data
      use b2mod_b2plot
      use b2mod_b2plot_values
      implicit none
      integer i,is,iCv
      character*256 stmp

      iCv=maxval(mpg%fcCv(imapfc(1),1:2))
      values(9)=wall_stack(1,0,iwwrk)/state%dv%ne(iCv)
      values(10)=wall_stack(1,0,iwwrk)/state%dv%ne(iCv)
      do is=0,nwdim(2,iwwrk)-1
        do i=1,nbnd
          iCv=maxval(mpg%fcCv(imapfc(i),1:2))
          wall_stack(i,is,iwwrk)=wall_stack(i,is,iwwrk)/
     .     state%dv%ne(iCv)
          values(9)=min(values(9),wall_stack(i,is,iwwrk))
          values(10)=max(values(10),wall_stack(i,is,iwwrk))
        enddo
      enddo
      stmp=trim(wall_label(iwwrk))//'/(m^{-3})'
      call customise_label(stmp,wall_label(iwwrk))
      return
      end subroutine pwne
