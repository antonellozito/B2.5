      subroutine b2pwrdld
c
c  version : 18.02.98 18:24
c
c======================================================================
c*** Calculates radiation load onto the specified wall sections
c*** The list of the wall segments is to be specified in lwldrad
c======================================================================
      use b2mod_types
      use b2mod_diag
      use b2mod_b2cmpa
      use b2mod_b2plot
      use b2mod_rdload
      use b2mod_neutrals_namelist
      use b2mod_b2plot_wall_loading
      use b2us_data
      use b2mod_dimensions
#ifdef B25_EIRENE
      use eirmod_wneutrals
#endif
      implicit none
      integer i,j,k,l,m,is,iCv
      real (kind=R8) :: u
      external b2pwrdsrc, xertst
#ifdef DBG
      integer lun_w
      parameter (lun_w=95)
      character tab
#endif
c======================================================================
#ifdef DBG
      tab=char(9)
#endif
c*** Fill in the input arrays for rlcomp.
c*** First, translate the information on the wall segments belonging
c*** to the "shadowing structure" and "plot zones" into the format used
c*** in rlcomp. That is, we make a "structure" polygon, which is closed,
c*** and the target list referring to the segments of this polygon.
c
      lrdlsp=.false.
      lrdlcr=.false.
      do i=1,nsgmx !{
        luwrk(i)=0
      end do !}
      nplgs=nwldlis(1)
      nbndry=nplgs
      do i=1,nplgs
        l=iabs(lwldlis(i,1))
        luwrk(l)=i
        xp(i)=pp1(1,l)
        yp(i)=pp1(2,l)
        lbnd(i)=i
        radld(i)=0.0_R8
      end do

#ifdef DBG
      print *,'b2pwrdld: nplgs,nbndry,nwldrad=',nplgs,nbndry,nwldrad
      if (nplgs.gt.0) then
        print *,'lwldlis(i,1):'
        print '(20i5)',(lwldlis(i,1),i=1,nplgs)
        print *,'lbnd(i):'
        print '(20i5)',(lbnd(i),i=1,nplgs)
      endif
      if (maxval(luwrk).gt.0) then
        print *,'luwrk:'
        print '(20i5)',(luwrk(i),i=1,nsgmx)
      endif
      if (maxval(lewrk).gt.0) then
        print *,'lewrk:'
        print '(20i5)',(lewrk(i),i=1,nsgmx)
      endif
#endif
      if (ntrg.eq.0) then
        k=0
        do i=1,nwldrad
          l=lwldrad(i)
          m=luwrk(l)
          if(m .eq. 0) then
#ifdef DBG
            write(*,'(2a,i4,a,i4)')
     ,             'b2pwrdld: radiation target must belong to ',
     ,             'the shadowing structure. Target ',i,',  segment ',l
#endif
            lwldrad(i)=0
          else
            k=k+1
            ltrg(k)=m
          end if
        end do
        ntrg=k
      end if
#ifdef DBG
      if (nplgs.gt.0) then
        print *,'Shadowing structure: ',
     ,           nplgs,nbndry,nwldlis(1),lwldlis(1,1)
        print '(20i5)',(lbnd(i),i=1,nplgs)
        print '(5(1x,2f9.4))',(xp(i),yp(i),i=1,nplgs)
        print '(a,12i5)','Radiation target ',ntrg,
     ,         (nwldlis(i),i=2,nwldprt)
        print '(20i5)',(ltrg(i),i=1,ntrg)

        open(lun_w,file='layout-w.dat')
        write(lun_w,'(2a)') '# Wall layout : ',trim(label)
        write(lun_w,'(3a)') '# R',tab,'Z'
        do i=1,nplgs
          write(lun_w,'(1p,2(1x,e15.8))') xp(i),yp(i)
        end do
        close(lun_w)
      endif
#endif
c
c*** Now specify the radiation sources - in W (remember, the additional
c*** sources are in MW !!!)
c
      call xertst (mpg%nCi+nadsrc.le.nsrcm,
     .  'b2pwrdld: more radiation sources than nsrcm')
      call b2pwrdsrc(0,ssrcm)
      do iCv=1,mpg%nCi
        xysrc(1,iCv)=geo%cvX(iCv)
        xysrc(2,iCv)=geo%cvY(iCv)
      end do
      k=mpg%nCi

#ifdef DBG
      print *,'nadsrc, add_rad = ',nadsrc,add_rad
#endif
      if(nadsrc.gt.0 .and. add_rad.gt.0.0_R8) then
        u=1.0e6_R8*add_rad/nadsrc
        do i=1,nadsrc
          k=k+1
          ssrcm(k)=u
          xysrc(1,k)=addsrc_xy(1,i)
          xysrc(2,k)=addsrc_xy(2,i)
        end do
      end if

      nsrc=k
#ifdef DBG
      u=0.0_R8
      do i=1,nsrc
        u=u+ssrcm(i)
      end do
      print *,'Sources: nsrc, power = ',nsrc,u
#endif
c
c*** Calculate the radiation load
c
      call rlcomp
c      write(*,*) 'b2pwrdld: after rlcomp'
      write(*,*)   '  Total radiation         : ',ttsrc
      if (add_rad.gt.0.0_R8)
     >  write(*,*) '(of which add. core rad.) : ',add_rad*1.0e6
      if (irdparm.eq.1) then
        if (grid_mode.eq.2) then
          write(*,*) '  Reaching the surfaces   : ',tttrg
        else
          write(*,*) '  Reaching the targets    : ',tttrg
        end if
      else if (irdparm.eq.2 .and. nwldprt.gt.1) then
        write(*,*) '  Reaching the plot zones : ',tttrg
      end if
c----------------------------------------------------------------------
c*** Extract the data and put them in the output module.
c*** After rlcomp, the radiation loads are stored in radld (b2mod_rdload)
c*** in the first ntrg elements, in the order given by the ltrg array.
c*** We need them put in wldrdld (b2mod_b2plot_wall_loading.F) in the
c*** positions just after lppp, in the following order:
c***  Western target, cells ny-1 to lytg;
c*** [ inner top target, cells lytg to ny-1;
c***   outer top target, cells ny-1 to lytg ];
c***  Eastern target, cells lytg to ny-1
c*** Normalisation per unit area is also done here.
c
      do i=1,ntrg
        l=ltrg(i)
        wldrdld(lwldlis(l,1))=radld(i)/segm_area(lwldlis(l,1))
      end do

      do i=1,nsgmx !{
        do j=1,nspecies !{
          wldrdlsp(i,j)=0.0_R8
        end do !}
      end do !}

      if(ildparm.eq.99) then !{

c*** Calculate the radiation load split by impurity species from B2

        do j=1,nspecies !{
          write(*,*) '  Radiation split: ',species_list(j)
          do i=1,nplgs !{
            radld(i)=0.0_R8
          end do !}
          call b2pwrdsrc(j,ssrcm)
          do iCv=1,mpg%nCi !{
            xysrc(1,iCv)=geo%cvX(iCv)
            xysrc(2,iCv)=geo%cvY(iCv)
          end do !}
          nsrc=mpg%nCi
          call rlcomp
          write(*,*)   '  Total radiation         : ',ttsrc
          if (irdparm.eq.1) then
            if (grid_mode.eq.2) then
              write(*,*) '  Reaching the surfaces   : ',tttrg
            else
              write(*,*) '  Reaching the targets    : ',tttrg
            end if
          else if (irdparm.eq.2 .and. nwldprt.gt.1) then
            write(*,*) '  Reaching the plot zones : ',tttrg
          end if
          do i=1,ntrg !{
            l=ltrg(i)
            wldrdlsp(lwldlis(l,1),j)=radld(i)/segm_area(lwldlis(l,1))
          end do !}
        end do !}
        lrdlsp=.true.
      end if !}

#ifdef DBG
      print *,'ildparm, nadsrc, add_rad: ',ildparm,nadsrc,add_rad
#endif
      if(ildparm.eq.98 .and. nadsrc.gt.0 .and. add_rad.gt.0.0_R8) then !{

c*** Calculate the radiation load from the core sources

        write(*,*) '  Extra radiation from the core       : ',
     .   1.0e6_R8*add_rad
        do i=1,nplgs !{
          radld(i)=0.0_R8
        end do !}
        k=0
        u=1.0e6_R8*add_rad/nadsrc
        do i=1,nadsrc !{
          k=k+1
          ssrcm(k)=u
          xysrc(1,k)=addsrc_xy(1,i)
          xysrc(2,k)=addsrc_xy(2,i)
        end do !}
        nsrc=k
        call rlcomp
#ifdef DBG
        write(*,*) '  Total additional core radiation     : ',ttsrc
#endif
        write(*,*) '  Core radiation reaching the targets : ',tttrg
        is=1
        do i=1,ntrg !{
          l=ltrg(i)
          wldrdlsp(lwldlis(l,1),is)=radld(i)/segm_area(lwldlis(l,1))
        end do !}
        lrdlcr=.true.
      end if !}
      return
c======================================================================
      end subroutine b2pwrdld
c
c======================================================================
c*** Power radiated in each B2.5 cell (W): line radiation and
c*** bremsstrahlung of the B2.5 species, radiation of the EIRENE atoms,
c*** molecules and test ions, and recombination radiation of the EIRENE
c*** bulk ions, the last only for the species whose line radiation does
c*** not contain it already (fix_recomb and zmax_recomb, as in b2stel).
c*** j = 0: all species, j > 0: those of the isonuclear sequence j.
c======================================================================
      subroutine b2pwrdsrc(j,src)
      use b2mod_types
      use b2mod_diag
      use b2mod_b2cmpa
      use b2mod_b2plot
      use b2mod_b2cmrc, only : fix_recomb, zmax_recomb
      use b2mod_neutrals_namelist
      use b2us_data
#ifdef B25_EIRENE
      use eirmod_wneutrals
#endif
      implicit none
      integer, intent(in) :: j
      real (kind=R8), intent(out) :: src(mpg%nCi)
      integer iCv, is, is0, is1, iatm, imol, iion, ipl
      logical keep
      intrinsic nint, size
c======================================================================
      if (j.eq.0) then
        is0=0
        is1=state%ns-1
      else
        is0=eb2spcr(j)
        is1=eb2spcr(j)+nfluids(j)
      end if
      do iCv=1,mpg%nCi
        src(iCv)=0.0_R8
        do is=is0,is1
          src(iCv)=src(iCv)+state%srw%rqrad(iCv,is)+
     .                      state%srw%rqbrm(iCv,is)
        end do
        do iatm=1,natmi
          if (j.gt.0) then
            if (latmscl(iatm).ne.j) cycle
          end if
          src(iCv)=src(iCv)-eneutrad(iCv,iatm,0)
        end do
        do imol=1,nmoli
          if (j.gt.0) then
            iatm=lmolscl(imol)
            if (iatm.eq.0) cycle
            if (latmscl(iatm).ne.j) cycle
          end if
          src(iCv)=src(iCv)-emolrad(iCv,imol,0)
        end do
        do iion=1,nioni
          if (j.gt.0) then
            iatm=lionscl(iion)
            if (iatm.eq.0) cycle
            if (latmscl(iatm).ne.j) cycle
          end if
          src(iCv)=src(iCv)-eionrad(iCv,iion,0)
        end do
        do is=is0,is1
          ipl=fluids_list(is)
          keep=ipl.ge.1 .and. ipl.le.size(erecrad,2)
          if (keep) keep=.not.(fix_recomb.ne.0 .and.
     .                         nint(zn(is)).le.zmax_recomb)
          if (keep) src(iCv)=src(iCv)-erecrad(iCv,ipl,0)
        end do
      end do
      return
      end subroutine b2pwrdsrc
c
c======================================================================
c*** Radiation load (W/m**2) on the elements of the wall layer, the
c*** segments from (r1,z1) to (r2,z2): power radiated by the B2.5 cells
c*** (b2pwrdsrc) and by the additional sources of param.dg, if present,
c*** in toroidal geometry (rad3d), over the area of each element.
c======================================================================
      subroutine b2pwrdwl(nb,r1,z1,r2,z2,q)
      use b2mod_types
      use b2mod_constants, only : pi
      use b2mod_b2plot_wall_loading
      use b2us_data
      implicit none
      integer, intent(in) :: nb
      real (kind=R8), intent(in) :: r1(nb), z1(nb), r2(nb), z2(nb)
      real (kind=R8), intent(out) :: q(nb)
      integer i, n, nadd
      integer, allocatable :: kt(:)
      real (kind=R8), allocatable :: xy(:,:), ps(:), es(:,:)
      real (kind=R8) :: a
      logical ex
      character*256 ufile
      real (kind=R8), allocatable, save :: qsave(:)
      intrinsic sqrt, sum
      external b2pwinp, b2pwrdsrc, rad3d, find_file
c======================================================================
c the plasma state does not change: computed once
      if (allocated(qsave)) then
        if (size(qsave).eq.nb) then
          q = qsave
          return
        end if
        deallocate(qsave)
      end if
c The additional sources need no shadowing structure in param.dg
      ufile='param.dg'
      call find_file(ufile,ex)
      nadd=0
      if (ex) then
        call b2pwinp
        if (nadsrc.gt.0 .and. add_rad.gt.0.0_R8) nadd=nadsrc
      end if
      n=mpg%nCi+nadd
      allocate(xy(2,n), ps(n), es(4,max(1,nb)), kt(max(1,nb)))
      call b2pwrdsrc(0,ps)
      do i=1,mpg%nCi
        xy(1,i)=geo%cvX(i)
        xy(2,i)=geo%cvY(i)
      end do
      do i=1,nadd
        ps(mpg%nCi+i)=1.0e6_R8*add_rad/nadsrc
        xy(1,mpg%nCi+i)=addsrc_xy(1,i)
        xy(2,mpg%nCi+i)=addsrc_xy(2,i)
      end do
      do i=1,nb
        es(1,i)=r1(i)
        es(2,i)=z1(i)
        es(3,i)=r2(i)
        es(4,i)=z2(i)
        kt(i)=i
      end do
      call rad3d(n,xy,ps,nb,es,nb,kt,q)
      write(*,*)   '  Total radiation         : ',sum(ps)
      if (nadd.gt.0)
     >  write(*,*) '(of which add. core rad.) : ',add_rad*1.0e6_R8
      write(*,*)   '  Reaching the wall       : ',sum(q(1:nb))
      do i=1,nb
        a=pi*(r1(i)+r2(i))*sqrt((r2(i)-r1(i))**2+(z2(i)-z1(i))**2)
        if (a.gt.0.0_R8) then
          q(i)=q(i)/a
        else
          q(i)=0.0_R8
        end if
      end do
      deallocate(xy, ps, es, kt)
      allocate(qsave(nb))
      qsave = q
      return
      end subroutine b2pwrdwl
