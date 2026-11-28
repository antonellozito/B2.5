C.....................................................................
c Program to solve the heat diffusion equation

c  Last modified by Manoj Warrier on (08-02-2004)

c Using the numerical method written by M. Kundu to
c SOLVE PARABOLIC EQUATIONS USING THE METHOD OF LINES

c.....................................................................

c This is the subroutine which solves for a heat diffusion equation
c of the form:
c              doe T         doe^2 T     doe D   doe T
c             ------- =  D  --------- + ------- ------- + sources/sinks
c              doe t         doe^2 x     doe x   doe x
c T is the temperature, t the time variable, x the space variable and
c D the heat diffusion coefficient.

c i.e., ....
c       DU = (1/rho/Cp) * (K*UXX + gradK*UX + S)
c       DU is the rate of change of temperature
c       UXX is "doe^2 T / doe^2 x"   and UX is "doe T / doe x"
c with Neumann boundary condition at the plasma facing surface side
c i.e.  UX is specified (the heat flux less radiation and evaporation),
c and the contact with the coolant or support at the back
c (heatdiff_contact: the heat flux h (T - Tcbc)). Without NAG
c (NO_NAG) the finite-volume form of lines2c, all depth points
c unknown, integrated implicitly (heatdiff_implicit); heatdiff1D_steady
c gives its steady state.

c Inputs:
c  Q = incident total heat flux (Watts/m^2/sec).
c  Tcbc = Temperature at the cooling end (K).
c  dt = Time step of simulation (sec)

c Outputs:
c  Tp = Target surface temperature (K).
      subroutine heatdiff1D(ibnd,ns0,Q,dt,Tp,output)
      use b2mod_types
      use b2mod_wall
      use b2mod_heatdiff_common
     &    , only: dx, x0, xn, x, Qtot, Tcbc, active, iwall, ns
      use b2mod_subsys
      implicit none
      external lines2c

      integer i
      integer ibnd, ns0, output

      real (kind=R8) :: Q, dt
      real (kind=R8) :: Tp(ndepth_nml)
      real (kind=R8) :: u(7*ndepth_nml), u1(7*ndepth_nml)
      real (kind=R8) :: t0, tn

#ifdef NO_NAG
      external heatdiff_implicit
#else
      integer :: itask, itrace, ind, ifail
      real (kind=R8) :: acc
      external pdedef, bndary, D03PCF
#endif
      integer, save :: ncall = 0

      integer, save :: npde, nw, niw
      integer, allocatable, save :: iw(:)
      real (kind=R8), allocatable, save :: U_nag(:,:), w(:)
#ifdef _OPENMP
!$OMP THREADPRIVATE(ncall, npde, nw, niw, U_nag, iw, w)
#endif
      external xertst

      call subini('heatdiff1D')
      if (ncall.eq.0) then
        ncall = ncall + 1
        npde = 1
        nw = (10+6*npde)*npde*ndepth_nml+
     .       (21+3*npde)*npde+7*ndepth_nml+56
        niw = npde*ndepth_nml+25
        allocate(U_nag(npde,ndepth_nml))
        allocate(iw(niw))
        allocate(w(nw))
      endif

      if (output.ge.10) then
        write(*,*)'Wall element number =',ibnd
        write(*,*)'Total energy flux =',Q,'Joules/m^2/s'
        write(*,*)'Temp at the cold end =',backplate_temp(ibnd),'K'
        write(*,*)'Temp at the warm end =',Tp(1),'K'
        write(*,*)'Thickness of the plate =',plate_thickness(ibnd),'m'
        if (inertial_cooling(ibnd)) then
          write(*,*)'Plate element is inertially cooled'
        else
          write(*,*)'Plate element is actively cooled'
        endif
        write(*,*)'Time step of simulation =',dt,'sec'
      endif

      Qtot = Q
      Tcbc = backplate_temp(ibnd)
      active = .not.inertial_cooling(ibnd)
      t0 = 0.0_R8
      tn = dt
      iwall = ibnd
      ns = ns0

c  Compute positions and intervals
      x(1:ndepth_nml) = target_depth(ibnd,1:ndepth_nml)
      x0 = target_depth(ibnd,1)
      xn = target_depth(ibnd,ndepth_nml)
      do i=1,ndepth_nml-1
        dx(i) = x(i+1)-x(i)
      enddo

c Checking if it is the first call and if it is, call fic(x,u), else
c use previous Tp profile as initial guess.
      if (dt.le.0.0_R8) then
c       Calling subroutine fic(ibnd,x,u) which supplies the
c       initial condition u(x) i.e.. the temperature "u" as a
c       function of the position "x" at ndepth_nml grid points.
        call fic(ibnd,x,u)
        do i=1,ndepth_nml
          Tp(i) = u(i)
        enddo
        call subend ()
        return
      else
        do i=1,ndepth_nml
          u1(i)= Tp(i)
        enddo
      endif

c No. of equations to be solved: one per depth point (the back point
c is in contact with the coolant, not at its temperature).
#ifdef NO_NAG
      call xertst(ndepth_nml.ge.2,'heatdiff1D needs NDEPTH_NML >= 2')
      call heatdiff_implicit(ibnd,ns,ndepth_nml,u1,dt,output)
      do i=1,ndepth_nml
        Tp(i) = min( u1(i), matsurf_vapor_temp(ibnd))
      enddo
#else
cxpb   New NAG solution to replace MSTEP calls
      ind = 0
      ifail = 1
      itask = 1
      itrace = -1
      do i = 1, ndepth_nml
        U_nag(1,i)=Tp(i)
      enddo
      acc = 1.0e-8_R8

      call D03PCF(npde, 0, t0, tn, pdedef,
     &   bndary, U_nag, ndepth_nml, x, acc, w, nw, iw, niw,
     &   itask, itrace, ind, ifail)
      if(ifail.ne.0.and.ifail.ne.10) write(*,*) 'D03PCF:ifail',ifail
      call xertst (ifail.eq.0.or.ifail.eq.10,
     & 'Abnormal exit from D03PCF')
      do i=1,ndepth_nml
        Tp(i) = min( U_nag(1,i), matsurf_vapor_temp(ibnd))
      enddo
#endif

      call subend ()
      return
      end subroutine heatdiff1D
c_______________________________________________________________________

c               SUBROUTINE FIC

c Subroutine to specify the initial condition: Temperature "uu" as a
c function of the spatial coordinate "xx".

c Both xx and uu are outputs from the subroutine.

c Thermal conductivity is assumed to vary as 1/(KA*T+KB) with
c temperature T (= uu in this routine).
c The above dependence follows the reference
      subroutine fic(ibnd,xx,uu)
      use b2mod_types
      use b2mod_wall
      use b2mod_layer
      use b2mod_heatdiff_common
     &    , only: Qtot, Tcbc, xn, x
      implicit none

      integer i, ibnd
      real (kind=R8) :: KAb, KBb, Tc, KAc, KBc, Tl, KAl, KBl
      real (kind=R8) :: xx(ndepth_nml), uu(ndepth_nml)


      KAb = matbulk_Kappa_KA(ibnd)
      KBb = matbulk_Kappa_KB(ibnd)
cxpb  Will need to update this later for a better computation
      if (coating_thickness(ibnd).gt.0.0_R8.and.
     &   (coating_thickness(ibnd)+layer_thickness(ibnd)).gt.0.0_R8)
     > then
        KAl = coating_Kappa_KA(ibnd)
        KBl = coating_Kappa_KB(ibnd)
      else
        KAl = matbulk_Kappa_KA(ibnd)
        KBl = matbulk_Kappa_KB(ibnd)
      endif
      if (coating_thickness(ibnd).gt.0.0_R8) then
        KAc = coating_Kappa_KA(ibnd)
        KBc = coating_Kappa_KB(ibnd)
        Tc = max(Tcbc, min( matsurf_vapor_temp(ibnd), (Tcbc + KBb/KAb) *
     .       exp(-KAb*Qtot*(coating_thickness(ibnd)-xn)) - KBb/KAb ))
      endif
      do i = ndepth_nml, 2, -1
        xx(i) = x(i)
        if (target_depth(ibnd,i).ge.coating_thickness(ibnd)) then
          uu(i) = max(Tcbc, min( matsurf_vapor_temp(ibnd),
     .     (Tcbc + KBb/KAb) * exp(-KAb*Qtot*(x(i)-xn)) - KBb/KAb ))
        else
          uu(i) = max(Tcbc, min( matsurf_vapor_temp(ibnd),
     .     (Tc + KBc/KAc) *
     .     exp(-KAc*Qtot*(x(i)-coating_thickness(ibnd))) - KBc/KAc ))
        endif
c        print*,x(i),' ',uu(i),' ',i,' ',xn,' ',xn
      enddo
      xx(1) = x(1)
      Tl = uu(2)
      uu(1) = max(Tl, min( matsurf_vapor_temp(ibnd),
     .     (Tl + Kbl/KAl) *
     .     exp(-KAl*Qtot*(x(1)-x(2))) - KBl/KAl ))

      return
      end subroutine fic

c_______________________________________________________________________
c Finite-volume form of the conduction through the plate (NO_NAG): the
c temperatures u(1:n) at the n = ndepth_nml depth points x(1:n), node i
c holding the volume between the midpoints to its neighbours; the
c fluxes of heatdiff_fluxes: the energy is conserved exactly. (The
c right-hand side for mstep; heatdiff1D integrates the same form with
c heatdiff_implicit.)
c_______________________________________________________________________
      subroutine lines2c(ibnd,ns,t,n,u,du)
      use b2mod_types
      use b2mod_wall
      use b2mod_heatdiff_common
     &    , only: x

      implicit none

      integer ns,n,i,ibnd
      real (kind=R8) :: u(n),du(n),t
      real (kind=R8) :: f(0:n), vol
      real (kind=R8) :: Cp, mass_density
      external Cp, mass_density, heatdiff_fluxes

      call heatdiff_fluxes(ibnd,ns,n,u,f)
      do i = 1, n
        if (i.eq.1) then
          vol = 0.5_R8*(x(2)-x(1))
        else if (i.eq.n) then
          vol = 0.5_R8*(x(n)-x(n-1))
        else
          vol = 0.5_R8*(x(i+1)-x(i-1))
        endif
        du(i) = (f(i-1)-f(i))/
     .   (mass_density(ibnd,i)*Cp(ibnd,i,u(i))*vol)
      enddo

      return
      end subroutine lines2c

c..................................................................
c Fluxes of the finite-volume form (W/m^2, positive towards the back):
c f(0) into the surface (the heat flux less radiation and evaporation),
c f(i) from node i to node i+1 (the harmonic mean of their
c conductivities), f(n) out through the back (heatdiff_contact).
      subroutine heatdiff_fluxes(ibnd,ns,n,u,f)
      use b2mod_types
      use b2mod_heatdiff_common
     &    , only: x, Qtot, Tcbc

      implicit none

      integer ibnd, ns, n, i
      real (kind=R8) :: u(n), f(0:n), k1, k2
      real (kind=R8) :: Kappa, qrad, qevap, heatdiff_contact
      external Kappa, qrad, qevap, heatdiff_contact

      f(0) = Qtot-qrad(ibnd,u(1),.true.)-qevap(ibnd,ns,u(1))
      do i = 1, n-1
        k1 = Kappa(ibnd,i,u(i))
        k2 = Kappa(ibnd,i+1,u(i+1))
        f(i) = 2.0_R8*k1*k2/((k1+k2)*(x(i+1)-x(i)))*(u(i)-u(i+1))
      enddo
      f(n) = heatdiff_contact(ibnd,u(n))*(u(n)-Tcbc)

      return
      end subroutine heatdiff_fluxes

c..................................................................
c Contact conductance (W m^-2 K^-1) between the back of element ibnd
c and its coolant or support at temperature T: BACKPLATE_CONDUCTANCE
c if given (zero: an insulated back), else cond_coef (active or
c inertial cooling, INERTIAL_COOLING).
      function heatdiff_contact(ibnd,T)
      use b2mod_types
      use b2mod_wall

      implicit none

      integer ibnd
      real (kind=R8) :: heatdiff_contact, T
      real (kind=R8) :: cond_coef
      external cond_coef

      if (backplate_conductance(ibnd).ge.0.0_R8) then
        heatdiff_contact = backplate_conductance(ibnd)
      else
        heatdiff_contact = cond_coef(inertial_cooling(ibnd),T)
      endif

      return
      end function heatdiff_contact

c..................................................................
c Steady state of the finite-volume form (heatdiff_fluxes): the depth
c temperatures at which the flux through every face equals the heat
c flux Q less radiation and evaporation at the surface and the contact
c flux at the back, by Newton iterations (tridiagonal Jacobian by
c differences, backtracking). From Tp if it holds a profile (Tp(1) > 0),
c else from a linear profile through the contact temperature. With an
c insulated back the plate settles where radiation and evaporation take
c Q away; if they cannot (and above the vaporisation temperature) the
c plate is at the vaporisation temperature (a message says so).
      subroutine heatdiff1D_steady(ibnd,ns0,Q,Tp,output)
      use b2mod_types
      use b2mod_wall
      use b2mod_heatdiff_common
     &    , only: dx, x0, xn, x, Qtot, Tcbc, active, iwall, ns
      use b2mod_subsys

      implicit none

      integer ibnd, ns0, output
      real (kind=R8) :: Q, Tp(ndepth_nml)

      integer itmax, nback
      parameter (itmax=200, nback=40)
      integer n, i, j, it, m, ib
      real (kind=R8) :: u(ndepth_nml), un(ndepth_nml), r(ndepth_nml),
     .  rp(ndepth_nml), du(ndepth_nml), sub(ndepth_nml),
     .  dia(ndepth_nml), sup(ndepth_nml), cc(ndepth_nml),
     .  f(0:ndepth_nml)
      real (kind=R8) :: hb, tb, tv, tlo, thi, tmid, k0, tol, rnorm,
     .  rnew, lam, step, del, denom
      logical ok
      real (kind=R8) :: Kappa, qrad, qevap, heatdiff_contact
      external Kappa, qrad, qevap, heatdiff_contact, heatdiff_fluxes

      call subini('heatdiff1D_steady')
      call xertst(ndepth_nml.ge.2,
     . 'heatdiff1D_steady needs NDEPTH_NML >= 2')
      n = ndepth_nml
      Qtot = Q
      Tcbc = backplate_temp(ibnd)
      active = .not.inertial_cooling(ibnd)
      iwall = ibnd
      ns = ns0
      x(1:n) = target_depth(ibnd,1:n)
      x0 = x(1)
      xn = x(n)
      do i=1,n-1
        dx(i) = x(i+1)-x(i)
      enddo
      tv = matsurf_vapor_temp(ibnd)

c Initial guess
      if (Tp(1).gt.0.0_R8) then
        u(1:n) = Tp(1:n)
      else
        hb = heatdiff_contact(ibnd,Tcbc)
        if (hb.gt.0.0_R8) then
          tb = Tcbc+max(Q,0.0_R8)/hb
          k0 = Kappa(ibnd,n,tb)
          do i=1,n
            u(i) = min(tb+max(Q,0.0_R8)*(xn-x(i))/k0, tv)
          enddo
        else
c         insulated: the temperature where radiation and evaporation
c         take Q away (bisection)
          tlo = Tcbc
          thi = tv
          do it=1,100
            tmid = 0.5_R8*(tlo+thi)
            if (qrad(ibnd,tmid,.true.)+qevap(ibnd,ns,tmid).lt.Q) then
              tlo = tmid
            else
              thi = tmid
            endif
          enddo
          u(1:n) = 0.5_R8*(tlo+thi)
        endif
      endif

c Newton iterations
      tol = 1.0e-10_R8*max(abs(Q),1.0_R8)
      call heatdiff_fluxes(ibnd,ns,n,u,f)
      do i=1,n
        r(i) = f(i-1)-f(i)
      enddo
      rnorm = maxval(abs(r(1:n)))
      ok = rnorm.le.tol
      it = 0
      do while (.not.ok .and. it.lt.itmax)
        it = it+1
c       tridiagonal Jacobian dr(i)/du(j), j = i-1, i, i+1, by
c       differences, three colours of nodes
        sub = 0.0_R8
        dia = 0.0_R8
        sup = 0.0_R8
        do m = 0, 2
          un(1:n) = u(1:n)
          do j = 1+m, n, 3
            un(j) = u(j)+1.0e-6_R8*max(abs(u(j)),1.0_R8)
          enddo
          call heatdiff_fluxes(ibnd,ns,n,un,f)
          do i=1,n
            rp(i) = f(i-1)-f(i)
          enddo
          do j = 1+m, n, 3
            del = un(j)-u(j)
            dia(j) = (rp(j)-r(j))/del
            if (j.gt.1) sup(j-1) = (rp(j-1)-r(j-1))/del
            if (j.lt.n) sub(j+1) = (rp(j+1)-r(j+1))/del
          enddo
        enddo
c       solve J du = -r (Thomas)
        cc(1) = sup(1)/dia(1)
        du(1) = -r(1)/dia(1)
        do i=2,n
          denom = dia(i)-sub(i)*cc(i-1)
          if (i.lt.n) cc(i) = sup(i)/denom
          du(i) = (-r(i)-sub(i)*du(i-1))/denom
        enddo
        do i=n-1,1,-1
          du(i) = du(i)-cc(i)*du(i+1)
        enddo
c       backtracking, steps of at most 500 K; a correction at round-off
c       level ends the iterations
        step = maxval(abs(du(1:n)))
        if (step.le.1.0e-11_R8*maxval(abs(u(1:n)))) then
          u(1:n) = u(1:n)+du(1:n)
          call heatdiff_fluxes(ibnd,ns,n,u,f)
          do i=1,n
            r(i) = f(i-1)-f(i)
          enddo
          rnorm = maxval(abs(r(1:n)))
          ok = .true.
          exit
        endif
        lam = min(1.0_R8, 500.0_R8/max(step,1.0e-300_R8))
        do ib = 1, nback
          un(1:n) = max(u(1:n)+lam*du(1:n), 1.0_R8)
          call heatdiff_fluxes(ibnd,ns,n,un,f)
          do i=1,n
            rp(i) = f(i-1)-f(i)
          enddo
          rnew = maxval(abs(rp(1:n)))
          if (rnew.lt.rnorm .or. rnew.le.tol) exit
          lam = 0.5_R8*lam
        enddo
        u(1:n) = un(1:n)
        r(1:n) = rp(1:n)
        rnorm = rnew
        ok = rnorm.le.tol .or. lam*step.le.1.0e-12_R8*maxval(u(1:n))
      enddo
      if (.not.ok) write(*,'(a,i5,a,1p,e12.4,a)')
     . ' heatdiff1D_steady: element ',ibnd,
     . ' not converged (flux residual ',rnorm,' W/m^2)'
      if (maxval(u(1:n)).gt.tv) write(*,'(a,i5,a)')
     . ' heatdiff1D_steady: element ',ibnd,
     . ' above the vaporisation temperature (capped)'
      do i=1,n
        Tp(i) = min(u(i), tv)
      enddo
      if (output.ge.10) write(*,*) 'heatdiff1D_steady: element ',ibnd,
     . ', iterations ',it,', surface ',Tp(1),' K, back ',Tp(n),' K'

      call subend()
      return
      end subroutine heatdiff1D_steady

c..................................................................
c Time integration of the finite-volume form (NO_NAG) over dt from the
c temperatures u: implicit (backward Euler) steps (heatdiff_be), the
c heat stored in enthalpy form. Each step h is taken whole and as two
c halves; their difference sets the next step (1e-2 K; the error at the
c end is about half of it), and the enthalpies are extrapolated,
c 2 H(two halves) - H(whole), which is second order and keeps the
c energy exact (the change of each node is a sum of face fluxes).
c Stable for any step: long discharges and cool-downs take few steps,
c and long steps tend to the steady state (heatdiff1D_steady).
      subroutine heatdiff_implicit(ibnd,ns,n,u,dt,output)
      use b2mod_types
      use b2mod_heatdiff_common
     &    , only: dx

      implicit none

      integer ibnd, ns, n, output
      real (kind=R8) :: u(n), dt

      real (kind=R8) :: tol
      parameter (tol=1.0e-2_R8)
      integer i, it, nstep, nrej
      real (kind=R8) :: uf(n), uh(n), u2(n), t, h, err, hexp, alpha,
     .  dh, ue, du
      logical ok1, ok2, ok3
      real (kind=R8) :: Kappa, Cp, mass_density, heatdiff_enthalpy
      external Kappa, Cp, mass_density, heatdiff_be, heatdiff_enthalpy

      if (dt.le.0.0_R8) return
c     first step: ten explicit steps of the finest cell
      hexp = dt
      do i=1,n-1
        alpha = Kappa(ibnd,i,u(i))/
     .   (mass_density(ibnd,i)*Cp(ibnd,i,u(i)))
        hexp = min(hexp, 0.25_R8*dx(i)**2/max(alpha,1.0e-30_R8))
      enddo
      h = min(dt, 10.0_R8*hexp)
      t = 0.0_R8
      nstep = 0
      nrej = 0
      do while (dt-t.gt.1.0e-12_R8*dt)
        h = min(h, dt-t)
        ok3 = .false.
        call heatdiff_be(ibnd,ns,n,u,h,uf,ok1)
        call heatdiff_be(ibnd,ns,n,u,0.5_R8*h,uh,ok2)
        if (ok2) call heatdiff_be(ibnd,ns,n,uh,0.5_R8*h,u2,ok3)
        if (.not.(ok1.and.ok2.and.ok3)) then
          h = 0.25_R8*h
          nrej = nrej+1
          call xertst(h.gt.1.0e-14_R8*dt .and. nrej.lt.100000,
     .     'heatdiff1D: no convergence')
          cycle
        endif
        err = maxval(abs(u2(1:n)-uf(1:n)))
        if (err.le.tol .or. h.le.1.0e-14_R8*dt) then
c         extrapolated enthalpy of each node, its temperature by Newton
          do i=1,n
            dh = 2.0_R8*heatdiff_enthalpy(ibnd,i,u(i),u2(i))
     .        -heatdiff_enthalpy(ibnd,i,u(i),uf(i))
            ue = 2.0_R8*u2(i)-uf(i)
            do it=1,30
              du = (dh-heatdiff_enthalpy(ibnd,i,u(i),ue))/
     .          Cp(ibnd,i,ue)
              ue = ue+du
              if (abs(du).le.1.0e-12_R8*abs(ue)) exit
            enddo
            if (ue.gt.0.0_R8 .and. abs(du).le.1.0e-9_R8*abs(ue)) then
              u(i) = ue
            else
              u(i) = u2(i)
            endif
          enddo
          t = t+h
          nstep = nstep+1
          h = h*min(4.0_R8, max(0.5_R8,
     .     0.9_R8*sqrt(tol/max(err,1.0e-300_R8))))
        else
          h = h*max(0.1_R8, 0.9_R8*sqrt(tol/err))
          nrej = nrej+1
        endif
      enddo
      if (output.ge.10) write(6,*) 'heatdiff1D: element ',ibnd,
     . ', steps ',nstep,', rejected ',nrej

      return
      end subroutine heatdiff_implicit

c..................................................................
c One implicit step h of the finite-volume form from u0: u solves
c rho V (H(u) - H(u0))/h = f(i-1) - f(i) at every node (H: the enthalpy,
c the integral of Cp; f: heatdiff_fluxes), by Newton iterations with a
c tridiagonal Jacobian (the fluxes by differences); ok false if they
c do not converge.
      subroutine heatdiff_be(ibnd,ns,n,u0,h,u,ok)
      use b2mod_types
      use b2mod_heatdiff_common
     &    , only: x

      implicit none

      integer ibnd, ns, n
      real (kind=R8) :: u0(n), h, u(n)
      logical ok

      integer i, j, m, it
      real (kind=R8) :: vol, rv(n), r(n), rf(n), rp(n), un(n), du(n),
     .  sub(n), dia(n), sup(n), cc(n), f(0:n), del, denom, step
      real (kind=R8) :: Cp, mass_density, heatdiff_enthalpy
      external Cp, mass_density, heatdiff_enthalpy, heatdiff_fluxes

      do i=1,n
        if (i.eq.1) then
          vol = 0.5_R8*(x(2)-x(1))
        else if (i.eq.n) then
          vol = 0.5_R8*(x(n)-x(n-1))
        else
          vol = 0.5_R8*(x(i+1)-x(i-1))
        endif
        rv(i) = mass_density(ibnd,i)*vol/h
      enddo
      u(1:n) = u0(1:n)
      ok = .false.
      do it = 1, 60
        call heatdiff_fluxes(ibnd,ns,n,u,f)
        do i=1,n
          rf(i) = -(f(i-1)-f(i))
          r(i) = rv(i)*heatdiff_enthalpy(ibnd,i,u0(i),u(i))+rf(i)
        enddo
c       Jacobian: the storage rho V Cp/h on the diagonal, the fluxes by
c       differences (three colours of nodes)
        sub = 0.0_R8
        sup = 0.0_R8
        do i=1,n
          dia(i) = rv(i)*Cp(ibnd,i,u(i))
        enddo
        do m = 0, 2
          un(1:n) = u(1:n)
          do j = 1+m, n, 3
            un(j) = u(j)+1.0e-6_R8*max(abs(u(j)),1.0_R8)
          enddo
          call heatdiff_fluxes(ibnd,ns,n,un,f)
          do i=1,n
            rp(i) = -(f(i-1)-f(i))
          enddo
          do j = 1+m, n, 3
            del = un(j)-u(j)
            dia(j) = dia(j)+(rp(j)-rf(j))/del
            if (j.gt.1) sup(j-1) = (rp(j-1)-rf(j-1))/del
            if (j.lt.n) sub(j+1) = (rp(j+1)-rf(j+1))/del
          enddo
        enddo
c       J du = -r (Thomas)
        cc(1) = sup(1)/dia(1)
        du(1) = -r(1)/dia(1)
        do i=2,n
          denom = dia(i)-sub(i)*cc(i-1)
          if (i.lt.n) cc(i) = sup(i)/denom
          du(i) = (-r(i)-sub(i)*du(i-1))/denom
        enddo
        do i=n-1,1,-1
          du(i) = du(i)-cc(i)*du(i+1)
        enddo
        step = maxval(abs(du(1:n)))
        if (.not.(step.lt.1.0e30_R8)) return
        if (step.gt.500.0_R8) du(1:n) = du(1:n)*500.0_R8/step
        u(1:n) = u(1:n)+du(1:n)
        if (minval(u(1:n)).le.0.0_R8) return
        if (step.le.1.0e-10_R8*max(maxval(abs(u(1:n))),1.0_R8)) then
          ok = .true.
          return
        endif
      enddo

      return
      end subroutine heatdiff_be

c..................................................................
c Enthalpy change (J/kg) of node i of element ibnd from temperature a to
c b: the integral of Cp, five-point Gauss-Legendre (exact for the cubic
c fits of Cp).
      function heatdiff_enthalpy(ibnd,i,a,b)
      use b2mod_types

      implicit none

      integer ibnd, i
      real (kind=R8) :: heatdiff_enthalpy, a, b

      integer k
      real (kind=R8) :: xg(5), wg(5), mid, half, sum
      real (kind=R8) :: Cp
      external Cp
      data xg /-0.906179845938664_R8, -0.538469310105683_R8, 0.0_R8,
     .  0.538469310105683_R8, 0.906179845938664_R8/
      data wg /0.236926885056189_R8, 0.478628670499366_R8,
     .  0.568888888888889_R8, 0.478628670499366_R8,
     .  0.236926885056189_R8/

      mid = 0.5_R8*(a+b)
      half = 0.5_R8*(b-a)
      sum = 0.0_R8
      do k=1,5
        sum = sum+wg(k)*Cp(ibnd,i,mid+half*xg(k))
      enddo
      heatdiff_enthalpy = half*sum

      return
      end function heatdiff_enthalpy

c.................................................................
c                      function Kappa
c  function Kappa(T) gives the coefficient of thermal conductivity
c  at temperature T.

      function Kappa(ibnd,i,T)
      use b2mod_types
      use b2mod_wall
      use b2mod_layer
      implicit none
      integer ibnd, i
      real (kind=R8) :: Kappa, T
      real (kind=R8) :: a(5), Tcels, Teff
      real (kind=R8) :: bulk_value, coating_value
      real (kind=R8) :: layer_Kappa
      external layer_Kappa, xerrab

c Units of Kappa are in "W/m/K"
      a(1) = matbulk_Kappa_A1(ibnd)
      a(2) = matbulk_Kappa_A2(ibnd)
      a(3) = matbulk_Kappa_A3(ibnd)
      a(4) = matbulk_Kappa_A4(ibnd)
      a(5) = matbulk_Kappa_A5(ibnd)
      Tcels = T - 273.15_R8
      if (bulk_temp_unit(ibnd).eq.'K') then
        Teff = T
      else if (bulk_temp_unit(ibnd).eq.'C') then
        Teff = Tcels
      else
        call xerrab ('Unrecognized temperature units (bulk Kappa) !')
      endif
      bulk_value = a(1) +
     . (a(2) + (a(3) + (a(4) +  a(5) * Teff) * Teff) * Teff) * Teff
      if (coating_thickness(ibnd).gt.0.0_R8) then
        a(1) = coating_Kappa_A1(ibnd)
        a(2) = coating_Kappa_A2(ibnd)
        a(3) = coating_Kappa_A3(ibnd)
        a(4) = coating_Kappa_A4(ibnd)
        a(5) = coating_Kappa_A5(ibnd)
        if (coating_temp_unit(ibnd).eq.'K') then
          Teff = T
        else if (coating_temp_unit(ibnd).eq.'C') then
          Teff = Tcels
        else
          call xerrab ('Unrecognized temperature units (coat Kappa) !')
        endif
        coating_value = a(1) +
     .   (a(2) + (a(3) + (a(4) +  a(5) * Teff) * Teff) * Teff) * Teff
      end if
      if (target_depth(ibnd,i).ge.coating_thickness(ibnd).and.
     .    target_depth(ibnd,i).ge.layer_depth(ibnd)) then
        Kappa = bulk_value
      else if (target_depth(ibnd,i).ge.layer_depth(ibnd)) then
        Kappa = coating_value
      else
        Kappa = layer_Kappa(ibnd,T)*
     .   (layer_depth(ibnd)-target_depth(ibnd,i))/
     .   (target_depth(ibnd,i+1)-target_depth(ibnd,i))
        if (layer_depth(ibnd).ge.coating_thickness(ibnd)) then
          Kappa = Kappa+bulk_value*
     .     (target_depth(ibnd,i+1)-layer_depth(ibnd))/
     .     (target_depth(ibnd,i+1)-target_depth(ibnd,i))
        else
          Kappa = Kappa+coating_value*
     .     (target_depth(ibnd,i+1)-layer_depth(ibnd))/
     .     (target_depth(ibnd,i+1)-target_depth(ibnd,i))
        end if
      end if

      return
      end function Kappa

c.................................................................
c                      function dKappabydT
c  function Kappa(T) gives the coefficient of thermal conductivity
c  at temperature T. We use the analytical expression resulting
c  from the first derivative of Kappa above.

      function dKappabydT(ibnd,i,T)
      use b2mod_types
      use b2mod_wall
      implicit none
      integer ibnd, i
      real (kind=R8) :: dKappabydT, T
      real (kind=R8) :: a(4), Tcels, Teff
      real (kind=R8) :: bulk_value, coating_value
      real (kind=R8) :: dlayer_KappabydT
      external dlayer_KappabydT, xerrab

      a(1) = matbulk_Kappa_A2(ibnd)
      a(2) = matbulk_Kappa_A3(ibnd)
      a(3) = matbulk_Kappa_A4(ibnd)
      a(4) = matbulk_Kappa_A5(ibnd)
      Tcels = T - 273.15_R8
      if (bulk_temp_unit(ibnd).eq.'K') then
        Teff = T
      else if (bulk_temp_unit(ibnd).eq.'C') then
        Teff = Tcels
      else
        call xerrab ('Unrecognized temperature units (bulk dKappadT) !')
      endif
      bulk_value = a(1) +
     . (2.0_R8*a(2) + (3.0_R8*a(3) + 4.0_R8*a(4) * Teff) * Teff) * Teff
      if (coating_thickness(ibnd).gt.0.0_R8) then
        a(1) = coating_Kappa_A2(ibnd)
        a(2) = coating_Kappa_A3(ibnd)
        a(3) = coating_Kappa_A4(ibnd)
        a(4) = coating_Kappa_A5(ibnd)
        if (coating_temp_unit(ibnd).eq.'K') then
          Teff = T
        else if (coating_temp_unit(ibnd).eq.'C') then
          Teff = Tcels
        else
          call xerrab
     .     ('Unrecognized temperature units (coat dKappadT) !')
        endif
        coating_value = a(1) +
     .   (2.0_R8*a(2) + (3.0_R8*a(3) +
     .                   4.0_R8*a(4) * Teff) * Teff) * Teff
      end if
      if (target_depth(ibnd,i).ge.coating_thickness(ibnd).and.
     .    target_depth(ibnd,i).ge.layer_depth(ibnd)) then
        dKappabydT = bulk_value
      else if (target_depth(ibnd,i).ge.layer_depth(ibnd)) then
        dKappabydT = coating_value
      else
        dKappabydT = dlayer_KappabydT(ibnd,T)*
     .   (layer_depth(ibnd)-target_depth(ibnd,i))/
     .   (target_depth(ibnd,i+1)-target_depth(ibnd,i))
        if (layer_depth(ibnd).ge.coating_thickness(ibnd)) then
          dKappabydT = dKappabydT+bulk_value*
     .     (target_depth(ibnd,i+1)-layer_depth(ibnd))/
     .     (target_depth(ibnd,i+1)-target_depth(ibnd,i))
        else
          dKappabydT = dKappabydT+coating_value*
     .     (target_depth(ibnd,i+1)-layer_depth(ibnd))/
     .     (target_depth(ibnd,i+1)-target_depth(ibnd,i))
        end if
      end if

      return
      end function dKappabydT

c.................................................................
c                      function Kappa_perp
c  function Kappa_perp(T) gives the coefficient of perpendicular
c  thermal conductivity at temperature T.

      function Kappa_perp(ibnd,i,T)
      use b2mod_types
      use b2mod_wall
      use b2mod_layer
      implicit none
      integer ibnd, i
      real (kind=R8) :: Kappa_perp, T
      real (kind=R8) :: a(5), Tcels, Teff
      real (kind=R8) :: bulk_value, coating_value
      real (kind=R8) :: layer_Kappa
      external layer_Kappa, xerrab

c Units are in "W/m/K"
      a(1) = matbulk_Kappa_perp_A1(ibnd)
      a(2) = matbulk_Kappa_perp_A2(ibnd)
      a(3) = matbulk_Kappa_perp_A3(ibnd)
      a(4) = matbulk_Kappa_perp_A4(ibnd)
      a(5) = matbulk_Kappa_perp_A5(ibnd)
      Tcels = T - 273.15_R8
      if (bulk_temp_unit(ibnd).eq.'K') then
        Teff = T
      else if (bulk_temp_unit(ibnd).eq.'C') then
        Teff = Tcels
      else
        call xerrab ('Unrecognized temperature units (bulk Kappa) !')
      endif
      bulk_value = a(1) +
     . (a(2) + (a(3) + (a(4) +  a(5) * Teff) * Teff) * Teff) * Teff
      if (coating_thickness(ibnd).gt.0.0_R8) then
        a(1) = coating_Kappa_perp_A1(ibnd)
        a(2) = coating_Kappa_perp_A2(ibnd)
        a(3) = coating_Kappa_perp_A3(ibnd)
        a(4) = coating_Kappa_perp_A4(ibnd)
        a(5) = coating_Kappa_perp_A5(ibnd)
        if (coating_temp_unit(ibnd).eq.'K') then
          Teff = T
        else if (coating_temp_unit(ibnd).eq.'C') then
          Teff = Tcels
        else
          call xerrab
     .     ('Unrecognized temperature units (coat Kappa) !')
        endif
        coating_value = a(1) +
     .   (a(2) + (a(3) + (a(4) +  a(5) * Teff) * Teff) * Teff) * Teff
      end if
      if (target_depth(ibnd,i).ge.coating_thickness(ibnd).and.
     .    target_depth(ibnd,i).ge.layer_depth(ibnd)) then
        Kappa_perp = bulk_value
      else if (target_depth(ibnd,i).ge.layer_depth(ibnd)) then
        Kappa_perp = coating_value
      else
cxpb  The layer is expected to be amorphous and therefore not to support
cxpb  the distinction between perp and parallel directions
        Kappa_perp = layer_Kappa(ibnd,T)*
     .   (layer_depth(ibnd)-target_depth(ibnd,i))/
     .   (target_depth(ibnd,i+1)-target_depth(ibnd,i))
        if (layer_depth(ibnd).ge.coating_thickness(ibnd)) then
          Kappa_perp = Kappa_perp+bulk_value*
     .     (target_depth(ibnd,i+1)-layer_depth(ibnd))/
     .     (target_depth(ibnd,i+1)-target_depth(ibnd,i))
        else
          Kappa_perp = Kappa_perp+coating_value*
     .     (target_depth(ibnd,i+1)-layer_depth(ibnd))/
     .     (target_depth(ibnd,i+1)-target_depth(ibnd,i))
        end if
      end if

      return
      end function Kappa_perp

c.................................................................
c                      function dKappa_perpbydT
c  function Kappa_perp(T) gives the coefficient of thermal conductivity
c  at temperature T in the perpendicular direction. We use the analytical
c  expression resulting from the first derivative of Kappa_perp above.

      function dKappa_perpbydT(ibnd,i,T)
      use b2mod_types
      use b2mod_wall
      use b2mod_layer
      implicit none
      integer ibnd, i
      real (kind=R8) :: dKappa_perpbydT, T
      real (kind=R8) :: a(4), Tcels, Teff
      real (kind=R8) :: bulk_value, coating_value
      real (kind=R8) :: dlayer_KappabydT
      external dlayer_KappabydT, xerrab

      a(1) = matbulk_Kappa_perp_A2(ibnd)
      a(2) = matbulk_Kappa_perp_A3(ibnd)
      a(3) = matbulk_Kappa_perp_A4(ibnd)
      a(4) = matbulk_Kappa_perp_A5(ibnd)
      Tcels = T - 273.15_R8
      if (bulk_temp_unit(ibnd).eq.'K') then
        Teff = T
      else if (bulk_temp_unit(ibnd).eq.'C') then
        Teff = Tcels
      else
        call xerrab ('Unrecognized temperature units (bulk dKappadT) !')
      endif
      bulk_value = a(1) +
     . (2.0_R8*a(2) + (3.0_R8*a(3) + 4.0_R8*a(4) * Teff) * Teff) * Teff
      if (coating_thickness(ibnd).gt.0.0_R8) then
        a(1) = coating_Kappa_perp_A2(ibnd)
        a(2) = coating_Kappa_perp_A3(ibnd)
        a(3) = coating_Kappa_perp_A4(ibnd)
        a(4) = coating_Kappa_perp_A5(ibnd)
        if (coating_temp_unit(ibnd).eq.'K') then
          Teff = T
        else if (coating_temp_unit(ibnd).eq.'C') then
          Teff = Tcels
        else
          call xerrab
     .     ('Unrecognized temperature units (coat dKappadT) !')
        endif
        coating_value = a(1) +
     .   (2.0_R8*a(2) + (3.0_R8*a(3) +
     .                   4.0_R8*a(4) * Teff) * Teff) * Teff
      end if
      if (target_depth(ibnd,i).ge.coating_thickness(ibnd).and.
     .    target_depth(ibnd,i).ge.layer_depth(ibnd)) then
        dKappa_perpbydT = bulk_value
      else if (target_depth(ibnd,i).ge.layer_depth(ibnd)) then
        dKappa_perpbydT = coating_value
      else
cxpb  The layer is expected to be amorphous and therefore not to support
cxpb  the distinction between perp and parallel directions
        dKappa_perpbydT = dlayer_Kappabydt(ibnd,T)*
     .   (layer_depth(ibnd)-target_depth(ibnd,i))/
     .   (target_depth(ibnd,i+1)-target_depth(ibnd,i))
        if (layer_depth(ibnd).ge.coating_thickness(ibnd)) then
          dKappa_perpbydT = dKappa_perpbydT+bulk_value*
     .     (target_depth(ibnd,i+1)-layer_depth(ibnd))/
     .     (target_depth(ibnd,i+1)-target_depth(ibnd,i))
        else
          dKappa_perpbydT = dKappa_perpbydT+coating_value*
     .     (target_depth(ibnd,i+1)-layer_depth(ibnd))/
     .     (target_depth(ibnd,i+1)-target_depth(ibnd,i))
        end if
      end if

      return
      end function dKappa_perpbydT

c....................................................................
c                      function Cp
c  function Cp(T) gives the specific heat capacity at const pressure
c  at temperature T.

      function Cp(ibnd,i,T)
      use b2mod_types
      use b2mod_wall
      use b2mod_layer
      implicit none
      integer ibnd, i
      real (kind=R8) :: Cp, T
      integer mf, cf
      real (kind=R8) :: b(4), Tcels, Teff
      real (kind=R8) :: Cp_bulk, Cp_coating
      real (kind=R8) :: layer_Cp
      external layer_Cp, xerrab

      Tcels = T - 273.15_R8
      b(1) = matbulk_Cp_B1(ibnd)
      b(2) = matbulk_Cp_B2(ibnd)
      b(3) = matbulk_Cp_B3(ibnd)
      b(4) = matbulk_Cp_B4(ibnd)
      mf = matbulk_format(ibnd)
      if (bulk_temp_unit(ibnd).eq.'K') then
        Teff = T
      else if (bulk_temp_unit(ibnd).eq.'C') then
        Teff = Tcels
      else
        call xerrab ('Unrecognized temperature units (bulk Cp) !')
      endif
      if (mf.eq.1) then
        Cp_bulk = b(1) + ( b(2) + ( b(3) + b(4) * Teff) * Teff ) * Teff
      else
        Cp_bulk = b(1) + b(2) * Teff - b(3)/Teff**2
      end if

      if (coating_thickness(ibnd).gt.0.0_R8) then
        b(1) = coating_Cp_B1(ibnd)
        b(2) = coating_Cp_B2(ibnd)
        b(3) = coating_Cp_B3(ibnd)
        b(4) = coating_Cp_B4(ibnd)
        cf = coating_format(ibnd)
        if (bulk_temp_unit(ibnd).eq.'K') then
          Teff = T
        else if (bulk_temp_unit(ibnd).eq.'C') then
          Teff = Tcels
        else
          call xerrab ('Unrecognized temperature units (coat Cp) !')
        endif
        if (cf.eq.1) then
          Cp_coating =
     .     b(1) + ( b(2) + ( b(3) + b(4) * Teff) * Teff ) * Teff
        else
          Cp_coating = b(1) + b(2) * Teff - b(3)/Teff**2
        end if
      end if
      if (target_depth(ibnd,i).ge.coating_thickness(ibnd).and.
     .    target_depth(ibnd,i).ge.layer_depth(ibnd)) then
        Cp = Cp_bulk
      else if (target_depth(ibnd,i).ge.layer_depth(ibnd)) then
        Cp = Cp_coating
      else
        Cp = (layer_Cp(ibnd,T))*
     .   (layer_depth(ibnd)-target_depth(ibnd,i))/
     .   (target_depth(ibnd,i+1)-target_depth(ibnd,i))
        if (layer_depth(ibnd).ge.coating_thickness(ibnd)) then
          Cp = Cp+Cp_bulk*
     .     (target_depth(ibnd,i+1)-layer_depth(ibnd))/
     .     (target_depth(ibnd,i+1)-target_depth(ibnd,i))
        else
          Cp = Cp+Cp_coating*
     .     (target_depth(ibnd,i+1)-layer_depth(ibnd))/
     .     (target_depth(ibnd,i+1)-target_depth(ibnd,i))
        end if
      end if

      return
      end function Cp

c Getting the heat loss due to evaporation and radiation from surface.
c Radiation: grey body, the Stefan-Boltzmann constant
c 2 pi^5 k^4/(15 h^3 c^2).
c....................................................................
      function qrad(ibnd,Ts,top_layer)
      use b2mod_types
      use b2mod_constants
      use b2mod_wall
      implicit none
      integer ibnd
#     include "emissivity.constants"
      real (kind=R8) :: Ts, Tw, SBconst, qrad
      logical top_layer
      parameter (SBconst =
     . pi**5 * (2.0_R8/15.0_R8) * (KBolt/const_h)**3
     . * KBolt / c**2)
      real (kind=R8) :: layer_emissivity
      external layer_emissivity

      Tw = backplate_temp(ibnd)      ! Temperature of surrounding environment
      if (layer_depth(ibnd).eq.0.0_R8 .or. .not.top_layer) then
        qrad = matsurf_blackbody_emiss(ibnd)*SBconst*(Ts**4-Tw**4)
      else
        qrad = layer_emissivity(ibnd)*SBconst*(Ts**4-Tw**4)
      endif

      return
      end function qrad

c....................................................................
c   Ref: 1) C. Garcia Rosales, Erosion processes in plasma-wall
c           interactions, J. Nucl. Mater. 211 (1994) Pg.202-214.
c        2) R.A. Langley, Evaporation, Data Compendium for Plasma
c           Surface Interactions, Special Issue, NUCLEAR FUSION,
c           IAEA, Vienna, (1984) 55-61.
c  Table of A and B for some first wall materials in Fusion research:
c  (Ref (2) above).
c  (note that A and B are provided such that the vappres is in Pascals).
c  element    |  A       |  B      |
c  Aluminum   | -16165   | 10.91   |
c  Beryllium  | -16720   | 11.61   |
c  Graphite   | -40181   | 14.80   |
c  Copper     | -17079   | 11.26   |
c  Molybdenum | -33069   | 12.02   |
c  Titanium   | -23340   | 9.79    |
c  Tungsten   | -44485   | 12.74   |
c  Titanium - |          |         |
c     Carbide | -39117.9 | 15.23   |
c  Stainless -|          |         |
c     Steel   | -14321.9 | 9.44    |

      function qevap(ibnd,ns,Ts)
      use b2mod_types
      use b2mod_wall
      use b2mod_b2cmpa
      use b2mod_sputter
      use b2mod_subsys
      implicit none
      integer ibnd, ns, is, itrack
      logical, save :: lfirst
      real (kind=R8) :: Ts, Tv, flxthev, qevap, Vapor_temp, qthev
      real (kind=R8), save :: therm_evap
      data lfirst /.true./, therm_evap /0.0_R8/
      external ipgetr
      intrinsic min

c  Ts       ->  temperature of the plate in Kelvin.
c  flxthev  ->  Flux of target atoms evaporated ( / m^2 /s )
c  qevap -> Energy lost (Joules / m^2 / s)

      call subini ('qevap')
      if (lfirst) then
#ifdef _OPENMP
!$OMP CRITICAL (qevap_heatdiff1D)
#endif
        call ipgetr('b2stbr_therm_evap', therm_evap)
        call xertst(0.0_R8.le.therm_evap.and.
     &              1.0_R8.ge.therm_evap,
     .    'faulty parameter therm_evap')
        lfirst = .false.
#ifdef _OPENMP
!$OMP END CRITICAL (qevap_heatdiff1D)
#endif
      endif

      qevap = 0.0_R8
      if (therm_evap.eq.0.0_R8) then
        call subend ()
        return
      endif
      do is = 0, ns-1
        itrack = track_index(is)
        if (is_neutral(is)) then
c  Checking inputs
          if (layer_depth(ibnd).eq.0.0_R8 .or. itrack.eq.0) then
            Vapor_temp = matsurf_vapor_temp(ibnd)
            Tv = min(Vapor_temp, Ts)
          else if (itrack.gt.0 .and. itrack.le.ntrack) then
            Vapor_temp = species_vapor_temp(itrack)
            Tv = min(Vapor_temp, Ts)
          else
            Tv = Ts
          endif
          if (Ts.gt.Vapor_temp) write(*,'(a,1p,2e12.5)')
     >     'Input error; Correct so that Ts <= Vapor_temp',
     .                                   Ts,   Vapor_temp
          call thermev(ibnd,is,Tv,flxthev,qthev)
          qevap = qevap + qthev * therm_evap
        endif
      enddo

      call subend ()
      return
      end function qevap

      subroutine pdedef(npde, T, X, U, UX, P, Q, R, ires)
      use b2mod_types
      use b2mod_wall
      use b2mod_heatdiff_common
     &    , only: iwall
      implicit none
      integer npde, ires
      real (kind=R8) :: T, X, U(NPDE), UX(NPDE),
     .  P(NPDE,NPDE), Q(NPDE), R(NPDE)
      integer i, ib, it, ih

      real (kind=R8) :: Cp, Kappa, mass_density
      external Cp, Kappa, mass_density

      if (X.le.target_depth(iwall,1)) then
        ih = 1
      else if (X.ge.target_depth(iwall,ndepth_nml)) then
        ih = ndepth_nml
      else
        ib = 1
        it = ndepth_nml
        do while (it-ib.gt.1)
          ih = (ib+it)/2
          if (X.lt.target_depth(iwall,ih)) then
            it = ih
          else if (X.gt.target_depth(iwall,ih)) then
            ib = ih
          else
            ib = ih
            it = ih
          endif
        end do
      end if

      do i = 1, npde
        P(i,1:npde) = mass_density(iwall, ih) * Cp(iwall, ih, U(i))
        Q(i) = 0.0_R8
        R(i) = Kappa(iwall, ih, U(i)) * UX(i)
      enddo

      return
      end subroutine pdedef

      subroutine bndary (npde, T, U, UX, iposb, beta, gamma, ires)
      use b2mod_types
      use b2mod_heatdiff_common
     &    , only: Qtot, Tcbc, iwall, ns
      implicit none
      integer npde, iposb, ires
      real (kind=R8) :: T, U(npde), UX(npde), beta(npde), gamma(npde)
      real (kind=R8) :: heatdiff_contact
      external heatdiff_contact

      real (kind=R8) :: qrad, qevap
      external qrad, qevap

      if (iposb.eq.0) then  ! LHS BC : x=0 : Heat flux imposed
        beta(1)=-1.0_R8
        gamma(1)=Qtot-qrad(iwall,U(1),.true.)-qevap(iwall,ns,U(1))
      else if (iposb.eq.1) then ! RHS BC : Cooling with some contact heat conductivity (low for inertial cooling, high for active cooling)
        beta(1)= 1.0_R8
        gamma(1)=-heatdiff_contact(iwall,U(1))*(U(1)-Tcbc)
!! FIXME: Could use additionally the code below once we have a way of applying
!!        the correct emissivity
cxpb        if (.not.active)
cxpb     >    gamma(1)=gamma(1)-qrad(iwall,U(1),.false.)
cxpb     .                     -qevap(iwall,ns,U(1))
c$$$ Old treatment where conductivity to coolant was assumed to be perfect
c$$$      else if (iposb.eq.1.and.active) then ! RHS BC : x=xn : Cooling temperature imposed
c$$$        beta(1)= 0.0_R8
c$$$        gamma(1)= U(1)-Tcbc
c$$$      else if (iposb.eq.1.and..not.active) then ! RHS BC : x=xn : Inertial Cooling
c$$$        beta(1)= 1.0_R8
c$$$        gamma(1)= -100.0_R8*(U(1)-Tcbc) ! we assume the contact heat conductivity of 100 W.m-2.K-1
      end if

      return
      end subroutine bndary
