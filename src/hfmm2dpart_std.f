cc Copyright (C) 2010-2011: Leslie Greengard and Zydrunas Gimbutas
cc SPDX-License-Identifier: GPL-2.0-or-later
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c     This file contains the main FMM routines and some related
c     subroutines for evaluating Helmholtz potentials and gradients due to
c     point charges and dipoles.  (FORTRAN 90 VERSION)
c
c     nd-vectorized version with tree-API support.
c
c hfmm2d: charge and dipstr are complex valued, x \in R^2
c
c \phi(x_i) = (ima/4) \sum_{j\ne i}   charge_j H^{(1)}_0(k |x_i-x_j| )
c             + dipstr_j (dipvec_j \dot \grad_j H^{(1)}_0(k |x_i-x_j| )
c
c or, more precisely,
c
c \phi(x_i) = (ima/4) \sum_{j\ne i}   charge_j H^{(1)}_0(k |x_i-x_j| )
c + dipstr_j (dipvec_j \dot (x_i-x_j) ) H^{(1)}_1(k |x_i-x_j|) *k/|x_i-x_j|
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     hfmm2dpart_std_nd - Helmholtz FMM in R^2: evaluate all pairwise particle
c         interactions (ignoring self-interaction)
c
c     hfmm2dpartself_std_nd - Helmholtz FMM in R^2: evaluate all pairwise particle
c         interactions (ignoring self-interaction)
c
c     hfmm2dparttarg_std_nd - Helmholtz FMM in R^2: evaluate all pairwise
c         particle interactions (ignoring self-interaction) +
c         interactions with targets
c
c     h2dpartdirect_std_nd - Helmholtz interactions in R^2:  evaluate all
c         pairwise particle interactions (ignoring self-interaction) +
c         interactions with targets via direct O(N^2) algorithm
c
c     hfmm2dpart_tree_build_std - Build FMM tree (geometry-dependent)
c     hfmm2dpart_tree_eval_std  - Evaluate FMM using pre-built tree
c     hfmm2dpart_tree_free_std  - Free FMM tree
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        this is the end of the debugging code and the beginning
c        of the Helmholtz particle FMM in R^2
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine hfmm2dpart_std_nd(ier,iprec,zk,nd,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
        implicit real *8 (a-h,o-z)
c
c       Helmholtz FMM in R^2: evaluate all pairwise particle
c       interactions (ignoring self-interaction).
c       We use H_0(kr)*(i/4) for the Green's function.
c       Self-interactions are not included.
c
c hfmm2d: charge and dipstr are complex valued, x \in R^2
c
c \phi(x_i) = (ima/4) \sum_{j\ne i}   charge_j H^{(1)}_0(k |x_i-x_j| )
c             + dipstr_j (dipvec_j \dot \grad_j H^{(1)}_0(k |x_i-x_j| )
c
c or, more precisely,
c
c \phi(x_i) = (ima/4) \sum_{j\ne i}   charge_j H^{(1)}_0(k |x_i-x_j| )
c + dipstr_j (dipvec_j \dot (x_i-x_j) ) H^{(1)}_1(k |x_i-x_j|) *k/|x_i-x_j|
c
c       The main FMM routine permits both evaluation at sources
c       and at a collection of targets.
c       This subroutine is used to simplify the user interface
c       (by setting the number of targets to zero) and calling the more
c       general FMM.
c
c       See below for explanation of calling sequence arguments.
c
        real *8 source(2,*)
        complex *16 charge(nd,*),zk
        complex *16 dipstr(nd,*)
        real *8 dipvec(2,*)
        complex *16 ima
        complex *16 pot(nd,*)
        complex *16 grad(nd,2,*)
        complex *16 hess(nd,3,*)
c
        real *8 target(2,1)
        complex *16 pottarg(nd,1)
        complex *16 gradtarg(nd,2,1)
        complex *16 hesstarg(nd,3,1)
c
        data ima/(0.0d0,1.0d0)/
c
        ntarget=0
        ifpottarg=0
        ifgradtarg=0
        ifhesstarg=0
c
        call hfmm2dparttarg_std_nd(ier,iprec,zk,nd,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     ntarget,target,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg)
c
        return
        end
c
c
c
c
c
        subroutine hfmm2dpartself_std_nd(ier,iprec,zk,nd,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
        implicit real *8 (a-h,o-z)
c
c       Helmholtz FMM in R^2: evaluate all pairwise particle
c       interactions (ignoring self-interaction).
c       We use H_0(kr)*(i/4) for the Green's function.
c       Self-interactions are not included.
c
c hfmm2d: charge and dipstr are complex valued, x \in R^2
c
c \phi(x_i) = (ima/4) \sum_{j\ne i}   charge_j H^{(1)}_0(k |x_i-x_j| )
c             + dipstr_j (dipvec_j \dot \grad_j H^{(1)}_0(k |x_i-x_j| )
c
c or, more precisely,
c
c \phi(x_i) = (ima/4) \sum_{j\ne i}   charge_j H^{(1)}_0(k |x_i-x_j| )
c + dipstr_j (dipvec_j \dot (x_i-x_j) ) H^{(1)}_1(k |x_i-x_j|) *k/|x_i-x_j|
c
c       The main FMM routine permits both evaluation at sources
c       and at a collection of targets.
c       This subroutine is used to simplify the user interface
c       (by setting the number of targets to zero) and calling the more
c       general FMM.
c
c       See below for explanation of calling sequence arguments.
c
        real *8 source(2,*)
        complex *16 charge(nd,*),zk
        complex *16 dipstr(nd,*)
        real *8 dipvec(2,*)
        complex *16 ima
        complex *16 pot(nd,*)
        complex *16 grad(nd,2,*)
        complex *16 hess(nd,3,*)
c
        real *8 target(2,1)
        complex *16 pottarg(nd,1)
        complex *16 gradtarg(nd,2,1)
        complex *16 hesstarg(nd,3,1)
c
        data ima/(0.0d0,1.0d0)/
c
        ntarget=0
        ifpottarg=0
        ifgradtarg=0
        ifhesstarg=0
c
        call hfmm2dparttarg_std_nd(ier,iprec,zk,nd,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     ntarget,target,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg)
c
        return
        end
c
c
c
c
c
        subroutine hfmm2dparttarg_std_nd(ier,iprec,zk,nd,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     ntarget,target,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg)
        implicit real *8 (a-h,o-z)
c
c       Helmholtz FMM in R^2: evaluate all pairwise particle
c       interactions (ignoring self-interaction)
c       and interactions with targets.
c
c       We use H_0(kr)*(i/4) for the Green's function.
c       Self-interactions are not included.
c
c hfmm2d: charge and dipstr are complex valued, x \in R^2
c
c \phi(x_i) = (ima/4) \sum_{j\ne i}   charge_j H^{(1)}_0(k |x_i-x_j| )
c             + dipstr_j (dipvec_j \dot \grad_j H^{(1)}_0(k |x_i-x_j| )
c
c or, more precisely,
c
c \phi(x_i) = (ima/4) \sum_{j\ne i}   charge_j H^{(1)}_0(k |x_i-x_j| )
c + dipstr_j (dipvec_j \dot (x_i-x_j) ) H^{(1)}_1(k |x_i-x_j|) *k/|x_i-x_j|
c
c       This is primarily a memory management code.
c       The actual work is carried out in subroutine hfmm2dparttargmain_std.
c
c       INPUT PARAMETERS:
c
c       iprec:  FMM precision flag
c
c                 -2 => tolerance =.5d0
c                 -1 => tolerance =.5d-1
c                  0 => tolerance =.5d-2
c                  1 => tolerance =.5d-3
c                  2 => tolerance =.5d-6
c                  3 => tolerance =.5d-9
c                  4 => tolerance =.5d-12
c                  5 => tolerance =.5d-15
c
c       zk: complex *16: Helmholtz parameter
c       nd: integer: number of densities
c       nsource: integer:  number of sources
c       source: real *8 (2,nsource):  source locations
c       ifcharge:  charge computation flag
c                  ifcharge = 1   =>  include charge contribution
c                                     otherwise do not
c       charge: complex *16 (nd,nsource): charge strengths
c       ifdipole:  dipole computation flag
c                  ifdipole = 1   =>  include dipole contribution
c                                     otherwise do not
c       dipstr: complex *16 (nd,nsource): dipole strengths
c       dipvec: real *8 (2,nsource): dipole orientation vectors.
c
c       ifpot:  potential flag (1=compute potential, otherwise no)
c       ifgrad:  gradient flag (1=compute gradient, otherwise no)
c       ifhess:  hessian flag (1=compute hessian, otherwise no)
c       ntarget: integer:  number of targets
c       target: real *8 (2,ntarget):  target locations
c       ifpottarg:  target potential flag
c                   (1=compute potential, otherwise no)
c       ifgradtarg:  target gradient flag
c                   (1=compute gradient, otherwise no)
c       ifhesstarg:  target hessian flag
c                   (1=compute hessian, otherwise no)
c
c       OUTPUT PARAMETERS:
c
c       ier   =  error return code
c                ier=0     =>  normal execution
c                ier=4     =>  cannot allocate tree workspace
c                ier=8     =>  cannot allocate bulk FMM  workspace
c                ier=16    =>  cannot allocate mpole expansion
c                              workspace in FMM
c
c       pot: complex *16 (nd,nsource): potential at source locations
c       grad: complex *16 (nd,2,nsource): gradient at source locations
c       hess: complex *16 (nd,3,nsource): hessian at source locations
c       pottarg: complex *16 (nd,ntarget): potential at target locations
c       gradtarg: complex *16 (nd,2,ntarget): gradient at target locations
c       hesstarg: complex *16 (nd,3,ntarget): hessian at target locations
c
c
c
        real *8 source(2,*)
        complex *16 charge(nd,*),zk
        complex *16 dipstr(nd,*)
        real *8 dipvec(2,*)
        complex *16 ima
        complex *16 pot(nd,*)
        complex *16 grad(nd,2,*)
        complex *16 hess(nd,3,*)
c
        real *8 target(2,*)
        complex *16 pottarg(nd,*)
        complex *16 gradtarg(nd,2,*)
        complex *16 hesstarg(nd,3,*)
c
        real *8 timeinfo(10)
c
        real *8 center(2)
c
c
c     Note: various arrays dimensioned here to 200.
c     That allows for 200 levels of refinement, which is
c     more than enough for any non-pathological case.
c
        integer laddr(2,200)
        real *8 scale(0:200)
        real *8 bsize(0:200)
        integer nterms(0:200)
c
        complex *16 ptemp,gtemp(2),htemp(3)
c
        integer box(15)
        real *8 center0(2),corners0(2,4)
c
        integer box1(15)
        real *8 center1(2),corners1(2,4)
c
        real *8, allocatable :: w(:)
        real *8, allocatable :: wlists(:)
        complex *16, allocatable :: wrmlexp(:)
c
        data ima/(0.0d0,1.0d0)/
c
        ier=0
        lused7=0
c
        done=1
        pi=4*atan(done)
c
c
c     ifprint is an internal information printing flag.
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c
        ifprint=1
c
c     set fmm tolerance based on iprec flag.
c
        if( iprec .eq. -2 ) epsfmm=.5d-0
        if( iprec .eq. -1 ) epsfmm=.5d-1
        if( iprec .eq. 0 ) epsfmm=.5d-2
        if( iprec .eq. 1 ) epsfmm=.5d-3
        if( iprec .eq. 2 ) epsfmm=.5d-6
        if( iprec .eq. 3 ) epsfmm=.5d-9
        if( iprec .eq. 4 ) epsfmm=.5d-12
        if( iprec .eq. 5 ) epsfmm=.5d-15
        if( iprec .eq. 6 ) epsfmm=0
c
        if(ifprint .eq. 1) call prin2('epsfmm=*',epsfmm,1)
c
c     set criterion for box subdivision (number of sources per box)
c
        if( iprec .eq. -2 ) nbox=3
        if( iprec .eq. -1 ) nbox=5
        if( iprec .eq. 0 ) nbox=8
        if( iprec .eq. 1 ) nbox=10
        if( iprec .eq. 2 ) nbox=15
        if( iprec .eq. 3 ) nbox=20
        if( iprec .eq. 4 ) nbox=25
        if( iprec .eq. 5 ) nbox=45
        if( iprec .eq. 6 ) nbox=nsource+ntarget
c
c
        if(ifprint .eq. 1) call prinf('nbox=*',nbox,1)
c
c
c     create quad-tree data structure
c
        t1=second()
C$        t1=omp_get_wtime()
        ntot = 100*(nsource+ntarget)+10000
        do ii = 1,10
           allocate (wlists(ntot))
           call hfmm2dparttree(ier,iprec,zk,
     $        nsource,source,ntarget,target,
     $        nbox,epsfmm,iisource,iitarget,iwlists,lwlists,
     $        nboxes,laddr,nlev,center,size,
     $        wlists,ntot,lused7)
           if (ier.ne.0) then
              deallocate(wlists)
              ntot = ntot*1.5
              call prinf(' increasing allocation, ntot is *',ntot,1)
           else
             goto 1200
           endif
        enddo
1200    continue
        if (ier.ne.0) then
           call prinf(' exceeded max allocation, ntot is *',ntot,1)
           ier = 4
           return
        endif
        t2=second()
C$        t2=omp_get_wtime()
        if( ifprint .eq. 1 ) call prin2('time in d2tstcr=*',t2-t1,1)
c
c     lused7 is counter that steps through workspace,
c     keeping track of total memory used.
c
        lused7=1
c
c
c       ... prepare data structures
c
        do i = 0,nlev
        scale(i) = 1.0d0
        boxsize = abs((size/2.0**i)*zk)
        if (boxsize .lt. 1) scale(i) = boxsize
        enddo
c
        if(ifprint .eq. 1) call prin2('scale=*',scale,nlev+1)
c
c
c       carve up workspace further
c
c     isourcesort is pointer for sorted source coordinates
c     itargetsort is pointer for sorted target locations
c     ichargesort is pointer for sorted charge densities
c     idipvecsort is pointer for sorted dipole orientation vectors
c     idipstrsort is pointer for sorted dipole densities
c
        isourcesort = lused7 + 5
        lsourcesort = 2*nsource
        itargetsort = isourcesort+lsourcesort
        ltargetsort = 2*ntarget
        ichargesort = itargetsort+ltargetsort
        lchargesort = 2*nsource*nd
        idipvecsort = ichargesort+lchargesort
        if (ifdipole.eq.1) then
          ldipvec = 2*nsource
          ldipstr = 2*nsource*nd
        else
          ldipvec = 2
          ldipstr = 2
        endif
        idipstrsort = idipvecsort + ldipvec
        lused7 = idipstrsort + ldipstr
c
c       ... allocate the potential and gradient arrays
c
        ipot = lused7
        lpot = 2*nsource*nd
        lused7=lused7+lpot
c
        igrad = lused7
        if( ifgrad .eq. 1) then
        lgrad = 2*2*nsource*nd
        else
        lgrad=4
        endif
        lused7=lused7+lgrad
c
        ihess = lused7
        if( ifhess .eq. 1) then
        lhess = 2*3*nsource*nd
        else
        lhess=6
        endif
        lused7=lused7+lhess
c
        ipottarg = lused7
        lpottarg = 2*ntarget*nd
        lused7=lused7+lpottarg
c
        igradtarg = lused7
        if( ifgradtarg .eq. 1) then
        lgradtarg = 2*2*ntarget*nd
        else
        lgradtarg=4
        endif
        lused7=lused7+lgradtarg
c
        ihesstarg = lused7
        if( ifhesstarg .eq. 1) then
        lhesstarg = 2*3*ntarget*nd
        else
        lhesstarg=6
        endif
        lused7=lused7+lhesstarg
c
        if(ifprint .eq. 1) call prinf(' lused7 is *',lused7,1)
c
c       based on FMM tolerance, compute expansion lengths nterms(i)
c
        nmax = 0

        do i = 0,nlev
           bsize(i)=size/2.0d0**i
           call h2dterms(bsize(i),zk,epsfmm, nterms(i), ier)
           if (nterms(i).gt. nmax .and. i.ge. 2) nmax = nterms(i)
        enddo
c
        if (ifprint.eq.1)
     $     call prin2('in hfmm2dpart_std_nd, bsize(0) zk/2 pi=*',
     $     abs(bsize(0)*zk)/2/pi,1)
c
        if (ifprint.eq.1) call prin2('zk=*',zk,2)
        if (ifprint.eq.1) call prin2('bsize=*',bsize,nlev+1)
        if (ifprint.eq.1) call prinf('nterms=*',nterms,nlev+1)
c
c
c     Multipole and local expansions will be held in workspace
c     in locations pointed to by array iaddr(2,nboxes).
c
c     iiaddr is pointer to iaddr array, itself contained in workspace.
c     imptemp is pointer for single expansion (dimensioned by nmax)
c
c       ... allocate iaddr and temporary arrays
c
        iiaddr = lused7
        imptemp = iiaddr + 2*nboxes
        lmptemp = (2*nmax+1)*2*nd
        lused7 = imptemp + lmptemp
        allocate(w(lused7),stat=ier)
        if (ier.ne.0) then
           call prinf(' cannot allocate bulk FMM workspace,
     1                  lused7 is *',lused7,1)
           ier = 8
           return
        endif
c
c     reorder sources, targets so that each box holds
c     contiguous list of source/target numbers.

c
        call h2dreorder(nd,nsource,source,ifcharge,charge,
     $     wlists(iisource),
     $     ifdipole,dipstr,dipvec,
     1     w(isourcesort),w(ichargesort),w(idipvecsort),
     $     w(idipstrsort))
c
        call h2dreordertarg(ntarget,target,wlists(iitarget),
     $     w(itargetsort))
c
        if(ifprint .eq. 1) then
        call prinf('finished reordering=*',ier,1)
        call prinf('ier=*',ier,1)
        call prinf('nboxes=*',nboxes,1)
        call prinf('nlev=*',nlev,1)
        call prinf('nboxes=*',nboxes,1)
        call prinf('lused7=*',lused7,1)
        endif
c
c
c     allocate memory need by multipole, local expansions at all
c     levels
c     irmlexp is pointer for workspace need by various fmm routines,
c
        call h2dmpalloc(nd,wlists(iwlists),w(iiaddr),nboxes,lmptot,
     $     nterms)
c
        if(ifprint .eq. 1) call prinf(' lmptot is *',lmptot,1)
c
        irmlexp = 1
        lused7 = irmlexp + lmptot
        if(ifprint .eq. 1) call prinf(' lused7 is *',lused7,1)
        allocate(wrmlexp(lused7),stat=ier)
        if (ier.ne.0) then
           call prinf(' cannot allocate mpole expansion workspace,
     1                  lused7 is *',lused7,1)
           ier = 16
           return
        endif
c
c
ccc        do i=lused7+1,lused7+1+100
ccc        w(i)=777
ccc        enddo
c
c     Memory allocation is complete.
c     Call main fmm routine. There are, unfortunately, a lot
c     of parameters here. ifevalfar and ifevalloc determine
c     whether far gradient and local gradients (respectively) are to
c     be evaluated. Setting both to 1 means that both will be
c     computed (which is the normal scenario).
c
        ifevalfar=1
        ifevalloc=1
c
        t1=second()
C$        t1=omp_get_wtime()
        call hfmm2dparttargmain_std(ier,iprec,zk,nd,
     $     ifevalfar,ifevalloc,
     $     nsource,w(isourcesort),wlists(iisource),
     $     ifcharge,w(ichargesort),
     $     ifdipole,w(idipstrsort),w(idipvecsort),
     $     ifpot,w(ipot),ifgrad,w(igrad),ifhess,w(ihess),
     $     ntarget,w(itargetsort),wlists(iitarget),
     $     ifpottarg,w(ipottarg),ifgradtarg,w(igradtarg),
     $     ifhesstarg,w(ihesstarg),
     $     epsfmm,w(iiaddr),wrmlexp(irmlexp),w(imptemp),lmptemp,
     $     nboxes,laddr,nlev,scale,bsize,nterms,
     $     wlists(iwlists),lwlists)
        t2=second()
C$        t2=omp_get_wtime()
        if( ifprint .eq. 1 ) call prin2('time in fmm main=*',t2-t1,1)
c
c       parameter ier from targmain routine is currently meaningless, reset to 0
        if( ier .ne. 0 ) ier = 0
c
        if(ifprint .eq. 1) call prinf('lwlists=*',lwlists,1)
        if(ifprint .eq. 1) then
        call prinf('lused total=*',lused7,1)
        call prinf('lused total(k)=*',lused7/1000,1)
        call prinf('lused total(M)=*',lused7/1000000,1)
        endif
c
        if(ifprint .eq. 1)
     $     call prin2('memory / point = *',(lused7)/dble(nsource),1)
c
ccc        call prin2('after w=*', w(1+lused7-100), 2*100)
c
        if(ifpot .eq. 1)
     $     call h2dpsort(nd,nsource,wlists(iisource),w(ipot),pot)
        if(ifgrad .eq. 1)
     $     call h2dgsort(nd,nsource,wlists(iisource),w(igrad),grad)
        if(ifhess .eq. 1)
     $     call h2dhsort(nd,nsource,wlists(iisource),w(ihess),hess)
c
        if(ifpottarg .eq. 1 )
     $     call h2dpsort(nd,ntarget,wlists(iitarget),w(ipottarg),
     $     pottarg)
        if(ifgradtarg .eq. 1)
     $     call h2dgsort(nd,ntarget,wlists(iitarget),w(igradtarg),
     $     gradtarg)
        if(ifhesstarg .eq. 1)
     $     call h2dhsort(nd,ntarget,wlists(iitarget),w(ihesstarg),
     $     hesstarg)
c
        return
        end
c
c
c
c
c
        subroutine hfmm2dparttargmain_std(ier,iprec,zk,nd,
     $     ifevalfar,ifevalloc,
     $     nsource,sourcesort,isource,
     $     ifcharge,chargesort,
     $     ifdipole,dipstrsort,dipvecsort,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,ntarget,
     $     targetsort,itarget,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg,
     $     epsfmm,iaddr,rmlexp,mptemp,lmptemp,
     $     nboxes,laddr,nlev,scale,bsize,nterms,
     $     wlists,lwlists)
        implicit real *8 (a-h,o-z)
        real *8 sourcesort(2,*)
        integer isource(*)
        complex *16 chargesort(nd,*),zk
        complex *16 dipstrsort(nd,*)
        real *8 dipvecsort(2,*)
        complex *16 ima
        complex *16 pot(nd,*)
        complex *16 grad(nd,2,*)
        complex *16 hess(nd,3,*)
        real *8 targetsort(2,*)
        integer itarget(*)
        complex *16 pottarg(nd,*)
        complex *16 gradtarg(nd,2,*)
        complex *16 hesstarg(nd,3,*)
        real *8 wlists(*)
        integer iaddr(2,nboxes)
        complex *16 rmlexp(*)
        complex *16 mptemp(lmptemp)
        real *8 timeinfo(10)
        real *8 center(3)
        integer laddr(2,200)
        real *8 scale(0:200)
        real *8 bsize(0:200)
        integer nterms(0:200)
        integer list(10 000)
        complex *16 ptemp,gtemp(2),htemp(3)
        integer box(15)
        real *8 center0(2),corners0(2,4)
        integer box1(15)
        real *8 center1(2),corners1(2,4)
        integer nterms_eval(4,0:200)
c
        data ima/(0.0d0,1.0d0)/
ccc        save
c
c
c     ifprint is an internal information printing flag.
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c     Prints timing breakdown, list information, and other things if ifprint=2.
c
        ifprint=1
c
c
c       ... set the potential and gradient to zero
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idim)
        do i=1,nsource
        do idim=1,nd
          if( ifpot .eq. 1) pot(idim,i)=0
          if( ifgrad .eq. 1) then
             grad(idim,1,i)=0
             grad(idim,2,i)=0
          endif
          if( ifhess .eq. 1) then
             hess(idim,1,i)=0
             hess(idim,2,i)=0
             hess(idim,3,i)=0
          endif
        enddo
        enddo
C$OMP END PARALLEL DO
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idim)
        do i=1,ntarget
        do idim=1,nd
          if( ifpottarg .eq. 1) pottarg(idim,i)=0
          if( ifgradtarg .eq. 1) then
             gradtarg(idim,1,i)=0
             gradtarg(idim,2,i)=0
          endif
          if( ifhesstarg .eq. 1) then
             hesstarg(idim,1,i)=0
             hesstarg(idim,2,i)=0
             hesstarg(idim,3,i)=0
          endif
        enddo
        enddo
C$OMP END PARALLEL DO
c
        do i=1,10
        timeinfo(i)=0
        enddo
c
c
        if( ifevalfar .eq. 0 ) goto 8000
c
c
        if (ifprint .ge. 2)
     $     call prinf('nterms_eval=*',nterms_eval,4*(nlev+1))
c
c       ... set all multipole and local expansions to zero
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level,ier)
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(4)
        do ibox = 1,nboxes
        call d2tgetb(ier,ibox,box,center0,corners0,wlists)
        level=box(1)
        if( level .ge. 2 ) then
        call h2dzero(nd,rmlexp(iaddr(1,ibox)),nterms(level))
        call h2dzero(nd,rmlexp(iaddr(2,ibox)),nterms(level))
        endif
        enddo
C$OMP END PARALLEL DO
c
c
        if(ifprint .ge. 1)
     $     call prinf('=== STEP 1 (form mp) ====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 1, locate all charges, assign them to boxes, and
c       form multipole expansions
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level,npts,nkids,radius)
C$OMP$PRIVATE(mptemp,lused,ier,i,j,ptemp,gtemp,htemp,cd)
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1)
        do 1200 ibox=1,nboxes
c
        call d2tgetb(ier,ibox,box,center0,corners0,wlists)
        call d2tnkids(box,nkids)
c
        level=box(1)
        if( level .lt. 2 ) goto 1200
c
c
        if (ifprint .ge. 2) then
           call prinf('ibox=*',ibox,1)
           call prinf('box=*',box,15)
           call prinf('nkids=*',nkids,1)
        endif
c
        if (nkids .eq. 0) then
c        ipts=box(9)
c        npts=box(10)
c        call prinf('ipts=*',ipts,1)
c        call prinf('npts=*',npts,1)
        npts=box(10)
        if (ifprint .ge. 2) then
           call prinf('npts=*',npts,1)
           call prinf('isource=*',isource(box(9)),box(10))
        endif
        endif
c
c       ... prune all sourceless boxes
c
        if( box(10) .eq. 0 ) goto 1200
c
        if (nkids .eq. 0) then
c
c       ... form multipole expansions
c
            call h2dzero(nd,rmlexp(iaddr(1,ibox)),nterms(level))
            if_use_trunc = 0

            if( ifcharge .eq. 1 ) then
            call h2dformmp_add(ier,nd,zk,scale(level),
     $          sourcesort(1,box(9)),
     1          chargesort(1,box(9)),npts,center0,nterms(level),
     2          rmlexp(iaddr(1,ibox)))
            endif
c
            if (ifdipole .eq. 1 ) then
            call h2dformmp_dp_add(ier,nd,zk,scale(level),
     $          sourcesort(1,box(9)),
     1          dipstrsort(1,box(9)),dipvecsort(1,box(9)),
     $          npts,center0,nterms(level),
     2          rmlexp(iaddr(1,ibox)))
            endif
         endif
c
 1200    continue
C$OMP END PARALLEL DO
 1300    continue
c
         t2=second()
C$        t2=omp_get_wtime()
ccc        call prin2('time=*',t2-t1,1)
         timeinfo(1)=t2-t1
c
         if(ifprint .ge. 1)
     $      call prinf('=== STEP 2 (form lo) ====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 2, adaptive part, form local expansions,
c           or evaluate the potentials and gradients directly
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level0,itype,list,nlist)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1,ifdirect3,radius)
C$OMP$PRIVATE(lused,ier,i,j,ptemp,gtemp,htemp,cd,ilist,npts)
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1)
         do 3251 ibox=1,nboxes
c
         call d2tgetb(ier,ibox,box,center0,corners0,wlists)
c
ccc         if( box(10) .eq. 0 ) goto 3251
c
         itype=4
         call d2tgetl(ier,ibox,itype,list,nlist,wlists)
         if (nlist .gt. 0) then
            if (ifprint .ge. 2) then
               call prinf('ibox=*',ibox,1)
               call prinf('list3=*',list,nlist)
            endif
         endif
c
c       ... prune all sourceless boxes
c
ccc         if( box(10) .eq. 0 ) nlist=0
c
c
c       ... note that lists 3 and 4 are dual
c
c       ... form local expansions for all boxes in list 3
c       ... if target is childless, evaluate directly (if cheaper)
c
         do 3250 ilist=1,nlist
            jbox=list(ilist)
            call d2tgetb(ier,jbox,box1,center1,corners1,wlists)
c
            npts=box1(10)
            if( npts .eq. 0 ) goto 3250
c
            level0=box(1)
            level1=box1(1)
c
c            ifdirect3 = 0
c
c            if( box1(10) .lt. (nterms(level1)+1)/2 .and.
c     $          box(10) .lt. (nterms(level1)+1)/2 ) ifdirect3 = 1
c
            ifdirect3 = 0
c
            if( ifdirect3 .eq. 0 ) then
               npts=box1(10)
               if_use_trunc = 0
c
               if( ifcharge .eq. 1 ) then
               call h2dformta_add(ier,nd,zk,scale(level0),
     $            sourcesort(1,box1(9)),
     1            chargesort(1,box1(9)),npts,center0,nterms(level0),
     2            rmlexp(iaddr(2,ibox)))
               endif
c
               if( ifdipole .eq. 1 ) then
               call h2dformta_dp_add(ier,nd,zk,scale(level0),
     1            sourcesort(1,box1(9)),dipstrsort(1,box1(9)),
     2            dipvecsort(1,box1(9)),npts,center0,
     3            nterms(level0),rmlexp(iaddr(2,ibox)))
               endif
c
            else

            call hfmm2dpart_direct_std(nd,zk,box1,box,sourcesort,
     $         ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort,
     $         ifpot,pot,ifgrad,grad,ifhess,hess,
     $         targetsort,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $         ifhesstarg,hesstarg)
c
            endif
 3250    continue
c
 3251    continue
C$OMP END PARALLEL DO
c
         t2=second()
C$        t2=omp_get_wtime()
ccc        call prin2('time=*',t2-t1,1)
         timeinfo(2)=t2-t1
c
c
        if(ifprint .ge. 1)
     $      call prinf('=== STEPS 3,4,5 ====*',i,0)
        ifprune_list2 = 1
        if (ifpot.eq.1) ifprune_list2 = 0
        if (ifgrad.eq.1) ifprune_list2 = 0
        if (ifhess.eq.1) ifprune_list2 = 0
        call hfmm2d_list2(nd,zk,
     $     bsize,nlev,laddr,scale,nterms,rmlexp,iaddr,epsfmm,
     $     timeinfo,wlists,mptemp,lmptemp,
     $     ifprune_list2)
c
        if(ifprint .ge. 1)
     $     call prinf('=== STEP 6 (eval mp) ====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 6, adaptive part, evaluate multipole expansions,
c           or evaluate the potentials and gradients directly
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,itype,list,nlist)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1,ifdirect4,radius)
C$OMP$PRIVATE(lused,ier,i,j,ptemp,gtemp,htemp,cd,ilist,level)
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1)
         do 3252 ibox=1,nboxes
         call d2tgetb(ier,ibox,box,center0,corners0,wlists)
c
c       ... prune all sourceless boxes
c
ccc         if( box(10) .eq. 0 ) goto 3252
c
         itype=3
         call d2tgetl(ier,ibox,itype,list,nlist,wlists)
         if (nlist .gt. 0) then
            if (ifprint .ge. 2) then
               call prinf('ibox=*',ibox,1)
               call prinf('list4=*',list,nlist)
            endif
         endif
c
c       ... prune all sourceless boxes
c
ccc         if( box(10) .eq. 0 ) nlist=0
c
c       ... note that lists 3 and 4 are dual
c
c       ... evaluate multipole expansions for all boxes in list 4
c       ... if source is childless, evaluate directly (if cheaper)
c
         do ilist=1,nlist
            jbox=list(ilist)
            call d2tgetb(ier,jbox,box1,center1,corners1,wlists)
c
            level=box1(1)
c
c            ifdirect4 = 0
c
c            if (box1(10) .lt. (nterms(level)+1)/2 .and.
c     $         box(10) .lt. (nterms(level)+1)/2 ) ifdirect4 = 1
c
            ifdirect4 = 0
c
            if (ifdirect4 .eq. 0) then

            if( box(10) .gt. 0 )
     $         call h2dmpevalall(nd,zk,scale(level),center1,
     $         rmlexp(iaddr(1,jbox)),nterms(level),
     $         sourcesort(1,box(9)),box(10),
     $         ifpot,pot(1,box(9)),
     $         ifgrad,grad(1,1,box(9)),
     $         ifhess,hess(1,1,box(9)))

            if( box(12) .gt. 0 )
     $         call h2dmpevalall(nd,zk,scale(level),center1,
     $         rmlexp(iaddr(1,jbox)),nterms(level),
     $         targetsort(1,box(11)),box(12),
     $         ifpottarg,pottarg(1,box(11)),
     $         ifgradtarg,gradtarg(1,1,box(11)),
     $         ifhesstarg,hesstarg(1,1,box(11)))

            else

            call hfmm2dpart_direct_std(nd,zk,box1,box,sourcesort,
     $         ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort,
     $         ifpot,pot,ifgrad,grad,ifhess,hess,
     $         targetsort,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $         ifhesstarg,hesstarg)

            endif
        enddo
c
 3252   continue
C$OMP END PARALLEL DO
c
        t2=second()
C$        t2=omp_get_wtime()
ccc     call prin2('time=*',t2-t1,1)
        timeinfo(6)=t2-t1
c

        if(ifprint .ge. 1)
     $     call prinf('=== STEP 7 (eval lo) ====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 7, evaluate local expansions
c       and all gradients directly
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level,npts,nkids,ier)
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1)
        do 6201 ibox=1,nboxes
c
        call d2tgetb(ier,ibox,box,center0,corners0,wlists)
        call d2tnkids(box,nkids)
c
        if (ifprint .ge. 2) then
           call prinf('ibox=*',ibox,1)
           call prinf('box=*',box,15)
           call prinf('nkids=*',nkids,1)
        endif
c
        if (nkids .eq. 0) then
            npts=box(10)
            if (ifprint .ge. 2) then
               call prinf('npts=*',npts,1)
               call prinf('isource=*',isource(box(9)),box(10))
            endif
        endif
c
        if (nkids .eq. 0) then
c
c       ... evaluate local expansions
c
        level=box(1)
        npts=box(10)
c
        if (level .ge. 2) then

        if( box(10) .gt. 0 )
     $     call h2dtaevalall(nd,zk,scale(level),center0,
     $     rmlexp(iaddr(2,ibox)),nterms(level),
     $     sourcesort(1,box(9)),box(10),
     $     ifpot,pot(1,box(9)),ifgrad,grad(1,1,box(9)),
     $     ifhess,hess(1,1,box(9)))

        if( box(12) .gt. 0 )
     $     call h2dtaevalall(nd,zk,scale(level),center0,
     $     rmlexp(iaddr(2,ibox)),nterms(level),
     $     targetsort(1,box(11)),box(12),
     $     ifpottarg,pottarg(1,box(11)),ifgradtarg,
     $     gradtarg(1,1,box(11)),
     $     ifhesstarg,hesstarg(1,1,box(11)))

        endif
c
        endif
c
 6201   continue
C$OMP END PARALLEL DO
        t2=second()
C$        t2=omp_get_wtime()
ccc     call prin2('time=*',t2-t1,1)
        timeinfo(7)=t2-t1
c
c
 8000   continue
c
c
        if( ifevalloc .eq. 0 ) goto 9000
c
        if(ifprint .ge. 1)
     $     call prinf('=== STEP 8 (direct) =====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 8, evaluate direct interactions
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,nkids,list,nlist,npts)
C$OMP$PRIVATE(jbox,box1,center1,corners1)
C$OMP$PRIVATE(ier,ilist,itype)
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1)
        do 6202 ibox=1,nboxes
c
        call d2tgetb(ier,ibox,box,center0,corners0,wlists)
        call d2tnkids(box,nkids)
c
        if (ifprint .ge. 2) then
           call prinf('ibox=*',ibox,1)
           call prinf('box=*',box,15)
           call prinf('nkids=*',nkids,1)
        endif
c
        if (nkids .eq. 0) then
            npts=box(10)
            if (ifprint .ge. 2) then
               call prinf('npts=*',npts,1)
               call prinf('isource=*',isource(box(9)),box(10))
            endif
        endif
c
c
        if (nkids .eq. 0) then
c
c       ... evaluate self interactions
c
        call hfmm2dpart_direct_self_std(nd,zk,box,sourcesort,
     $     ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     targetsort,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg)
c
c
c       ... retrieve list #1
c
c       ... evaluate interactions with the nearest neighbours
c
        itype=1
        call d2tgetl(ier,ibox,itype,list,nlist,wlists)
        if (ifprint .ge. 2) call prinf('list1=*',list,nlist)
c
c       ... for all pairs in list #1,
c       evaluate the potentials and gradients directly
c
            do 6203 ilist=1,nlist
               jbox=list(ilist)
               call d2tgetb(ier,jbox,box1,center1,corners1,wlists)
c
c       ... prune all sourceless boxes
c
         if( box1(10) .eq. 0 ) goto 6203
c
            call hfmm2dpart_direct_std(nd,zk,box1,box,sourcesort,
     $         ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort,
     $         ifpot,pot,ifgrad,grad,ifhess,hess,
     $         targetsort,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $         ifhesstarg,hesstarg)
c
 6203       continue
        endif
c
 6202   continue
C$OMP END PARALLEL DO
c
ccc        call prin2('inside fmm, pot=*',pot,2*nd*nsource)
ccc        call prin2('inside fmm, grad=*',grad,2*2*nd*nsource)
ccc        call prin2('inside fmm, hess=*',hess,2*3*nd*nsource)
c
c
        t2=second()
C$        t2=omp_get_wtime()
ccc     call prin2('time=*',t2-t1,1)
        timeinfo(8)=t2-t1
c
 9000   continue
c
ccc        call prinf('=== DOWNWARD PASS COMPLETE ===*',i,0)
c
        if (ifprint .ge. 1) call prin2('timeinfo=*',timeinfo,8)
c
        d=0
        do i=1,8
        d=d+timeinfo(i)
        enddo
c
        if (ifprint .ge. 1) call prin2('sum(timeinfo)=*',d,1)
c
        if (ifprint .ge. 1) then
        call prinf('nboxes=*',nboxes,1)
        call prinf('nsource=*',nsource,1)
        call prinf('ntarget=*',ntarget,1)
        endif
c
        return
        end
c
c
c
c
c
        subroutine hfmm2dpart_direct_self_std(nd,zk,box,
     $     source,ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     target,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg)
        implicit real *8 (a-h,o-z)
c
        integer box(15),box1(15)
c
        real *8 source(2,*),dipvec(2,*)
        complex *16 charge(nd,*),dipstr(nd,*),zk
        real *8 target(2,*)
c
        complex *16 pot(nd,*),grad(nd,2,*),hess(nd,3,*)
        complex *16 pottarg(nd,*),gradtarg(nd,2,*),hesstarg(nd,3,*)
        complex *16 ptemp,gtemp(2),htemp(3)
c
c
        do 6160 j=box(9),box(9)+box(10)-1
        do 6150 i=box(9),box(9)+box(10)-1
            if (i .eq. j) goto 6150
            do idim=1,nd
            call hpotgrad2d_sdp(zk,source(1,i),
     $         ifcharge,charge(idim,i),ifdipole,dipstr(idim,i),
     $         dipvec(1,i),
     1         source(1,j),ifpot,ptemp,ifgrad,gtemp,ifhess,htemp)
            if (ifpot .eq. 1) pot(idim,j)=pot(idim,j)+ptemp
            if (ifgrad .eq. 1) then
               grad(idim,1,j)=grad(idim,1,j)+gtemp(1)
               grad(idim,2,j)=grad(idim,2,j)+gtemp(2)
            endif
            if (ifhess .eq. 1) then
               hess(idim,1,j)=hess(idim,1,j)+htemp(1)
               hess(idim,2,j)=hess(idim,2,j)+htemp(2)
               hess(idim,3,j)=hess(idim,3,j)+htemp(3)
            endif
            enddo
 6150   continue
 6160   continue
        do j=box(11),box(11)+box(12)-1
        do i=box( 9),box( 9)+box(10)-1
            do idim=1,nd
            call hpotgrad2d_sdp(zk,source(1,i),
     $         ifcharge,charge(idim,i),ifdipole,dipstr(idim,i),
     $         dipvec(1,i),
     1         target(1,j),ifpottarg,ptemp,ifgradtarg,gtemp,
     $         ifhesstarg,htemp)
            if (ifpottarg .eq. 1) pottarg(idim,j)=pottarg(idim,j)+ptemp
            if (ifgradtarg .eq. 1) then
               gradtarg(idim,1,j)=gradtarg(idim,1,j)+gtemp(1)
               gradtarg(idim,2,j)=gradtarg(idim,2,j)+gtemp(2)
            endif
            if (ifhesstarg .eq. 1) then
               hesstarg(idim,1,j)=hesstarg(idim,1,j)+htemp(1)
               hesstarg(idim,2,j)=hesstarg(idim,2,j)+htemp(2)
               hesstarg(idim,3,j)=hesstarg(idim,3,j)+htemp(3)
            endif
            enddo
        enddo
        enddo
c
        return
        end
c
c
c
c
c
        subroutine hfmm2dpart_direct_std(nd,zk,box,box1,
     $     source,ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     target,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg)
        implicit real *8 (a-h,o-z)
c
        integer box(15),box1(15)
c
        real *8 source(2,*),dipvec(2,*)
        complex *16 charge(nd,*),dipstr(nd,*),zk
        real *8 target(2,*)
c
        complex *16 pot(nd,*),grad(nd,2,*),hess(nd,3,*)
        complex *16 pottarg(nd,*),gradtarg(nd,2,*),hesstarg(nd,3,*)
        complex *16 ptemp,gtemp(2),htemp(3)
c
c
        do j=box1(9),box1(9)+box1(10)-1
          do i=box(9),box(9)+box(10)-1
            do idim=1,nd
            call hpotgrad2d_sdp(zk,source(1,i),
     $         ifcharge,charge(idim,i),ifdipole,dipstr(idim,i),
     $         dipvec(1,i),
     1         source(1,j),ifpot,ptemp,ifgrad,gtemp,ifhess,htemp)
            if (ifpot .eq. 1) pot(idim,j)=pot(idim,j)+ptemp
            if (ifgrad .eq. 1) then
               grad(idim,1,j)=grad(idim,1,j)+gtemp(1)
               grad(idim,2,j)=grad(idim,2,j)+gtemp(2)
            endif
            if (ifhess .eq. 1) then
               hess(idim,1,j)=hess(idim,1,j)+htemp(1)
               hess(idim,2,j)=hess(idim,2,j)+htemp(2)
               hess(idim,3,j)=hess(idim,3,j)+htemp(3)
            endif
            enddo
          enddo
        enddo
c
        do j=box1(11),box1(11)+box1(12)-1
          do i=box(9),box(9)+box(10)-1
            do idim=1,nd
            call hpotgrad2d_sdp(zk,source(1,i),
     $         ifcharge,charge(idim,i),ifdipole,dipstr(idim,i),
     $         dipvec(1,i),
     1         target(1,j),ifpottarg,ptemp,ifgradtarg,gtemp,
     $         ifhesstarg,htemp)
            if (ifpottarg .eq. 1) pottarg(idim,j)=pottarg(idim,j)+ptemp
            if (ifgradtarg .eq. 1) then
               gradtarg(idim,1,j)=gradtarg(idim,1,j)+gtemp(1)
               gradtarg(idim,2,j)=gradtarg(idim,2,j)+gtemp(2)
            endif
            if (ifhesstarg .eq. 1) then
               hesstarg(idim,1,j)=hesstarg(idim,1,j)+htemp(1)
               hesstarg(idim,2,j)=hesstarg(idim,2,j)+htemp(2)
               hesstarg(idim,3,j)=hesstarg(idim,3,j)+htemp(3)
            endif
            enddo
          enddo
        enddo
c
        return
        end
c
c
c
c
c
        subroutine h2dpartdirect_std_nd(nd,zk,nsource,
     $     source,ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     ntarget,target,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg)
        implicit real *8 (a-h,o-z)
c
c       Helmholtz interactions in R^2: evaluate all pairwise particle
c       interactions (ignoring self-interaction)
c       and interactions with targets via direct O(N^2) algorithm.
c
c hfmm2d: charge and dipstr are complex valued, x \in R^2
c
c \phi(x_i) = \sum_{j\ne i}   charge_j H^{(1)}_0(k |x_i-x_j| ) * (i/4)
c             + dipstr_j (dipvec_j \dot \grad_j H^{(1)}_0(k |x_i-x_j| )*(i/4))
c
c or, more precisely,
c
c \phi(x_i) = \sum_{j\ne i}   charge_j H^{(1)}_0(k |x_i-x_j| ) * (i/4)
c +dipstr_j (dipvec_j \dot (x_i-x_j) ) H^{(1)}_1(k |x_i-x_j|) *k/|x_i-x_j|*(i/4)
c
c       We use H_0(kr)*(i/4) for the Green's function.
c       Self-interactions are not included.
c
c       INPUT PARAMETERS:
c
c       nd: integer: number of densities
c       zk: complex *16: Helmholtz parameter
c       nsource: integer:  number of sources
c       source: real *8 (2,nsource):  source locations
c       ifcharge:  charge computation flag
c                  ifcharge = 1   =>  include charge contribution
c                                     otherwise do not
c       charge: complex *16 (nd,nsource): charge strengths
c       ifdipole:  dipole computation flag
c                  ifdipole = 1   =>  include dipole contribution
c                                     otherwise do not
c       dipstr: complex *16 (nd,nsource): dipole strengths
c       dipvec: real *8 (2,nsource): dipole orientation vectors.
c
c       ifpot:  potential flag (1=compute potential, otherwise no)
c       ifgrad:  gradient flag (1=compute gradient, otherwise no)
c       ifhess:  hessian flag (1=compute hessian, otherwise no)
c       ntarget: integer:  number of targets
c       target: real *8 (2,ntarget):  target locations
c       ifpottarg:  target potential flag
c                   (1=compute potential, otherwise no)
c       ifgradtarg:  target gradient flag
c                   (1=compute gradient, otherwise no)
c       ifhesstarg:  target hessian flag
c                   (1=compute hessian, otherwise no)
c
c       OUTPUT PARAMETERS:
c
c       pot: complex *16 (nd,nsource): potential at source locations
c       grad: complex *16 (nd,2,nsource): gradient at source locations
c       hess: complex *16 (nd,3,nsource): hessian at source locations
c       pottarg: complex *16 (nd,ntarget): potential at target locations
c       gradtarg: complex *16 (nd,2,ntarget): gradient at target locations
c       hesstarg: complex *16 (nd,3,ntarget): hessian at target locations
c
c
c
        real *8 source(2,*),dipvec(2,*)
        complex *16 charge(nd,*),dipstr(nd,*),zk
        real *8 target(2,*)
c
        complex *16 pot(nd,*),grad(nd,2,*),hess(nd,3,*)
        complex *16 pottarg(nd,*),gradtarg(nd,2,*),hesstarg(nd,3,*)
        complex *16 ptemp,gtemp(2),htemp(3)
c
c
        do i=1,nsource
        do idim=1,nd
        if( ifpot .eq. 1) pot(idim,i)=0
        if( ifgrad .eq. 1) then
           grad(idim,1,i)=0
           grad(idim,2,i)=0
        endif
        if( ifhess .eq. 1) then
           hess(idim,1,i)=0
           hess(idim,2,i)=0
           hess(idim,3,i)=0
        endif
        enddo
        enddo
c
        do i=1,ntarget
        do idim=1,nd
        if( ifpottarg .eq. 1) pottarg(idim,i)=0
        if( ifgradtarg .eq. 1) then
           gradtarg(idim,1,i)=0
           gradtarg(idim,2,i)=0
        endif
        if( ifhesstarg .eq. 1) then
           hesstarg(idim,1,i)=0
           hesstarg(idim,2,i)=0
           hesstarg(idim,3,i)=0
        endif
        enddo
        enddo
c
        if( ifpot .eq. 1 .or. ifgrad .eq. 1 .or. ifhess .eq. 1) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(idim,i,j,ptemp,gtemp,htemp)
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(4)
        do 6160 j=1,nsource
        do 6150 i=1,nsource
            if (i .eq. j) goto 6150
            do idim=1,nd
            call hpotgrad2d_sdp(zk,source(1,i),
     $         ifcharge,charge(idim,i),ifdipole,dipstr(idim,i),
     $         dipvec(1,i),
     1         source(1,j),ifpot,ptemp,ifgrad,gtemp,ifhess,htemp)
            if (ifpot .eq. 1) pot(idim,j)=pot(idim,j)+ptemp
            if (ifgrad .eq. 1) then
               grad(idim,1,j)=grad(idim,1,j)+gtemp(1)
               grad(idim,2,j)=grad(idim,2,j)+gtemp(2)
            endif
            if (ifhess .eq. 1) then
               hess(idim,1,j)=hess(idim,1,j)+htemp(1)
               hess(idim,2,j)=hess(idim,2,j)+htemp(2)
               hess(idim,3,j)=hess(idim,3,j)+htemp(3)
            endif
            enddo
 6150   continue
 6160   continue
C$OMP END PARALLEL DO
        endif
c
        if( ifpottarg .eq. 1 .or. ifgradtarg .eq. 1
     $      .or. ifhesstarg .eq. 1) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(idim,i,j,ptemp,gtemp,htemp)
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(4)
        do j=1,ntarget
        do i=1,nsource
            do idim=1,nd
            call hpotgrad2d_sdp(zk,source(1,i),
     $         ifcharge,charge(idim,i),ifdipole,dipstr(idim,i),
     $         dipvec(1,i),
     1         target(1,j),ifpottarg,ptemp,ifgradtarg,gtemp,
     $         ifhesstarg,htemp)
            if (ifpottarg .eq. 1) pottarg(idim,j)=pottarg(idim,j)+ptemp
            if (ifgradtarg .eq. 1) then
               gradtarg(idim,1,j)=gradtarg(idim,1,j)+gtemp(1)
               gradtarg(idim,2,j)=gradtarg(idim,2,j)+gtemp(2)
            endif
            if (ifhesstarg .eq. 1) then
               hesstarg(idim,1,j)=hesstarg(idim,1,j)+htemp(1)
               hesstarg(idim,2,j)=hesstarg(idim,2,j)+htemp(2)
               hesstarg(idim,3,j)=hesstarg(idim,3,j)+htemp(3)
            endif
            enddo
        enddo
        enddo
C$OMP END PARALLEL DO
        endif
c
        return
        end
c
c
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Tree reuse API - build tree once, evaluate multiple times
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine hfmm2dpart_tree_build_std(ier,iprec,zk,
     $     nsource,source,ntarget,target,tree_handle)
        use fmm2d_tree_module
        implicit real *8 (a-h,o-z)
        integer tree_handle
c
c       Build FMM tree structure (geometry-dependent only, NO nd).
c       Tree can be reused for multiple evaluations with different
c       nd values and charge distributions.
c
c       INPUT:
c         iprec - FMM precision flag (see hfmm2dparttarg_std_nd docs)
c         zk - complex *16: Helmholtz parameter
c         nsource - number of sources
c         source(2,nsource) - source locations
c         ntarget - number of targets
c         target(2,ntarget) - target locations
c
c       OUTPUT:
c         ier - error code (0=success, 1=out of tree handles, 4=tree alloc failed)
c         tree_handle - integer handle to store tree (pass to tree_eval and tree_free)
c
        real *8 source(2,*)
        real *8 target(2,*)
        complex *16 zk
        type(fmm2d_tree_data), pointer :: tree
        integer laddr_tmp(2,0:200)
        real *8 scale_tmp(0:200)
        real *8 bsize_tmp(0:200)
        integer nterms_tmp(0:200)
c
        ier = 0
c
c       Allocate tree handle
c
        call fmm2d_tree_alloc(tree_handle,ier)
        if (ier .ne. 0) then
          call prinf(' cannot allocate tree handle *',ier,1)
          return
        endif
c
c       Get pointer to tree data
c
        call fmm2d_tree_get(tree_handle,tree,ier)
        if (ier .ne. 0) then
          return
        endif
c
c       Store geometry parameters
c
        tree%nsource = nsource
        tree%ntarget = ntarget
        tree%iprec = iprec
        tree%zk = zk
c
c       Set FMM tolerance based on iprec
c
        if( iprec .eq. -2 ) tree%epsfmm=.5d-0
        if( iprec .eq. -1 ) tree%epsfmm=.5d-1
        if( iprec .eq. 0 ) tree%epsfmm=.5d-2
        if( iprec .eq. 1 ) tree%epsfmm=.5d-3
        if( iprec .eq. 2 ) tree%epsfmm=.5d-6
        if( iprec .eq. 3 ) tree%epsfmm=.5d-9
        if( iprec .eq. 4 ) tree%epsfmm=.5d-12
        if( iprec .eq. 5 ) tree%epsfmm=.5d-15
        if( iprec .eq. 6 ) tree%epsfmm=0
c
c       Set box subdivision criterion
c
        if( iprec .eq. -2 ) tree%nbox=3
        if( iprec .eq. -1 ) tree%nbox=5
        if( iprec .eq. 0 ) tree%nbox=8
        if( iprec .eq. 1 ) tree%nbox=10
        if( iprec .eq. 2 ) tree%nbox=15
        if( iprec .eq. 3 ) tree%nbox=20
        if( iprec .eq. 4 ) tree%nbox=25
        if( iprec .eq. 5 ) tree%nbox=45
        if( iprec .eq. 6 ) tree%nbox=nsource+ntarget
c
c       Build quad-tree structure
c
        ntot = 100*(nsource+ntarget)+10000
        do ii = 1,10
           allocate (tree%wlists(ntot),stat=ier)
           if (ier .ne. 0) then
             ier = 4
             return
           endif
           call hfmm2dparttree(ier,iprec,zk,
     $        nsource,source,ntarget,target,
     $        tree%nbox,tree%epsfmm,tree%iisource,tree%iitarget,
     $        tree%iwlists,tree%lwlists,
     $        tree%nboxes,laddr_tmp,tree%nlev,tree%center,tree%size,
     $        tree%wlists,ntot,lused7)
           if (ier.ne.0) then
              deallocate(tree%wlists)
              ntot = ntot*1.5
           else
             goto 1200
           endif
        enddo
1200    continue
        if (ier.ne.0) then
           ier = 4
           call fmm2d_tree_dealloc(tree_handle,ier2)
           return
        endif
c
c       Allocate and copy laddr
c
        allocate(tree%laddr(2,0:tree%nlev),stat=ier)
        if (ier .ne. 0) then
          ier = 4
          call fmm2d_tree_dealloc(tree_handle,ier2)
          return
        endif
        do i = 0,tree%nlev
          tree%laddr(1,i) = laddr_tmp(1,i)
          tree%laddr(2,i) = laddr_tmp(2,i)
        enddo
c
c       Compute scale, bsize, nterms arrays (geometry-dependent)
c
        allocate(tree%scale(0:tree%nlev),stat=ier)
        if (ier .ne. 0) then
          ier = 4
          call fmm2d_tree_dealloc(tree_handle,ier2)
          return
        endif
        allocate(tree%bsize(0:tree%nlev),stat=ier)
        if (ier .ne. 0) then
          ier = 4
          call fmm2d_tree_dealloc(tree_handle,ier2)
          return
        endif
        allocate(tree%nterms(0:tree%nlev),stat=ier)
        if (ier .ne. 0) then
          ier = 4
          call fmm2d_tree_dealloc(tree_handle,ier2)
          return
        endif
c
        do i = 0,tree%nlev
          tree%bsize(i) = tree%size/2.0d0**i
          boxsize = abs((tree%size/2.0**i)*zk)
          tree%scale(i) = 1.0d0
          if (boxsize .lt. 1) tree%scale(i) = boxsize
          call h2dterms(tree%bsize(i),zk,tree%epsfmm,
     $       tree%nterms(i),ier)
        enddo
c
        return
        end
c
c
c
c
c
        subroutine hfmm2dpart_tree_free_std(tree_handle)
        use fmm2d_tree_module
        implicit real *8 (a-h,o-z)
        integer tree_handle
c
c       Free FMM tree structure and deallocate all associated memory.
c
c       INPUT:
c         tree_handle - handle returned by hfmm2dpart_tree_build_std
c
        call fmm2d_tree_dealloc(tree_handle,ier)
        return
        end
c
c
c
c
c
        subroutine hfmm2dpart_tree_eval_std(ier,tree_handle,nd,
     $     source,ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     target,ifpottarg,pottarg,ifgradtarg,gradtarg,ifhesstarg,
     $     hesstarg)
        use fmm2d_tree_module
        implicit real *8 (a-h,o-z)
        integer tree_handle
c
c       Evaluate FMM using pre-built tree with arbitrary nd and charges.
c       Tree must have been built with hfmm2dpart_tree_build_std.
c
c       IMPORTANT: source and target locations must be IDENTICAL to
c       those used in hfmm2dpart_tree_build_std. Only charges can change.
c
c       INPUT:
c         tree_handle - handle from hfmm2dpart_tree_build_std
c         nd - number of densities
c         source(2,nsource) - source locations (must match tree_build)
c         target(2,ntarget) - target locations (must match tree_build)
c         ifcharge,ifdipole - flags for charge/dipole evaluation
c         charge(nd,nsource) - charge strengths
c         dipstr(nd,nsource) - dipole strengths
c         dipvec(2,nsource) - dipole orientation vectors
c         ifpot,ifgrad,ifhess - flags for output at sources
c         ifpottarg,ifgradtarg,ifhesstarg - flags for output at targets
c
c       OUTPUT:
c         ier - error code (0=success)
c         pot(nd,nsource) - potential at sources
c         grad(nd,2,nsource) - gradient at sources
c         hess(nd,3,nsource) - hessian at sources
c         pottarg(nd,ntarget) - potential at targets
c         gradtarg(nd,2,ntarget) - gradient at targets
c         hesstarg(nd,3,ntarget) - hessian at targets
c
        real *8 source(2,*)
        real *8 target(2,*)
        complex *16 charge(nd,*)
        complex *16 dipstr(nd,*)
        real *8 dipvec(2,*)
        complex *16 pot(nd,*)
        complex *16 grad(nd,2,*)
        complex *16 hess(nd,3,*)
        complex *16 pottarg(nd,*)
        complex *16 gradtarg(nd,2,*)
        complex *16 hesstarg(nd,3,*)
c
        type(fmm2d_tree_data), pointer :: tree
        real *8, allocatable :: w(:)
        complex *16, allocatable :: wrmlexp(:)
c
        ier = 0
c
c       Get tree from handle
c
        call fmm2d_tree_get(tree_handle,tree,ier)
        if (ier .ne. 0) then
          call prinf(' invalid tree handle *',tree_handle,1)
          return
        endif
c
c       Allocate workspace for evaluation (nd-dependent)
c
        lused7 = 1
c
        isourcesort = lused7 + 5
        lsourcesort = 2*tree%nsource
        itargetsort = isourcesort + lsourcesort
        ltargetsort = 2*tree%ntarget
        ichargesort = itargetsort + ltargetsort
        lchargesort = 2*tree%nsource*nd
        idipvecsort = ichargesort + lchargesort
        if (ifdipole.eq.1) then
          ldipvec = 2*tree%nsource
          ldipstr = 2*tree%nsource*nd
        else
          ldipvec = 2
          ldipstr = 2
        endif
        idipstrsort = idipvecsort + ldipvec
        lused7 = idipstrsort + ldipstr
c
        ipot = lused7
        lpot = 2*tree%nsource*nd
        lused7 = lused7 + lpot
c
        igrad = lused7
        if (ifgrad .eq. 1) then
          lgrad = 2*2*tree%nsource*nd
        else
          lgrad = 4
        endif
        lused7 = lused7 + lgrad
c
        ihess = lused7
        if (ifhess .eq. 1) then
          lhess = 2*3*tree%nsource*nd
        else
          lhess = 6
        endif
        lused7 = lused7 + lhess
c
        ipottarg = lused7
        lpottarg = 2*tree%ntarget*nd
        lused7 = lused7 + lpottarg
c
        igradtarg = lused7
        if (ifgradtarg .eq. 1) then
          lgradtarg = 2*2*tree%ntarget*nd
        else
          lgradtarg = 4
        endif
        lused7 = lused7 + lgradtarg
c
        ihesstarg = lused7
        if (ifhesstarg .eq. 1) then
          lhesstarg = 2*3*tree%ntarget*nd
        else
          lhesstarg = 6
        endif
        lused7 = lused7 + lhesstarg
c
c       Compute nmax for mptemp allocation
c
        nmax = 0
        do i = 0,tree%nlev
          if (tree%nterms(i).gt. nmax .and. i.ge. 2) then
            nmax = tree%nterms(i)
          endif
        enddo
c
        iiaddr = lused7
        imptemp = iiaddr + 2*tree%nboxes
        lmptemp = (2*nmax+1)*2*nd
        lused7 = imptemp + lmptemp
c
        allocate(w(lused7),stat=ier)
        if (ier.ne.0) then
          ier = 8
          return
        endif
c
c       Reorder sources and charges using tree's sorting indices
c
        call h2dreorder(nd,tree%nsource,source,ifcharge,charge,
     $     tree%wlists(tree%iisource),ifdipole,dipstr,dipvec,
     1     w(isourcesort),w(ichargesort),w(idipvecsort),
     $     w(idipstrsort))
c
        if (tree%ntarget .gt. 0) then
          call h2dreordertarg(tree%ntarget,target,
     $       tree%wlists(tree%iitarget),w(itargetsort))
        endif
c
c       Allocate multipole/local expansion workspace (nd-dependent)
c
        call h2dmpalloc(nd,tree%wlists(tree%iwlists),w(iiaddr),
     $     tree%nboxes,lmptot,tree%nterms)
c
        irmlexp = 1
        lused7_rml = irmlexp + lmptot
        allocate(wrmlexp(lused7_rml),stat=ier)
        if (ier.ne.0) then
          ier = 16
          deallocate(w)
          return
        endif
c
c       Call main FMM evaluation routine
c
        ifevalfar = 1
        ifevalloc = 1
c
        call hfmm2dparttargmain_std(ier,tree%iprec,tree%zk,nd,
     $     ifevalfar,ifevalloc,
     $     tree%nsource,w(isourcesort),tree%wlists(tree%iisource),
     $     ifcharge,w(ichargesort),
     $     ifdipole,w(idipstrsort),w(idipvecsort),
     $     ifpot,w(ipot),ifgrad,w(igrad),ifhess,w(ihess),
     $     tree%ntarget,w(itargetsort),tree%wlists(tree%iitarget),
     $     ifpottarg,w(ipottarg),ifgradtarg,w(igradtarg),
     $     ifhesstarg,w(ihesstarg),
     $     tree%epsfmm,w(iiaddr),wrmlexp(irmlexp),w(imptemp),lmptemp,
     $     tree%nboxes,tree%laddr,tree%nlev,tree%scale,tree%bsize,
     $     tree%nterms,tree%wlists(tree%iwlists),tree%lwlists)
c
c       Unsort results back to original ordering
c
        if (ifpot .eq. 1)
     $     call h2dpsort(nd,tree%nsource,tree%wlists(tree%iisource),
     $     w(ipot),pot)
        if (ifgrad .eq. 1)
     $     call h2dgsort(nd,tree%nsource,tree%wlists(tree%iisource),
     $     w(igrad),grad)
        if (ifhess .eq. 1)
     $     call h2dhsort(nd,tree%nsource,tree%wlists(tree%iisource),
     $     w(ihess),hess)
c
        if (ifpottarg .eq. 1)
     $     call h2dpsort(nd,tree%ntarget,tree%wlists(tree%iitarget),
     $     w(ipottarg),pottarg)
        if (ifgradtarg .eq. 1)
     $     call h2dgsort(nd,tree%ntarget,tree%wlists(tree%iitarget),
     $     w(igradtarg),gradtarg)
        if (ifhesstarg .eq. 1)
     $     call h2dhsort(nd,tree%ntarget,tree%wlists(tree%iitarget),
     $     w(ihesstarg),hesstarg)
c
c       Deallocate temporary arrays
c
        deallocate(wrmlexp)
        deallocate(w)
c
        if (ier .ne. 0) ier = 0
c
        return
        end
c
c
c
c
c
