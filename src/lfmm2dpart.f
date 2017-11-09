cc Copyright (C) 2009-2012: Leslie Greengard and Zydrunas Gimbutas
cc Contact: greengard@cims.nyu.edu
cc 
cc This software is being released under a modified FreeBSD license
cc (see COPYING in home directory). 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date$
c    $Revision$
c
c       
c     This file contains the main FMM routines and some related
c     subroutines for evaluating complex-valued Laplace potentials and
c     gradients due to point charges and dipoles.  (FORTRAN 90 VERSION)
c
c lfmm2d: charge and dipstr are complex valued, x \in R^2
c
c \phi(x_i) = \sum_{j\ne i}   charge_j \log |x_i-x_j|  
c                    + dipstr_j (dipvec_j \dot \grad_j log|x_i-x_j|)
c
c or 
c
c \phi(x_i) = \sum_{j\ne i}   charge_j \log |x_i-x_j|  
c                   + dipstr_j (dipvec_j \dot (x_i-x_j)) * (-1/|x_i-x_j|^2)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     lfmm2dpart - Laplace FMM in R^2: evaluate all pairwise particle
c         interactions (ignoring self-interaction)
c
c     lfmm2dpartself - Laplace FMM in R^2: evaluate all pairwise particle
c         interactions (ignoring self-interaction)
c
c     lfmm2dparttarg - Laplace FMM in R^2: evaluate all pairwise
c         particle interactions (ignoring self-interaction) +
c         interactions with targets
c
c     l2dpartdirect - Laplace interactions in R^2:  evaluate all
c         pairwise particle interactions (ignoring self-interaction) +
c         interactions with targets via direct O(N^2) algorithm
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        this is the end of the debugging code and the beginning 
c        of the Laplace particle FMM in R^2
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine lfmm2dpart(ier,iprec,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
        implicit real *8 (a-h,o-z)
c              
c       Laplace FMM in R^2: evaluate all pairwise particle
c       interactions (ignoring self-interaction). 
c       We use log(r) for the Green's function.
c       Self-interactions are not included.
c   
c lfmm2d: charge and dipstr are complex valued, x \in R^2
c
c \phi(x_i) = \sum_{j\ne i}   charge_j \log |x_i-x_j|  
c                    + dipstr_j (dipvec_j \dot \grad_j log|x_i-x_j|)
c
c or 
c
c \phi(x_i) = \sum_{j\ne i}   charge_j \log |x_i-x_j|  
c                   + dipstr_j (dipvec_j \dot (x_i-x_j)) * (-1/|x_i-x_j|^2)
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
        complex *16 charge(*)
        complex *16 dipstr(*)
        real *8 dipvec(2,*)
        complex *16 ima
        complex *16 pot(*)
        complex *16 grad(2,*)
        complex *16 hess(3,*)
c
        real *8 target(2,1)
        complex *16 pottarg(1)
        complex *16 gradtarg(2,1)
        complex *16 hesstarg(3,1)        
c
        data ima/(0.0d0,1.0d0)/
c       
        ntarget=0
        ifpottarg=0
        ifgradtarg=0
        ifhesstarg=0
c
        call lfmm2dparttarg(ier,iprec,nsource,source,
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
        subroutine lfmm2dpartself(ier,iprec,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
        implicit real *8 (a-h,o-z)
c              
c       Laplace FMM in R^2: evaluate all pairwise particle
c       interactions (ignoring self-interaction). 
c       We use log(r) for the Green's function.
c       Self-interactions are not included.
c   
c lfmm2d: charge and dipstr are complex valued, x \in R^2
c
c \phi(x_i) = \sum_{j\ne i}   charge_j \log |x_i-x_j|  
c                    + dipstr_j (dipvec_j \dot \grad_j log|x_i-x_j|)
c
c or 
c
c \phi(x_i) = \sum_{j\ne i}   charge_j \log |x_i-x_j|  
c                   + dipstr_j (dipvec_j \dot (x_i-x_j)) * (-1/|x_i-x_j|^2)
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
        complex *16 charge(*)
        complex *16 dipstr(*)
        real *8 dipvec(2,*)
        complex *16 ima
        complex *16 pot(*)
        complex *16 grad(2,*)
        complex *16 hess(3,*)
c
        real *8 target(2,1)
        complex *16 pottarg(1)
        complex *16 gradtarg(2,1)
        complex *16 hesstarg(3,1)        
c
        data ima/(0.0d0,1.0d0)/
c       
        ntarget=0
        ifpottarg=0
        ifgradtarg=0
        ifhesstarg=0
c
        call lfmm2dparttarg(ier,iprec,nsource,source,
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
        subroutine lfmm2dparttarg(ier,iprec,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     ntarget,target,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg)
        implicit real *8 (a-h,o-z)
c       
c       Laplace FMM in R^2: evaluate all pairwise particle
c       interactions (ignoring self-interaction) 
c       and interactions with targets.
c
c       We use log(r) for the Green's function.
c       Self-interactions are not included.
c   
c       This is primarily a memory management code. 
c       The actual work is carried out in subroutine cfmm2dparttarg
c
c lfmm2d: charge and dipstr are complex valued, x \in R^2
c
c \phi(x_i) = \sum_{j\ne i}   charge_j \log |x_i-x_j|  
c                    + dipstr_j (dipvec_j \dot \grad_j log|x_i-x_j|)
c
c or 
c
c \phi(x_i) = \sum_{j\ne i}   charge_j \log |x_i-x_j|  
c                   + dipstr_j (dipvec_j \dot (x_i-x_j)) * (-1/|x_i-x_j|^2)
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
c       nsource: integer:  number of sources
c       source: real *8 (2,nsource):  source locations
c       ifcharge:  charge computation flag
c                  ifcharge = 1   =>  include charge contribution
c                                     otherwise do not
c       charge: complex *16 (nsource): charge strengths
c       ifdipole:  dipole computation flag
c                  ifdipole = 1   =>  include dipole contribution
c                                     otherwise do not
c       dipstr: complex *16 (nsource): dipole strengths
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
c       ihesstarg:  target hessian flag 
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
c       pot: complex *16 (nsource): potential at source locations
c       grad: complex *16 (2,nsource): gradient  at source locations
c       hess: complex *16 (3,nsource): hessian at source locations
c       pottarg: complex *16 (ntarget): potential at target locations 
c       gradtarg: complex *16 (2,ntarget): gradient  at target locations 
c       hesstarg: complex *16 (3,ntarget): hessian at target locations
c
c
c
        real *8 source(2,*)
        complex *16 charge(*)
        complex *16 dipstr(*)
        real *8 dipvec(2,*)
        complex *16 ima
        complex *16 pot(*)
        complex *16 grad(2,*)
        complex *16 hess(3,*)
c
        real *8 target(2,*)
        complex *16 pottarg(*)
        complex *16 gradtarg(2,*)        
        complex *16 hesstarg(3,*)
c
        complex *16, allocatable :: charge1(:)
        complex *16, allocatable :: dipstr1(:)
        complex *16, allocatable :: pot1(:)
        complex *16, allocatable :: grad1(:)
        complex *16, allocatable :: hess1(:)
        complex *16, allocatable :: pottarg1(:)
        complex *16, allocatable :: gradtarg1(:)
        complex *16, allocatable :: hesstarg1(:)
c
        data ima/(0.0d0,1.0d0)/
c
c
        if( ifcharge .eq. 1 ) then
        allocate(charge1(nsource))
        else
        allocate(charge1(1))
        endif
c
        if( ifdipole .eq. 1 ) then
        allocate(dipstr1(nsource))
        else
        allocate(dipstr1(1))
        endif
c
        if( ifpot .eq. 1 ) then
        allocate(pot1(nsource))
        else
        allocate(pot1(1))
        endif
c
        if( ifgrad .eq. 1 ) then
        allocate(grad1(nsource))
        else
        allocate(grad1(1))
        endif
c
        if( ifhess .eq. 1 ) then
        allocate(hess1(nsource))
        else
        allocate(hess1(1))
        endif
c
        if( ifpottarg .eq. 1 ) then
        allocate(pottarg1(ntarget))
        else
        allocate(pottarg1(1))
        endif
c
        if( ifgradtarg .eq. 1 ) then
        allocate(gradtarg1(ntarget))
        else
        allocate(gradtarg1(1))
        endif
c
        if( ifhesstarg .eq. 1 ) then
        allocate(hesstarg1(ntarget))
        else
        allocate(hesstarg1(1))
        endif
c
        if( ifcharge .eq. 1 ) then
        do i=1,nsource
        charge1(i)=dble(charge(i))
        enddo
        endif
c
        if( ifdipole .eq. 1 ) then
        do i=1,nsource
        dipstr1(i)=-dble(dipstr(i))*(dipvec(1,i)+ima*dipvec(2,i))
        enddo
        endif
c
        call cfmm2dparttarg(ier,iprec,
     $     nsource,source,
     $     ifcharge,charge1,ifdipole,dipstr1,
     $     ifpot,pot1,ifgrad,grad1,ifhess,hess1,
     $     ntarget,target,ifpottarg,pottarg1,ifgradtarg,gradtarg1,
     $     ifhesstarg,hesstarg1)
c
        do i=1,nsource
        if( ifpot .eq. 1 ) then
        pot(i)=dble(pot1(i))
        endif
        if( ifgrad .eq. 1 ) then
        grad(1,i)=dble(grad1(i))
        grad(2,i)=-imag(grad1(i))
        endif
        if( ifhess .eq. 1 ) then
        hess(1,i)=dble(hess1(i))
        hess(2,i)=-imag(hess1(i))
        hess(3,i)=-dble(hess1(i))
        endif
        enddo
c
        do i=1,ntarget
        if( ifpottarg .eq. 1 ) then
        pottarg(i)=dble(pottarg1(i))
        endif
        if( ifgradtarg .eq. 1 ) then
        gradtarg(1,i)=dble(gradtarg1(i))
        gradtarg(2,i)=-imag(gradtarg1(i))
        endif
        if( ifhesstarg .eq. 1 ) then
        hesstarg(1,i)=dble(hesstarg1(i))
        hesstarg(2,i)=-imag(hesstarg1(i))
        hesstarg(3,i)=-dble(hesstarg1(i))
        endif
        enddo
c
c
        if( ifcharge .eq. 1 ) then
        do i=1,nsource
        charge1(i)=imag(charge(i))
        enddo
        endif
c
        if( ifdipole .eq. 1 ) then
        do i=1,nsource
        dipstr1(i)=-imag(dipstr(i))*(dipvec(1,i)+ima*dipvec(2,i))
        enddo
        endif
c
        call cfmm2dparttarg(ier,iprec,
     $     nsource,source,
     $     ifcharge,charge1,ifdipole,dipstr1,
     $     ifpot,pot1,ifgrad,grad1,ifhess,hess1,
     $     ntarget,target,ifpottarg,pottarg1,ifgradtarg,gradtarg1,
     $     ifhesstarg,hesstarg1)
c
        do i=1,nsource
        if( ifpot .eq. 1 ) then
        pot(i)=pot(i)+ima*dble(pot1(i))
        endif
        if( ifgrad .eq. 1 ) then
        grad(1,i)=grad(1,i)+ima*dble(grad1(i))
        grad(2,i)=grad(2,i)-ima*imag(grad1(i))
        endif
        if( ifhess .eq. 1 ) then
        hess(1,i)=hess(1,i)+ima*dble(hess1(i))
        hess(2,i)=hess(2,i)-ima*imag(hess1(i))
        hess(3,i)=hess(3,i)-ima*dble(hess1(i))
        endif
        enddo
c
        do i=1,ntarget
        if( ifpottarg .eq. 1 ) then
        pottarg(i)=pottarg(i)+ima*dble(pottarg1(i))
        endif
        if( ifgradtarg .eq. 1 ) then
        gradtarg(1,i)=gradtarg(1,i)+ima*dble(gradtarg1(i))
        gradtarg(2,i)=gradtarg(2,i)-ima*imag(gradtarg1(i))
        endif
        if( ifhesstarg .eq. 1 ) then
        hesstarg(1,i)=hesstarg(1,i)+ima*dble(hesstarg1(i))
        hesstarg(2,i)=hesstarg(2,i)-ima*imag(hesstarg1(i))
        hesstarg(3,i)=hesstarg(3,i)-ima*dble(hesstarg1(i))
        endif
        enddo
c
        return
        end
c
c
c
c
c
        subroutine l2dpartdirect(nsource,
     $     source,ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     ntarget,target,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg)
        implicit real *8 (a-h,o-z)
c
c       Laplace interactions in R^2: evaluate all pairwise particle
c       interactions (ignoring self-interaction) 
c       and interactions with targets via direct O(N^2) algorithm.
c
c       We use log(r) for the Green's function.
c       Self-interactions are not included.
c   
c l2d: charge and dipstr are complex valued, x \in R^2
c
c \phi(x_i) = \sum_{j\ne i}   charge_j \log |x_i-x_j|  
c                    + dipstr_j (dipvec_j \dot \grad_j log|x_i-x_j|)
c
c or 
c
c \phi(x_i) = \sum_{j\ne i}   charge_j \log |x_i-x_j|  
c                   + dipstr_j (dipvec_j \dot (x_i-x_j)) * (-1/|x_i-x_j|^2)
c
c       INPUT PARAMETERS:
c
c       nsource: integer:  number of sources
c       source: real *8 (2,nsource):  source locations
c       ifcharge:  charge computation flag
c                  ifcharge = 1   =>  include charge contribution
c                                     otherwise do not
c       charge: complex *16 (nsource): charge strengths
c       ifdipole:  dipole computation flag
c                  ifdipole = 1   =>  include dipole contribution
c                                     otherwise do not
c       dipstr: complex *16 (nsource): dipole strengths
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
c       ihesstarg:  target hessian flag 
c                   (1=compute hessian, otherwise no)
c
c       OUTPUT PARAMETERS:
c
c       pot: complex *16 (nsource): potential at source locations
c       grad: complex *16 (2,nsource): gradient  at source locations
c       hess: complex *16 (3,nsource): hessian at source locations
c       pottarg: complex *16 (ntarget): potential at target locations 
c       gradtarg: complex *16 (2,ntarget): gradient  at target locations 
c       hesstarg: complex *16 (3,ntarget): hessian at target locations
c
c
c
        real *8 source(2,*),dipvec(2,*)
        complex *16 charge(*),dipstr(*)
        real *8 target(2,*)
c
        complex *16 pot(*),grad(2,*),hess(3,*)
        complex *16 pottarg(*),gradtarg(2,*),hesstarg(3,*)
        complex *16 ptemp,gtemp(2),htemp(3)
c
c
        do i=1,nsource
        if( ifpot .eq. 1) pot(i)=0
        if( ifgrad .eq. 1) then
           grad(1,i)=0
           grad(2,i)=0
        endif
        if( ifhess .eq. 1) then
           hess(1,i)=0
           hess(2,i)=0
           hess(3,i)=0
        endif
        enddo
c       
        do i=1,ntarget
        if( ifpottarg .eq. 1) pottarg(i)=0
        if( ifgradtarg .eq. 1) then
           gradtarg(1,i)=0
           gradtarg(2,i)=0
        endif
        if( ifhesstarg .eq. 1) then
           hesstarg(1,i)=0
           hesstarg(2,i)=0
           hesstarg(3,i)=0
        endif
        enddo
c
        if( ifpot .eq. 1 .or. ifgrad .eq. 1 .or. ifhess .eq. 1) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,ptemp,gtemp,htemp) 
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do j=1,nsource
        do i=1,nsource
            if (i .eq. j) cycle
            call rcpotgrad2d_sdp(source(1,i),
     $         ifcharge,charge(i),ifdipole,dipstr(i),dipvec(1,i),
     1         source(1,j),ifpot,ptemp,ifgrad,gtemp,ifhess,htemp)
            if (ifpot .eq. 1) pot(j)=pot(j)+ptemp
            if (ifgrad .eq. 1) then
               grad(1,j)=grad(1,j)+gtemp(1)
               grad(2,j)=grad(2,j)+gtemp(2)
            endif
            if (ifhess .eq. 1) then
               hess(1,j)=hess(1,j)+htemp(1)
               hess(2,j)=hess(2,j)+htemp(2)
               hess(3,j)=hess(3,j)+htemp(3)
            endif
        enddo
        enddo
C$OMP END PARALLEL DO
        endif
c
        if( ifpottarg .eq. 1 .or. ifgradtarg .eq. 1 
     $      .or. ifhesstarg .eq. 1) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,ptemp,gtemp,htemp) 
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do j=1,ntarget
        do i=1,nsource
            call rcpotgrad2d_sdp(source(1,i),
     $         ifcharge,charge(i),ifdipole,dipstr(i),dipvec(1,i),
     1         target(1,j),ifpottarg,ptemp,ifgradtarg,gtemp,
     $         ifhesstarg,htemp)
            if (ifpottarg .eq. 1) pottarg(j)=pottarg(j)+ptemp
            if (ifgradtarg .eq. 1) then
               gradtarg(1,j)=gradtarg(1,j)+gtemp(1)
               gradtarg(2,j)=gradtarg(2,j)+gtemp(2)
            endif
            if (ifhesstarg .eq. 1) then
               hesstarg(1,j)=hesstarg(1,j)+htemp(1)
               hesstarg(2,j)=hesstarg(2,j)+htemp(2)
               hesstarg(3,j)=hesstarg(3,j)+htemp(3)
            endif
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
