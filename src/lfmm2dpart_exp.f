cc SPDX-License-Identifier: BSD-3-Clause
cc Copyright (C) 2010-2011: Leslie Greengard and Zydrunas Gimbutas
cc See LICENSE file for full license text
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
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
c     lfmm2dpart_exp_nd - Laplace FMM in R^2: evaluate all pairwise particle
c         interactions (ignoring self-interaction)
c
c     lfmm2dpartself_exp_nd - Laplace FMM in R^2: evaluate all pairwise particle
c         interactions (ignoring self-interaction)
c
c     lfmm2dparttarg_exp_nd - Laplace FMM in R^2: evaluate all pairwise
c         particle interactions (ignoring self-interaction) +
c         interactions with targets
c
c     l2dpartdirect_exp_nd - Laplace interactions in R^2:  evaluate all
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
        subroutine lfmm2dpart_exp_nd(ier,iprec,nd,nsource,source,
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
        complex *16 charge(nd,*)
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
        call lfmm2dparttarg_exp_nd(ier,iprec,nd,nsource,source,
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
        subroutine lfmm2dpartself_exp_nd(ier,iprec,nd,nsource,source,
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
        complex *16 charge(nd,*)
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
        call lfmm2dparttarg_exp_nd(ier,iprec,nd,nsource,source,
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
        subroutine lfmm2dparttarg_exp_nd(ier,iprec,nd,nsource,source,
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
c       The actual work is carried out in subroutine cfmm2dparttarg_exp_nd
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
        complex *16 charge(nd,*)
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
        complex *16, allocatable :: charge1(:,:)
        complex *16, allocatable :: dipstr1(:,:)
        complex *16, allocatable :: pot1(:,:)
        complex *16, allocatable :: grad1(:,:)
        complex *16, allocatable :: hess1(:,:)
        complex *16, allocatable :: pottarg1(:,:)
        complex *16, allocatable :: gradtarg1(:,:)
        complex *16, allocatable :: hesstarg1(:,:)
c
        data ima/(0.0d0,1.0d0)/
c
c
        if( ifcharge .eq. 1 ) then
        allocate(charge1(nd,nsource))
        else
        allocate(charge1(1,1))
        endif
c
        if( ifdipole .eq. 1 ) then
        allocate(dipstr1(nd,nsource))
        else
        allocate(dipstr1(1,1))
        endif
c
        if( ifpot .eq. 1 ) then
        allocate(pot1(nd,nsource))
        else
        allocate(pot1(1,1))
        endif
c
        if( ifgrad .eq. 1 ) then
        allocate(grad1(nd,nsource))
        else
        allocate(grad1(1,1))
        endif
c
        if( ifhess .eq. 1 ) then
        allocate(hess1(nd,nsource))
        else
        allocate(hess1(1,1))
        endif
c
        if( ifpottarg .eq. 1 ) then
        allocate(pottarg1(nd,ntarget))
        else
        allocate(pottarg1(1,1))
        endif
c
        if( ifgradtarg .eq. 1 ) then
        allocate(gradtarg1(nd,ntarget))
        else
        allocate(gradtarg1(1,1))
        endif
c
        if( ifhesstarg .eq. 1 ) then
        allocate(hesstarg1(nd,ntarget))
        else
        allocate(hesstarg1(1,1))
        endif
c
        if( ifcharge .eq. 1 ) then
        do i=1,nsource
           do idim=1,nd
           charge1(idim,i)=dble(charge(idim,i))
           enddo
        enddo
        endif
c
        if( ifdipole .eq. 1 ) then
        do i=1,nsource
           do idim=1,nd
           dipstr1(idim,i)=-dble(dipstr(idim,i))
     $        *(dipvec(1,i)+ima*dipvec(2,i))
           enddo
        enddo
        endif
c
        call cfmm2dparttarg_exp_nd(ier,iprec,nd,
     $     nsource,source,
     $     ifcharge,charge1,ifdipole,dipstr1,
     $     ifpot,pot1,ifgrad,grad1,ifhess,hess1,
     $     ntarget,target,ifpottarg,pottarg1,ifgradtarg,gradtarg1,
     $     ifhesstarg,hesstarg1)
c
        do i=1,nsource
           do idim=1,nd
           if( ifpot .eq. 1 ) then
           pot(idim,i)=dble(pot1(idim,i))
           endif
           if( ifgrad .eq. 1 ) then
           grad(idim,1,i)=dble(grad1(idim,i))
           grad(idim,2,i)=-imag(grad1(idim,i))
           endif
           if( ifhess .eq. 1 ) then
           hess(idim,1,i)=dble(hess1(idim,i))
           hess(idim,2,i)=-imag(hess1(idim,i))
           hess(idim,3,i)=-dble(hess1(idim,i))
           endif
           enddo
        enddo
c
        do i=1,ntarget
           do idim=1,nd
           if( ifpottarg .eq. 1 ) then
           pottarg(idim,i)=dble(pottarg1(idim,i))
           endif
           if( ifgradtarg .eq. 1 ) then
           gradtarg(idim,1,i)=dble(gradtarg1(idim,i))
           gradtarg(idim,2,i)=-imag(gradtarg1(idim,i))
           endif
           if( ifhesstarg .eq. 1 ) then
           hesstarg(idim,1,i)=dble(hesstarg1(idim,i))
           hesstarg(idim,2,i)=-imag(hesstarg1(idim,i))
           hesstarg(idim,3,i)=-dble(hesstarg1(idim,i))
           endif
           enddo
        enddo
c
c
        if( ifcharge .eq. 1 ) then
        do i=1,nsource
           do idim=1,nd
           charge1(idim,i)=imag(charge(idim,i))
           enddo
        enddo
        endif
c
        if( ifdipole .eq. 1 ) then
        do i=1,nsource
           do idim=1,nd
           dipstr1(idim,i)=-imag(dipstr(idim,i))
     $        *(dipvec(1,i)+ima*dipvec(2,i))
           enddo
        enddo
        endif
c
        call cfmm2dparttarg_exp_nd(ier,iprec,nd,
     $     nsource,source,
     $     ifcharge,charge1,ifdipole,dipstr1,
     $     ifpot,pot1,ifgrad,grad1,ifhess,hess1,
     $     ntarget,target,ifpottarg,pottarg1,ifgradtarg,gradtarg1,
     $     ifhesstarg,hesstarg1)
c
        do i=1,nsource
           do idim=1,nd
           if( ifpot .eq. 1 ) then
           pot(idim,i)=pot(idim,i)+ima*dble(pot1(idim,i))
           endif
           if( ifgrad .eq. 1 ) then
           grad(idim,1,i)=grad(idim,1,i)+ima*dble(grad1(idim,i))
           grad(idim,2,i)=grad(idim,2,i)-ima*imag(grad1(idim,i))
           endif
           if( ifhess .eq. 1 ) then
           hess(idim,1,i)=hess(idim,1,i)+ima*dble(hess1(idim,i))
           hess(idim,2,i)=hess(idim,2,i)-ima*imag(hess1(idim,i))
           hess(idim,3,i)=hess(idim,3,i)-ima*dble(hess1(idim,i))
           endif
           enddo
        enddo
c
        do i=1,ntarget
           do idim=1,nd
           if( ifpottarg .eq. 1 ) then
           pottarg(idim,i)=pottarg(idim,i)
     $        +ima*dble(pottarg1(idim,i))
           endif
           if( ifgradtarg .eq. 1 ) then
           gradtarg(idim,1,i)=gradtarg(idim,1,i)
     $        +ima*dble(gradtarg1(idim,i))
           gradtarg(idim,2,i)=gradtarg(idim,2,i)
     $        -ima*imag(gradtarg1(idim,i))
           endif
           if( ifhesstarg .eq. 1 ) then
           hesstarg(idim,1,i)=hesstarg(idim,1,i)
     $        +ima*dble(hesstarg1(idim,i))
           hesstarg(idim,2,i)=hesstarg(idim,2,i)
     $        -ima*imag(hesstarg1(idim,i))
           hesstarg(idim,3,i)=hesstarg(idim,3,i)
     $        -ima*dble(hesstarg1(idim,i))
           endif
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
        subroutine l2dpartdirect_exp_nd(nd,nsource,
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
c       ihesstarg:  target hessian flag
c                   (1=compute hessian, otherwise no)
c
c       OUTPUT PARAMETERS:
c
c       pot: complex *16 (nd,nsource): potential at source locations
c       grad: complex *16 (nd,nsource): gradient  at source locations
c       hess: complex *16 (nd,nsource): hessian at source locations
c       pottarg: complex *16 (nd,ntarget): potential at target locations
c       gradtarg: complex *16 (nd,ntarget): gradient  at target locations
c       hesstarg: complex *16 (nd,ntarget): hessian at target locations
c
c
c
        real *8 source(2,*),dipvec(2,*)
        complex *16 charge(nd,*),dipstr(nd,*)
        real *8 target(2,*)
c
        complex *16 pot(nd,*),grad(nd,2,*),hess(nd,3,*)
        complex *16 pottarg(nd,*),gradtarg(nd,2,*),hesstarg(nd,3,*)
        complex *16 ptemp,gtemp(2),htemp(3)
c
c
        do idim=1,nd
           do i=1,nsource
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
        do idim=1,nd
           do i=1,ntarget
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
cccC$OMP$NUM_THREADS(1)
        do j=1,nsource
        do i=1,nsource
            if (i .eq. j) cycle
            do idim=1,nd
            call rcpotgrad2d_sdp(source(1,i),
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
C$OMP END PARALLEL DO
        endif
c
        if( ifpottarg .eq. 1 .or. ifgradtarg .eq. 1
     $      .or. ifhesstarg .eq. 1) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(idim,i,j,ptemp,gtemp,htemp)
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1)
        do j=1,ntarget
        do i=1,nsource
            do idim=1,nd
            call rcpotgrad2d_sdp(source(1,i),
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Tree reuse API for Laplace FMM
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c

        subroutine lfmm2dpart_tree_build_exp(ier,iprec,
     $     nsource,source,ntarget,target,tree_handle)
        use fmm2d_tree_module
        implicit real *8 (a-h,o-z)
c
c       Build geometry-only tree for lfmm2dpart_exp_nd (Complex Laplace FMM).
c       Tree is independent of nd and charge/dipole values.
c
c       INPUT PARAMETERS:
c       iprec:  FMM precision flag
c       nsource: number of sources
c       source: real *8 (2,nsource): source locations
c       ntarget: number of targets
c       target: real *8 (2,ntarget): target locations
c
c       OUTPUT PARAMETERS:
c       ier: error code (0=success)
c       tree_handle: integer handle for tree reuse
c
        integer tree_handle
        real *8 source(2,*)
        real *8 target(2,*)
c
c       Just call cfmm tree build routine
        call cfmm2dpart_tree_build_exp(ier,iprec,
     $     nsource,source,ntarget,target,tree_handle)
c
        return
        end
c
c
c
c
c
        subroutine lfmm2dpart_tree_eval_exp(ier,tree_handle,nd,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     target,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg)
        use fmm2d_tree_module
        implicit real *8 (a-h,o-z)
c
c       Evaluate lfmm2dpart_exp_nd (Complex Laplace FMM) using pre-built tree.
c       This routine processes real and imaginary parts separately by
c       calling cfmm2dpart_tree_eval_exp twice.
c
c       INPUT PARAMETERS:
c       tree_handle: integer handle from lfmm2dpart_tree_build_exp
c       nd: number of densities
c       source: real *8 (2,nsource): source locations
c       ifcharge: charge flag
c       charge: complex *16 (nd,nsource): charge strengths
c       ifdipole: dipole flag
c       dipstr: complex *16 (nd,nsource): dipole strengths
c       dipvec: real *8 (2,nsource): dipole orientations
c       ifpot, ifgrad, ifhess: source output flags
c       target: real *8 (2,ntarget): target locations
c       ifpottarg, ifgradtarg, ifhesstarg: target output flags
c
c       OUTPUT PARAMETERS:
c       ier: error code (0=success)
c       pot: complex *16 (nd,nsource): potential at sources
c       grad: complex *16 (nd,2,nsource): gradient at sources
c       hess: complex *16 (nd,3,nsource): hessian at sources
c       pottarg: complex *16 (nd,ntarget): potential at targets
c       gradtarg: complex *16 (nd,2,ntarget): gradient at targets
c       hesstarg: complex *16 (nd,3,ntarget): hessian at targets
c
        integer tree_handle
        type(fmm2d_tree_data), pointer :: tree
c
        real *8 source(2,*)
        complex *16 charge(nd,*)
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
        complex *16, allocatable :: charge1(:,:)
        complex *16, allocatable :: dipstr1(:,:)
        complex *16, allocatable :: pot1(:,:)
        complex *16, allocatable :: grad1(:,:)
        complex *16, allocatable :: hess1(:,:)
        complex *16, allocatable :: pottarg1(:,:)
        complex *16, allocatable :: gradtarg1(:,:)
        complex *16, allocatable :: hesstarg1(:,:)
c
        data ima/(0.0d0,1.0d0)/
c
c       Get tree info to know nsource and ntarget
        call fmm2d_tree_get(tree_handle, tree, ier)
        nsource = tree%nsource
        ntarget = tree%ntarget
c
c       Allocate temporary arrays
        if( ifcharge .eq. 1 ) then
        allocate(charge1(nd,nsource))
        else
        allocate(charge1(1,1))
        endif
c
        if( ifdipole .eq. 1 ) then
        allocate(dipstr1(nd,nsource))
        else
        allocate(dipstr1(1,1))
        endif
c
        if( ifpot .eq. 1 ) then
        allocate(pot1(nd,nsource))
        else
        allocate(pot1(1,1))
        endif
c
        if( ifgrad .eq. 1 ) then
        allocate(grad1(nd,nsource))
        else
        allocate(grad1(1,1))
        endif
c
        if( ifhess .eq. 1 ) then
        allocate(hess1(nd,nsource))
        else
        allocate(hess1(1,1))
        endif
c
        if( ifpottarg .eq. 1 ) then
        allocate(pottarg1(nd,ntarget))
        else
        allocate(pottarg1(1,1))
        endif
c
        if( ifgradtarg .eq. 1 ) then
        allocate(gradtarg1(nd,ntarget))
        else
        allocate(gradtarg1(1,1))
        endif
c
        if( ifhesstarg .eq. 1 ) then
        allocate(hesstarg1(nd,ntarget))
        else
        allocate(hesstarg1(1,1))
        endif
c
c       Process real part
        if( ifcharge .eq. 1 ) then
        do i=1,nsource
           do idim=1,nd
           charge1(idim,i)=dble(charge(idim,i))
           enddo
        enddo
        endif
c
        if( ifdipole .eq. 1 ) then
        do i=1,nsource
           do idim=1,nd
           dipstr1(idim,i)=-dble(dipstr(idim,i))
     $        *(dipvec(1,i)+ima*dipvec(2,i))
           enddo
        enddo
        endif
c
        call cfmm2dpart_tree_eval_exp(ier,tree_handle,nd,source,
     $     ifcharge,charge1,ifdipole,dipstr1,
     $     ifpot,pot1,ifgrad,grad1,ifhess,hess1,
     $     target,ifpottarg,pottarg1,ifgradtarg,gradtarg1,
     $     ifhesstarg,hesstarg1)
c
        do i=1,nsource
           do idim=1,nd
           if( ifpot .eq. 1 ) then
           pot(idim,i)=dble(pot1(idim,i))
           endif
           if( ifgrad .eq. 1 ) then
           grad(idim,1,i)=dble(grad1(idim,i))
           grad(idim,2,i)=-imag(grad1(idim,i))
           endif
           if( ifhess .eq. 1 ) then
           hess(idim,1,i)=dble(hess1(idim,i))
           hess(idim,2,i)=-imag(hess1(idim,i))
           hess(idim,3,i)=-dble(hess1(idim,i))
           endif
           enddo
        enddo
c
        do i=1,ntarget
           do idim=1,nd
           if( ifpottarg .eq. 1 ) then
           pottarg(idim,i)=dble(pottarg1(idim,i))
           endif
           if( ifgradtarg .eq. 1 ) then
           gradtarg(idim,1,i)=dble(gradtarg1(idim,i))
           gradtarg(idim,2,i)=-imag(gradtarg1(idim,i))
           endif
           if( ifhesstarg .eq. 1 ) then
           hesstarg(idim,1,i)=dble(hesstarg1(idim,i))
           hesstarg(idim,2,i)=-imag(hesstarg1(idim,i))
           hesstarg(idim,3,i)=-dble(hesstarg1(idim,i))
           endif
           enddo
        enddo
c
c       Process imaginary part
        if( ifcharge .eq. 1 ) then
        do i=1,nsource
           do idim=1,nd
           charge1(idim,i)=imag(charge(idim,i))
           enddo
        enddo
        endif
c
        if( ifdipole .eq. 1 ) then
        do i=1,nsource
           do idim=1,nd
           dipstr1(idim,i)=-imag(dipstr(idim,i))
     $        *(dipvec(1,i)+ima*dipvec(2,i))
           enddo
        enddo
        endif
c
        call cfmm2dpart_tree_eval_exp(ier,tree_handle,nd,source,
     $     ifcharge,charge1,ifdipole,dipstr1,
     $     ifpot,pot1,ifgrad,grad1,ifhess,hess1,
     $     target,ifpottarg,pottarg1,ifgradtarg,gradtarg1,
     $     ifhesstarg,hesstarg1)
c
        do i=1,nsource
           do idim=1,nd
           if( ifpot .eq. 1 ) then
           pot(idim,i)=pot(idim,i)+ima*dble(pot1(idim,i))
           endif
           if( ifgrad .eq. 1 ) then
           grad(idim,1,i)=grad(idim,1,i)+ima*dble(grad1(idim,i))
           grad(idim,2,i)=grad(idim,2,i)-ima*imag(grad1(idim,i))
           endif
           if( ifhess .eq. 1 ) then
           hess(idim,1,i)=hess(idim,1,i)+ima*dble(hess1(idim,i))
           hess(idim,2,i)=hess(idim,2,i)-ima*imag(hess1(idim,i))
           hess(idim,3,i)=hess(idim,3,i)-ima*dble(hess1(idim,i))
           endif
           enddo
        enddo
c
        do i=1,ntarget
           do idim=1,nd
           if( ifpottarg .eq. 1 ) then
           pottarg(idim,i)=pottarg(idim,i)
     $        +ima*dble(pottarg1(idim,i))
           endif
           if( ifgradtarg .eq. 1 ) then
           gradtarg(idim,1,i)=gradtarg(idim,1,i)
     $        +ima*dble(gradtarg1(idim,i))
           gradtarg(idim,2,i)=gradtarg(idim,2,i)
     $        -ima*imag(gradtarg1(idim,i))
           endif
           if( ifhesstarg .eq. 1 ) then
           hesstarg(idim,1,i)=hesstarg(idim,1,i)
     $        +ima*dble(hesstarg1(idim,i))
           hesstarg(idim,2,i)=hesstarg(idim,2,i)
     $        -ima*imag(hesstarg1(idim,i))
           hesstarg(idim,3,i)=hesstarg(idim,3,i)
     $        -ima*dble(hesstarg1(idim,i))
           endif
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
        subroutine lfmm2dpart_tree_free_exp(tree_handle)
        use fmm2d_tree_module
        implicit real *8 (a-h,o-z)
c
c       Free tree resources for lfmm2dpart_exp_nd.
c
c       INPUT PARAMETERS:
c       tree_handle: integer handle from lfmm2dpart_tree_build_exp
c
        integer tree_handle
c
c       Just call cfmm tree free routine
        call cfmm2dpart_tree_free_exp(tree_handle)
c
        return
        end
c
c
c
c
c
