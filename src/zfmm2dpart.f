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
c     subroutines for evaluating complex-valued Cauchy sums and
c     gradients due to point dipoles.  (FORTRAN 90 VERSION)
c
c zfmm2d: dipstr are complex valued, z are complex numbers
c
c \phi(z_i) = \sum_{j\ne i} dipstr_j / (z_i-z_j)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     zfmm2dpart - Cauchy FMM in R^2: evaluate all pairwise particle
c         interactions (ignoring self-interaction)
c
c     zfmm2dpartself - Cauchy FMM in R^2: evaluate all pairwise particle
c         interactions (ignoring self-interaction)
c
c     zfmm2dparttarg - Cauchy FMM in R^2: evaluate all pairwise
c         particle interactions (ignoring self-interaction) +
c         interactions with targets
c
c     z2dpartdirect - Cauchy sums in R^2:  evaluate all
c         pairwise particle interactions (ignoring self-interaction) +
c         interactions with targets via direct O(N^2) algorithm
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        this is the end of the debugging code and the beginning 
c        of the Cauchy FMM in R^2
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine zfmm2dpart(ier,iprec,nsource,source,dipstr,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
        implicit real *8 (a-h,o-z)
c              
c       Cauchy FMM in R^2: evaluate all pairwise particle
c       interactions (ignoring self-interaction). 
c
c       Self-interactions are not included.
c   
c zfmm2d: dipstr are complex valued, z are complex numbers
c
c \phi(z_i) = \sum_{j\ne i} dipstr_j / (z_i-z_j)
c
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
        complex *16 charge0,dipstr(*)
        complex *16 ima
        complex *16 pot(*)
        complex *16 grad(*)
        complex *16 hess(*)
c
        real *8 target(2,1)
        complex *16 pottarg(1)
        complex *16 gradtarg(1)
        complex *16 hesstarg(1)        
c
        data ima/(0.0d0,1.0d0)/
c       
        ntarget=0
        ifpottarg=0
        ifgradtarg=0
        ifhesstarg=0
c
        call zfmm2dparttarg(ier,iprec,nsource,source,dipstr,
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
        subroutine zfmm2dpartself(ier,iprec,nsource,source,dipstr,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
        implicit real *8 (a-h,o-z)
c              
c       Cauchy FMM in R^2: evaluate all pairwise particle
c       interactions (ignoring self-interaction). 
c
c       Self-interactions are not included.
c   
c zfmm2d: dipstr are complex valued, z are complex numbers
c
c \phi(z_i) = \sum_{j\ne i} dipstr_j / (z_i-z_j)
c
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
        complex *16 charge0,dipstr(*)
        complex *16 ima
        complex *16 pot(*)
        complex *16 grad(*)
        complex *16 hess(*)
c
        real *8 target(2,1)
        complex *16 pottarg(1)
        complex *16 gradtarg(1)
        complex *16 hesstarg(1)        
c
        data ima/(0.0d0,1.0d0)/
c       
        ntarget=0
        ifpottarg=0
        ifgradtarg=0
        ifhesstarg=0
c
        call zfmm2dparttarg(ier,iprec,nsource,source,dipstr,
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
        subroutine zfmm2dparttarg(ier,iprec,nsource,source,dipstr,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     ntarget,target,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg)
        implicit real *8 (a-h,o-z)
c       
c       Cauchy FMM in R^2: evaluate all pairwise particle
c       interactions (ignoring self-interaction) 
c       and interactions with targets.
c
c       Self-interactions are not included.
c   
c       This is primarily a memory management code. 
c       The actual work is carried out in subroutine cfmm2dparttarg
c
c zfmm2d: dipstr are complex valued, z are complex numbers
c
c \phi(z_i) = \sum_{j\ne i} dipstr_j / (z_i-z_j)
c
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
c       dipstr: complex *16 (nsource): dipole strengths
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
        complex *16 charge0,dipstr(*)
        complex *16 ima
        complex *16 pot(*)
        complex *16 grad(*)
        complex *16 hess(*)
c
        real *8 target(2,*)
        complex *16 pottarg(*)
        complex *16 gradtarg(*)        
        complex *16 hesstarg(*)
c
        complex *16, allocatable :: charge1(:)
c
        data ima/(0.0d0,1.0d0)/
c
        ifcharge=0
        ifdipole=1
c
        call cfmm2dparttarg(ier,iprec,
     $     nsource,source,
     $     ifcharge,charge0,ifdipole,dipstr,
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
        subroutine z2dpartdirect(nsource,source,dipstr,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     ntarget,target,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg)
        implicit real *8 (a-h,o-z)
c
c       Cauchy sums in R^2: evaluate all pairwise particle
c       interactions (ignoring self-interaction) 
c       and interactions with targets via direct O(N^2) algorithm.
c
c       Self-interactions are not included.
c   
c z2d: dipstr are complex valued, z are complex numbers.
c
c \phi(z_i) = \sum_{j\ne i} dipstr_j / (z_i-z_j)
c
c        In this routine, we define the gradient as the first
c        derivative with respect to z, and the hessian as the second
c        derivative with respect to z.
c
c \gradient \phi(z_i) = \frac{\partial \phi(z_i)}{\partial z}
c \hessian  \phi(z_i) = \frac{\partial^2 \phi(z_i)}{\partial z^2}
c
c       INPUT PARAMETERS:
c
c       nsource: integer:  number of sources
c       source: real *8 (2,nsource):  source locations
c       dipstr: complex *16 (nsource): dipole strengths
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
c       grad: complex *16 (nsource): gradient  at source locations
c       hess: complex *16 (nsource): hessian at source locations
c       pottarg: complex *16 (ntarget): potential at target locations 
c       gradtarg: complex *16 (ntarget): gradient  at target locations 
c       hesstarg: complex *16 (ntarget): hessian at target locations
c
c
c
        real *8 source(2,*)
        complex *16 charge0,dipstr(*)
        real *8 target(2,*)
c
        complex *16 pot(*),grad(*),hess(*)
        complex *16 pottarg(*),gradtarg(*),hesstarg(*)
        complex *16 ptemp,gtemp,htemp
c
        ifcharge=0
        ifdipole=1
c
        charge0=0
c
        do i=1,nsource
        if( ifpot .eq. 1) pot(i)=0
        if( ifgrad .eq. 1) then
           grad(i)=0
           grad(i)=0
        endif
        if( ifhess .eq. 1) then
           hess(i)=0
           hess(i)=0
           hess(i)=0
        endif
        enddo
c       
        do i=1,ntarget
        if( ifpottarg .eq. 1) pottarg(i)=0
        if( ifgradtarg .eq. 1) then
           gradtarg(i)=0
           gradtarg(i)=0
        endif
        if( ifhesstarg .eq. 1) then
           hesstarg(i)=0
           hesstarg(i)=0
           hesstarg(i)=0
        endif
        enddo
c
        if( ifpot .eq. 1 .or. ifgrad .eq. 1 .or. ifhess .eq. 1) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,ptemp,gtemp,htemp) 
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(4) 
        do 6160 j=1,nsource
        do 6150 i=1,nsource
            if (i .eq. j) goto 6150
            call cpotgrad2d_sdp(source(1,i),
     $         ifcharge,charge0,ifdipole,dipstr(i),
     1         source(1,j),ifpot,ptemp,ifgrad,gtemp,ifhess,htemp)
            if (ifpot .eq. 1) pot(j)=pot(j)+ptemp
            if (ifgrad .eq. 1) then
               grad(j)=grad(j)+gtemp
            endif
            if (ifhess .eq. 1) then
               hess(j)=hess(j)+htemp
            endif
 6150   continue
 6160   continue
C$OMP END PARALLEL DO
        endif
c
        if( ifpottarg .eq. 1 .or. ifgradtarg .eq. 1 
     $      .or. ifhesstarg .eq. 1) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,ptemp,gtemp,htemp) 
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(4) 
        do j=1,ntarget
        do i=1,nsource
            call cpotgrad2d_sdp(source(1,i),
     $         ifcharge,charge0,ifdipole,dipstr(i),
     1         target(1,j),ifpottarg,ptemp,ifgradtarg,gtemp,
     $         ifhesstarg,htemp)
            if (ifpottarg .eq. 1) pottarg(j)=pottarg(j)+ptemp
            if (ifgradtarg .eq. 1) then
               gradtarg(j)=gradtarg(j)+gtemp
            endif
            if (ifhesstarg .eq. 1) then
               hesstarg(j)=hesstarg(j)+htemp
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
