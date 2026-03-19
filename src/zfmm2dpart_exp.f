cc SPDX-License-Identifier: BSD-3-Clause
cc Copyright (C) 2010-2011: Leslie Greengard and Zydrunas Gimbutas
cc See LICENSE file for full license text
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
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
c     zfmm2dpart_exp_nd - Cauchy FMM in R^2: evaluate all pairwise particle
c         interactions (ignoring self-interaction)
c
c     zfmm2dpartself_exp_nd - Cauchy FMM in R^2: evaluate all pairwise particle
c         interactions (ignoring self-interaction)
c
c     zfmm2dparttarg_exp_nd - Cauchy FMM in R^2: evaluate all pairwise
c         particle interactions (ignoring self-interaction) +
c         interactions with targets
c
c     z2dpartdirect_exp_nd - Cauchy sums in R^2:  evaluate all
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
        subroutine zfmm2dpart_exp_nd(ier,iprec,nd,
     $     nsource,source,dipstr,
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
        complex *16 charge0,dipstr(nd,*)
        complex *16 ima
        complex *16 pot(nd,*)
        complex *16 grad(nd,*)
        complex *16 hess(nd,*)
c
        real *8 target(2,1)
        complex *16 pottarg(nd,1)
        complex *16 gradtarg(nd,1)
        complex *16 hesstarg(nd,1)
c
        data ima/(0.0d0,1.0d0)/
c
        ntarget=0
        ifpottarg=0
        ifgradtarg=0
        ifhesstarg=0
c
        call zfmm2dparttarg_exp_nd(ier,iprec,nd,
     $     nsource,source,dipstr,
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
        subroutine zfmm2dpartself_exp_nd(ier,iprec,nd,
     $     nsource,source,dipstr,
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
        complex *16 charge0,dipstr(nd,*)
        complex *16 ima
        complex *16 pot(nd,*)
        complex *16 grad(nd,*)
        complex *16 hess(nd,*)
c
        real *8 target(2,1)
        complex *16 pottarg(nd,1)
        complex *16 gradtarg(nd,1)
        complex *16 hesstarg(nd,1)
c
        data ima/(0.0d0,1.0d0)/
c
        ntarget=0
        ifpottarg=0
        ifgradtarg=0
        ifhesstarg=0
c
        call zfmm2dparttarg_exp_nd(ier,iprec,nd,
     $     nsource,source,dipstr,
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
        subroutine zfmm2dparttarg_exp_nd(ier,iprec,nd,
     $     nsource,source,dipstr,
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
c       The actual work is carried out in subroutine cfmm2dparttarg_exp_nd
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
c       nd: integer: number of densities
c       nsource: integer:  number of sources
c       source: real *8 (2,nsource):  source locations
c       dipstr: complex *16 (nd,nsource): dipole strengths
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
c       pot: complex *16 (nd,nsource): potential at source locations
c       grad: complex *16 (nd,nsource): gradient  at source locations
c       hess: complex *16 (nd,nsource): hessian at source locations
c       pottarg: complex *16 (nd,ntarget): potential at target locations
c       gradtarg: complex *16 (nd,ntarget): gradient  at target locations
c       hesstarg: complex *16 (nd,ntarget): hessian at target locations
c
c
c
        real *8 source(2,*)
        complex *16 charge0,dipstr(nd,*)
        complex *16 ima
        complex *16 pot(nd,*)
        complex *16 grad(nd,*)
        complex *16 hess(nd,*)
c
        real *8 target(2,*)
        complex *16 pottarg(nd,*)
        complex *16 gradtarg(nd,*)
        complex *16 hesstarg(nd,*)
c
        complex *16, allocatable :: charge1(:,:)
c
        data ima/(0.0d0,1.0d0)/
c
        ifcharge=0
        ifdipole=1
c
        call cfmm2dparttarg_exp_nd(ier,iprec,nd,
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
        subroutine z2dpartdirect_exp_nd(nd,nsource,source,dipstr,
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
c       nd: integer: number of densities
c       nsource: integer:  number of sources
c       source: real *8 (2,nsource):  source locations
c       dipstr: complex *16 (nd,nsource): dipole strengths
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
        real *8 source(2,*)
        complex *16 charge0,dipstr(nd,*)
        real *8 target(2,*)
c
        complex *16 pot(nd,*),grad(nd,*),hess(nd,*)
        complex *16 pottarg(nd,*),gradtarg(nd,*),hesstarg(nd,*)
        complex *16 ptemp,gtemp,htemp
c
        ifcharge=0
        ifdipole=1
c
        charge0=0
c
        do idim=1,nd
           do i=1,nsource
           if( ifpot .eq. 1) pot(idim,i)=0
           if( ifgrad .eq. 1) grad(idim,i)=0
           if( ifhess .eq. 1) hess(idim,i)=0
           enddo
        enddo
c
        do idim=1,nd
           do i=1,ntarget
           if( ifpottarg .eq. 1) pottarg(idim,i)=0
           if( ifgradtarg .eq. 1) gradtarg(idim,i)=0
           if( ifhesstarg .eq. 1) hesstarg(idim,i)=0
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
            call cpotgrad2d_sdp(source(1,i),
     $         ifcharge,charge0,ifdipole,dipstr(idim,i),
     1         source(1,j),ifpot,ptemp,ifgrad,gtemp,ifhess,htemp)
            if (ifpot .eq. 1) pot(idim,j)=pot(idim,j)+ptemp
            if (ifgrad .eq. 1) grad(idim,j)=grad(idim,j)+gtemp
            if (ifhess .eq. 1) hess(idim,j)=hess(idim,j)+htemp
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
            call cpotgrad2d_sdp(source(1,i),
     $         ifcharge,charge0,ifdipole,dipstr(idim,i),
     1         target(1,j),ifpottarg,ptemp,ifgradtarg,gtemp,
     $         ifhesstarg,htemp)
            if (ifpottarg .eq. 1) pottarg(idim,j)=pottarg(idim,j)+ptemp
            if (ifgradtarg .eq. 1) gradtarg(idim,j)=gradtarg(idim,j)+gtemp
            if (ifhesstarg .eq. 1) hesstarg(idim,j)=hesstarg(idim,j)+htemp
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
c     Tree reuse API for Cauchy FMM
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine zfmm2dpart_tree_build_exp(ier,iprec,
     $     nsource,source,ntarget,target,tree_handle)
        use fmm2d_tree_module
        implicit real *8 (a-h,o-z)
c
c       Cauchy FMM in R^2: build tree structure for reuse
c
c       This routine builds the FMM tree structure from geometry only.
c       The tree can be reused for multiple evaluations with different
c       nd values or dipole strengths.
c
c       INPUT PARAMETERS:
c
c       iprec:  FMM precision flag
c       nsource: integer:  number of sources
c       source: real *8 (2,nsource):  source locations
c       ntarget: integer:  number of targets
c       target: real *8 (2,ntarget):  target locations
c
c       OUTPUT PARAMETERS:
c
c       ier   =  error return code
c                ier=0     =>  normal execution
c                ier=1     =>  no free tree handles (max 100)
c                ier=4     =>  cannot allocate tree workspace
c
c       tree_handle: integer: handle to tree structure (1-100)
c
c
        real *8 source(2,*)
        real *8 target(2,*)
        integer tree_handle
c
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
        subroutine zfmm2dpart_tree_eval_exp(ier,tree_handle,nd,source,
     $     dipstr,ifpot,pot,ifgrad,grad,ifhess,hess,
     $     target,ifpottarg,pottarg,ifgradtarg,gradtarg,ifhesstarg,
     $     hesstarg)
        use fmm2d_tree_module
        implicit real *8 (a-h,o-z)
c
c       Cauchy FMM in R^2: evaluate using pre-built tree
c
c       This routine evaluates the FMM using a pre-built tree structure.
c       Source and target locations must match those used in tree_build.
c
c       INPUT PARAMETERS:
c
c       tree_handle: integer: handle to tree structure (from tree_build)
c       nd: integer: number of densities
c       source: real *8 (2,nsource):  source locations (must match tree)
c       dipstr: complex *16 (nd,nsource): dipole strengths
c       ifpot:  potential flag (1=compute potential, otherwise no)
c       ifgrad:  gradient flag (1=compute gradient, otherwise no)
c       ifhess:  hessian flag (1=compute hessian, otherwise no)
c       target: real *8 (2,ntarget):  target locations (must match tree)
c       ifpottarg:  target potential flag
c       ifgradtarg:  target gradient flag
c       ifhesstarg:  target hessian flag
c
c       OUTPUT PARAMETERS:
c
c       ier   =  error return code
c                ier=0     =>  normal execution
c                ier=1     =>  invalid tree handle
c                ier=2     =>  tree not allocated
c                ier=8     =>  cannot allocate bulk FMM  workspace
c                ier=16    =>  cannot allocate mpole expansion workspace
c
c       pot: complex *16 (nd,nsource): potential at source locations
c       grad: complex *16 (nd,nsource): gradient  at source locations
c       hess: complex *16 (nd,nsource): hessian at source locations
c       pottarg: complex *16 (nd,ntarget): potential at target locations
c       gradtarg: complex *16 (nd,ntarget): gradient  at target locations
c       hesstarg: complex *16 (nd,ntarget): hessian at target locations
c
c
        real *8 source(2,*)
        complex *16 dipstr(nd,*)
        complex *16 ima
        complex *16 pot(nd,*)
        complex *16 grad(nd,*)
        complex *16 hess(nd,*)
c
        real *8 target(2,*)
        complex *16 pottarg(nd,*)
        complex *16 gradtarg(nd,*)
        complex *16 hesstarg(nd,*)
c
        integer tree_handle
c
        complex *16, allocatable :: charge0(:,:)
c
        data ima/(0.0d0,1.0d0)/
c
c       Allocate dummy charge array (not used, but required by cfmm API)
        allocate(charge0(1,1))
c
        ifcharge=0
        ifdipole=1
c
        call cfmm2dpart_tree_eval_exp(ier,tree_handle,nd,source,
     $     ifcharge,charge0,ifdipole,dipstr,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     target,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg)
c
        deallocate(charge0)
c
        return
        end
c
c
c
c
c
        subroutine zfmm2dpart_tree_free_exp(tree_handle)
        use fmm2d_tree_module
        implicit real *8 (a-h,o-z)
c
c       Free tree handle and deallocate tree resources
c
c       INPUT PARAMETERS:
c
c       tree_handle: integer: handle to tree structure
c
c
        integer tree_handle
c
        call cfmm2dpart_tree_free_exp(tree_handle)
c
        return
        end
c
c
c
c
c
