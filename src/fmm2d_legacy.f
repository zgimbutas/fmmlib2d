cc Copyright (C) 2010-2024: Leslie Greengard and Zydrunas Gimbutas
cc SPDX-License-Identifier: BSD-3-Clause
c
c     Legacy wrappers: preserve old nd=1 calling conventions
c     by forwarding to new nd-vectorized routines with nd=1.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Generalized Cauchy FMM wrappers
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        subroutine cfmm2dpart(ier,iprec,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
        implicit real *8 (a-h,o-z)
c
        real *8 source(2,nsource)
        complex *16 charge(nsource)
        complex *16 dipstr(nsource)
        complex *16 pot(nsource)
        complex *16 grad(nsource)
        complex *16 hess(nsource)
c
        call cfmm2dpart_nd(ier,iprec,1,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
c
        return
        end
c
c
c
        subroutine cfmm2dpartself(ier,iprec,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
        implicit real *8 (a-h,o-z)
c
        real *8 source(2,nsource)
        complex *16 charge(nsource)
        complex *16 dipstr(nsource)
        complex *16 pot(nsource)
        complex *16 grad(nsource)
        complex *16 hess(nsource)
c
        call cfmm2dpartself_nd(ier,iprec,1,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
c
        return
        end
c
c
c
        subroutine cfmm2dparttarg(ier,iprec,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     ntarget,target,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg)
        implicit real *8 (a-h,o-z)
c
        real *8 source(2,nsource)
        complex *16 charge(nsource)
        complex *16 dipstr(nsource)
        complex *16 pot(nsource)
        complex *16 grad(nsource)
        complex *16 hess(nsource)
        real *8 target(2,ntarget)
        complex *16 pottarg(ntarget)
        complex *16 gradtarg(ntarget)
        complex *16 hesstarg(ntarget)
c
        call cfmm2dparttarg_nd(ier,iprec,1,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     ntarget,target,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg)
c
        return
        end
c
c
c
        subroutine c2dpartdirect(nsource,
     $     source,ifcharge,charge,ifdipole,dipstr,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     ntarget,target,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg)
        implicit real *8 (a-h,o-z)
c
        real *8 source(2,nsource)
        complex *16 charge(nsource)
        complex *16 dipstr(nsource)
        complex *16 pot(nsource)
        complex *16 grad(nsource)
        complex *16 hess(nsource)
        real *8 target(2,ntarget)
        complex *16 pottarg(ntarget)
        complex *16 gradtarg(ntarget)
        complex *16 hesstarg(ntarget)
c
        call c2dpartdirect_nd(1,nsource,
     $     source,ifcharge,charge,ifdipole,dipstr,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     ntarget,target,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg)
c
        return
        end
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Laplace FMM (complex) wrappers
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        subroutine lfmm2dpart(ier,iprec,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
        implicit real *8 (a-h,o-z)
c
        real *8 source(2,nsource)
        complex *16 charge(nsource)
        complex *16 dipstr(nsource)
        real *8 dipvec(2,nsource)
        complex *16 pot(nsource)
        complex *16 grad(nsource)
        complex *16 hess(nsource)
c
        call lfmm2dpart_nd(ier,iprec,1,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
c
        return
        end
c
c
c
        subroutine lfmm2dpartself(ier,iprec,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
        implicit real *8 (a-h,o-z)
c
        real *8 source(2,nsource)
        complex *16 charge(nsource)
        complex *16 dipstr(nsource)
        real *8 dipvec(2,nsource)
        complex *16 pot(nsource)
        complex *16 grad(nsource)
        complex *16 hess(nsource)
c
        call lfmm2dpartself_nd(ier,iprec,1,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
c
        return
        end
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
        real *8 source(2,nsource)
        complex *16 charge(nsource)
        complex *16 dipstr(nsource)
        real *8 dipvec(2,nsource)
        complex *16 pot(nsource)
        complex *16 grad(nsource)
        complex *16 hess(nsource)
        real *8 target(2,ntarget)
        complex *16 pottarg(ntarget)
        complex *16 gradtarg(ntarget)
        complex *16 hesstarg(ntarget)
c
        call lfmm2dparttarg_nd(ier,iprec,1,nsource,source,
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
        subroutine l2dpartdirect(nsource,
     $     source,ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     ntarget,target,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg)
        implicit real *8 (a-h,o-z)
c
        real *8 source(2,nsource)
        complex *16 charge(nsource)
        complex *16 dipstr(nsource)
        real *8 dipvec(2,nsource)
        complex *16 pot(nsource)
        complex *16 grad(nsource)
        complex *16 hess(nsource)
        real *8 target(2,ntarget)
        complex *16 pottarg(ntarget)
        complex *16 gradtarg(ntarget)
        complex *16 hesstarg(ntarget)
c
        call l2dpartdirect_nd(1,nsource,
     $     source,ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     ntarget,target,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg)
c
        return
        end
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Laplace FMM (real) wrappers
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        subroutine rfmm2dpart(ier,iprec,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
        implicit real *8 (a-h,o-z)
c
        real *8 source(2,nsource)
        real *8 charge(nsource)
        real *8 dipstr(nsource)
        real *8 dipvec(2,nsource)
        real *8 pot(nsource)
        real *8 grad(2,nsource)
        real *8 hess(3,nsource)
c
        call rfmm2dpart_nd(ier,iprec,1,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
c
        return
        end
c
c
c
        subroutine rfmm2dpartself(ier,iprec,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
        implicit real *8 (a-h,o-z)
c
        real *8 source(2,nsource)
        real *8 charge(nsource)
        real *8 dipstr(nsource)
        real *8 dipvec(2,nsource)
        real *8 pot(nsource)
        real *8 grad(2,nsource)
        real *8 hess(3,nsource)
c
        call rfmm2dpartself_nd(ier,iprec,1,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
c
        return
        end
c
c
c
        subroutine rfmm2dparttarg(ier,iprec,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     ntarget,target,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg)
        implicit real *8 (a-h,o-z)
c
        real *8 source(2,nsource)
        real *8 charge(nsource)
        real *8 dipstr(nsource)
        real *8 dipvec(2,nsource)
        real *8 pot(nsource)
        real *8 grad(2,nsource)
        real *8 hess(3,nsource)
        real *8 target(2,ntarget)
        real *8 pottarg(ntarget)
        real *8 gradtarg(2,ntarget)
        real *8 hesstarg(3,ntarget)
c
        call rfmm2dparttarg_nd(ier,iprec,1,nsource,source,
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
        subroutine r2dpartdirect(nsource,
     $     source,ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     ntarget,target,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg)
        implicit real *8 (a-h,o-z)
c
        real *8 source(2,nsource)
        real *8 charge(nsource)
        real *8 dipstr(nsource)
        real *8 dipvec(2,nsource)
        real *8 pot(nsource)
        real *8 grad(2,nsource)
        real *8 hess(3,nsource)
        real *8 target(2,ntarget)
        real *8 pottarg(ntarget)
        real *8 gradtarg(2,ntarget)
        real *8 hesstarg(3,ntarget)
c
        call r2dpartdirect_nd(1,nsource,
     $     source,ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     ntarget,target,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg)
c
        return
        end
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Cauchy FMM wrappers
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        subroutine zfmm2dpart(ier,iprec,nsource,source,dipstr,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
        implicit real *8 (a-h,o-z)
c
        real *8 source(2,nsource)
        complex *16 dipstr(nsource)
        complex *16 pot(nsource)
        complex *16 grad(nsource)
        complex *16 hess(nsource)
c
        call zfmm2dpart_nd(ier,iprec,1,nsource,source,dipstr,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
c
        return
        end
c
c
c
        subroutine zfmm2dpartself(ier,iprec,nsource,source,dipstr,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
        implicit real *8 (a-h,o-z)
c
        real *8 source(2,nsource)
        complex *16 dipstr(nsource)
        complex *16 pot(nsource)
        complex *16 grad(nsource)
        complex *16 hess(nsource)
c
        call zfmm2dpartself_nd(ier,iprec,1,nsource,source,dipstr,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
c
        return
        end
c
c
c
        subroutine zfmm2dparttarg(ier,iprec,nsource,source,dipstr,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     ntarget,target,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg)
        implicit real *8 (a-h,o-z)
c
        real *8 source(2,nsource)
        complex *16 dipstr(nsource)
        complex *16 pot(nsource)
        complex *16 grad(nsource)
        complex *16 hess(nsource)
        real *8 target(2,ntarget)
        complex *16 pottarg(ntarget)
        complex *16 gradtarg(ntarget)
        complex *16 hesstarg(ntarget)
c
        call zfmm2dparttarg_nd(ier,iprec,1,nsource,source,dipstr,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     ntarget,target,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg)
c
        return
        end
c
c
c
        subroutine z2dpartdirect(nsource,source,dipstr,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     ntarget,target,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg)
        implicit real *8 (a-h,o-z)
c
        real *8 source(2,nsource)
        complex *16 dipstr(nsource)
        complex *16 pot(nsource)
        complex *16 grad(nsource)
        complex *16 hess(nsource)
        real *8 target(2,ntarget)
        complex *16 pottarg(ntarget)
        complex *16 gradtarg(ntarget)
        complex *16 hesstarg(ntarget)
c
        call z2dpartdirect_nd(1,nsource,source,dipstr,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     ntarget,target,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg)
c
        return
        end
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Helmholtz FMM wrappers
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        subroutine hfmm2dpart(ier,iprec,zk,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
        implicit real *8 (a-h,o-z)
c
        complex *16 zk
        real *8 source(2,nsource)
        complex *16 charge(nsource)
        complex *16 dipstr(nsource)
        real *8 dipvec(2,nsource)
        complex *16 pot(nsource)
        complex *16 grad(2,nsource)
        complex *16 hess(3,nsource)
c
        call hfmm2dpart_nd(ier,iprec,zk,1,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
c
        return
        end
c
c
c
        subroutine hfmm2dpartself(ier,iprec,zk,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
        implicit real *8 (a-h,o-z)
c
        complex *16 zk
        real *8 source(2,nsource)
        complex *16 charge(nsource)
        complex *16 dipstr(nsource)
        real *8 dipvec(2,nsource)
        complex *16 pot(nsource)
        complex *16 grad(2,nsource)
        complex *16 hess(3,nsource)
c
        call hfmm2dpartself_nd(ier,iprec,zk,1,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
c
        return
        end
c
c
c
        subroutine hfmm2dparttarg(ier,iprec,zk,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     ntarget,target,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg)
        implicit real *8 (a-h,o-z)
c
        complex *16 zk
        real *8 source(2,nsource)
        complex *16 charge(nsource)
        complex *16 dipstr(nsource)
        real *8 dipvec(2,nsource)
        complex *16 pot(nsource)
        complex *16 grad(2,nsource)
        complex *16 hess(3,nsource)
        real *8 target(2,ntarget)
        complex *16 pottarg(ntarget)
        complex *16 gradtarg(2,ntarget)
        complex *16 hesstarg(3,ntarget)
c
        call hfmm2dparttarg_nd(ier,iprec,zk,1,nsource,source,
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
        subroutine h2dpartdirect(zk,nsource,
     $     source,ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     ntarget,target,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg)
        implicit real *8 (a-h,o-z)
c
        complex *16 zk
        real *8 source(2,nsource)
        complex *16 charge(nsource)
        complex *16 dipstr(nsource)
        real *8 dipvec(2,nsource)
        complex *16 pot(nsource)
        complex *16 grad(2,nsource)
        complex *16 hess(3,nsource)
        real *8 target(2,ntarget)
        complex *16 pottarg(ntarget)
        complex *16 gradtarg(2,ntarget)
        complex *16 hesstarg(3,ntarget)
c
        call h2dpartdirect_nd(1,zk,nsource,
     $     source,ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     ntarget,target,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg)
c
        return
        end
c
