cc Copyright (C) 2010-2024: Leslie Greengard and Zydrunas Gimbutas
cc SPDX-License-Identifier: BSD-3-Clause
c
c     Scalar (nd=1) wrappers for FMM variant routines.
c
c     Naming convention:
c       *_VARIANT_nd  = nd-vectorized (takes nd parameter)
c       *_VARIANT     = scalar wrapper (nd=1, no nd parameter)
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Generalized Cauchy FMM exponential wrappers
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        subroutine cfmm2dpart_exp(ier,iprec,nsource,source,
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
        call cfmm2dpart_exp_nd(ier,iprec,1,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
c
        return
        end
c
c
c
        subroutine cfmm2dpartself_exp(ier,iprec,nsource,source,
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
        call cfmm2dpartself_exp_nd(ier,iprec,1,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
c
        return
        end
c
c
c
        subroutine cfmm2dparttarg_exp(ier,iprec,nsource,source,
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
        call cfmm2dparttarg_exp_nd(ier,iprec,1,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,
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
c     Laplace FMM (complex) exponential wrappers
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        subroutine lfmm2dpart_exp(ier,iprec,nsource,source,
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
        call lfmm2dpart_exp_nd(ier,iprec,1,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
c
        return
        end
c
c
c
        subroutine lfmm2dpartself_exp(ier,iprec,nsource,source,
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
        call lfmm2dpartself_exp_nd(ier,iprec,1,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
c
        return
        end
c
c
c
        subroutine lfmm2dparttarg_exp(ier,iprec,nsource,source,
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
        call lfmm2dparttarg_exp_nd(ier,iprec,1,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
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
c     Laplace FMM (complex) optimized wrappers
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        subroutine lfmm2dpart_opt(ier,iprec,nsource,source,
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
        call lfmm2dpart_opt_nd(ier,iprec,1,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
c
        return
        end
c
c
c
        subroutine lfmm2dpartself_opt(ier,iprec,nsource,source,
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
        call lfmm2dpartself_opt_nd(ier,iprec,1,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
c
        return
        end
c
c
c
        subroutine lfmm2dparttarg_opt(ier,iprec,nsource,source,
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
        call lfmm2dparttarg_opt_nd(ier,iprec,1,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
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
c     Laplace FMM (real) exponential wrappers
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        subroutine rfmm2dpart_exp(ier,iprec,nsource,source,
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
        call rfmm2dpart_exp_nd(ier,iprec,1,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
c
        return
        end
c
c
c
        subroutine rfmm2dpartself_exp(ier,iprec,nsource,source,
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
        call rfmm2dpartself_exp_nd(ier,iprec,1,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
c
        return
        end
c
c
c
        subroutine rfmm2dparttarg_exp(ier,iprec,nsource,source,
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
        call rfmm2dparttarg_exp_nd(ier,iprec,1,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
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
c     Cauchy FMM exponential wrappers
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        subroutine zfmm2dpart_exp(ier,iprec,nsource,source,dipstr,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
        implicit real *8 (a-h,o-z)
c
        real *8 source(2,nsource)
        complex *16 dipstr(nsource)
        complex *16 pot(nsource)
        complex *16 grad(nsource)
        complex *16 hess(nsource)
c
        call zfmm2dpart_exp_nd(ier,iprec,1,nsource,source,dipstr,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
c
        return
        end
c
c
c
        subroutine zfmm2dpartself_exp(ier,iprec,nsource,source,dipstr,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
        implicit real *8 (a-h,o-z)
c
        real *8 source(2,nsource)
        complex *16 dipstr(nsource)
        complex *16 pot(nsource)
        complex *16 grad(nsource)
        complex *16 hess(nsource)
c
        call zfmm2dpartself_exp_nd(ier,iprec,1,nsource,source,dipstr,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
c
        return
        end
c
c
c
        subroutine zfmm2dparttarg_exp(ier,iprec,nsource,source,dipstr,
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
        call zfmm2dparttarg_exp_nd(ier,iprec,1,nsource,source,dipstr,
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
c     Helmholtz FMM standard wrappers
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        subroutine hfmm2dpart_std(ier,iprec,zk,nsource,source,
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
        call hfmm2dpart_std_nd(ier,iprec,zk,1,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
c
        return
        end
c
c
c
        subroutine hfmm2dpartself_std(ier,iprec,zk,nsource,source,
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
        call hfmm2dpartself_std_nd(ier,iprec,zk,1,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
c
        return
        end
c
c
c
        subroutine hfmm2dparttarg_std(ier,iprec,zk,nsource,source,
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
        call hfmm2dparttarg_std_nd(ier,iprec,zk,1,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
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
c     Helmholtz FMM exponential wrappers
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        subroutine hfmm2dpart_exp(ier,iprec,zk,nsource,source,
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
        call hfmm2dpart_exp_nd(ier,iprec,zk,1,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
c
        return
        end
c
c
c
        subroutine hfmm2dpartself_exp(ier,iprec,zk,nsource,source,
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
        call hfmm2dpartself_exp_nd(ier,iprec,zk,1,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
c
        return
        end
c
c
c
        subroutine hfmm2dparttarg_exp(ier,iprec,zk,nsource,source,
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
        call hfmm2dparttarg_exp_nd(ier,iprec,zk,1,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
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
c     Helmholtz FMM far-field wrappers
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        subroutine hfmm2dpart_far(ier,iprec,zk,nsource,source,
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
        call hfmm2dpart_far_nd(ier,iprec,zk,1,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
c
        return
        end
c
c
c
        subroutine hfmm2dpartself_far(ier,iprec,zk,nsource,source,
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
        call hfmm2dpartself_far_nd(ier,iprec,zk,1,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
c
        return
        end
c
c
c
        subroutine hfmm2dparttarg_far(ier,iprec,zk,nsource,source,
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
        call hfmm2dparttarg_far_nd(ier,iprec,zk,1,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     ntarget,target,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg)
c
        return
        end
c
c
