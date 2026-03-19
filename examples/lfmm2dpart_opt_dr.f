c
c     testing code for FMM - tests charges and dipoles against
c     O(N^2) direct method 
c
c
        implicit real *8 (a-h,o-z)
        real *8, allocatable :: source(:,:)
        complex *16, allocatable :: charge(:,:)
        complex *16, allocatable :: dipstr(:,:)
        real *8, allocatable :: dipvec(:,:)
        complex *16, allocatable :: pot(:,:)
        complex *16, allocatable :: grad(:,:,:)
        complex *16, allocatable :: hess(:,:,:)
c
        complex *16, allocatable :: pot2(:,:)
        complex *16, allocatable :: grad2(:,:,:)
        complex *16, allocatable :: hess2(:,:,:)
c
        real *8, allocatable :: target(:,:)
        complex *16, allocatable :: pottarg(:,:)
        complex *16, allocatable :: gradtarg(:,:,:)
        complex *16, allocatable :: hesstarg(:,:,:)
c
        complex *16 ptemp,gtemp(2),htemp(3)
c
c       Tree reuse variables
        integer tree_handle
        real*8 ttree, teval
c
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c
c
        done=1
        pi=4*atan(done)
c
c
c       SET ALL PARAMETERS
c
        call prini(6,13)
ccc        call prini(0,13)
c
        print *, 'ENTER n, nd'
        read *, nsource, nd
c
        allocate(source(2,nsource))
        allocate(charge(nd,nsource))
        allocate(dipstr(nd,nsource))
        allocate(dipvec(2,nsource))
        allocate(pot(nd,nsource))
        allocate(grad(nd,2,nsource))
        allocate(hess(nd,3,nsource))
c
        allocate(pot2(nd,nsource))
        allocate(grad2(nd,2,nsource))
        allocate(hess2(nd,3,nsource))
c
        allocate(target(2,nsource))
        allocate(pottarg(nd,nsource))
        allocate(gradtarg(nd,2,nsource))
        allocate(hesstarg(nd,3,nsource))
c
        call prinf('nsource=*',nsource,1)
        call prinf('nd=*',nd,1)
c
c
c
c       ... construct randomly located charge distribution on a unit circle
c 
        done=1
        pi=4*atan(done)
c
        d=hkrand(0)
        do i=1,nsource
c
        phi=hkrand(0)*2*pi
        source(1,i)=cos(phi)
        source(2,i)=sin(phi)
ccc        source(1,i)=hkrand(0)
ccc        source(2,i)=hkrand(0)
        enddo
c
c
c       ... construct target distribution on a unit circle
c
        ntarget=nsource*1
        do i=1,ntarget
        phi=hkrand(0)*2*pi
        target(1,i)=cos(phi) + 3
        target(2,i)=sin(phi)
ccc        target(1,i)=hkrand(0)
ccc        target(2,i)=hkrand(0)
        enddo
c
        call prinf('ntarget=*',ntarget,1)
c       
c
        scale=1d+0
        do i=1,nsource
        source(1,i)=source(1,i)*scale
        source(2,i)=source(2,i)*scale
        enddo
        do i=1,ntarget
        target(1,i)=target(1,i)*scale
        target(2,i)=target(2,i)*scale
        enddo


        iprec=4
c
        call prinf('iprec=*',iprec,1)
c
        ifpot=1
        ifgrad=1
c
        ifcharge=1
        ifdipole=1
c
        ifpottarg=1
        ifgradtarg=1
c
        ifhess=1
        ifhesstarg=1
c
        if( ifpottarg .eq. 0 .and. ifgradtarg .eq. 0
     $     .and. ifhesstarg .eq. 0 ) ntarget = 0
c
        if (ifcharge .eq. 1 ) then
c
        do idim=1,nd
           do i=1,nsource
           charge(idim,i)=idim*10 + i + i*ima/2
           enddo
        enddo
c
        endif
c
        if (ifdipole .eq. 1) then
c
        do idim=1,nd
           do i=1,nsource
           dipstr(idim,i)=idim*10 + i + i*ima/2
           dipvec(1,i)=1
           dipvec(2,i)=2
           enddo
        enddo
c
        endif
c
c
c
c       Build tree (one-time cost)
        t1=second()
C$        t1=omp_get_wtime()
c
        call lfmm2dpart_tree_build(ier,iprec,
     $     nsource,source,ntarget,target,tree_handle)
c
        t2=second()
C$        t2=omp_get_wtime()
        ttree = t2-t1
c
c       Evaluate FMM using pre-built tree
        t1=second()
C$        t1=omp_get_wtime()
c
        call lfmm2dpart_tree_eval(ier,tree_handle,nd,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     target,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg)
c
        t2=second()
C$        t2=omp_get_wtime()
        teval = t2-t1
c
c
        call prinf('nsource=*',nsource,1)
        call prinf('ntarget=*',ntarget,1)
        call prin2('tree build time (sec)=*',ttree,1)
        call prin2('fmm eval time (sec)=*',teval,1)
        call prin2('total time (sec)=*',ttree+teval,1)
ccc        call prin2('after fmm, speed (points/sec)=*',nsource/(t2-t1),1)
        call prin2('after fmm, speed (points+targets/sec)=*',
     $     (nsource+ntarget)/(ttree+teval),1)
c       
c
ccc        m=nsource
        m=min(nsource,12)
c
c
        ifprint=0
        if (ifprint .eq. 1) then
        call prin2('source=*',source,3*nsource)
        endif

        ifprint=1
        if (ifprint .eq. 1) then
        call prin2('after fmm, pot=*',pot,2*m)
        if( ifgrad.eq.1 ) call prin2('after fmm, grad=*',grad,2*2*m)
        if( ifhess.eq.1 ) call prin2('after fmm, hess=*',hess,2*3*m)
        endif
c
c
c
        do idim=1,nd
           do i=1,nsource
           if (ifpot .eq. 1) pot2(idim,i)=0
           if (ifgrad .eq. 1) then
              grad2(idim,1,i)=0
              grad2(idim,2,i)=0
           endif
           if (ifhess .eq. 1) then
              hess2(idim,1,i)=0
              hess2(idim,2,i)=0
              hess2(idim,3,i)=0
           endif
           enddo
        enddo
c
        t1=second()
C$        t1=omp_get_wtime()
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(idim,i,j,ptemp,gtemp,htemp)
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(4)
        do j=1,m
        do i=1,nsource
        if( i .eq. j ) cycle
        do idim=1,nd
        call rcpotgrad2d_sdp(source(1,i),
     $     ifcharge,charge(idim,i),ifdipole,dipstr(idim,i),
     $     dipvec(1,i),
     1     source(1,j),ifpot,ptemp,ifgrad,gtemp,ifhess,htemp)
        if (ifpot .eq. 1) pot2(idim,j)=pot2(idim,j)+ptemp
        if (ifgrad .eq. 1) then
           grad2(idim,1,j)=grad2(idim,1,j)+gtemp(1)
           grad2(idim,2,j)=grad2(idim,2,j)+gtemp(2)
        endif
        if (ifhess .eq. 1) then
           hess2(idim,1,j)=hess2(idim,1,j)+htemp(1)
           hess2(idim,2,j)=hess2(idim,2,j)+htemp(2)
           hess2(idim,3,j)=hess2(idim,3,j)+htemp(3)
        endif
        enddo
        enddo
        enddo
C$OMP END PARALLEL DO
c
        t2=second()
C$        t2=omp_get_wtime()
c
        if (ifprint .eq. 1) then
        if( ifpot.eq.1 ) call prin2('directly, pot=*',pot2,2*m)
        if( ifgrad.eq.1 ) call prin2('directly, grad=*',grad2,2*2*m)
        if( ifhess.eq.1 ) call prin2('directly, hess=*',hess2,2*3*m)
        endif
c
        call prin2('directly, estimated time (sec)=*',
     $     (t2-t1)*dble(nsource)/dble(m),1)
        call prin2('directly, estimated speed (points/sec)=*',
     $     m/(t2-t1),1)
c       
        if (ifpot .eq. 1)  then
        do idim=1,nd
        call l2derror(pot(idim,1),pot2(idim,1),m,aerr,rerr)
        call prinf('  density idim=*',idim,1)
        call prin2('  relative L2 error in potential=*',rerr,1)
        call l2derror2(pot(idim,1),pot2(idim,1),m,aerr,rerr)
        call prin2('  relative L2 error in re(potential)=*',rerr,1)
        enddo
        endif
c
        if (ifgrad .eq. 1) then
ccc        call prin2('after fmm, grad=*',grad2,2*m*2)
ccc        call prin2('directly, grad=*',grad2,2*m*2)
        call l2derror(grad,grad2,2*m,aerr,rerr)
ccc         call prin2('absolute L2 error in gradient=*',aerr,1)
        call prin2('relative L2 error in gradient=*',rerr,1)
        endif
c       
        if (ifhess .eq. 1) then
        call l2derror(hess,hess2,3*m,aerr,rerr)
ccc         call prin2('absolute L2 error in hessian=*',aerr,1)
        call prin2('relative L2 error in hessian=*',rerr,1)
        endif
c
c
        if( ntarget .eq. 0 ) goto 9999
c
        do idim=1,nd
           do i=1,ntarget
           if (ifpottarg .eq. 1) pot2(idim,i)=0
           if (ifgradtarg .eq. 1) then
              grad2(idim,1,i)=0
              grad2(idim,2,i)=0
           endif
           if (ifhesstarg .eq. 1) then
              hess2(idim,1,i)=0
              hess2(idim,2,i)=0
              hess2(idim,3,i)=0
           endif
           enddo
        enddo
c
        t1=second()
C$        t1=omp_get_wtime()
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(idim,i,j,ptemp,gtemp,htemp)
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(4)
        do 8160 j=1,m
        do 8150 i=1,nsource
        do idim=1,nd
        call rcpotgrad2d_sdp(source(1,i),
     $     ifcharge,charge(idim,i),ifdipole,dipstr(idim,i),
     $     dipvec(1,i),
     1     target(1,j),
     $     ifpottarg,ptemp,ifgradtarg,gtemp,ifhesstarg,htemp)
        if (ifpottarg .eq. 1) pot2(idim,j)=pot2(idim,j)+ptemp
        if (ifgradtarg .eq. 1) then
           grad2(idim,1,j)=grad2(idim,1,j)+gtemp(1)
           grad2(idim,2,j)=grad2(idim,2,j)+gtemp(2)
        endif
        if (ifhesstarg .eq. 1) then
           hess2(idim,1,j)=hess2(idim,1,j)+htemp(1)
           hess2(idim,2,j)=hess2(idim,2,j)+htemp(2)
           hess2(idim,3,j)=hess2(idim,3,j)+htemp(3)
        endif
        enddo
 8150   continue
 8160   continue
C$OMP END PARALLEL DO
c
        t2=second()
C$        t2=omp_get_wtime()
c
        if (ifprint .eq. 1) then
        if (ifpottarg .eq. 1) 
     $     call prin2('after fmm, pottarg=*',pottarg,2*m)
        if( ifgradtarg.eq.1 ) 
     $     call prin2('after fmm, gradtarg=*',gradtarg,2*2*m)
        if( ifhesstarg.eq.1 ) 
     $     call prin2('after fmm, hesstarg=*',hesstarg,2*3*m)
        if (ifpottarg .eq. 1) 
     $     call prin2('directly, pottarg=*',pot2,2*m)
        if( ifgradtarg.eq.1 ) 
     $     call prin2('directly, gradtarg=*',grad2,2*2*m)
        if( ifhesstarg.eq.1 ) 
     $     call prin2('directly, hesstarg=*',hess2,2*3*m)
        endif
c
        call prin2('directly, estimated time (sec)=*',
     $     (t2-t1)*dble(ntarget)/dble(m),1)
        call prin2('directly, estimated speed (targets/sec)=*',
     $     m/(t2-t1),1)
c       
        if (ifpottarg .eq. 1) then
        do idim=1,nd
        call l2derror(pottarg(idim,1),pot2(idim,1),m,aerr,rerr)
        call prinf('  target density idim=*',idim,1)
        call prin2('  relative L2 error in target potential=*',rerr,1)
        call l2derror2(pottarg(idim,1),pot2(idim,1),m,aerr,rerr)
        call prin2
     $     ('  relative L2 error in target re(potential)=*',rerr,1)
        enddo
        endif
c
        if (ifgradtarg .eq. 1) then
        call l2derror(gradtarg,grad2,2*m,aerr,rerr)
ccc         call prin2('absolute L2 error in gradient=*',aerr,1)
        call prin2('relative L2 error in target gradient=*',rerr,1)
        endif
c       
        if (ifhesstarg .eq. 1) then
        call l2derror(hesstarg,hess2,3*m,aerr,rerr)
ccc         call prin2('absolute L2 error in hessian=*',aerr,1)
        call prin2('relative L2 error in target hessian=*',rerr,1)
        endif
c
 9999   continue
c
c       Free tree handle
        call lfmm2dpart_tree_free(tree_handle)
c
        stop
        end
c
c
c
c
c
        subroutine l2dmperr(mpole1,mpole2,nterms,d)
        implicit real *8 (a-h,o-z)
c       
        complex *16 mpole1(-nterms:nterms)
        complex *16 mpole2(-nterms:nterms)
c       
        d=0
c       
        do n=-nterms,nterms
        d=d+abs(mpole1(n)-mpole2(n))**2
        enddo
c       
        d=d/(2*nterms+1)
        d=sqrt(d)
c       
        return
        end
c
c
c
c
c
        subroutine l2dmpnorm(mpole,nterms,d)
        implicit real *8 (a-h,o-z)
c
        complex *16 mpole(-nterms:nterms)
c
        d=0
c
        do n=-nterms,nterms
        d=d+abs(mpole(n))**2
        enddo
c
        d=d/(2*nterms+1)
        d=sqrt(d)
c
        return
        end
c
c
c
c
c
        subroutine l2derror(pot1,pot2,n,ae,re)
        implicit real *8 (a-h,o-z)
c
c       evaluate absolute and relative errors
c
        complex *16 pot1(n),pot2(n)
c
        d=0
        a=0
c       
        do i=1,n
        d=d+abs(pot1(i)-pot2(i))**2
        a=a+abs(pot1(i))**2
        enddo
c       
        d=d/n
        d=sqrt(d)
        a=a/n
        a=sqrt(a)
c       
        ae=d
        re=d/a
c       
        return
        end
c
c
c
c
c
        subroutine l2derror2(pot1,pot2,n,ae,re)
        implicit real *8 (a-h,o-z)
c
c       evaluate absolute and relative errors of the real part
c
        complex *16 pot1(n),pot2(n)
c
        d=0
        a=0
c       
        do i=1,n
        d=d+dble(pot1(i)-pot2(i))**2
        a=a+dble(pot1(i))**2
        enddo
c       
        d=d/n
        d=sqrt(d)
        a=a/n
        a=sqrt(a)
c       
        ae=d
        re=d/a
c       
        return
        end
c
c
c
