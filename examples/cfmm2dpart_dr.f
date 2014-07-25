c
c     testing code for FMM - tests charges and dipoles against
c     O(N^2) direct method 
c
c
        implicit real *8 (a-h,o-z)
        real *8 source(2,2 000 000)
        complex *16 charge(2 000 000)
        complex *16 dipstr(2 000 000)
        real *8 dipvec(2,2 000 000)
        complex *16 pot(2 000 000)
        complex *16 grad(2,2 000 000)
        complex *16 hess(3,2 000 000)
c       
        complex *16 pot2(2 000 000)
        complex *16 grad2(2 000 000)
        complex *16 hess2(2 000 000)
c       
        real *8 target(2,2 000 000)
        complex *16 pottarg(2 000 000)
        complex *16 gradtarg(2 000 000)
        complex *16 hesstarg(2 000 000)
c
        complex *16 ptemp,gtemp,htemp
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
        print *, 'ENTER n'
        read *, nsource
c
c
        call prinf('nsource=*',nsource,1)
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
ccc        target(1,i)=hkrand(0) + 3
ccc        target(2,i)=hkrand(0)
c
        enddo
c
        call prinf('ntarget=*',ntarget,1)
c       
c
        scale=1
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
        do i=1,nsource
        charge(i)=i+i*ima/2*0
ccc        charge(i)=1+ima
        enddo
c
        endif
c       
        if (ifdipole .eq. 1) then
c       
        do i=1,nsource
           dipstr(i)=i+i*ima/2
ccc           dipstr(i)=1+ima
        enddo
c
        endif
c
c
c
        t1=second()
C$        t1=omp_get_wtime()
c       
        call cfmm2dparttarg(ier,iprec,
     $     nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     ntarget,target,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg)
c       
        t2=second()
C$        t2=omp_get_wtime()
c       
c       
        call prinf('nsource=*',nsource,1)
        call prinf('ntarget=*',ntarget,1)
        call prin2('after fmm, time (sec)=*',t2-t1,1)
ccc        call prin2('after fmm, speed (points/sec)=*',nsource/(t2-t1),1)
        call prin2('after fmm, speed (points+targets/sec)=*',
     $     (nsource+ntarget)/(t2-t1),1)
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
        if( ifpot.eq.1 ) call prin2('after fmm, pot=*',pot,2*m)
        if( ifgrad.eq.1 ) call prin2('after fmm, grad=*',grad,2*m)
        if( ifhess.eq.1 ) call prin2('after fmm, hess=*',hess,2*m)
        endif
c
c
c
        do i=1,nsource
        if (ifpot .eq. 1) pot2(i)=0
        if (ifgrad .eq. 1) then
           grad2(i)=0
        endif
        if (ifhess .eq. 1) then
           hess2(i)=0
        endif
        enddo
c        
        t1=second()
C$        t1=omp_get_wtime()
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,ptemp,gtemp,htemp) 
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(4) 
        do 7160 j=1,m
        do 7150 i=1,nsource       
        if( i .eq. j ) goto 7150
        call cpotgrad2d_sdp(source(1,i),
     $     ifcharge,charge(i),ifdipole,dipstr(i),
     1     source(1,j),ifpot,ptemp,ifgrad,gtemp,ifhess,htemp)
        if (ifpot .eq. 1) pot2(j)=pot2(j)+ptemp
        if (ifgrad .eq. 1) then
           grad2(j)=grad2(j)+gtemp
        endif
        if (ifhess .eq. 1) then
           hess2(j)=hess2(j)+htemp
        endif
 7150   continue
 7160   continue
C$OMP END PARALLEL DO
c
        t2=second()
C$        t2=omp_get_wtime()
c
        if (ifprint .eq. 1) then
        if( ifpot.eq.1 ) call prin2('directly, pot=*',pot2,2*m)
        if( ifgrad.eq.1 ) call prin2('directly, grad=*',grad2,2*m)
        if( ifhess.eq.1 ) call prin2('directly, hess=*',hess2,2*m)
        endif
c
        call prin2('directly, estimated time (sec)=*',
     $     (t2-t1)*dble(nsource)/dble(m),1)
        call prin2('directly, estimated speed (points/sec)=*',
     $     m/(t2-t1),1)
c       
        if (ifpot .eq. 1)  then
        call l2derror(pot,pot2,m,aerr,rerr)
ccc        call prin2('absolute L2 error in potential=*',aerr,1)
        call prin2('relative L2 error in potential=*',rerr,1)
        call l2derror2(pot,pot2,m,aerr,rerr)
ccc        call prin2('absolute L2 error in re(potential)=*',aerr,1)
        call prin2('relative L2 error in re(potential)=*',rerr,1)
        endif
c
        if (ifgrad .eq. 1) then
        call l2derror(grad,grad2,m,aerr,rerr)
ccc         call prin2('absolute L2 error in gradient=*',aerr,1)
        call prin2('relative L2 error in gradient=*',rerr,1)
        endif
c       
        if (ifhess .eq. 1) then
        call l2derror(hess,hess2,m,aerr,rerr)
ccc         call prin2('absolute L2 error in hessian=*',aerr,1)
        call prin2('relative L2 error in hessian=*',rerr,1)
        endif
c       
c
        if( ntarget .eq. 0 ) stop
c
        do i=1,ntarget
        if (ifpottarg .eq. 1) pot2(i)=0
        if (ifgradtarg .eq. 1) then
           grad2(i)=0
        endif
        if (ifhesstarg .eq. 1) then
           hess2(i)=0
        endif
        enddo
c        
        t1=second()
C$        t1=omp_get_wtime()
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,ptemp,gtemp,htemp) 
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(4) 
        do 8160 j=1,m
        do 8150 i=1,nsource        
        call cpotgrad2d_sdp(source(1,i),
     $     ifcharge,charge(i),ifdipole,dipstr(i),
     1     target(1,j),
     $     ifpottarg,ptemp,ifgradtarg,gtemp,ifhesstarg,htemp)
        if (ifpottarg .eq. 1) pot2(j)=pot2(j)+ptemp
        if (ifgradtarg .eq. 1) then
           grad2(j)=grad2(j)+gtemp
        endif
        if (ifhesstarg .eq. 1) then
           hess2(j)=hess2(j)+htemp
        endif
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
     $     call prin2('after fmm, gradtarg=*',gradtarg,2*m)
        if( ifhesstarg.eq.1 ) 
     $     call prin2('after fmm, hesstarg=*',hesstarg,2*m)
        if (ifpottarg .eq. 1) 
     $     call prin2('directly, pottarg=*',pot2,2*m)
        if( ifgradtarg.eq.1 ) 
     $     call prin2('directly, gradtarg=*',grad2,2*m)
        if( ifhesstarg.eq.1 ) 
     $     call prin2('directly, hesstarg=*',hess2,2*m)
        endif
c
        call prin2('directly, estimated time (sec)=*',
     $     (t2-t1)*dble(ntarget)/dble(m),1)
        call prin2('directly, estimated speed (targets/sec)=*',
     $     m/(t2-t1),1)
c       
        if (ifpottarg .eq. 1) then
        call l2derror(pottarg,pot2,m,aerr,rerr)
ccc        call prin2('absolute L2 error in potential=*',aerr,1)
        call prin2('relative L2 error in target potential=*',rerr,1)
        call l2derror2(pottarg,pot2,m,aerr,rerr)
ccc        call prin2('absolute L2 error in re(potential)=*',aerr,1)
        call prin2
     $     ('relative L2 error in target re(potential)=*',rerr,1)
        endif
c
        if (ifgradtarg .eq. 1) then
        call l2derror(gradtarg,grad2,m,aerr,rerr)
ccc         call prin2('absolute L2 error in gradient=*',aerr,1)
        call prin2('relative L2 error in target gradient=*',rerr,1)
        endif
c       
        if (ifhesstarg .eq. 1) then
        call l2derror(hesstarg,hess2,m,aerr,rerr)
ccc         call prin2('absolute L2 error in hessian=*',aerr,1)
        call prin2('relative L2 error in target hessian=*',rerr,1)
        endif
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
