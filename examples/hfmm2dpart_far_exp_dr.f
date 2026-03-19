c
c     testing code for Helmholtz FMM - tests charges and dipoles against
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
        complex *16 zk
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c
c       Tree reuse variables
        integer tree_handle
        real*8 ttree, teval
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
c       Set Helmholtz wavenumber (moderate wavenumber for testing)
        zk = dcmplx(32*2*pi, 0)
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
        call prin2('zk=*',zk,2)
c
c
c
c       ... construct randomly located charge distribution (circle)
c
        done=1
        pi=4*atan(done)
c
        call random_seed()
        do i=1,nsource
        call random_number(phi)
        phi=phi*2*pi
        source(1,i)=cos(phi)
        source(2,i)=sin(phi)
        enddo
c
c
c       ... construct target distribution (circle, shifted)
c
        ntarget=nsource*1
        do i=1,ntarget
        call random_number(phi)
        phi=phi*2*pi
        target(1,i)=cos(phi) + 3
        target(2,i)=sin(phi)
        enddo
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
        ifhess=0
c
        ifcharge=1
        ifdipole=1
c
        ifpottarg=1
        ifgradtarg=0
        ifhesstarg=0
c
        if( ifpottarg .eq. 0 .and. ifgradtarg .eq. 0
     $     .and. ifhesstarg .eq. 0 ) ntarget = 0
c
        call prinf('ntarget=*',ntarget,1)
c
c
        if (ifcharge .eq. 1 ) then
c
        do idim=1,nd
        do i=1,nsource
        charge(idim,i)=dcmplx(i+idim, 0.0d0)
ccc        charge(idim,i)=dcmplx(1.0d0, 0.0d0)
        enddo
        enddo
c
        endif
c
        if (ifdipole .eq. 1) then
c
        do idim=1,nd
        do i=1,nsource
           dipstr(idim,i)=dcmplx(dble(i+idim), 0.0d0)
ccc           dipstr(idim,i)=i+i*ima/2
ccc           dipstr(idim,i)=1+ima
        enddo
        enddo
c
c       ... set dipole orientation vectors
        do i=1,nsource
        call random_number(phi)
        phi=phi*2*pi
        dipvec(1,i)=cos(phi)
        dipvec(2,i)=sin(phi)
        enddo
c
        endif
c
c
c
c       Call FMM directly (not tree API) to activate far-field
        t1=second()
C$        t1=omp_get_wtime()
c
        ttree = 0
c
        call hfmm2dparttarg(ier,iprec,zk,nd,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     ntarget,target,ifpottarg,pottarg,ifgradtarg,gradtarg,
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
        m=min(nsource,50)
c
c
        ifprint=0
        if (ifprint .eq. 1) then
        call prin2('source=*',source,3*nsource)
        endif

        ifprint=1
        if (ifprint .eq. 1) then
        if( ifpot.eq.1 ) call prin2('after fmm, pot=*',pot,2*nd*m)
        if( ifgrad.eq.1 )
     $     call prin2('after fmm, grad=*',grad,2*nd*2*m)
        if( ifhess.eq.1 )
     $     call prin2('after fmm, hess=*',hess,2*nd*3*m)
        endif
c
c
c
        do i=1,nsource
        do idim=1,nd
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
        do 7160 j=1,m
        do 7150 i=1,nsource
        if( i .eq. j ) goto 7150
        do idim=1,nd
        call hpotgrad2d_sdp(zk,source(1,i),
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
 7150   continue
 7160   continue
C$OMP END PARALLEL DO
c
        t2=second()
C$        t2=omp_get_wtime()
c
        if (ifprint .eq. 1) then
        if( ifpot.eq.1 ) call prin2('directly, pot=*',pot2,2*nd*m)
        if( ifgrad.eq.1 )
     $     call prin2('directly, grad=*',grad2,2*nd*2*m)
        if( ifhess.eq.1 )
     $     call prin2('directly, hess=*',hess2,2*nd*3*m)
        endif
c
        call prin2('directly, estimated time (sec)=*',
     $     (t2-t1)*dble(nsource)/dble(m),1)
        call prin2('directly, estimated speed (points/sec)=*',
     $     m/(t2-t1),1)
c
        if (ifpot .eq. 1)  then
        call l2derror_nd(pot,pot2,nd,m,aerr,rerr)
ccc        call prin2('absolute L2 error in potential=*',aerr,1)
        call prin2('relative L2 error in potential=*',rerr,1)
        write(*,'(a,e12.5)') 'Relative L2 error in potential: ',rerr
        call l2derror2_nd(pot,pot2,nd,m,aerr,rerr)
ccc        call prin2('absolute L2 error in re(potential)=*',aerr,1)
        call prin2('relative L2 error in re(potential)=*',rerr,1)
        write(*,'(a,e12.5)') 'Relative L2 error in re(pot): ',rerr
        endif
c
        if (ifgrad .eq. 1) then
        d=0
        a=0
        do i=1,m
        do idim=1,nd
        d=d+abs(grad(idim,1,i)-grad2(idim,1,i))**2
        d=d+abs(grad(idim,2,i)-grad2(idim,2,i))**2
        a=a+abs(grad(idim,1,i))**2+abs(grad(idim,2,i))**2
        enddo
        enddo
        rerr=sqrt(d)/sqrt(a)
        call prin2('relative L2 error in gradient=*',rerr,1)
        write(*,'(a,e12.5)') 'Relative L2 error in gradient: ',rerr
        endif
c
        if (ifhess .eq. 1) then
        d=0
        a=0
        do i=1,m
        do idim=1,nd
        d=d+abs(hess(idim,1,i)-hess2(idim,1,i))**2
        d=d+abs(hess(idim,2,i)-hess2(idim,2,i))**2
        d=d+abs(hess(idim,3,i)-hess2(idim,3,i))**2
        a=a+abs(hess(idim,1,i))**2+abs(hess(idim,2,i))**2
     $     +abs(hess(idim,3,i))**2
        enddo
        enddo
        rerr=sqrt(d)/sqrt(a)
        call prin2('relative L2 error in hessian=*',rerr,1)
        write(*,'(a,e12.5)') 'Relative L2 error in hessian: ',rerr
        endif
c
c
        if( ntarget .eq. 0 ) stop
c
        do i=1,ntarget
        do idim=1,nd
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
        call hpotgrad2d_sdp(zk,source(1,i),
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
     $     call prin2('after fmm, pottarg=*',pottarg,2*nd*m)
        if( ifgradtarg.eq.1 )
     $     call prin2('after fmm, gradtarg=*',gradtarg,2*nd*2*m)
        if( ifhesstarg.eq.1 )
     $     call prin2('after fmm, hesstarg=*',hesstarg,2*nd*3*m)
        if (ifpottarg .eq. 1)
     $     call prin2('directly, pottarg=*',pot2,2*nd*m)
        if( ifgradtarg.eq.1 )
     $     call prin2('directly, gradtarg=*',grad2,2*nd*2*m)
        if( ifhesstarg.eq.1 )
     $     call prin2('directly, hesstarg=*',hess2,2*nd*3*m)
        endif
c
        call prin2('directly, estimated time (sec)=*',
     $     (t2-t1)*dble(ntarget)/dble(m),1)
        call prin2('directly, estimated speed (targets/sec)=*',
     $     m/(t2-t1),1)
c
        if (ifpottarg .eq. 1) then
        call l2derror_nd(pottarg,pot2,nd,m,aerr,rerr)
ccc        call prin2('absolute L2 error in potential=*',aerr,1)
        call prin2('relative L2 error in target potential=*',rerr,1)
        write(*,'(a,e12.5)') 'Relative L2 error in targ pot: ',rerr
        call l2derror2_nd(pottarg,pot2,nd,m,aerr,rerr)
ccc        call prin2('absolute L2 error in re(potential)=*',aerr,1)
        call prin2
     $     ('relative L2 error in target re(potential)=*',rerr,1)
        write(*,'(a,e12.5)') 'Relative L2 error in targ re(pot): ',rerr
        endif
c
        if (ifgradtarg .eq. 1) then
        d=0
        a=0
        do i=1,m
        do idim=1,nd
        d=d+abs(gradtarg(idim,1,i)-grad2(idim,1,i))**2
        d=d+abs(gradtarg(idim,2,i)-grad2(idim,2,i))**2
        a=a+abs(gradtarg(idim,1,i))**2+abs(gradtarg(idim,2,i))**2
        enddo
        enddo
        rerr=sqrt(d)/sqrt(a)
        call prin2('relative L2 error in target gradient=*',rerr,1)
        write(*,'(a,e12.5)')
     $     'Relative L2 error in targ gradient: ',rerr
        endif
c
        if (ifhesstarg .eq. 1) then
        d=0
        a=0
        do i=1,m
        do idim=1,nd
        d=d+abs(hesstarg(idim,1,i)-hess2(idim,1,i))**2
        d=d+abs(hesstarg(idim,2,i)-hess2(idim,2,i))**2
        d=d+abs(hesstarg(idim,3,i)-hess2(idim,3,i))**2
        a=a+abs(hesstarg(idim,1,i))**2+abs(hesstarg(idim,2,i))**2
     $     +abs(hesstarg(idim,3,i))**2
        enddo
        enddo
        rerr=sqrt(d)/sqrt(a)
        call prin2('relative L2 error in target hessian=*',rerr,1)
        write(*,'(a,e12.5)')
     $     'Relative L2 error in targ hessian: ',rerr
        endif
c
c       Free tree handle
ccc        call hfmm2dpart_tree_free(tree_handle)
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
        subroutine l2derror_nd(pot1,pot2,nd,n,ae,re)
        implicit real *8 (a-h,o-z)
c
c       evaluate absolute and relative errors for nd densities
c
        complex *16 pot1(nd,n),pot2(nd,n)
c
        d=0
        a=0
c
        do i=1,n
        do idim=1,nd
        d=d+abs(pot1(idim,i)-pot2(idim,i))**2
        a=a+abs(pot1(idim,i))**2
        enddo
        enddo
c
        d=d/(n*nd)
        d=sqrt(d)
        a=a/(n*nd)
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
        subroutine l2derror2_nd(pot1,pot2,nd,n,ae,re)
        implicit real *8 (a-h,o-z)
c
c       evaluate absolute and relative errors of the real part for nd densities
c
        complex *16 pot1(nd,n),pot2(nd,n)
c
        d=0
        a=0
c
        do i=1,n
        do idim=1,nd
        d=d+dble(pot1(idim,i)-pot2(idim,i))**2
        a=a+dble(pot1(idim,i))**2
        enddo
        enddo
c
        d=d/(n*nd)
        d=sqrt(d)
        a=a/(n*nd)
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
