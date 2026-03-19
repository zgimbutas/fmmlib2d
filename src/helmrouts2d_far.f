ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        this is the end of the debugging code and the beginning 
c        of the routines for high-frequency Helmholtz FMM in R^2
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
        subroutine hfar2dterms(zk, c1, c2, tol, nterms, ier)
        implicit real *8 (a-h,o-z)
        complex *16 zk
        dimension c1(2),c2(2)
c
c       Determine number of terms in far to local signatures
c       for box of size "size" with Helmholtz parameter zk.
c
c     The method is based on examining the decay of H_n.
c
c     Maximum number of terms is 20000000, which 
c     works for boxes up to 1600000 wavelengths in size     
c
c
	r = sqrt( (c2(1)-c1(1))**2 + (c2(2)-c1(2))**2)

        call hfar2dterms1(zk, r, tol, nterms, ier)
        return
        end


        
        subroutine hfar2dterms1(zk, r, tol, nterms, ier)
        implicit real *8 (a-h,o-z)
c       
c       Determine number of terms in mpole expansions for box of size
c       "size" with Helmholtz parameter zk.
c       
c       The method is based on examining the decay of H_n.
c       
c       Maximum number of terms is 10000000, which 
c       works for boxes up to 1600000 wavelengths in size     
c       
c       
c-----------------------------------------------------------------------------
c       
        complex *16  zk, zksize, z1, z2, z3, ht0,
     1     ht1, ht2, fjder(0:1), ztmp,fhder(0:1)
c
        complex *16, allocatable :: hfun(:)
c        
        done = 1.0d0
        pi = 4*atan(done)
c       
        ier = 0
c        
        zksize = zk*r
c       
c       
c       the code will run out memory if frequency is too small 
c       set frequency to something more reasonable, nterms is 
c       approximately the same for all small frequencies
c       
        if( abs(zksize) .lt. 1d-8 ) zksize = 1d-8
c
        z1 = abs(dble(zksize))
c       
        ntmax = abs(dble(zk*r))/(2*pi)*10
        if( ntmax .lt. 1000 ) ntmax = 1000
        if( ntmax .gt. 20000000 ) ntmax = 20000000
c
        allocate(hfun(0:ntmax))
c
        ifder = 0
        rscale = 1.0d0
        if ( abs(zksize) .lt. 1.0d0) rscale = abs(zksize)
        call h2dall(ntmax,z1,rscale,hfun,ifder,fhder)
c
ccc        call prin2('bsize=*',r,1)
ccc        call prin2('zksize=*',zksize,2)
ccc        call prin2('hfuns=*',hfun(0),ntmax*2-2)
c       
        done = 1
        pi = 4*atan(done)
        
        nterms = 1
        do j = 1, ntmax
        xtemp = abs(hfun(j)) * sqrt(pi/2) * sqrt(zksize)
ccc        write(*,*) j, xtemp
        if(xtemp .gt. 1/tol) then
        nterms = j - 1
        return
        endif
c       
        enddo
c
c        
c       ... computational box is too big, set nterms to 20000001
c       
        ier = 13
        nterms=20000001
c       
        return
        end
c       
c       
c       
c       
c       
        subroutine h2dfardir(ndir,dirs)
        implicit real *8 (a-h,o-z)
        dimension dirs(2,ndir)
        
        done=1.0d0
        pi=4*atan(done)

        do i = 1,ndir
        theta = 2*pi*(i-1)/dble(ndir)
        x = cos(theta)
        y = sin(theta)
        dirs(1,i) = x
        dirs(2,i) = y
        enddo

        return
        end
c       
c       
c       
c       
c       
        subroutine h2dmpfar(mpole,nterms,mpfar,nfar)
        implicit real *8 (a-h,o-z)
        complex *16 mpole(-nterms:nterms)
        complex *16 mpfar(1:nfar)
        complex *16 ima
ccc        complex *16 wsave(100000)
        complex *16, allocatable :: wsave(:)
        data ima/(0.0d0,1.0d0)/

        allocate(wsave(4*nfar+15),stat=ier)
        
        do m = 1,nfar
        mpfar(m) = 0
        enddo

        mpfar(1) = mpole(0)
        do m = 1,nterms
        mpfar(m+1) = mpole(m) * (-1.0*ima)**m
        mpfar(nfar-m+1) = mpole(-m) * (-1.0*ima)**(-m)
        enddo

        call zffti(nfar, wsave)
        call zfftb(nfar, mpfar, wsave)
        do m = 1,nfar
        mpfar(m) = mpfar(m)
        enddo

        return
        end
c
c
c
c
c       
        subroutine h2dfarfar(zk,center1,mpfar1,nfar1,
     $                          center2,mpfar2,nfar2)
        implicit real *8 (a-h,o-z)
        complex *16 zk,trans
        dimension center1(2),center2(2)
        complex *16 mpfar1(1:nfar1),mpfar2(1:nfar2)
        complex *16 ima
ccc        complex *16 wsave1(100000), wsave2(100000)
ccc        complex *16 mptmp(100000),mpint(100000)
        complex *16, allocatable :: wsave1(:), wsave2(:)
        complex *16, allocatable :: mptmp(:),mpint(:)
        data ima/(0.0d0,1.0d0)/

        allocate(wsave1(4*nfar1+15),stat=ier)
        allocate(wsave2(4*nfar2+15),stat=ier)
        allocate(mptmp(nfar1),stat=ier)
        allocate(mpint(nfar2),stat=ier)
        
        do m = 1,nfar1
        mptmp(m) = mpfar1(m)
        enddo

        call zffti(nfar1, wsave1)
        call zfftf(nfar1, mptmp, wsave1)

        do m = 1,nfar2
        mpint(m) = 0
        enddo

        do m = 0,(nfar1-1)/2
        mpint(m+1) = mptmp(m+1)
        enddo

        do m = 1,nfar1/2
        mpint(nfar2-m+1) = mptmp(nfar1-m+1)
        enddo

        call zffti(nfar2, wsave2)
        call zfftb(nfar2, mpint, wsave2)

        dx = center2(1)-center1(1)
        dy = center2(2)-center1(2)

        done=1.0d0
        pi=4*atan(done)

        do m = 1,nfar2
        theta = 2*pi*(m-1)/dble(nfar2)
        dir1 = cos(theta)
        dir2 = sin(theta)
        proj = dx*dir1 + dy*dir2
        trans = exp(ima*zk*proj)
        mpfar2(m) = mpfar2(m) + mpint(m) * trans / dble(nfar1)
        enddo

        return
        end
c
c
c
c
c
        subroutine h2dfarpwl_init(zk,center1,center2,
     $     nterms,nfar,trans)
        implicit real *8 (a-h,o-z)
        complex *16 zk
        dimension center1(2),center2(2),zdiff(2)
        complex *16 ima, z, zmul, zinv, ztemp1, ztemp2
        complex *16, allocatable :: hval(:), hder(:), htemp(:)
        complex *16 trans(*)
        data ima/(0.0d0,1.0d0)/
c
        done=1.0d0
        pi=4*atan(done)
c
        allocate(hval(0:nterms+5), stat=ier)
        if (ier .eq. 1) then
        return
        endif

        allocate(hder(0:nterms+5), stat=ier)
        if (ier .eq. 1) then
        return
        endif

        allocate(htemp(-nterms:nterms), stat=ier)
        if (ier .eq. 1) then
        return
        endif

        zdiff(1)=center2(1)-center1(1)
        zdiff(2)=center2(2)-center1(2)
        call h2cart2polar(zdiff,r,theta)

        theta=theta-pi

        z=zk*r
        ifder=0
        rscale=1.0d0

        call h2dall(nterms+1,z,rscale,hval,ifder,hder)

        htemp(0) = hval(0)
        zmul=exp(-ima*theta)
        zinv=dconjg(zmul)
        ztemp1= zmul
        ztemp2=-zinv
        do j = 1,nterms
        htemp(+j) = ztemp1*hval(j)
        htemp(-j) = ztemp2*hval(j)
        ztemp1= ztemp1*zmul
        ztemp2=-ztemp2*zinv
        enddo

        call h2dmpfar(htemp,nterms,trans,nfar)

        return
        end
c
c
c
c
c       
        subroutine h2dfarpwl(zk,center1,mpfar1,nfar1,
     $                          center2,pwloc2,nloc2,trans)
        implicit real *8 (a-h,o-z)
        complex *16 zk
        dimension center1(2),center2(2)
        complex *16 mpfar1(1:nfar1),pwloc2(1:nloc2),trans(*)
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
    
        do j = 1,nfar1
        pwloc2(j) = pwloc2(j) + mpfar1(j)*trans(j)
        enddo

        return
        end
c
c
c
c
c
        subroutine h2dpwlpwl(zk,center1,pwloc1,nloc1,
     $                          center2,pwloc2,nloc2)
        implicit real *8 (a-h,o-z)
        complex *16 zk,trans
        dimension center1(2),center2(2)
        complex *16 pwloc1(1:nloc1),pwloc2(1:nloc2)
        complex *16 ima
ccc        complex *16 wsave1(100000), wsave2(100000)
ccc        complex *16 mptmp(100000),mpint(100000)
        complex *16, allocatable :: wsave1(:), wsave2(:)
        complex *16, allocatable :: mptmp(:),mpint(:)
        data ima/(0.0d0,1.0d0)/

        allocate(wsave1(4*nloc1+15),stat=ier)
        allocate(wsave2(4*nloc2+15),stat=ier)
        allocate(mptmp(nloc1),stat=ier)
        allocate(mpint(nloc2),stat=ier)
        
        do m = 1,nloc1
        mptmp(m) = pwloc1(m)
        enddo

        call zffti(nloc1, wsave1)
        call zfftf(nloc1, mptmp, wsave1)

        do m = 1,nloc2
        mpint(m) = 0
        enddo

        do m = 0,(nloc2-1)/2
        mpint(m+1) = mptmp(m+1)
        enddo

        do m = 1,nloc2/2
        mpint(nloc2-m+1) = mptmp(nloc1-m+1)
        enddo

        call zffti(nloc2, wsave2)
        call zfftb(nloc2, mpint, wsave2)

        dx = center2(1)-center1(1)
        dy = center2(2)-center1(2)

        done=1.0d0
        pi=4*atan(done)

        do m = 1,nloc2
        theta = 2*pi*(m-1)/dble(nloc2)
        dir1 = cos(theta)
        dir2 = sin(theta)
        proj = dx*dir1 + dy*dir2
        trans = exp(ima*zk*proj)
        pwloc2(m) = pwloc2(m) + mpint(m) * trans / dble(nloc1)
        enddo

        return
        end
c
c
c
c
c       
        subroutine h2dpwlta(mploc,nterms,pwloc,nloc)
        implicit real *8 (a-h,o-z)
        complex *16 mploc(-nterms:nterms)
        complex *16 pwloc(1:nloc)
        complex *16 ima
ccc        complex *16 wsave(100000)
        complex *16, allocatable :: wsave(:)
        data ima/(0.0d0,1.0d0)/

        allocate(wsave(4*nloc+15),stat=ier)
        
        call zffti(nloc, wsave)
        call zfftf(nloc, pwloc, wsave)

        do m = 1,nloc
        pwloc(m) = pwloc(m) / dble(nloc)
        enddo
        
        mploc(0) = mploc(0) + pwloc(1)
        do m = 1,nterms
        mploc(+m) = mploc(+m) + pwloc(m+1) / (-1.0*ima)**m
        mploc(-m) = mploc(-m) + pwloc(nloc-m+1) / (-1.0*ima)**(-m)
        enddo

        return
        end
c
c
c
c
c
        subroutine charge2far(zk,xy,charge,dir,cfar)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates the far field
c       signatures in the direction dir due to the monochromatic
c       charge located at the at the location xy
c
c          Input parameters:
c
c       zk (complex *16)  - the frequency parameter
c       xy (real *8 ) - the location of the charge in R^2
c       dir (real *8 ) - the far field direction vector in R^2 
c       charge (complex *16) - the strength of the charge
c
c          Output parameters:
c
c       cfar (complex*16) - the far field signature of the potential
c
        dimension xy(2),dir(2),far(2)
        complex *16 zk,charge,cfar,ima
        data ima/(0.0d0,1.0d0)/
c
c       ... first, find the point on the unit circle
c        
        d=sqrt(dir(1)**2+dir(2)**2)
        far(1)=dir(1)/d
        far(2)=dir(2)/d
c
c       ... find the far-field signature of the unit charge located at
c       the location xy
c
        d=xy(1)*far(1)+xy(2)*far(2)
        cfar=exp(-ima*zk*d)
c
        cfar=cfar*charge
c
        return
        end
c
c
c
c
c
        subroutine dipole2far(zk,xy,dipstr,dipvec,dir,cfar)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates the far field
c       signatures in the direction dir due to the monochromatic
c       dipole located at the at the location xy
c
c          Input parameters:
c
c       zk (complex *16)  - the frequency parameter
c       xy (real *8 ) - the location of the charge in R^2
c       dir (real *8 ) - the far field direction vector in R^2 
c       dipstr (complex *16) - the strength of the dipole
c       dipvec (real *8) - the orientation vector of the dipole
c
c          Output parameters:
c
c       cfar (complex*16) - the far field signature of the potential
c
        dimension xy(2),dir(2),far(2),dipvec(2)
        complex *16 zk,dipstr,cfar,ima
        data ima/(0.0d0,1.0d0)/
c
c       ... first, find the point on the unit circle
c        
        d=sqrt(dir(1)**2+dir(2)**2)
        far(1)=dir(1)/d
        far(2)=dir(2)/d
c
c       ... find the far-field signature of the unit dipole located at
c       the location xy
c
        d=xy(1)*far(1)+xy(2)*far(2)
        f=dipvec(1)*far(1)+dipvec(2)*far(1)
c
        cfar=-ima*zk*f*exp(-ima*zk*d)
c
        cfar=cfar*dipstr
c
        return
        end
c
c
c
c
c
C***********************************************************************
        subroutine h2dformfar(ier,zk,rscale,source,charge,ns,center,
     1                       nfar,mpfar)
        implicit real *8 (a-h,o-z)
C***********************************************************************
c
c     This subroutine constructs a far field signature about CENTER due
c     to NS sources located at SOURCES(2,*).
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     zk              : Helmholtz parameter 
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     charge(ns)      : source strengths
c     ns              : number of sources
c     center(2)       : expansion center
c     nfar            : length of the far field signature
c
c     OUTPUT:
c
c     ier       : error return code
c                 ier=0 returned successfully;
c
c     mpfar     : far field signature
c
        complex *16 zk,mpfar(1:nfar),charge(*),cfar
        real *8 center(2),source(2,*),zdiff(2),dir(2)
c
        data ima/(0.0d0,1.0d0)/
c
c
        do n=1,nfar
        mpfar(n)=0
        enddo

        done = 1.0d0
        pi = 4*atan(done)
        
        do n=1,nfar
        theta = 2*pi*(n-1)/dble(nfar)
        x = cos(theta)
        y = sin(theta)
        dir(1) = x
        dir(2) = y
        do j=1,ns
c
        zdiff(1)=source(1,j)-center(1)
        zdiff(2)=source(2,j)-center(2)
c
        call charge2far(zk,zdiff,charge(j),dir,cfar)
        mpfar(n) = mpfar(n) + cfar
        
        enddo
        enddo

        return
        end
c
c
c
c
c
        subroutine h2dpwleval(zk,center,pwloc,nloc,target,ntarget,pot)
        implicit real *8 (a-h,o-z)
        complex *16 zk,trans
        dimension center(2),target(2,ntarget),dir(2)
        complex *16 ima,ima4,ima4inv
        complex *16 pwloc(1:nloc),pot(ntarget)
        data ima/(0.0d0,1.0d0)/

        done=1.0d0
        pi=4*atan(done)
        
        ima4=-4*ima
        ima4inv=ima/4

        do k = 1,ntarget

        dx = target(1,k)-center(1)
        dy = target(2,k)-center(2)

        pot(k) = dcmplx(0.0d0,0.0d0)

        do i = 1,nloc
        theta = 2*pi*(i-1)/dble(nloc)
        x = cos(theta)
        y = sin(theta)
        dir(1) = x
        dir(2) = y
        
        proj = dx*dir(1) + dy*dir(2)
        trans = exp(ima*zk*proj)
        pot(k) = pot(k) + pwloc(i) * trans
        enddo
        
        pot(k) = pot(k) / dble(nloc) * ima4inv
        enddo
        
        return
        end
