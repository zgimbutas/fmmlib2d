ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        this is the end of the debugging code and the beginning 
c        of the exponential representation routines
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
        subroutine hwts2e(ier,iprec,rk,cxs,cws,zcc,n,err)
        implicit real *8 (a-h,o-z)
        complex *16 cxs(1),cws(1),zcc
        complex *16 rk,ima,rk_image
        complex *16 xk,wk,z,h0,h1
        dimension xs(1 000), ws(1 000)
        data ima/(0.0d0,1.0d0)/
c
c       
c       This subroutine returns the nodes and weights for 
c       the exponential representation 
c
c       H_0(rk*sqrt(u^2+v^2)) = 2 / (ima pi)
c       \int_0^{+\infty} exp(-sqrt(t^2-rk^2)*u) cos(tv) /sqrt(t^2-rk^2) dt
c
c       u in [1,4], v in [0, 4]
c
c       All tables are valid in 
c
c       rk range [-16*pi, 16*pi] \times [0, 12*pi]
c
c
c
c       Input parameter:
c       
c       iprec - precision required:
c           iprec = 0 - 2 digits
c           iprec = 1 - 3 digits
c           iprec = 2 - 6 digits
c           iprec = 3 - 9 digits
c           iprec = 4 - 12 digits
c
c       Output parameters:
c
c       n - the number of nodes
c       cxs - the nodes: complex *16 (n)
c       cws - the weights: complex *16 (n)
c
c        
        if( dble(rk) .gt. 0 ) then
        call hwts2ep(ier,iprec,rk,cxs,cws,zcc,n,err)
        else
        rk_image = dcmplx(-dble(rk),imag(rk))
        call hwts2ep(ier,iprec,rk_image,cxs,cws,zcc,n,err)
        do i = 1,n
        cxs(i) = dconjg(cxs(i))
        cws(i) = dconjg(cws(i))
        enddo
        zcc = dcmplx(-dble(zcc),imag(zcc))
        endif
c
        return
        end
c
c        
c
c
c
        subroutine hwts2ep(ier,iprec,rk,cxs,cws,zcc,n,err)
        implicit real *8 (a-h,o-z)
        complex *16 cxs(1),cws(1),zcc
        complex *16 rk,ima
        complex *16 xk,wk,z,h0,h1
        dimension xs(1 000), ws(1 000)
        data ima/(0.0d0,1.0d0)/
c
c       
c       This subroutine returns the nodes and weights for 
c       the exponential representation 
c
c       H_0(rk*sqrt(u^2+v^2)) = 2 / (ima pi)
c       \int_0^{+\infty} exp(-sqrt(t^2-rk^2)*u) cos(tv) /sqrt(t^2-rk^2) dt
c
c       u in [1,4], v in [0, 4]
c
c       All tables are valid in 
c
c       rk range [0, 16*pi] \times [0, 12*pi]
c
c
c
c       Input parameter:
c       
c       iprec - precision required:
c           iprec = 0 - 2 digits
c           iprec = 1 - 3 digits
c           iprec = 2 - 6 digits
c           iprec = 3 - 9 digits
c           iprec = 4 - 12 digits
c
c       Output parameters:
c
c       n - the number of nodes
c       cxs - the nodes: complex *16 (n)
c       cws - the weights: complex *16 (n)
c
        done = 1
        pi = 4*atan(done)
c
        n = 0
c
        ier = 0
c       
        if( iprec .eq. 0 ) errmax = .1d-2
        if( iprec .eq. 1 ) errmax = .1d-3
        if( iprec .eq. 2 ) errmax = .1d-6
        if( iprec .eq. 3 ) errmax = .1d-9
        if( iprec .eq. 4 ) errmax = .1d-12
        if( iprec .eq. 5 ) errmax = .1d-15
        
        zcc = 0
        err = 0

        if( dble(rk) .lt. 0 ) ier = 1
        if( imag(rk) .lt. 0 ) ier = 1
        if( ier .ne. 0 ) return

        call hwts2getd(ier,rk,idomain)
c        call prinf('ier=*',ier,1)
ccc        write(*,*) 'idomain=*', idomain
c        call prinf('idomain=*',idomain,1)
c
        if( idomain .gt. 1 ) then
c
        if( iprec .eq. 0 ) call hwts2p0(idomain,xs,ws,n,iquadtype,err)
        if( iprec .eq. 1 ) call hwts2p1(idomain,xs,ws,n,iquadtype,err)
        if( iprec .eq. 2 ) call hwts2p2(idomain,xs,ws,n,iquadtype,err)
        if( iprec .eq. 3 ) call hwts2p3(idomain,xs,ws,n,iquadtype,err)
        if( iprec .eq. 4 ) call hwts2p4(idomain,xs,ws,n,iquadtype,err)
        if( iprec .eq. 5 ) call hwts2p5(idomain,xs,ws,n,iquadtype,err)
        
        else

        call hwts2getdz(ier,rk,idomain)
c        call prinf('after getdz, ier=*',ier,1)
ccc        write(*,*) 'idomain=*', idomain
c        call prinf('after getdz, idomain=*',idomain,1)
        
        if( iprec .eq. 0 ) call hwts2zp0(idomain,xs,ws,n,iquadtype,err)
        if( iprec .eq. 1 ) call hwts2zp1(idomain,xs,ws,n,iquadtype,err)
        if( iprec .eq. 2 ) call hwts2zp2(idomain,xs,ws,n,iquadtype,err)
        if( iprec .eq. 3 ) call hwts2zp3(idomain,xs,ws,n,iquadtype,err)
        if( iprec .eq. 4 ) call hwts2zp4(idomain,xs,ws,n,iquadtype,err)
        if( iprec .eq. 5 ) call hwts2zp5(idomain,xs,ws,n,iquadtype,err)

        endif
        
 1100   continue
        
c
c        call prinf('n=*',n,1)
c        call prin2('xs=*',xs,n)
c        call prin2('ws=*',ws,n)
c        call prinf('iquadtype=*',iquadtype,1)
c
c        
        do 1200 k=1,n
c
        if( iquadtype .eq. 1 ) then
        cxs(k) = xs(k)
        cws(k) = ws(k)
        endif
c        
        if( iquadtype .eq. 3 ) then
        cxs(k) = xs(k)
        cws(k) = ws(k)
        xk = cxs(k)
        wk = cws(k)
        zcc = zcc - wk*exp(-1.0d0*sqrt(xk**2-rk**2))*2*cos(0.0d0*xk)/
     $     sqrt(xk**2-rk**2) /pi/ima       
        endif
c        
        if( iquadtype .eq. 5 ) then
        cxs(k) = xs(k)-ima*atan(xs(k))
        cws(k) = ws(k)*(1-ima/(1+xs(k)**2))
        endif
c
        if( iquadtype .eq. 6 ) then
        cxs(k) = xs(k)-ima*atan(xs(k))
        cws(k) = ws(k)*(1-ima/(1+xs(k)**2))
        xk = cxs(k)
        wk = cws(k)
        zcc = zcc - wk*exp(-1.0d0*sqrt(xk**2-rk**2))*2*cos(0.0d0*xk)/
     $     sqrt(xk**2-rk**2) /pi/ima       
        endif
c
        if( iquadtype .eq. 7 ) then
        cxs(k) = xs(k)-ima*xs(k)/(1+xs(k)**2)
        cws(k) = ws(k)*(1-ima*(1-xs(k)**2)/(1+xs(k)**2)**2)
        endif
c
        if( iquadtype .eq. 8 ) then
        cxs(k) = xs(k)-ima*xs(k)/(1+xs(k)**2)
        cws(k) = ws(k)*(1-ima*(1-xs(k)**2)/(1+xs(k)**2)**2)
        xk = cxs(k)
        wk = cws(k)
        zcc = zcc - wk*exp(-1.0d0*sqrt(xk**2-rk**2))*2*cos(0.0d0*xk)/
     $     sqrt(xk**2-rk**2) /pi/ima       
        endif
c
 1200   continue
c
        if( iquadtype .eq. 3 ) then
        rho = sqrt(1.0d0**2+0.0d0**2)
        z = rk*rho
        theta = atan2(0.0d0,1.0d0)
        ifexpon=1
        call hank103(z,h0,h1,ifexpon)  
        zcc = zcc + h0
        endif
        
        if( iquadtype .eq. 6 ) then
        rho = sqrt(1.0d0**2+0.0d0**2)
        z = rk*rho
        theta = atan2(0.0d0,1.0d0)
        ifexpon=1
        call hank103(z,h0,h1,ifexpon)
        zcc = zcc + h0
        endif
        
        if( iquadtype .eq. 8 ) then
        rho = sqrt(1.0d0**2+0.0d0**2)
        z = rk*rho
        theta = atan2(0.0d0,1.0d0)
        ifexpon=1
        call hank103(z,h0,h1,ifexpon)
        zcc = zcc + h0
        endif
        
        return
        end
c
c
c
c
        subroutine hwts2getd(ier,rk,idomain)
        implicit real *8 (a-h,o-z)
        complex *16 ima,rk
        dimension cx(100),cy(100)
        data ima/(0.0d0,1.0d0)/
c
        done=1
        pi=4*atan(done)
c
c       rk range [0, 16*pi] \times [0, 12*pi]
c
        nx=23
        ny=19
c
        cx(1)=0.00d0
        cx(2)=0.02d0
        cx(3)=0.04d0
        cx(4)=0.10d0
        cx(5)=0.20d0
        cx(6)=0.40d0
        cx(7)=1.00d0
        cx(8)=pi/2
        cx(9)=pi
        cx(10)=2*pi
        cx(11)=3*pi
        cx(12)=4*pi
        cx(13)=5*pi
        cx(14)=6*pi
        cx(15)=7*pi
        cx(16)=8*pi
        cx(17)=9*pi
        cx(18)=10*pi
        cx(19)=11*pi
        cx(20)=12*pi
        cx(21)=13*pi
        cx(22)=14*pi
        cx(23)=15*pi
        cx(24)=16*pi
c
        cy(1)=0.00d0
        cy(2)=0.02d0
        cy(3)=0.04d0
        cy(4)=0.10d0
        cy(5)=0.20d0
        cy(6)=0.40d0
        cy(7)=1.00d0
        cy(8)=pi/2
        cy(9)=pi
        cy(10)=2*pi
        cy(11)=3*pi
        cy(12)=4*pi
        cy(13)=5*pi
        cy(14)=6*pi
        cy(15)=7*pi
        cy(16)=8*pi
        cy(17)=9*pi
        cy(18)=10*pi
        cy(19)=11*pi
        cy(20)=12*pi
        cy(21)=13*pi
        cy(22)=14*pi
        cy(23)=15*pi
        cy(24)=16*pi
c
        rkrea=dble(rk)
        rkima=imag(rk)
c
        idomain = 0
        ier = 0
c
        ix=0
        do 1200 i=1,nx
        if( rkrea .ge. cx(i+1) ) ix=i+1
 1200   continue
c
 1300   continue
c
        iy=0
        do 2200 i=1,ny
        if( rkima .ge. cy(i+1) ) iy=i+1
 2200   continue
c
 2300   continue
c
ccc        write(*,*) rk
ccc        write(*,*) ix,iy
c
        if( ix .gt. nx ) ier = 4
        if( iy .gt. ny ) ier = 4
        if( ix .gt. nx ) return
        if( iy .gt. ny ) return
c       
        if( ix .eq. 0 ) ix = 1
        if( iy .eq. 0 ) iy = 1
c
        if( ix .gt. nx ) ix = nx
        if( iy .gt. ny ) iy = ny
c
        idomain = ix + (iy-1) * nx
c
ccc        write(*,*) ix,iy
ccc        write(*,*) idomain
c
ccc        stop
c
        return
        end
c
c
c
c
        subroutine hwts2getdz(ier,rk,idomain)
        implicit real *8 (a-h,o-z)
        complex *16 ima,rk,rkmin,rkmax
        dimension cx(100),cy(100)
        data ima/(0.0d0,1.0d0)/
c
        done=1
        pi=4*atan(done)
c
c       rk range [0, 0.02] \times [0, 0.02]
c
        ier = 0
        ndomains = 76
        do inum = 1,ndomains
        call getrkinfo0(inum,ndomains,rkmin,rkmax)
        if( (dble(rk) .ge. dble(rkmin)) .and.
     $      (imag(rk) .ge. imag(rkmin)) .and.
     $      (dble(rk) .le. dble(rkmax)) .and.
     $      (imag(rk) .le. imag(rkmax))) then
        idomain = inum
        return
        endif
        enddo
        
        if( dble(rk) .lt. 0 ) ier = 1
        if( imag(rk) .lt. 0 ) ier = 1
        if( dble(rk) .gt. 0.02d0 ) ier = 4
        if( imag(rk) .gt. 0.02d0 ) ier = 4
        idomain = 0
        
        return
        end
c
c        
c        
c
        subroutine getrkinfo0(inum,ndomains,rkmin,rkmax)
        implicit real *8 (a-h,o-z)
        complex *16 rkmin,rkmax,ima
        complex *16 cx(100),cy(100)
        data ima/(0.0d0,1.0d0)/
c
        done=1
        pi=4*atan(done)
c
c       rk range [0, 0.02] \times [0, 0.02]
c
        cx(1)=0d0
        cx(2)=1d-10
        cx(3)=2d-10
        cx(4)=4d-10
        cx(5)=1d-9
        cx(6)=2d-9
        cx(7)=4d-9
        cx(8)=1d-8
        cx(9)=2d-8
        cx(10)=4d-8
        cx(11)=1d-7
        cx(12)=2d-7
        cx(13)=4d-7
        cx(14)=1d-6
        cx(15)=2d-6
        cx(16)=4d-6
        cx(17)=1d-5
        cx(18)=2d-5
        cx(19)=4d-5
        cx(20)=1d-4
        cx(21)=2d-4
        cx(22)=4d-4
        cx(23)=1d-3
        cx(24)=2d-3
        cx(25)=4d-3
        cx(26)=1d-2
        cx(27)=2d-2
c
        ndomains=25*3+1
c
        rkmin=0
        rkmax=0
c
        if( inum .eq. 1 ) then
        rkmin=cx(1)+ima*cx(1)
        rkmax=cx(2)+ima*cx(2)
        endif

        if( inum .gt. 1 .and. inum .le. ndomains ) then
c
        iy=mod(inum-2,3)+1
        ix=(inum-2)/3+1
c
ccc        write(*,*) ix, iy
        if( iy .eq. 1 ) then
        rkmin=0+ima*cx(ix+1)
        rkmax=cx(ix+1)+ima*cx(ix+2)
        endif
        if( iy .eq. 2 ) then
        rkmin=cx(ix+1)+ima*0
        rkmax=cx(ix+2)+ima*cx(ix+1)
        endif
        if( iy .eq. 3 ) then
        rkmin=cx(ix+1)+ima*cx(ix+1)
        rkmax=cx(ix+2)+ima*cx(ix+2)
        endif
c
        
        endif
        
        return
        end
c
c
c
c
c
        subroutine getrkinfo2(inum,ndomains,rkmin,rkmax)
        implicit real *8 (a-h,o-z)
        complex *16 rkmin,rkmax,ima
        complex *16 cx(100),cy(100)
        data ima/(0.0d0,1.0d0)/
c
        done=1
        pi=4*atan(done)
c
c       rk range [0, 16*pi] \times [0, 12*pi]
c
        nx=23
        ny=19
c
        cx(1)=0.00d0
        cx(2)=0.02d0
        cx(3)=0.04d0
        cx(4)=0.10d0
        cx(5)=0.20d0
        cx(6)=0.40d0
        cx(7)=1.00d0
        cx(8)=pi/2
        cx(9)=pi
        cx(10)=2*pi
        cx(11)=3*pi
        cx(12)=4*pi
        cx(13)=5*pi
        cx(14)=6*pi
        cx(15)=7*pi
        cx(16)=8*pi
        cx(17)=9*pi
        cx(18)=10*pi
        cx(19)=11*pi
        cx(20)=12*pi
        cx(21)=13*pi
        cx(22)=14*pi
        cx(23)=15*pi
        cx(24)=16*pi
c
        cy(1)=0.00d0
        cy(2)=0.02d0
        cy(3)=0.04d0
        cy(4)=0.10d0
        cy(5)=0.20d0
        cy(6)=0.40d0
        cy(7)=1.00d0
        cy(8)=pi/2
        cy(9)=pi
        cy(10)=2*pi
        cy(11)=3*pi
        cy(12)=4*pi
        cy(13)=5*pi
        cy(14)=6*pi
        cy(15)=7*pi
        cy(16)=8*pi
        cy(17)=9*pi
        cy(18)=10*pi
        cy(19)=11*pi
        cy(20)=12*pi
        cy(21)=13*pi
        cy(22)=14*pi
        cy(23)=15*pi
        cy(24)=16*pi
c
        ndomains=nx*ny
c
        ix=mod(inum-1,nx)+1
        iy=(inum-1)/nx+1
c
        rkmin=cx(ix)+ima*cy(iy)
        rkmax=cx(ix+1)+ima*cy(iy+1)
c
        return
        end
c
c
c
c
c
        subroutine hwts2p0(idomain,xs,ws,n,iquadtype,err)
        implicit real *8 (a-h,o-z)
        dimension xs(*), ws(*)
c
        n = 0
        iquadtype = 0
c
        INCLUDE 'hwts_2d_iprec0.txt'
c
        return
        end
c        
c
c
c
c
        subroutine hwts2p1(idomain,xs,ws,n,iquadtype,err)
        implicit real *8 (a-h,o-z)
        dimension xs(*), ws(*)
c
        n = 0
        iquadtype = 0
c
        INCLUDE 'hwts_2d_iprec1.txt'
c
        return
        end
c        
c
c
c
c
        subroutine hwts2p2(idomain,xs,ws,n,iquadtype,err)
        implicit real *8 (a-h,o-z)
        dimension xs(*), ws(*)
c
        n = 0
        iquadtype = 0
c
        INCLUDE 'hwts_2d_iprec2.txt'
c
        return
        end        
c        
c
c
c
c
        subroutine hwts2p3(idomain,xs,ws,n,iquadtype,err)
        implicit real *8 (a-h,o-z)
        dimension xs(*), ws(*)
c
        n = 0
        iquadtype = 0
c
        INCLUDE 'hwts_2d_iprec3.txt'
c
        return
        end
c        
c
c
c
c
        subroutine hwts2p4(idomain,xs,ws,n,iquadtype,err)
        implicit real *8 (a-h,o-z)
        dimension xs(*), ws(*)
c
        n = 0
        iquadtype = 0
c
        INCLUDE 'hwts_2d_iprec4.txt'
c
        return
        end
c
c
c
c
c
        subroutine hwts2p5(idomain,xs,ws,n,iquadtype,err)
        implicit real *8 (a-h,o-z)
        dimension xs(*), ws(*)
c
        n = 0
        iquadtype = 0
c
        INCLUDE 'hwts_2d_iprec5.txt'
c
        return
        end
c
c
c
c
c
        subroutine hwts2zp0(idomain,xs,ws,n,iquadtype,err)
        implicit real *8 (a-h,o-z)
        dimension xs(*), ws(*)
c
        n = 0
        iquadtype = 0
c
        INCLUDE 'hwts_2dz_iprec0.txt'
c
        return
        end
c        
c
c
c
c
        subroutine hwts2zp1(idomain,xs,ws,n,iquadtype,err)
        implicit real *8 (a-h,o-z)
        dimension xs(*), ws(*)
c
        n = 0
        iquadtype = 0
c
        INCLUDE 'hwts_2dz_iprec1.txt'
c
        return
        end
c        
c
c
c
c
        subroutine hwts2zp2(idomain,xs,ws,n,iquadtype,err)
        implicit real *8 (a-h,o-z)
        dimension xs(*), ws(*)
c
        n = 0
        iquadtype = 0
c
        INCLUDE 'hwts_2dz_iprec2.txt'
c
        return
        end        
c        
c
c
c
c
        subroutine hwts2zp3(idomain,xs,ws,n,iquadtype,err)
        implicit real *8 (a-h,o-z)
        dimension xs(*), ws(*)
c
        n = 0
        iquadtype = 0
c
        INCLUDE 'hwts_2dz_iprec3.txt'
c
        return
        end
c        
c
c
c
c
        subroutine hwts2zp4(idomain,xs,ws,n,iquadtype,err)
        implicit real *8 (a-h,o-z)
        dimension xs(*), ws(*)
c
        n = 0
        iquadtype = 0
c
        INCLUDE 'hwts_2dz_iprec4.txt'
c
        return
        end
c        
c
c
c
c
        subroutine hwts2zp5(idomain,xs,ws,n,iquadtype,err)
        implicit real *8 (a-h,o-z)
        dimension xs(*), ws(*)
c
        n = 0
        iquadtype = 0
c
        INCLUDE 'hwts_2dz_iprec5.txt'
c
        return
        end
