
        subroutine hwtsexp2(ier,iprec,rk,xs,ws,zcc,nexp,err)
        implicit real *8 (a-h,o-z)
        complex *16 rk,xs(*),ws(*),zcc
c
        call hwts2e(ier,iprec,rk,xs,ws,zcc,nexp,err)
        do i = 1,nexp
        xs(i+nexp) = xs(i)
        ws(i+nexp) = ws(i)
        xs(i) = -xs(i)
        enddo
        nexp = 2*nexp
c
        return
        end
c
c
c
c
c
        subroutine h2expini(rk,rscale,nterms,xs,ws,zcc,nexp,
     $     cscale,cshift,zmuls1f,zmuls2f,zmuls1b,zmuls2b)
        implicit real *8 (a-h,o-z)
        complex *16 cshift(nexp,-3:3,-3:3),rk,cosbeta,sinbeta
        complex *16 xs(nexp),ws(nexp),zcc,ima,khat
        complex *16 cscale(nexp)
        complex *16 zmul1,zmul2,ztemp1,ztemp2
        complex *16 zmuls1f(nexp,nterms),zmuls2f(nexp,nterms)
        complex *16 zmuls1b(nexp,nterms),zmuls2b(nexp,nterms)
        data ima/(0.0d0,1.0d0)/
c
        done=1.0d0
        pi=4*atan(done)
c
c       ... precompute exponential translation tables
c
        do k = 1,nexp
        do i = -3,3
        do j = -3,3
        cshift(k,i,j) = exp(-sqrt(xs(k)**2-rk**2)*i + ima*xs(k)*j)
        enddo
        enddo
        enddo

        do k = 1,nexp
            cscale(k) = ws(k)/sqrt(xs(k)**2-rk**2)/pi/ima
        enddo

        do k = 1,nexp
        khat = (sqrt(xs(k)**2-rk**2)+xs(k))/rk
        if( k .le. nexp/2) khat = -rk/(sqrt(xs(k)**2-rk**2)-xs(k))
        ztemp1 = khat*rscale
        ztemp2 = rscale/khat
        zmul1 = ztemp1
        zmul2 = ztemp2
        do n = 1,nterms-3,4
           zmuls1f(k,n+0)=zmul1
           zmul1 = zmul1*ztemp1
           zmuls1f(k,n+1)=zmul1
           zmul1 = zmul1*ztemp1
           zmuls1f(k,n+2)=zmul1
           zmul1 = zmul1*ztemp1
           zmuls1f(k,n+3)=zmul1
           zmul1 = zmul1*ztemp1
        enddo
        do n = 1,nterms-3,4
           zmuls2f(k,n+0)=zmul2
           zmul2 = zmul2*ztemp2
           zmuls2f(k,n+1)=zmul2
           zmul2 = zmul2*ztemp2
           zmuls2f(k,n+2)=zmul2
           zmul2 = zmul2*ztemp2
           zmuls2f(k,n+3)=zmul2
           zmul2 = zmul2*ztemp2
        enddo
        enddo

        do k = 1,nexp
        khat = (sqrt(xs(k)**2-rk**2)-xs(k))/rk
        if( k .gt. nexp/2) khat = -rk/(sqrt(xs(k)**2-rk**2)+xs(k))
        ztemp1 = khat*rscale
        ztemp2 = rscale/khat
        zmul1 = ztemp1
        zmul2 = ztemp2
        do n = 1,nterms-3,4
           zmuls1b(k,n+0)=zmul1
           zmul1 = zmul1*ztemp1
           zmuls1b(k,n+1)=zmul1
           zmul1 = zmul1*ztemp1
           zmuls1b(k,n+2)=zmul1
           zmul1 = zmul1*ztemp1
           zmuls1b(k,n+3)=zmul1
           zmul1 = zmul1*ztemp1
        enddo
        do n = 1,nterms-3,4
           zmuls2b(k,n+0)=zmul2
           zmul2 = zmul2*ztemp2
           zmuls2b(k,n+1)=zmul2
           zmul2 = zmul2*ztemp2
           zmuls2b(k,n+2)=zmul2
           zmul2 = zmul2*ztemp2
           zmuls2b(k,n+3)=zmul2
           zmul2 = zmul2*ztemp2
        enddo
        enddo

        return
        end
c
c
c
c
c
        subroutine h2dfft4f(r1,r2,r3,r4,s1,s2,s3,s4)
        implicit real *8 (a-h,o-z)
        complex *16 r1,r2,r3,r4,s1,s2,s3,s4
        complex *16 t1,t2,t3,t4
c
c       s1 = (r1 + r2     + r3 + r4    )
c       s2 = (r1 - r2*1im - r3 + r4*1im)
c       s3 = (r1 - r2     + r3 - r4    )
c       s4 = (r1 + r2*1im - r3 - r4*1im)
c
        t1 = r1 + r3
        t2 = r1 - r3
        t3 = r2 + r4
        t4 = r2 - r4
        s1 = t1 + t3
        s2 = t2 - dcmplx(-imag(t4), dble(t4))
        s3 = t1 - t3
        s4 = t2 + dcmplx(-imag(t4), dble(t4))
c
        return
        end
c
c
c
c
c
        subroutine h2dfft4b(r1,r2,r3,r4,s1,s2,s3,s4)
        implicit real *8 (a-h,o-z)
        complex *16 r1,r2,r3,r4,s1,s2,s3,s4
        complex *16 t1,t2,t3,t4
c
c       s1 = (r1 + r2     + r3 + r4    )
c       s2 = (r1 + r2*1im - r3 - r4*1im)
c       s3 = (r1 - r2     + r3 - r4    )
c       s4 = (r1 - r2*1im - r3 + r4*1im)
c
        t1 = r1 + r3
        t2 = r1 - r3
        t3 = r2 + r4
        t4 = r2 - r4
        s1 = t1 + t3
        s2 = t2 - dcmplx(imag(t4), -dble(t4))
        s3 = t1 - t3
        s4 = t2 + dcmplx(imag(t4), -dble(t4))
c
        return
        end
c
c
c
c
c
        subroutine h2dmpexp_fast(nd,rk,rscale,center,mpole,nterms,
     $     expout,nexp,xs,ws,zcc,cscale,zmuls1,zmuls2)
        implicit real *8 (a-h,o-z)
        dimension center(2)
        complex *16 mpole(nd,-nterms:nterms),rk
        complex *16 expout(nd,0:nexp,4)
        complex *16 xs(nexp),ws(nexp),zcc
        complex *16 cscale(nexp)
        complex *16 khat,zmul1,zmul2,ztemp1,ztemp2,cmul
        complex *16 s1,s2,s3,s4,r1,r2,r3,r4
        complex *16 zmuls1(nexp,nterms),zmuls2(nexp,nterms)
        complex *16 ima
        data ima/(0.0d0,1.0d0)/

        do k = 1,nexp
        do idim=1,nd
        expout(idim,k,1) = expout(idim,k,1) + mpole(idim,0)
        expout(idim,k,2) = expout(idim,k,2) + mpole(idim,0)
        expout(idim,k,3) = expout(idim,k,3) + mpole(idim,0)
        expout(idim,k,4) = expout(idim,k,4) + mpole(idim,0)
        enddo
        do idim=1,nd
        s1 = dcmplx(0.0d0,0.0d0)
        s2 = dcmplx(0.0d0,0.0d0)
        s3 = dcmplx(0.0d0,0.0d0)
        s4 = dcmplx(0.0d0,0.0d0)
        do n = 1,nterms-3,4
           s2 = s2 + zmuls1(k,n+0) * mpole(idim,n+0)
           s3 = s3 + zmuls1(k,n+1) * mpole(idim,n+1)
           s4 = s4 + zmuls1(k,n+2) * mpole(idim,n+2)
           s1 = s1 + zmuls1(k,n+3) * mpole(idim,n+3)
        enddo
        call h2dfft4b(s1, s2, s3, s4, r1, r2, r3, r4)
        expout(idim,k,1) = expout(idim,k,1) + r1
        expout(idim,k,2) = expout(idim,k,2) + r2
        expout(idim,k,3) = expout(idim,k,3) + r3
        expout(idim,k,4) = expout(idim,k,4) + r4
        s1 = dcmplx(0.0d0,0.0d0)
        s2 = dcmplx(0.0d0,0.0d0)
        s3 = dcmplx(0.0d0,0.0d0)
        s4 = dcmplx(0.0d0,0.0d0)
        do n = 1,nterms-3,4
           s2 = s2 + zmuls2(k,n+0) * mpole(idim,-n-0)
           s3 = s3 + zmuls2(k,n+1) * mpole(idim,-n-1)
           s4 = s4 + zmuls2(k,n+2) * mpole(idim,-n-2)
           s1 = s1 + zmuls2(k,n+3) * mpole(idim,-n-3)
        enddo
        call h2dfft4f(s1, s2, s3, s4, r1, r2, r3, r4)
        expout(idim,k,1) = expout(idim,k,1) + r1
        expout(idim,k,2) = expout(idim,k,2) + r2
        expout(idim,k,3) = expout(idim,k,3) + r3
        expout(idim,k,4) = expout(idim,k,4) + r4
        enddo
        enddo


        do k = 1,nexp
        do idim=1,nd
        expout(idim,k,1) = expout(idim,k,1) * cscale(k)
        expout(idim,k,2) = expout(idim,k,2) * cscale(k)
        expout(idim,k,3) = expout(idim,k,3) * cscale(k)
        expout(idim,k,4) = expout(idim,k,4) * cscale(k)
        enddo
        enddo

c
c       ... constant terms
c
        do idim=1,nd
        expout(idim,0,1) = mpole(idim,0) * zcc
        expout(idim,0,2) = mpole(idim,0) * zcc
        expout(idim,0,3) = mpole(idim,0) * zcc
        expout(idim,0,4) = mpole(idim,0) * zcc
        enddo

        return
        end
c
c
c
c
c
        subroutine h2dmpexp(nd,rk,rscale,center,mpole,nterms,
     $     expout,nexp,xs,ws,zcc,cscale)
        implicit real *8 (a-h,o-z)
        dimension center(2)
        complex *16 mpole(nd,-nterms:nterms),rk
        complex *16 expout(nd,0:nexp,4)
        complex *16 xs(nexp),ws(nexp),zcc
        complex *16 cscale(nexp)
        complex *16 khat,zmul1,zmul2,ztemp1,ztemp2,cmul
        complex *16 s1,s2,s3,s4,r1,r2,r3,r4
        complex *16 zmuls1(nexp,nterms),zmuls2(nexp,nterms)
        complex *16 ima
        data ima/(0.0d0,1.0d0)/

        if( 2 .eq. 2 ) then

        do k = 1,nexp

        khat = (sqrt(xs(k)**2-rk**2)+xs(k))/rk
        if( k .le. nexp/2) khat = -rk/(sqrt(xs(k)**2-rk**2)-xs(k))

        ztemp1 = khat*rscale
        ztemp2 = rscale/khat
        zmul1 = ztemp1
        zmul2 = ztemp2
        do n = 1,nterms-3,4
           zmuls1(k,n+0)=zmul1
           zmul1 = zmul1*ztemp1
           zmuls1(k,n+1)=zmul1
           zmul1 = zmul1*ztemp1
           zmuls1(k,n+2)=zmul1
           zmul1 = zmul1*ztemp1
           zmuls1(k,n+3)=zmul1
           zmul1 = zmul1*ztemp1
        enddo
        do n = 1,nterms-3,4
           zmuls2(k,n+0)=zmul2
           zmul2 = zmul2*ztemp2
           zmuls2(k,n+1)=zmul2
           zmul2 = zmul2*ztemp2
           zmuls2(k,n+2)=zmul2
           zmul2 = zmul2*ztemp2
           zmuls2(k,n+3)=zmul2
           zmul2 = zmul2*ztemp2
        enddo

        enddo


        do k = 1,nexp
        do idim=1,nd
        expout(idim,k,1) = expout(idim,k,1) + mpole(idim,0)
        expout(idim,k,2) = expout(idim,k,2) + mpole(idim,0)
        expout(idim,k,3) = expout(idim,k,3) + mpole(idim,0)
        expout(idim,k,4) = expout(idim,k,4) + mpole(idim,0)
        enddo
        do idim=1,nd
        s1 = dcmplx(0.0d0,0.0d0)
        s2 = dcmplx(0.0d0,0.0d0)
        s3 = dcmplx(0.0d0,0.0d0)
        s4 = dcmplx(0.0d0,0.0d0)
        do n = 1,nterms-3,4
           s2 = s2 + zmuls1(k,n+0) * mpole(idim,n+0)
           s3 = s3 + zmuls1(k,n+1) * mpole(idim,n+1)
           s4 = s4 + zmuls1(k,n+2) * mpole(idim,n+2)
           s1 = s1 + zmuls1(k,n+3) * mpole(idim,n+3)
        enddo
        call h2dfft4b(s1, s2, s3, s4, r1, r2, r3, r4)
        expout(idim,k,1) = expout(idim,k,1) + r1
        expout(idim,k,2) = expout(idim,k,2) + r2
        expout(idim,k,3) = expout(idim,k,3) + r3
        expout(idim,k,4) = expout(idim,k,4) + r4
        s1 = dcmplx(0.0d0,0.0d0)
        s2 = dcmplx(0.0d0,0.0d0)
        s3 = dcmplx(0.0d0,0.0d0)
        s4 = dcmplx(0.0d0,0.0d0)
        do n = 1,nterms-3,4
           s2 = s2 + zmuls2(k,n+0) * mpole(idim,-n-0)
           s3 = s3 + zmuls2(k,n+1) * mpole(idim,-n-1)
           s4 = s4 + zmuls2(k,n+2) * mpole(idim,-n-2)
           s1 = s1 + zmuls2(k,n+3) * mpole(idim,-n-3)
        enddo
        call h2dfft4f(s1, s2, s3, s4, r1, r2, r3, r4)
        expout(idim,k,1) = expout(idim,k,1) + r1
        expout(idim,k,2) = expout(idim,k,2) + r2
        expout(idim,k,3) = expout(idim,k,3) + r3
        expout(idim,k,4) = expout(idim,k,4) + r4
        enddo
        enddo


        do k = 1,nexp
        do idim=1,nd
        expout(idim,k,1) = expout(idim,k,1) * cscale(k)
        expout(idim,k,2) = expout(idim,k,2) * cscale(k)
        expout(idim,k,3) = expout(idim,k,3) * cscale(k)
        expout(idim,k,4) = expout(idim,k,4) * cscale(k)
        enddo
        enddo

c
c       ... constant terms
c
        do idim=1,nd
        expout(idim,0,1) = mpole(idim,0) * zcc
        expout(idim,0,2) = mpole(idim,0) * zcc
        expout(idim,0,3) = mpole(idim,0) * zcc
        expout(idim,0,4) = mpole(idim,0) * zcc
        enddo

        endif


        if( 1 .eq. 2 ) then

        do k = 1,nexp

        khat = (sqrt(xs(k)**2-rk**2)+xs(k))/rk
        if( k .le. nexp/2) khat = -rk/(sqrt(xs(k)**2-rk**2)-xs(k))

        do idim=1,nd
        expout(idim,k,1) = expout(idim,k,1) + mpole(idim,0)
        expout(idim,k,2) = expout(idim,k,2) + mpole(idim,0)
        expout(idim,k,3) = expout(idim,k,3) + mpole(idim,0)
        expout(idim,k,4) = expout(idim,k,4) + mpole(idim,0)
        enddo
        ztemp1 = khat*rscale
        ztemp2 = rscale/khat
        zmul1 = ztemp1
        zmul2 = ztemp2
        do idim=1,nd
        s1 = dcmplx(0.0d0,0.0d0)
        s2 = dcmplx(0.0d0,0.0d0)
        s3 = dcmplx(0.0d0,0.0d0)
        s4 = dcmplx(0.0d0,0.0d0)
        do n = 1,nterms-3,4
           s2 = s2 + zmul1 * mpole(idim,n+0)
           zmul1 = zmul1*ztemp1
           s3 = s3 + zmul1 * mpole(idim,n+1)
           zmul1 = zmul1*ztemp1
           s4 = s4 + zmul1 * mpole(idim,n+2)
           zmul1 = zmul1*ztemp1
           s1 = s1 + zmul1 * mpole(idim,n+3)
           zmul1 = zmul1*ztemp1
        enddo
        call h2dfft4b(s1, s2, s3, s4, r1, r2, r3, r4)
        expout(idim,k,1) = expout(idim,k,1) + r1
        expout(idim,k,2) = expout(idim,k,2) + r2
        expout(idim,k,3) = expout(idim,k,3) + r3
        expout(idim,k,4) = expout(idim,k,4) + r4
        s1 = dcmplx(0.0d0,0.0d0)
        s2 = dcmplx(0.0d0,0.0d0)
        s3 = dcmplx(0.0d0,0.0d0)
        s4 = dcmplx(0.0d0,0.0d0)
        do n = 1,nterms-3,4
           s2 = s2 + zmul2 * mpole(idim,-n-0)
           zmul2 = zmul2*ztemp2
           s3 = s3 + zmul2 * mpole(idim,-n-1)
           zmul2 = zmul2*ztemp2
           s4 = s4 + zmul2 * mpole(idim,-n-2)
           zmul2 = zmul2*ztemp2
           s1 = s1 + zmul2 * mpole(idim,-n-3)
           zmul2 = zmul2*ztemp2
        enddo
        call h2dfft4f(s1, s2, s3, s4, r1, r2, r3, r4)
        expout(idim,k,1) = expout(idim,k,1) + r1
        expout(idim,k,2) = expout(idim,k,2) + r2
        expout(idim,k,3) = expout(idim,k,3) + r3
        expout(idim,k,4) = expout(idim,k,4) + r4
        enddo
        enddo


ccc        done = 1.0d0
ccc        pi = 4*atan(done)
        do k = 1,nexp
ccc        cmul = ws(k) /sqrt(xs(k)**2-rk**2) /pi/ima
        do idim=1,nd
        expout(idim,k,1) = expout(idim,k,1) * cscale(k)
        expout(idim,k,2) = expout(idim,k,2) * cscale(k)
        expout(idim,k,3) = expout(idim,k,3) * cscale(k)
        expout(idim,k,4) = expout(idim,k,4) * cscale(k)
        enddo
        enddo

c
c       ... constant terms
c
        do idim=1,nd
        expout(idim,0,1) = mpole(idim,0) * zcc
        expout(idim,0,2) = mpole(idim,0) * zcc
        expout(idim,0,3) = mpole(idim,0) * zcc
        expout(idim,0,4) = mpole(idim,0) * zcc
        enddo

        endif

        return
        end
c
c
c
c
c
        subroutine h2dexpta_fast(nd,rk,rscale,center,mploc,nterms,
     $     expinc,nexp,xs,ws,zmuls1,zmuls2)
        implicit real *8 (a-h,o-z)
        dimension center(2)
        complex *16 mploc(nd,-nterms:nterms),rk
        complex *16 expinc(nd,0:nexp,4)
        complex *16 xs(nexp),ws(nexp)
        complex *16 khat,zmul1,zmul2,ztemp1,ztemp2
        complex *16 s1,s2,s3,s4,r1,r2,r3,r4
        complex *16 zmuls1(nexp,nterms),zmuls2(nexp,nterms)
        complex *16 ima
        data ima/(0.0d0,1.0d0)/

        do k = 1,nexp
        do idim=1,nd
        mploc(idim,0) = mploc(idim,0) + expinc(idim,k,1)
     $     + expinc(idim,k,2) + expinc(idim,k,3) + expinc(idim,k,4)
        enddo
        do idim=1,nd
        call h2dfft4f(expinc(idim,k,3),expinc(idim,k,4),
     $     expinc(idim,k,1),expinc(idim,k,2),
     $     s1, s2, s3, s4)
        do n = 1,nterms-3,4
            mploc(idim,n+0) = mploc(idim,n+0) + zmuls1(k,n+0) * s2
            mploc(idim,n+1) = mploc(idim,n+1) + zmuls1(k,n+1) * s3
            mploc(idim,n+2) = mploc(idim,n+2) + zmuls1(k,n+2) * s4
            mploc(idim,n+3) = mploc(idim,n+3) + zmuls1(k,n+3) * s1
        enddo
        call h2dfft4b(expinc(idim,k,3),expinc(idim,k,4),
     $     expinc(idim,k,1),expinc(idim,k,2),
     $     s1, s2, s3, s4)
        do n = 1,nterms-3,4
            mploc(idim,-n-0) = mploc(idim,-n-0) + zmuls2(k,n+0) * s2
            mploc(idim,-n-1) = mploc(idim,-n-1) + zmuls2(k,n+1) * s3
            mploc(idim,-n-2) = mploc(idim,-n-2) + zmuls2(k,n+2) * s4
            mploc(idim,-n-3) = mploc(idim,-n-3) + zmuls2(k,n+3) * s1
        enddo
        enddo
        enddo

c
c       ... constant terms
c

        do idim=1,nd
        mploc(idim,0) = mploc(idim,0) + expinc(idim,0,1)
     $     + expinc(idim,0,2) + expinc(idim,0,3) + expinc(idim,0,4)
        enddo

        return
        end
c
c
c
c
c
        subroutine h2dexpta(nd,rk,rscale,center,mploc,nterms,
     $     expinc,nexp,xs,ws)
        implicit real *8 (a-h,o-z)
        dimension center(2)
        complex *16 mploc(nd,-nterms:nterms),rk
        complex *16 expinc(nd,0:nexp,4)
        complex *16 xs(nexp),ws(nexp)
        complex *16 khat,zmul1,zmul2,ztemp1,ztemp2
        complex *16 s1,s2,s3,s4,r1,r2,r3,r4
        complex *16 zmuls1(nexp,nterms),zmuls2(nexp,nterms)
        complex *16 ima
        data ima/(0.0d0,1.0d0)/


        if( 2 .eq. 2 ) then

        do k = 1,nexp

        khat = (sqrt(xs(k)**2-rk**2)-xs(k))/rk
        if( k .gt. nexp/2) khat = -rk/(sqrt(xs(k)**2-rk**2)+xs(k))

        ztemp1 = khat*rscale
        ztemp2 = rscale/khat
        zmul1 = ztemp1
        zmul2 = ztemp2
        do n = 1,nterms-3,4
           zmuls1(k,n+0)=zmul1
           zmul1 = zmul1*ztemp1
           zmuls1(k,n+1)=zmul1
           zmul1 = zmul1*ztemp1
           zmuls1(k,n+2)=zmul1
           zmul1 = zmul1*ztemp1
           zmuls1(k,n+3)=zmul1
           zmul1 = zmul1*ztemp1
        enddo
        do n = 1,nterms-3,4
           zmuls2(k,n+0)=zmul2
           zmul2 = zmul2*ztemp2
           zmuls2(k,n+1)=zmul2
           zmul2 = zmul2*ztemp2
           zmuls2(k,n+2)=zmul2
           zmul2 = zmul2*ztemp2
           zmuls2(k,n+3)=zmul2
           zmul2 = zmul2*ztemp2
        enddo

        enddo


        do k = 1,nexp
        do idim=1,nd
        mploc(idim,0) = mploc(idim,0) + expinc(idim,k,1)
     $     + expinc(idim,k,2) + expinc(idim,k,3) + expinc(idim,k,4)
        enddo
        do idim=1,nd
        call h2dfft4f(expinc(idim,k,3),expinc(idim,k,4),
     $     expinc(idim,k,1),expinc(idim,k,2),
     $     s1, s2, s3, s4)
        do n = 1,nterms-3,4
            mploc(idim,n+0) = mploc(idim,n+0) + zmuls1(k,n+0) * s2
            mploc(idim,n+1) = mploc(idim,n+1) + zmuls1(k,n+1) * s3
            mploc(idim,n+2) = mploc(idim,n+2) + zmuls1(k,n+2) * s4
            mploc(idim,n+3) = mploc(idim,n+3) + zmuls1(k,n+3) * s1
        enddo
        call h2dfft4b(expinc(idim,k,3),expinc(idim,k,4),
     $     expinc(idim,k,1),expinc(idim,k,2),
     $     s1, s2, s3, s4)
        do n = 1,nterms-3,4
            mploc(idim,-n-0) = mploc(idim,-n-0) + zmuls2(k,n+0) * s2
            mploc(idim,-n-1) = mploc(idim,-n-1) + zmuls2(k,n+1) * s3
            mploc(idim,-n-2) = mploc(idim,-n-2) + zmuls2(k,n+2) * s4
            mploc(idim,-n-3) = mploc(idim,-n-3) + zmuls2(k,n+3) * s1
        enddo
        enddo
        enddo

c
c       ... constant terms
c

        do idim=1,nd
        mploc(idim,0) = mploc(idim,0) + expinc(idim,0,1)
     $     + expinc(idim,0,2) + expinc(idim,0,3) + expinc(idim,0,4)
        enddo

        endif


        if( 1 .eq. 2 ) then

        do k = 1,nexp

        khat = (sqrt(xs(k)**2-rk**2)-xs(k))/rk
        if( k .gt. nexp/2) khat = -rk/(sqrt(xs(k)**2-rk**2)+xs(k))

        do idim=1,nd
        mploc(idim,0) = mploc(idim,0) + expinc(idim,k,1)
     $     + expinc(idim,k,2) + expinc(idim,k,3) + expinc(idim,k,4)
        enddo
        ztemp1 = khat*rscale
        ztemp2 = rscale/khat
        zmul1 = ztemp1
        do idim=1,nd
        call h2dfft4f(expinc(idim,k,3),expinc(idim,k,4),
     $     expinc(idim,k,1),expinc(idim,k,2),
     $     s1, s2, s3, s4)
        do n = 1,nterms-3,4
            mploc(idim,n+0) = mploc(idim,n+0) + zmul1 * s2
            zmul1 = zmul1*ztemp1
            mploc(idim,n+1) = mploc(idim,n+1) + zmul1 * s3
            zmul1 = zmul1*ztemp1
            mploc(idim,n+2) = mploc(idim,n+2) + zmul1 * s4
            zmul1 = zmul1*ztemp1
            mploc(idim,n+3) = mploc(idim,n+3) + zmul1 * s1
            zmul1 = zmul1*ztemp1
        enddo
        zmul2 = ztemp2
        call h2dfft4b(expinc(idim,k,3),expinc(idim,k,4),
     $     expinc(idim,k,1),expinc(idim,k,2),
     $     s1, s2, s3, s4)
        do n = 1,nterms-3,4
            mploc(idim,-n-0) = mploc(idim,-n-0) + zmul2 * s2
            zmul2 = zmul2*ztemp2
            mploc(idim,-n-1) = mploc(idim,-n-1) + zmul2 * s3
            zmul2 = zmul2*ztemp2
            mploc(idim,-n-2) = mploc(idim,-n-2) + zmul2 * s4
            zmul2 = zmul2*ztemp2
            mploc(idim,-n-3) = mploc(idim,-n-3) + zmul2 * s1
            zmul2 = zmul2*ztemp2
        enddo
        enddo
        enddo

c
c       ... constant terms
c

        do idim=1,nd
        mploc(idim,0) = mploc(idim,0) + expinc(idim,0,1)
     $     + expinc(idim,0,2) + expinc(idim,0,3) + expinc(idim,0,4)
        enddo

        endif

        return
        end
c
c
c
c
c
        subroutine h2dexpexp_e(nd,expout,expinc,nexp,cshift,m,n)
        implicit real *8 (a-h,o-z)
        complex *16 expout(nd,0:nexp,4)
        complex *16 expinc(nd,0:nexp,4)
        complex *16 cshift(1:nexp,-3:3,-3:3)

        do idim=1,nd
        expinc(idim,0,1) = expinc(idim,0,1) + expout(idim,0,1)
        enddo
        do i = 1,nexp
        do idim=1,nd
        expinc(idim,i,1) = expinc(idim,i,1)
     $     + expout(idim,i,1) * cshift(i,m,n)
        enddo
        enddo

        return
        end
c
c
c
c
c
        subroutine h2dexpexp_n(nd,expout,expinc,nexp,cshift,m,n)
        implicit real *8 (a-h,o-z)
        complex *16 expout(nd,0:nexp,4)
        complex *16 expinc(nd,0:nexp,4)
        complex *16 cshift(1:nexp,-3:3,-3:3)

        do idim=1,nd
        expinc(idim,0,2) = expinc(idim,0,2) + expout(idim,0,2)
        enddo
        do i = 1,nexp
        do idim=1,nd
        expinc(idim,i,2) = expinc(idim,i,2)
     $     + expout(idim,i,2) * cshift(i,m,n)
        enddo
        enddo

        return
        end
c
c
c
c
c
        subroutine h2dexpexp_w(nd,expout,expinc,nexp,cshift,m,n)
        implicit real *8 (a-h,o-z)
        complex *16 expout(nd,0:nexp,4)
        complex *16 expinc(nd,0:nexp,4)
        complex *16 cshift(1:nexp,-3:3,-3:3)

        do idim=1,nd
        expinc(idim,0,3) = expinc(idim,0,3) + expout(idim,0,3)
        enddo
        do i = 1,nexp
        do idim=1,nd
        expinc(idim,i,3) = expinc(idim,i,3)
     $     + expout(idim,i,3) * cshift(i,m,n)
        enddo
        enddo

        return
        end
c
c
c
c
c
        subroutine h2dexpexp_s(nd,expout,expinc,nexp,cshift,m,n)
        implicit real *8 (a-h,o-z)
        complex *16 expout(nd,0:nexp,4)
        complex *16 expinc(nd,0:nexp,4)
        complex *16 cshift(1:nexp,-3:3,-3:3)

        do idim=1,nd
        expinc(idim,0,4) = expinc(idim,0,4) + expout(idim,0,4)
        enddo
        do i = 1,nexp
        do idim=1,nd
        expinc(idim,i,4) = expinc(idim,i,4)
     $     + expout(idim,i,4) * cshift(i,m,n)
        enddo
        enddo

        return
        end
c
c
c
c
c
