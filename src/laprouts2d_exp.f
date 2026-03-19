        subroutine l2expini(nterms,xs,ws,nexp,xsum,comp,cshift)
        implicit real *8 (a-h,o-z)
        complex *16 cshift(nexp,-3:3,-3:3)
        dimension xs(nexp),ws(nexp)
        dimension comp(0:nterms,nexp)
c       
c       ... precompute exponential translation tables
c
        do k = 1,nexp
        do i = -3,3
        do j = -3,3
        cshift(k,i,j) = exp(-xs(k)*dcmplx(1.0d0*i,1.0d0*j))
        enddo
        enddo
        enddo

        xsum = 0.0d0
        do i = 1,nexp
        xsum = xsum + exp(-xs(i))*ws(i)/xs(i)
        enddo

        do i = 1,nexp
        comp(0,i) = 1.0d0
        do j = 1,nterms
            comp(j,i) = comp(j-1,i) * xs(i)/j
        enddo
        enddo

        return
        end
c
c        
c
c
c
        subroutine l2dfft4f(r1,r2,r3,r4,s1,s2,s3,s4)
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
        subroutine l2dmpexp(rscale,center,hexp,nterms,
     $     expout,nexp,xs,ws,xsum,comp)
        implicit real *8 (a-h,o-z)
        dimension center(3)
        complex *16 hexp(0:nterms)
        complex *16 expout(0:nexp,4)
        dimension xs(nexp),ws(nexp),comp(0:nterms,nexp)
        complex *16 r1,r2,r3,r4,s1,s2,s3,s4
        complex *16 cd, ima
        data ima/(0.0d0,1.0d0)/
        
        if( 2 .eq. 2 ) then
c
c       ... 1/z, 1/z^2, 1/z^3, ..., 1/z^n terms
c
        do i = 1,nexp
        s1 = dcmplx(0.0d0,0.0d0)
        s2 = dcmplx(0.0d0,0.0d0)
        s3 = dcmplx(0.0d0,0.0d0)
        s4 = dcmplx(0.0d0,0.0d0)
        do j = 1,nterms-3,4
            s1 = s1 + comp(j-1+3,i) * hexp(j+3)
            s2 = s2 + comp(j-1+0,i) * hexp(j+0)
            s3 = s3 + comp(j-1+1,i) * hexp(j+1)
            s4 = s4 + comp(j-1+2,i) * hexp(j+2)
        enddo
        call l2dfft4f(s1, s2, s3, s4, r1, r2, r3, r4)
        expout(i,1) = expout(i,1) + ws(i) * r1
        expout(i,2) = expout(i,2) + ws(i) * r2
        expout(i,3) = expout(i,3) + ws(i) * r3
        expout(i,4) = expout(i,4) + ws(i) * r4
        enddo
c
c       ... log(z) term
c
        do i = 1,nexp
        cd = ws(i) / xs(i) * hexp(0)
        expout(i,1) = expout(i,1) - cd
        expout(i,2) = expout(i,2) - cd
        expout(i,3) = expout(i,3) - cd
        expout(i,4) = expout(i,4) - cd
        enddo

        cd = hexp(0) * (xsum + log(rscale))
        expout(0,1) = expout(0,1) + cd
        expout(0,2) = expout(0,2) + cd
        expout(0,3) = expout(0,3) + cd
        expout(0,4) = expout(0,4) + cd

        endif

        if( 1 .eq. 2 ) then
c
c       ... 1/z, 1/z^2, 1/z^3, ..., 1/z^n terms
c
        do j = 1,nterms
        do i = 1,nexp
        expout(i,1) = expout(i,1) +
     $     ws(i) * comp(j-1,i) * hexp(j) * (+1)**j
        expout(i,2) = expout(i,2) +
     $     ws(i) * comp(j-1,i) * hexp(j) * (-ima)**j
        expout(i,3) = expout(i,3) +
     $     ws(i) * comp(j-1,i) * hexp(j) * (-1)**j
        expout(i,4) = expout(i,4) +
     $     ws(i) * comp(j-1,i) * hexp(j) * (+ima)**j
        enddo
        enddo
c
c       ... log(z) term
c
        do i = 1,nexp
        expout(i,1) = expout(i,1) - ws(i) / xs(i) * hexp(0)
        expout(i,2) = expout(i,2) - ws(i) / xs(i) * hexp(0)
        expout(i,3) = expout(i,3) - ws(i) / xs(i) * hexp(0)
        expout(i,4) = expout(i,4) - ws(i) / xs(i) * hexp(0)
        enddo

        expout(0,1) = expout(0,1) + hexp(0) * (xsum + log(rscale))
        expout(0,2) = expout(0,2) + hexp(0) * (xsum + log(rscale))
        expout(0,3) = expout(0,3) + hexp(0) * (xsum + log(rscale))
        expout(0,4) = expout(0,4) + hexp(0) * (xsum + log(rscale))

        endif

        if ( 1 .eq. 2 ) then
c
c       ... 1/z, 1/z^2, 1/z^3, ..., 1/z^n terms
c
        do i = 1,nexp
        do j = 1,nterms-3,4

        call l2dfft4f(
     $     comp(j-1+3,i) * hexp(j+3),comp(j-1+0,i) * hexp(j+0),
     $     comp(j-1+1,i) * hexp(j+1),comp(j-1+2,i) * hexp(j+2),
     $     s1, s2, s3, s4)
        expout(i,1) = expout(i,1) + ws(i) * s1
        expout(i,2) = expout(i,2) + ws(i) * s2
        expout(i,3) = expout(i,3) + ws(i) * s3
        expout(i,4) = expout(i,4) + ws(i) * s4

        enddo
        enddo
c
c       ... log(z) term
c
        do i = 1,nexp
        cd = ws(i) / xs(i) * hexp(0)
        expout(i,1) = expout(i,1) - cd
        expout(i,2) = expout(i,2) - cd
        expout(i,3) = expout(i,3) - cd
        expout(i,4) = expout(i,4) - cd
        enddo

        cd = hexp(0) * (xsum + log(rscale))
        expout(0,1) = expout(0,1) + cd
        expout(0,2) = expout(0,2) + cd
        expout(0,3) = expout(0,3) + cd
        expout(0,4) = expout(0,4) + cd
        
        endif

        if( 1 .eq. 2 ) then 
c
c       ... 1/z, 1/z^2, 1/z^3, ..., 1/z^n terms
c
        do i = 1,nexp
        s1 = dcmplx(0.0d0,0.0d0)
        s2 = dcmplx(0.0d0,0.0d0)
        s3 = dcmplx(0.0d0,0.0d0)
        s4 = dcmplx(0.0d0,0.0d0)
        do j = 1,nterms-3,4
        call l2dfft4f(
     $     comp(j-1+3,i) * hexp(j+3), comp(j-1+0,i) * hexp(j+0),
     $     comp(j-1+1,i) * hexp(j+1), comp(j-1+2,i) * hexp(j+2),
     $     r1, r2, r3, r4)
        s1 = s1 + r1
        s2 = s2 + r2
        s3 = s3 + r3
        s4 = s4 + r4
        enddo
        expout(i,1) = expout(i,1) + ws(i) * s1
        expout(i,2) = expout(i,2) + ws(i) * s2
        expout(i,3) = expout(i,3) + ws(i) * s3
        expout(i,4) = expout(i,4) + ws(i) * s4
        enddo
c
c       ... log(z) term
c
        do i = 1,nexp
        cd = ws(i) / xs(i) * hexp(0)
        expout(i,1) = expout(i,1) - cd
        expout(i,2) = expout(i,2) - cd
        expout(i,3) = expout(i,3) - cd
        expout(i,4) = expout(i,4) - cd
        enddo

        cd = hexp(0) * (xsum + log(rscale))
        expout(0,1) = expout(0,1) + cd
        expout(0,2) = expout(0,2) + cd
        expout(0,3) = expout(0,3) + cd
        expout(0,4) = expout(0,4) + cd

        endif

        return
        end
c
c
c
c
c
        subroutine l2dexpta(rscale,center,jexp,nterms,
     $     expinc,nexp,comp)
        implicit real *8 (a-h,o-z)
        dimension center(3)
        complex *16 jexp(0:nterms)
        complex *16 expinc(0:nexp,4)
        dimension comp(0:nterms,nexp)
        complex *16 s1,s2,s3,s4
        complex *16 cd, ima
        data ima/(0.0d0,1.0d0)/

        if( 2 .eq. 2 ) then
c
c       ... exp(spin*s*xk(i)) Taylor expansions
c
        do i = 1,nexp
        call l2dfft4f(
     $     expinc(i,3), expinc(i,4), expinc(i,1), expinc(i,2),
     $     s1, s2, s3, s4)
        do j = 0,nterms-3,4
            jexp(j+0) = jexp(j+0) + comp(j+0,i) * s1
            jexp(j+1) = jexp(j+1) + comp(j+1,i) * s2
            jexp(j+2) = jexp(j+2) + comp(j+2,i) * s3
            jexp(j+3) = jexp(j+3) + comp(j+3,i) * s4
        enddo
        enddo
c
c       ... constant terms
c
        jexp(0) = jexp(0) + expinc(0,1)
        jexp(0) = jexp(0) + expinc(0,2)
        jexp(0) = jexp(0) + expinc(0,3)
        jexp(0) = jexp(0) + expinc(0,4)
        
        endif

        if( 1 .eq. 2 ) then
c
c       ... exp(spin*s*xk(i)) Taylor expansions
c
        do j = 0,nterms
        do i = 1,nexp
        jexp(j) = jexp(j) + expinc(i,1) * comp(j,i) * (-1)**j
        jexp(j) = jexp(j) + expinc(i,2) * comp(j,i) * (+ima)**j
        jexp(j) = jexp(j) + expinc(i,3) * comp(j,i) * (+1)**j
        jexp(j) = jexp(j) + expinc(i,4) * comp(j,i) * (-ima)**j
        enddo
        enddo
c
c       ... constant terms
c
        jexp(0) = jexp(0) + expinc(0,1)
        jexp(0) = jexp(0) + expinc(0,2)
        jexp(0) = jexp(0) + expinc(0,3)
        jexp(0) = jexp(0) + expinc(0,4)
        
        endif

        return
        end
c
c
c
c
c
        subroutine l2dexpexp_e(expout,expinc,nexp,cshift,m,n)
        implicit real *8 (a-h,o-z)
        complex *16 expout(0:nexp,4)
        complex *16 expinc(0:nexp,4)
        complex *16 cshift(1:nexp,-3:3,-3:3)
        
        expinc(0,1) = expinc(0,1) + expout(0,1)
        do i = 1,nexp
        expinc(i,1) = expinc(i,1) + expout(i,1) * cshift(i,m,n)
        enddo

        return
        end
c
c
c
c
c
        subroutine l2dexpexp_n(expout,expinc,nexp,cshift,m,n)
        implicit real *8 (a-h,o-z)
        complex *16 expout(0:nexp,4)
        complex *16 expinc(0:nexp,4)
        complex *16 cshift(1:nexp,-3:3,-3:3)
        
        expinc(0,2) = expinc(0,2) + expout(0,2)
        do i = 1,nexp
        expinc(i,2) = expinc(i,2) + expout(i,2) * cshift(i,m,n)
        enddo

        return
        end
c
c
c
c
c
        subroutine l2dexpexp_w(expout,expinc,nexp,cshift,m,n)
        implicit real *8 (a-h,o-z)
        complex *16 expout(0:nexp,4)
        complex *16 expinc(0:nexp,4)
        complex *16 cshift(1:nexp,-3:3,-3:3)
        
        expinc(0,3) = expinc(0,3) + expout(0,3)
        do i = 1,nexp
        expinc(i,3) = expinc(i,3) + expout(i,3) * cshift(i,m,n)
        enddo

        return
        end
c
c
c
c
c
        subroutine l2dexpexp_s(expout,expinc,nexp,cshift,m,n)
        implicit real *8 (a-h,o-z)
        complex *16 expout(0:nexp,4)
        complex *16 expinc(0:nexp,4)
        complex *16 cshift(1:nexp,-3:3,-3:3)
        
        expinc(0,4) = expinc(0,4) + expout(0,4)
        do i = 1,nexp
        expinc(i,4) = expinc(i,4) + expout(i,4) * cshift(i,m,n)
        enddo

        return
        end
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine l2expini_sep2(nterms,xs,ws,nexp,xsum,comp,cshift)
        implicit real *8 (a-h,o-z)
        complex *16 cshift(nexp,-5:5,-5:5)
        dimension xs(nexp),ws(nexp)
        dimension comp(0:nterms,nexp)
c       
c       ... precompute exponential translation tables
c
        do k = 1,nexp
        do i = -5,5
        do j = -5,5
        cshift(k,i,j) = exp(-xs(k)*dcmplx(1.0d0*i,1.0d0*j))
        enddo
        enddo
        enddo

        xsum = 0.0d0
        do i = 1,nexp
        xsum = xsum + exp(-xs(i))*ws(i)/xs(i)
        enddo

        do i = 1,nexp
        comp(0,i) = 1.0d0
        do j = 1,nterms
            comp(j,i) = comp(j-1,i) * xs(i)/j
        enddo
        enddo

        return
        end
c
c        
c
c
c
        subroutine l2dexpexp_sep2_e(expout,expinc,nexp,cshift,m,n)
        implicit real *8 (a-h,o-z)
        complex *16 expout(0:nexp,4)
        complex *16 expinc(0:nexp,4)
        complex *16 cshift(1:nexp,-5:5,-5:5)
        
        expinc(0,1) = expinc(0,1) + expout(0,1)
        do i = 1,nexp
        expinc(i,1) = expinc(i,1) + expout(i,1) * cshift(i,m,n)
        enddo

        return
        end
c
c
c
c
c
        subroutine l2dexpexp_sep2_n(expout,expinc,nexp,cshift,m,n)
        implicit real *8 (a-h,o-z)
        complex *16 expout(0:nexp,4)
        complex *16 expinc(0:nexp,4)
        complex *16 cshift(1:nexp,-5:5,-5:5)
        
        expinc(0,2) = expinc(0,2) + expout(0,2)
        do i = 1,nexp
        expinc(i,2) = expinc(i,2) + expout(i,2) * cshift(i,m,n)
        enddo

        return
        end
c
c
c
c
c
        subroutine l2dexpexp_sep2_w(expout,expinc,nexp,cshift,m,n)
        implicit real *8 (a-h,o-z)
        complex *16 expout(0:nexp,4)
        complex *16 expinc(0:nexp,4)
        complex *16 cshift(1:nexp,-5:5,-5:5)
        
        expinc(0,3) = expinc(0,3) + expout(0,3)
        do i = 1,nexp
        expinc(i,3) = expinc(i,3) + expout(i,3) * cshift(i,m,n)
        enddo

        return
        end
c
c
c
c
c
        subroutine l2dexpexp_sep2_s(expout,expinc,nexp,cshift,m,n)
        implicit real *8 (a-h,o-z)
        complex *16 expout(0:nexp,4)
        complex *16 expinc(0:nexp,4)
        complex *16 cshift(1:nexp,-5:5,-5:5)
        
        expinc(0,4) = expinc(0,4) + expout(0,4)
        do i = 1,nexp
        expinc(i,4) = expinc(i,4) + expout(i,4) * cshift(i,m,n)
        enddo

        return
        end
