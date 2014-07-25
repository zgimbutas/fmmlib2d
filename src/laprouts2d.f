cc Copyright (C) 2010-2011: Leslie Greengard and Zydrunas Gimbutas
cc Contact: greengard@cims.nyu.edu
cc 
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date$
c    $Revision$
c
c
c      This file contains the basic subroutines for 
c      forming and evaluating multipole (partial wave) expansions
c      in two dimensions.
c
c      Since log(z) is a multivalued complex function, we use
c      the real part Re(log(z)) = log(abs(z)) in all computations.
c
c      All multipole and local expansions are properly scaled 
c
c-----------------------------------------------------------------------
c
c      L2DFORMMP: creates multipole expansion (outgoing) due to 
c                 a collection of charge sources.
c
c      L2DFORMTA: creates local expansion due to 
c                 a collection of charge sources.
c
c      L2DFORMMP_DP: creates multipole expansion (outgoing) due to 
c                 a collection of dipole sources.
c
c      L2DFORMTA_DP: creates local expansion due to 
c                 a collection of dipole sources.
c
c-----------------------------------------------------------------------
c
c      Multipole and local translation operators
c
c      L2DMPMP:     converts multipole expansion to a multipole expansion.
c      L2DMPLOC:     converts multipole expansion to a local expansion.
c      L2DLOCLOC:     converts local expansion to a local expansion.
c
c-----------------------------------------------------------------------
c
c      Multipole and local translation operators with precomputed storage
c
c      L2DMPMP_CARRAY:   converts multipole expansion to a multipole expansion.
c      L2DMPLOC_CARRAY:   converts multipole expansion to a local expansion.
c      L2DLOCLOC_CARRAY:   converts local expansion to a local expansion.
c
c-----------------------------------------------------------------------
c       
c      Complex valued Cauchy sums
c       __ depreciated functions __
c
c      L2DMPEVALALL: computes potential and grad(potential)
c                 due to a multipole expansion for a collection of targets
c      L2DMPEVAL: computes potential and grad(potential)
c                 due to a multipole expansion.
c
c      L2DTAEVALALL: computes potential and grad(potential) 
c                  due to local expansion for a collection of targets
c      L2DTAEVAL: computes potential and grad(potential) 
c                  due to local expansion.
c
c      LPOTGRAD2DALL:  direct calculation for a collection of charge sources
c      LPOTGRAD2D : direct calculation for a single charge source
c
c      LPOTGRAD2DALL_DP:  direct calculation for a collection of dipoles
c      LPOTGRAD2D_DP : direct calculation for a single dipole
c
c
c      LPOTGRAD2DALL_SDP:  direct calculation for 
c                 a collection of charges and dipoles
c      LPOTGRAD2D_SDP : direct calculation for a single charge and a dipole
c
c-----------------------------------------------------------------------
c
c      Complex valued Cauchy sums
c
c      C2DMPEVALALL: computes potential and grad(potential)
c                 due to a multipole expansion for a collection of targets
c      C2DMPEVAL: computes potential and grad(potential)
c                 due to a multipole expansion.
c
c      C2DTAEVALALL: computes potential and grad(potential) 
c                  due to local expansion for a collection of targets
c      C2DTAEVAL: computes potential and grad(potential) 
c                  due to local expansion.
c
c      CPOTGRAD2DALL_SDP:  direct calculation for 
c                 a collection of charges and dipoles
c      CPOTGRAD2D_SDP : direct calculation for a single charge and a dipole
c
c      CPOTGRAD2D_SDP_SYM : direct calculation for
c      a pair of charges and a dipoles, uses symmetries
c
c-----------------------------------------------------------------------
c
c      Complex valued Laplace potentials
c
c      RCPOTGRAD2DALL_SDP:  direct calculation for 
c                 a collection of charges and dipoles
c      RCPOTGRAD2D_SDP : direct calculation for a single charge and a dipole
c
c-----------------------------------------------------------------------
c
c      Real valued Laplace potentials
c
c      RPOTGRAD2DALL_SDP:  direct calculation for 
c                 a collection of charges and dipoles
c      RPOTGRAD2D_SDP : direct calculation for a single charge and a dipole
c
c
c
c**********************************************************************
      subroutine l2dmpeval(rscale,center,mpole,nterms,ztarg,
     1      pot,ifgrad,grad,ifhess,hess)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient of the 
c     potential due to an outgoing partial wave expansion.
c               +nterms
c     pot  =      sum   mpole_n / z^n  + mpole_0 log(abs(z))
c                n=1  
c     grad  = gradient(pot) if ifgrad = 1.
c     hess = hessian if ifhess = 1.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale :    scaling parameter 
c     center :    expansion center
c     mpole  :    multipole expansion 
c     nterms :    order of the multipole expansion
c     ztarg  :    target location
c     ifgrad :   flag controlling evaluation of gradient:
c                   ifgrad = 0, do not compute gradient.
c                   ifgrad = 1, compute gradient.
c     ifhess :   flag for computing Hessian:
c	            ifhess = 0 -> don't compute 
c		    ifhess = 1 -> do compute 
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    potential at ztarg
c     grad   :    gradient at ztarg (if requested)
c     hess   :    hessian at ztarg (if requested)
c
c-----------------------------------------------------------------------
c
        complex *16 pot,grad(2),hess(3),mpole(0:nterms)
        real *8 center(2),ztarg(2),zdiff(2)
        complex *16 ztemp1, ztemp2
ccc        complex *16, allocatable :: z0pow(:), z0powm(:)
        complex *16 z0pow(0:1000)
c
        complex *16 ima,z,cd
        data ima/(0.0d0,1.0d0)/
c
c
        zdiff(1)=ztarg(1)-center(1)
        zdiff(2)=ztarg(2)-center(2)
ccc        call h2cart2polar(zdiff,r,theta)
        z=dcmplx(zdiff(1),zdiff(2))
c
c
        nmax = nterms + 3
ccc        allocate( z0pow(0:nmax) )
        ztemp1=rscale/z
        ztemp2=ztemp1
        z0pow(0)=1
        do i=1,nmax
        z0pow(i)=ztemp1
        ztemp1=ztemp1*ztemp2
        enddo
c
c
        pot=mpole(0)*log(abs(z))
        do n=1,nterms
        cd=mpole(n)*z0pow(n)
        pot=pot+cd
        enddo
c
c
        if( ifgrad .eq. 1 ) then

        rinv=1/rscale
        cd=mpole(0)*z0pow(1)
        grad(1)=cd
        grad(2)=cd*ima
        do n=1,nterms
        cd=-mpole(n)*z0pow(n+1)*n
        grad(1)=grad(1)+cd
        grad(2)=grad(2)+cd*ima 
        enddo
        grad(1)=grad(1)*rinv
        grad(2)=grad(2)*rinv

        endif
c
c
        if( ifhess .eq. 1 ) then

        rinv2=1/rscale**2
        cd=-mpole(0)*z0pow(2)
        hess(1)=cd*(1)
        hess(2)=cd*(ima)
        hess(3)=cd*(-1)
        do n=1,nterms
        cd=mpole(n)*z0pow(n+2)*n*(n+1)
        hess(1)=hess(1)+cd*(1) 
        hess(2)=hess(2)+cd*(ima)
        hess(3)=hess(3)+cd*(-1) 
        enddo
        hess(1)=hess(1)*rinv2
        hess(2)=hess(2)*rinv2
        hess(3)=hess(3)*rinv2

        endif
c
c
        return
        end
c
c
c
c
c
c**********************************************************************
      subroutine l2dtaeval(rscale,center,mpole,nterms,ztarg,
     1      pot,ifgrad,grad,ifhess,hess)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient of the 
c     potential due to an incoming partial wave expansion.
c               +nterms
c     pot  =      sum   mpole_n * z^n  + mpole_0 
c                n=1  
c     grad  = gradient(pot) if ifgrad = 1.
c     hess = hessian if ifhess = 1.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale :    scaling parameter 
c     center :    expansion center
c     mpole  :    multipole expansion 
c     nterms :    order of the local expansion
c     ztarg  :    target location
c     ifgrad :   flag controlling evaluation of gradient:
c                   ifgrad = 0, do not compute gradient.
c                   ifgrad = 1, compute gradient.
c     ifhess :   flag for computing Hessian:
c	            ifhess = 0 -> don't compute 
c		    ifhess = 1 -> do compute 
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    potential at ztarg
c     grad   :    gradient at ztarg (if requested)
c     hess   :    hessian at ztarg (if requested)
c
c-----------------------------------------------------------------------
c
        complex *16 pot,grad(2),hess(3),mpole(0:nterms)
        real *8 center(2),ztarg(2),zdiff(2)
        complex *16 ztemp1, ztemp2
ccc        complex *16, allocatable :: z0pow(:), z0powm(:)
        complex *16 z0pow(0:1000)
c
        complex *16 ima,z,cd
        data ima/(0.0d0,1.0d0)/
c
c
        zdiff(1)=ztarg(1)-center(1)
        zdiff(2)=ztarg(2)-center(2)
ccc        call h2cart2polar(zdiff,r,theta)
        z=dcmplx(zdiff(1),zdiff(2))
c
c
c
        nmax = nterms
ccc        allocate( z0pow(0:nmax) )
        ztemp1=z/rscale
        ztemp2=ztemp1
        z0pow(0)=1
        do i=1,nmax
        z0pow(i)=ztemp1
        ztemp1=ztemp1*ztemp2
        enddo
c
c
c
        pot=mpole(0)
        do n=1,nterms
        cd=mpole(n)*z0pow(n)
        pot=pot+cd
        enddo
c
c
        if( ifgrad .eq. 1 ) then

        grad(1)=0
        grad(2)=0

        rinv=1/rscale
        do n=1,nterms
        cd=mpole(n)*z0pow(n-1)*n
        grad(1)=grad(1)+cd
        grad(2)=grad(2)+cd*ima
        enddo
        grad(1)=grad(1)*rinv
        grad(2)=grad(2)*rinv

        endif
c
c
        if( ifhess .eq. 1 ) then

        hess(1)=0
        hess(2)=0
        hess(3)=0

        rinv2=1/rscale**2
        do n=2,nterms
        cd=mpole(n)*z0pow(n-2)*n*(n-1)
        hess(1)=hess(1)+cd*(1) 
        hess(2)=hess(2)+cd*(ima)
        hess(3)=hess(3)+cd*(-1)
        enddo
        hess(1)=hess(1)*rinv2
        hess(2)=hess(2)*rinv2
        hess(3)=hess(3)*rinv2

        endif
c
c
        return
        end
c
c
c
c
c
c**********************************************************************
      subroutine l2dmpevalall(rscale,center,mpole,nterms,ztarg,ntarg,
     1      ifpot,pot,ifgrad,grad,ifhess,hess)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient of the 
c     potential due to an outgoing partial wave expansion.
c               +nterms
c     pot  =      sum   mpole_n / z^n  + mpole_0 log(abs(z))
c                n=1  
c     grad  = gradient(pot) if ifgrad = 1.
c     hess = hessian if ifhess = 1.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale :    scaling parameter 
c     center :    expansion center
c     mpole  :    multipole expansion 
c     nterms :    order of the multipole expansion
c     ztarg  :    target location
c     ifgrad :   flag controlling evaluation of gradient:
c                   ifgrad = 0, do not compute gradient.
c                   ifgrad = 1, compute gradient.
c     ifhess :   flag for computing Hessian:
c	            ifhess = 0 -> don't compute 
c		    ifhess = 1 -> do compute 
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    potential at ztarg
c     grad   :    gradient at ztarg (if requested)
c     hess   :    hessian at ztarg (if requested)
c
c-----------------------------------------------------------------------
c
        complex *16 pot(1),grad(2,1),hess(3,1),mpole(0:nterms)
        complex *16 potloc,gradloc(2),hessloc(3)
        real *8 center(2),ztarg(2,ntarg),zdiff(2)
c
        complex *16 ima,z
        data ima/(0.0d0,1.0d0)/
c
        do i = 1,ntarg
        call l2dmpeval(rscale,center,mpole,nterms,ztarg(1,i),
     1     potloc,ifgrad,gradloc,ifhess,hessloc)
        if (ifpot.eq.1) pot(i) = pot(i) + potloc
        if (ifgrad.eq.1) then
        grad(1,i) = grad(1,i) + gradloc(1)
        grad(2,i) = grad(2,i) + gradloc(2)
        endif
        if (ifhess.eq.1) then
        hess(1,i) = hess(1,i) + hessloc(1)
        hess(2,i) = hess(2,i) + hessloc(2)
        hess(3,i) = hess(3,i) + hessloc(3)
        endif
        enddo
c
        return
        end
c
c
c
c
c
c**********************************************************************
      subroutine l2dtaevalall(rscale,center,mpole,nterms,ztarg,ntarg,
     1      ifpot,pot,ifgrad,grad,ifhess,hess)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient of the 
c     potential due to an incoming partial wave expansion.
c               +nterms
c     pot  =      sum   mpole_n * z^n  + mpole_0 
c                n=1  
c     grad  = gradient(pot) if ifgrad = 1.
c     hess = hessian if ifhess = 1.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale :    scaling parameter 
c     center :    expansion center
c     mpole  :    multipole expansion 
c     nterms :    order of the local expansion
c     ztarg  :    target location
c     ifgrad :   flag controlling evaluation of gradient:
c                   ifgrad = 0, do not compute gradient.
c                   ifgrad = 1, compute gradient.
c     ifhess :   flag for computing Hessian:
c	            ifhess = 0 -> don't compute 
c		    ifhess = 1 -> do compute 
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    potential at ztarg
c     grad   :    gradient at ztarg (if requested)
c     hess   :    hessian at ztarg (if requested)
c
c-----------------------------------------------------------------------
c
        complex *16 pot(1),grad(2,1),hess(3,1),mpole(0:nterms)
        complex *16 potloc,gradloc(2),hessloc(3)
        real *8 center(2),ztarg(2,ntarg),zdiff(2)
c
        complex *16 ima,z
        data ima/(0.0d0,1.0d0)/
c
        do i = 1,ntarg
        call l2dtaeval(rscale,center,mpole,nterms,ztarg(1,i),
     1     potloc,ifgrad,gradloc,ifhess,hessloc)
        if (ifpot.eq.1) pot(i) = pot(i) + potloc
        if (ifgrad.eq.1) then
        grad(1,i) = grad(1,i) + gradloc(1)
        grad(2,i) = grad(2,i) + gradloc(2)
        endif
        if (ifhess.eq.1) then
        hess(1,i) = hess(1,i) + hessloc(1)
        hess(2,i) = hess(2,i) + hessloc(2)
        hess(3,i) = hess(3,i) + hessloc(3)
        endif
        enddo
c
        return
        end
c
c
c
c
c
C***********************************************************************
        subroutine l2dformmp(ier,rscale,source,charge,ns,center,
     1                       nterms,mpole)
        implicit real *8 (a-h,o-z)
C***********************************************************************
c
c     This subroutine constructs a multipole expansion about CENTER due
c     to NS sources located at SOURCES(2,*).
c
c     mpole_0  =  sum charge_j 
c                  j  
c
c     mpole_n  = -sum charge_j 1/n (z)^n /rscale^n
c                  j  
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     charge(ns)      : source strengths
c     ns              : number of sources
c     center(2)       : expansion center
c     nterms          : order of multipole expansion
c
c     OUTPUT:
c
c     ier       : error return code
c                 ier=0 returned successfully;
c
c     mpole     : coeffs for the multipole-expansion
c
        complex *16 mpole(0:nterms),charge(ns)
        real *8 center(2),source(2,ns),zdiff(2)

        complex *16 zmul,zinv,ztemp1,ztemp2
c
        complex *16 ima,z,z0
        data ima/(0.0d0,1.0d0)/
c
c
        do n=0,nterms
        mpole(n)=0
        enddo

        do j=1,ns
c
        zdiff(1)=source(1,j)-center(1)
        zdiff(2)=source(2,j)-center(2)
ccc        call h2cart2polar(zdiff,r,theta)
        z0=dcmplx(zdiff(1),zdiff(2))
c
        mpole(0)=mpole(0)+charge(j)
c
        ztemp1=z0/rscale
        ztemp2=ztemp1
        do n=1,nterms
        mpole(n)=mpole(n)-charge(j)*ztemp1/n
        ztemp1=ztemp1*ztemp2
        enddo

        enddo

        return
        end
c
c
c
c
c
C***********************************************************************
        subroutine l2dformta(ier,rscale,source,charge,ns,center,
     1                       nterms,mpole)
        implicit real *8 (a-h,o-z)
C***********************************************************************
c
c     This subroutine constructs a local expansion about CENTER due
c     to NS sources located at SOURCES(2,*).
c
c     mpole_0  =  sum charge_j log(abs(z)) 
c                  j  
c
c     mpole_n  = -sum charge_j 1/n (1/z)^n *rscale^n
c                  j  
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     charge(ns)      : source strengths
c     ns              : number of sources
c     center(2)       : expansion center
c     nterms          : order of local expansion
c
c     OUTPUT:
c
c     ier       : error return code
c                 ier=0 returned successfully;
c
c     mpole     : coeffs for the multipole-expansion
c
        complex *16 mpole(0:nterms),charge(ns)
        real *8 center(2),source(2,ns),zdiff(2)

        complex *16 zmul,zinv,ztemp1,ztemp2
c
        complex *16 ima,z,z0
        data ima/(0.0d0,1.0d0)/
c
c
        do n=0,nterms
        mpole(n)=0
        enddo

        do j=1,ns
c
        zdiff(1)=source(1,j)-center(1)
        zdiff(2)=source(2,j)-center(2)
ccc        call h2cart2polar(zdiff,r,theta)
        z0=dcmplx(zdiff(1),zdiff(2))
c
        mpole(0)=mpole(0)+charge(j)*log(abs(-z0))
c
        zinv=rscale/z0
        ztemp1=zinv
        do n=1,nterms
        mpole(n)=mpole(n)-charge(j)*ztemp1/n
        ztemp1=ztemp1*zinv
        enddo

        enddo

        return
        end
c
c
c
c
c
C***********************************************************************
        subroutine l2dformta_add(ier,rscale,source,charge,ns,center,
     1                       nterms,mpole)
        implicit real *8 (a-h,o-z)
C***********************************************************************
c
c     This subroutine constructs a local expansion about CENTER due
c     to NS sources located at SOURCES(2,*).
c
c     mpole_0  =  sum charge_j log(abs(z)) 
c                  j  
c
c     mpole_n  = -sum charge_j 1/n (1/z)^n *rscale^n
c                  j  
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     charge(ns)      : source strengths
c     ns              : number of sources
c     center(2)       : expansion center
c     nterms          : order of local expansion
c
c     OUTPUT:
c
c     ier       : error return code
c                 ier=0 returned successfully;
c
c     mpole     : coeffs for the multipole-expansion
c
        complex *16 mpole(0:nterms),charge(ns)
        real *8 center(2),source(2,ns),zdiff(2)

        complex *16 zmul,zinv,ztemp1,ztemp2
c
        complex *16 ima,z,z0
        data ima/(0.0d0,1.0d0)/
c
c
        do j=1,ns
c
        zdiff(1)=source(1,j)-center(1)
        zdiff(2)=source(2,j)-center(2)
ccc        call h2cart2polar(zdiff,r,theta)
        z0=dcmplx(zdiff(1),zdiff(2))
c
        mpole(0)=mpole(0)+charge(j)*log(abs(-z0))
c
        zinv=rscale/z0
        ztemp1=zinv
        do n=1,nterms
        mpole(n)=mpole(n)-charge(j)*ztemp1/n
        ztemp1=ztemp1*zinv
        enddo

        enddo

        return
        end
c
c
c
c
c
        subroutine l2dmpmp(
     $     rscale1,center1,hexp,nterms1,
     $     rscale2,center2,jexp,nterms2)
        implicit real *8 (a-h,o-z)
C
C     Usage:
C
C           Converts multipole expansion to a multipole expansion.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           rscale1 = scaling parameter for original multipole expansion
C           center1 = center of original multiple expansion
C           hexp    = coefficients of original multiple expansion
C           nterms1 = order of original multipole expansion
C           rscale2 = scaling parameter for shifted multipole expansion
C           center2 = center of shifted multipole expansion
C           nterms2 = order of shifted multipole expansion
C
C     OUTPUT:
C
C           jexp    = coefficients of shifted multipole expansion
c
        complex *16 hexp(0:nterms1),jexp(0:nterms2)
        real *8 center1(2),center2(2),zdiff(2)
        complex *16 z,z0,ima, zmul,zinv,ztemp1,ztemp2
        real *8, allocatable :: carray(:,:)
        complex *16, allocatable :: z0pow(:), z0powm(:)
ccc        complex *16, allocatable :: z0pow1(:), z0pow2(:)
        complex *16 z0pow1(0:1000), z0pow2(0:1000)
c
        data ima/(0.0d0,1.0d0)/
c
cc        done=1
cc        pi=4*atan(done)
c
        nterms = nterms1+nterms2
c
        nmax = max(nterms1,nterms2)
        allocate( carray(0:nmax,0:nmax) )
c
        do l = 0,nmax
        carray(l,0) = 1.0d0
        enddo
        do m=1,nmax
        carray(m,m) = 1.0d0
        do l=m+1,nmax
            carray(l,m)=carray(l-1,m)+carray(l-1,m-1)
        enddo
        enddo
c
c
        zdiff(1)=center2(1)-center1(1)
        zdiff(2)=center2(2)-center1(2)
ccc        call h2cart2polar(zdiff,r,theta)
        z0=dcmplx(-zdiff(1),-zdiff(2))
c
c
ccc        allocate( z0pow1(0:nmax) )
        ztemp1=z0/rscale1
        ztemp2=ztemp1
        z0pow1(0)=1
        do i=1,nmax
        z0pow1(i)=ztemp1
        ztemp1=ztemp1*ztemp2
        enddo
c
ccc        allocate( z0pow2(0:nmax) )
        ztemp1=z0/rscale2
        ztemp2=ztemp1
        z0pow2(0)=1
        do i=1,nmax
        z0pow2(i)=ztemp1
        ztemp1=ztemp1*ztemp2
        enddo
c
        do i = 0,nterms2
        jexp(i) = 0
        enddo
c
        jexp(0) = hexp(0)
        do i = 1,nterms2
        jexp(i) = jexp(i) - hexp(0)*z0pow2(i)/i
        do j = 1,min(i,nterms1)
        jexp(i) = jexp(i) +
     $     hexp(j)*z0pow2(i)/z0pow1(j)*carray(i-1,j-1)
        enddo
        enddo
c
ccc        call prin2('jexp=*',jexp,2*(nterms2+1))
c
        return
        end
c
c
c
c
c
        subroutine l2dmploc(
     $     rscale1,center1,hexp,nterms1,
     $     rscale2,center2,jexp,nterms2)
        implicit real *8 (a-h,o-z)
C
C     Usage:
C
C           Converts multipole expansion to a local expansion.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           rscale1 = scaling parameter for original multipole expansion
C           center1 = center of original multiple expansion
C           hexp    = coefficients of original multiple expansion
C           nterms1 = order of original multipole expansion
C           rscale2 = scaling parameter for shifted local expansion
C           center2 = center of shifted local expansion
C           nterms2 = order of shifted local expansion
C
C     OUTPUT:
C
C           jexp    = coefficients of shifted local expansion
c
        complex *16 hexp(0:nterms1),jexp(0:nterms2)
        real *8 center1(2),center2(2),zdiff(2)
        complex *16 z,z0,ima, zmul,zinv,ztemp1,ztemp2
        real *8, allocatable :: carray(:,:)
ccc        complex *16, allocatable :: z0pow(:), z0powm(:)
        complex *16 z0pow(0:1000), z0powm(0:1000)
c
        data ima/(0.0d0,1.0d0)/
c
        done=1
        pi=4*atan(done)
c
        nterms = nterms1+nterms2
c
        nmax = nterms1+nterms2
        allocate( carray(0:nmax,0:nmax) )
c
        do l = 0,nmax
        carray(l,0) = 1.0d0
        enddo
        do m=1,nmax
        carray(m,m) = 1.0d0
        do l=m+1,nmax
            carray(l,m)=carray(l-1,m)+carray(l-1,m-1)
        enddo
        enddo
c
        zdiff(1)=center2(1)-center1(1)
        zdiff(2)=center2(2)-center1(2)
ccc        call h2cart2polar(zdiff,r,theta)
        z0=dcmplx(-zdiff(1),-zdiff(2))
c
c
ccc        allocate( z0pow(0:nmax) )
        ztemp1=1/z0
        ztemp2=ztemp1
        z0pow(0)=1
        do i=1,nmax
        z0pow(i)=ztemp1
        ztemp1=ztemp1*ztemp2        
        enddo
c
ccc        allocate( z0powm(0:nmax) )
        ztemp1=1/(-z0)
        ztemp2=ztemp1
        z0powm(0)=1
        do i=1,nmax
        z0powm(i)=ztemp1
        ztemp1=ztemp1*ztemp2
        enddo
c
c
        do i = 0,nterms2
        jexp(i) = 0
        enddo
c
        jexp(0) = hexp(0)*log(abs(-z0))
        do j = 1,nterms1
        jexp(0) = jexp(0) + hexp(j)*rscale1**j*z0powm(j)
        enddo

        do i = 1,nterms2
        jexp(i) = jexp(i) - hexp(0)/i
        do j = 1,nterms1
        jexp(i) = jexp(i) + hexp(j)*rscale1**j
     $     *z0powm(j)*carray(i+j-1,j-1)
        enddo
        jexp(i) = jexp(i)*z0pow(i)
        enddo
c
        do i = 1,nterms2
        jexp(i) = jexp(i) *rscale2**i
        enddo
c
        return
        end
c
c
c
c
c
        subroutine l2dmploc_add(
     $     rscale1,center1,hexp,nterms1,
     $     rscale2,center2,jexp,nterms2)
        implicit real *8 (a-h,o-z)
C
C     Usage:
C
C           Converts multipole expansion to a local expansion.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           rscale1 = scaling parameter for original multipole expansion
C           center1 = center of original multiple expansion
C           hexp    = coefficients of original multiple expansion
C           nterms1 = order of original multipole expansion
C           rscale2 = scaling parameter for shifted local expansion
C           center2 = center of shifted local expansion
C           nterms2 = order of shifted local expansion
C
C     OUTPUT:
C
C           jexp    = coefficients of shifted local expansion
c
        complex *16 hexp(0:nterms1),jexp(0:nterms2)
        real *8 center1(2),center2(2),zdiff(2)
        complex *16 z,cd,z0,ima, zmul,zinv,ztemp1,ztemp2
        real *8, allocatable :: carray(:,:)
        complex *16, allocatable :: z0pow(:), z0powm(:)
c
        data ima/(0.0d0,1.0d0)/
c
        done=1
        pi=4*atan(done)
c
        nterms = nterms1+nterms2
c
        nmax = nterms1+nterms2
        allocate( carray(0:nmax,0:nmax) )
c
        do l = 0,nmax
        carray(l,0) = 1.0d0
        enddo
        do m=1,nmax
        carray(m,m) = 1.0d0
        do l=m+1,nmax
            carray(l,m)=carray(l-1,m)+carray(l-1,m-1)
        enddo
        enddo
c
        zdiff(1)=center2(1)-center1(1)
        zdiff(2)=center2(2)-center1(2)
ccc        call h2cart2polar(zdiff,r,theta)
        z0=dcmplx(-zdiff(1),-zdiff(2))
c
c
        allocate( z0pow(0:nmax) )
        ztemp1=1/z0
        ztemp2=ztemp1
        z0pow(0)=1
        do i=1,nmax
        z0pow(i)=ztemp1
        ztemp1=ztemp1*ztemp2
        enddo
c
        allocate( z0powm(0:nmax) )
        ztemp1=1/(-z0)
        ztemp2=ztemp1
        z0powm(0)=1
        do i=1,nmax
        z0powm(i)=ztemp1
        ztemp1=ztemp1*ztemp2
        enddo
c
c
ccc        do i = 0,nterms2
ccc        jexp(i) = 0
ccc        enddo
c
        jexp(0) = jexp(0) + hexp(0)*log(abs(-z0))
        do j = 1,nterms1
        jexp(0) = jexp(0) + hexp(j)*rscale1**j*z0powm(j)
        enddo

        do i = 1,nterms2
        cd = - hexp(0)/i
        do j = 1,nterms1
        cd = cd + hexp(j)*rscale1**j*z0powm(j)*carray(i+j-1,j-1)
        enddo
        jexp(i) = jexp(i) + cd*z0pow(i)*rscale2**i
        enddo
c
        return
        end
c
c
c
c
c
        subroutine l2dlocloc(
     $     rscale1,center1,hexp,nterms1,
     $     rscale2,center2,jexp,nterms2)
        implicit real *8 (a-h,o-z)
C
C     Usage:
C
C           Converts local expansion to a local expansion.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           rscale1 = scaling parameter for original local expansion
C           center1 = center of original multiple expansion
C           hexp    = coefficients of original multiple expansion
C           nterms1 = order of original local expansion
C           rscale2 = scaling parameter for shifted local expansion
C           center2 = center of shifted local expansion
C           nterms2 = order of shifted local expansion
C
C     OUTPUT:
C
C           jexp    = coefficients of shifted local expansion
c
        complex *16 hexp(0:nterms1),jexp(0:nterms2)
        real *8 center1(2),center2(2),zdiff(2)
        complex *16 z,z0,ima, zmul,zinv,ztemp1,ztemp2
        real *8, allocatable :: carray(:,:)
        complex *16, allocatable :: z0pow(:), z0powm(:)
ccc        complex *16, allocatable :: z0powm1(:), z0powm2(:)
        complex *16 z0powm1(0:1000), z0powm2(0:1000)
c
        data ima/(0.0d0,1.0d0)/
c
        done=1
        pi=4*atan(done)
c
        nterms = nterms1+nterms2
c
        nmax = max(nterms1,nterms2)
        allocate( carray(0:nmax,0:nmax) )
c
        do l = 0,nmax
        carray(l,0) = 1.0d0
        enddo
        do m=1,nmax
        carray(m,m) = 1.0d0
        do l=m+1,nmax
            carray(l,m)=carray(l-1,m)+carray(l-1,m-1)
        enddo
        enddo
c
c
        zdiff(1)=center2(1)-center1(1)
        zdiff(2)=center2(2)-center1(2)
ccc        call h2cart2polar(zdiff,r,theta)
        z0=dcmplx(-zdiff(1),-zdiff(2))
c
c
ccc        allocate( z0powm1(0:nmax) )
        ztemp1=(-z0)/rscale1
        ztemp2=ztemp1
        z0powm1(0)=1
        do i=1,nmax
        z0powm1(i)=ztemp1
        ztemp1=ztemp1*ztemp2
        enddo
c
ccc        allocate( z0powm2(0:nmax) )
        ztemp1=(-z0)/rscale2
        ztemp2=ztemp1
        z0powm2(0)=1
        do i=1,nmax
        z0powm2(i)=ztemp1
        ztemp1=ztemp1*ztemp2
        enddo
c
        do i = 0,nterms2
        jexp(i) = 0
        enddo
c
        do i = 0,nterms2
        do j = i,nterms1
        jexp(i) = jexp(i) + hexp(j)
     $     *z0powm1(j)/z0powm2(i)*carray(j,i)
        enddo
        enddo
c
        return
        end
c
c
c
c
c
        subroutine l2dmpmp_carray(
     $     rscale1,center1,hexp,nterms1,
     $     rscale2,center2,jexp,nterms2,
     $     carray,ldc)
        implicit real *8 (a-h,o-z)
C
C     Usage:
C
C           Converts multipole expansion to a multipole expansion.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           rscale1 = scaling parameter for original multipole expansion
C           center1 = center of original multiple expansion
C           hexp    = coefficients of original multiple expansion
C           nterms1 = order of original multipole expansion
C           rscale2 = scaling parameter for shifted multipole expansion
C           center2 = center of shifted multipole expansion
C           nterms2 = order of shifted multipole expansion
C
C     OUTPUT:
C
C           jexp    = coefficients of shifted multipole expansion
c
        complex *16 hexp(0:nterms1),jexp(0:nterms2)
        real *8 center1(2),center2(2),zdiff(2)
        complex *16 z,z0,ima, zmul,zinv,ztemp1,ztemp2
        real *8 carray(0:ldc,0:ldc)
        complex *16, allocatable :: z0pow(:), z0powm(:)
ccc        complex *16, allocatable :: z0pow1(:), z0pow2(:)
        complex *16 z0pow1(0:1000), z0pow2(0:1000)
        real *8 rfactors(0:1000)
        complex *16 hexp1(0:1000)
c
        data ima/(0.0d0,1.0d0)/
c
cc        done=1
cc        pi=4*atan(done)
c
        nterms = nterms1+nterms2
c
        nmax = max(nterms1,nterms2)
c
        zdiff(1)=center2(1)-center1(1)
        zdiff(2)=center2(2)-center1(2)
ccc        call h2cart2polar(zdiff,r,theta)
        z0=dcmplx(-zdiff(1),-zdiff(2))
c
c
ccc        allocate( z0pow1(0:nmax) )
        ztemp1=1/(z0/rscale1)
        ztemp2=ztemp1
        z0pow1(0)=1
        do i=1,nmax
        z0pow1(i)=ztemp1
        ztemp1=ztemp1*ztemp2
        enddo
c
ccc        allocate( z0pow2(0:nmax) )
        ztemp1=z0/rscale2
        ztemp2=ztemp1
        z0pow2(0)=1
        do i=1,nmax
        z0pow2(i)=ztemp1
        ztemp1=ztemp1*ztemp2
        enddo
c
        rfactors(0)=1
        do i=1,max(nterms1,nterms2)
        rfactors(i)=rfactors(i-1)*(rscale1/rscale2)
        enddo
c
        do i = 0,nterms2
        jexp(i) = 0
        enddo
c
        do i = 0,nterms1
        hexp1(i) = hexp(i)*z0pow1(i)
        enddo
c
        jexp(0) = hexp(0)
        do i = 1,nterms2
        jexp(i) = jexp(i) - hexp1(0)/i
        do j = 1,min(i,nterms1)
        jexp(i) = jexp(i) + hexp1(j)*carray(i-1,j-1)
        enddo
        jexp(i)=jexp(i)*z0pow2(i)
        enddo
c
ccc        call prin2('jexp=*',jexp,2*(nterms2+1))
c
        return
        end
c
c
c
c
c
        subroutine l2dmploc_carray(
     $     rscale1,center1,hexp,nterms1,
     $     rscale2,center2,jexp,nterms2,
     $     carray,ldc)
        implicit real *8 (a-h,o-z)
C
C     Usage:
C
C           Converts multipole expansion to a local expansion.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           rscale1 = scaling parameter for original multipole expansion
C           center1 = center of original multiple expansion
C           hexp    = coefficients of original multiple expansion
C           nterms1 = order of original multipole expansion
C           rscale2 = scaling parameter for shifted local expansion
C           center2 = center of shifted local expansion
C           nterms2 = order of shifted local expansion
C
C     OUTPUT:
C
C           jexp    = coefficients of shifted local expansion
c
        complex *16 hexp(0:nterms1),jexp(0:nterms2)
        real *8 center1(2),center2(2),zdiff(2)
        complex *16 z,z0,ima, zmul,zinv,ztemp1,ztemp2,ztemp3
        real *8 carray(0:ldc,0:ldc)
cccc        complex *16, allocatable :: z0pow(:), z0powm(:)
        complex *16 z0pow(0:1000), z0powm(0:1000)
        complex *16 hexp1(0:1000)
c
        data ima/(0.0d0,1.0d0)/
c
c        done=1
c        pi=4*atan(done)
c
ccc        nmax = nterms1+nterms2
        nmax = max(nterms1,nterms2)
c
c
        zdiff(1)=center2(1)-center1(1)
        zdiff(2)=center2(2)-center1(2)
ccc        call h2cart2polar(zdiff,r,theta)
        z0=dcmplx(-zdiff(1),-zdiff(2))
c
c
ccc        allocate( z0pow(0:nmax) )
ccc        allocate( z0powm(0:nmax) )
        ztemp1=1/z0
        z0pow(0)=1
        z0powm(0)=1
        ztemp2=+ztemp1*rscale2
        ztemp3=-ztemp1*rscale1
        do i=1,nmax
        z0pow(i)=ztemp2
        z0powm(i)=ztemp3
        ztemp2=+ztemp2*ztemp1*rscale2
        ztemp3=-ztemp3*ztemp1*rscale1
        enddo

c
        do i = 0,nterms2
        jexp(i) = 0
        enddo
c
        do i = 0,nterms1
        hexp1(i) = hexp(i)*z0powm(i)
        enddo
c
        jexp(0) = hexp1(0)*log(abs(-z0))
        do j = 1,nterms1
        jexp(0) = jexp(0) + hexp1(j)
        enddo

        do i = 1,nterms2
        jexp(i) = jexp(i) - hexp1(0)/i
        do j = 1,nterms1
        jexp(i) = jexp(i) + hexp1(j)*carray(i+j-1,j-1)
        enddo
        jexp(i) = jexp(i)*z0pow(i) 
        enddo
c
ccc        call prin2('jexp=*',jexp,2*(nterms2+1))
c
        return
        end
c
c
c
c
c
        subroutine l2dlocloc_carray(
     $     rscale1,center1,hexp,nterms1,
     $     rscale2,center2,jexp,nterms2,
     $     carray,ldc)
        implicit real *8 (a-h,o-z)
C
C     Usage:
C
C           Converts local expansion to a local expansion.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           rscale1 = scaling parameter for original local expansion
C           center1 = center of original multiple expansion
C           hexp    = coefficients of original multiple expansion
C           nterms1 = order of original local expansion
C           rscale2 = scaling parameter for shifted local expansion
C           center2 = center of shifted local expansion
C           nterms2 = order of shifted local expansion
C
C     OUTPUT:
C
C           jexp    = coefficients of shifted local expansion
c
        complex *16 hexp(0:nterms1),jexp(0:nterms2)
        real *8 center1(2),center2(2),zdiff(2)
        complex *16 z,z0,ima, zmul,zinv,ztemp1,ztemp2
        real *8 carray(0:ldc,0:ldc)
ccc        complex *16, allocatable :: z0powm1(:), z0powm2(:)
        complex *16 z0powm1(0:ldc), z0powm2(0:ldc)
        complex *16 hexp1(0:1000)
c
        data ima/(0.0d0,1.0d0)/
c
        done=1
        pi=4*atan(done)
c
        nterms = nterms1+nterms2
c
        nmax = max(nterms1,nterms2)
c
        zdiff(1)=center2(1)-center1(1)
        zdiff(2)=center2(2)-center1(2)
ccc        call h2cart2polar(zdiff,r,theta)
        z0=dcmplx(-zdiff(1),-zdiff(2))
c
c
ccc        allocate( z0powm1(0:nmax) )
        ztemp1=(-z0)/rscale1
        ztemp2=ztemp1
        z0powm1(0)=1
        do i=1,nmax
        z0powm1(i)=ztemp1
        ztemp1=ztemp1*ztemp2
        enddo
c
ccc        allocate( z0powm2(0:nmax) )
        ztemp1=1/((-z0)/rscale2)
        ztemp2=ztemp1
        z0powm2(0)=1
        do i=1,nmax
        z0powm2(i)=ztemp1
        ztemp1=ztemp1*ztemp2
        enddo
c
        do i = 0,nterms2
        jexp(i) = 0
        enddo
c
        do i = 0,nterms1
        hexp1(i) = hexp(i)*z0powm1(i)
        enddo
c
        do i = 0,nterms2
        do j = i,nterms1
        jexp(i) = jexp(i) + hexp1(j)*carray(j,i)
        enddo
        jexp(i) = jexp(i)*z0powm2(i)
        enddo
c
        return
        end
c
c
c
c
c
C***********************************************************************
        subroutine l2dformmp_dp(ier,rscale,source,dipstr,ns,center,
     1                       nterms,mpole)
        implicit real *8 (a-h,o-z)
C***********************************************************************
c
c     This subroutine constructs a multipole expansion about CENTER due
c     to NS sources located at SOURCES(2,*).
c
c     mpole_0  =  0
c                 
c
c     mpole_n  =  sum dipstr_j (z_0)^(n-1)/z^n /rscale^n
c                  j  
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     dipstr(ns)      : source strengths
c     ns              : number of sources
c     center(2)       : expansion center
c     nterms          : order of multipole expansion
c
c     OUTPUT:
c
c     ier       : error return code
c                 ier=0 returned successfully;
c
c     mpole     : coeffs for the multipole-expansion
c
        complex *16 mpole(0:nterms),dipstr(ns)
        real *8 center(2),source(2,ns),zdiff(2)

        complex *16 zmul,zinv,ztemp1,ztemp2
c
        complex *16 ima,z,z0
        data ima/(0.0d0,1.0d0)/
c
c
        do n=0,nterms
        mpole(n)=0
        enddo

        do j=1,ns
c
        zdiff(1)=source(1,j)-center(1)
        zdiff(2)=source(2,j)-center(2)
ccc        call h2cart2polar(zdiff,r,theta)
        z0=dcmplx(zdiff(1),zdiff(2))
c
c
        zmul=z0/rscale
        ztemp1=1/rscale
        do n=1,nterms
        mpole(n)=mpole(n)+dipstr(j)*ztemp1
        ztemp1=ztemp1*zmul
        enddo

        enddo

        return
        end
c
c
c
c
c
C***********************************************************************
        subroutine l2dformta_dp(ier,rscale,source,dipstr,ns,center,
     1                       nterms,mpole)
        implicit real *8 (a-h,o-z)
C***********************************************************************
c
c     This subroutine constructs a local expansion about CENTER due
c     to NS sources located at SOURCES(2,*).
c
c     mpole_n  = -sum dipstr_j (z_0)^(n-1)/z^n /rscale^n
c                  j  
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     dipstr(ns)      : source strengths
c     ns              : number of sources
c     center(2)       : expansion center
c     nterms          : order of local expansion
c
c     OUTPUT:
c
c     ier       : error return code
c                 ier=0 returned successfully;
c
c     mpole     : coeffs for the multipole-expansion
c
        complex *16 mpole(0:nterms),dipstr(ns)
        real *8 center(2),source(2,ns),zdiff(2)

        complex *16 zmul,zinv,ztemp1,ztemp2
c
        complex *16 ima,z,z0
        data ima/(0.0d0,1.0d0)/
c
c
        do n=0,nterms
        mpole(n)=0
        enddo

        do j=1,ns
c
        zdiff(1)=source(1,j)-center(1)
        zdiff(2)=source(2,j)-center(2)
ccc        call h2cart2polar(zdiff,r,theta)
        z0=dcmplx(zdiff(1),zdiff(2))
c
        zinv=rscale/z0
        ztemp1=zinv/rscale
        do n=0,nterms
        mpole(n)=mpole(n)-dipstr(j)*ztemp1
        ztemp1=ztemp1*zinv
        enddo

        enddo

        return
        end
c
c
c
c
c
C***********************************************************************
        subroutine l2dformta_dp_add(ier,rscale,source,dipstr,ns,center,
     1                       nterms,mpole)
        implicit real *8 (a-h,o-z)
C***********************************************************************
c
c     This subroutine constructs a local expansion about CENTER due
c     to NS sources located at SOURCES(2,*).
c
c     mpole_n  = -sum dipstr_j (z_0)^(n-1)/z^n /rscale^n
c                  j  
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     dipstr(ns)      : source strengths
c     ns              : number of sources
c     center(2)       : expansion center
c     nterms          : order of local expansion
c
c     OUTPUT:
c
c     ier       : error return code
c                 ier=0 returned successfully;
c
c     mpole     : coeffs for the multipole-expansion
c
        complex *16 mpole(0:nterms),dipstr(ns)
        real *8 center(2),source(2,ns),zdiff(2)

        complex *16 zmul,zinv,ztemp1,ztemp2
c
        complex *16 ima,z,z0
        data ima/(0.0d0,1.0d0)/
c
c
        do j=1,ns
c
        zdiff(1)=source(1,j)-center(1)
        zdiff(2)=source(2,j)-center(2)
ccc        call h2cart2polar(zdiff,r,theta)
        z0=dcmplx(zdiff(1),zdiff(2))
c
        zinv=rscale/z0
        ztemp1=zinv/rscale
        do n=0,nterms
        mpole(n)=mpole(n)-dipstr(j)*ztemp1
        ztemp1=ztemp1*zinv
        enddo

        enddo

        return
        end
c
c
c
c
c
c**********************************************************************
c
c       Direct evaluation of Cauchy type sums, obsoleted version
c
c**********************************************************************
      subroutine lpotgrad2dall(ifgrad,ifhess,sources,charge,ns,
     1                   target,pot,grad,hess)
c**********************************************************************
c
c     This subroutine calculates the potential POT, gradient GRAD, and
c     Hessian at the target point TARGET, due to a collection of 
c     charges at SOURCE(2,ns). 
c     The scaling is that required of the delta function
c     response: i.e.,
c     
c              	pot = log(abs(z))
c		grad = gradient = (d/dx, d/dy)
c		hess = Hessian = (d^2/dx^2, d^2/dxdy, d^2/dy^2)
c
c---------------------------------------------------------------------
c     INPUT:
c
c     ifgrad        : flag for computing gradient
c	                 	   ifgrad = 0 -> don't compute 
c		                   ifgrad = 1 -> do compute 
c     ifhess        : flag for computing Hessian
c	                 	   ifhess = 0 -> don't compute 
c		                   ifhess = 1 -> do compute 
c     sources(2,*)  : location of the sources
c     charge        : charge strengths
c     ns            : number of sources
c     target        : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot   (complex *16)      : calculated potential
c     grad  (complex *16)      : calculated gradient
c     hess  (complex *16)      : calculated Hessian
c
c---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 sources(2,ns),target(2)
      complex *16 pot,grad(2),hess(3),potloc,gradloc(2),hessloc(3)
      complex *16 h0,h1,cd,eye,z
      complex *16 charge(ns)
c
      data eye/(0.0d0,1.0d0)/
c
      pot = 0.0d0
      if (ifgrad.eq.1) then
         grad(1) = 0.0d0
         grad(2) = 0.0d0
      endif
      if (ifhess.eq.1) then
         hess(1) = 0.0d0
         hess(2) = 0.0d0
         hess(3) = 0.0d0
      endif
c
      do i = 1,ns
         call lpotgrad2d(ifgrad,ifhess,sources(1,i),charge(i),target,
     1        potloc,gradloc,hessloc)
         pot = pot + potloc
         if (ifgrad.eq.1) then
         grad(1) = grad(1) + gradloc(1)
         grad(2) = grad(2) + gradloc(2)
         endif
         if (ifhess.eq.1) then
         hess(1) = hess(1) + hessloc(1)
         hess(2) = hess(2) + hessloc(2)
         hess(3) = hess(3) + hessloc(3)
         endif
      enddo
      return
      end
c
c
c
c
c**********************************************************************
      subroutine lpotgrad2d(ifgrad,ifhess,source,charge,target,
     1                pot,grad,hess)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine calculates the potential POT, gradient GRAD
c     and Hessian HESS at the target point TARGET, due to a charge at 
c     SOURCE. The scaling is that required of the delta function
c     response: i.e.,
c     
c              	pot = log(abs(z))
c		grad = gradient = (d/dx, d/dy)
c		hess = Hessian = (d^2/dx^2, d^2/dxdy, d^2/dy^2)
c
c---------------------------------------------------------------------
c     INPUT:
c
c     ifgrad     : flag for computing gradient
c	                 	ifgrad = 0 -> don't compute 
c		                ifgrad = 1 -> do compute 
c     ifhess     : flag for computing hessian
c	                 	ifhess = 0 -> don't compute 
c		                ifhess = 1 -> do compute 
c     source    : location of the source 
c     charge    : charge strength
c     target    : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot       : calculated potential
c     grad      : calculated gradient
c     hess      : calculated hessian
c
c---------------------------------------------------------------------
      real *8 source(2),target(2)
      complex *16 pot,grad(2),hess(3)
      complex *16 charge
      complex *16 z, cd, zk, ima
c
      data ima/(0.0d0,1.0d0)/
c
c ... Calculate offsets and distance
c
      xdiff=target(1)-source(1)
      ydiff=target(2)-source(2)
ccc      rr=xdiff*xdiff+ydiff*ydiff
ccc      r=sqrt(rr)
c
      z=dcmplx(xdiff,ydiff)
c
      pot=charge*log(abs(z))
c
      if (ifgrad.eq.1) then
         cd = charge/z
         grad(1) = cd
         grad(2) = cd*ima
      endif
c
      if (ifhess.eq.1) then
         cd = -charge/z**2
         hess(1) = cd
         hess(2) = cd*ima
         hess(3) = -cd
      endif
c
      return
      end
c
c
c
c**********************************************************************
      subroutine lpotgrad2dall_dp(ifgrad,ifhess,sources,dipstr,ns,
     1                   target,pot,grad,hess)
c**********************************************************************
c
c     This subroutine calculates the potential POT, gradient GRAD, and
c     Hessian at the target point TARGET, due to a collection of dipoles
c     at SOURCE(2,ns).  The scaling is that required of the delta
c     function response: i.e.,
c     
c              	pot = 1/z
c		grad = gradient = (d/dx, d/dy)
c		hess = Hessian = (d^2/dx^2, d^2/dxdy, d^2/dy^2)
c
c---------------------------------------------------------------------
c     INPUT:
c
c     ifgrad        : flag for computing gradient
c	                 	   ifgrad = 0 -> don't compute 
c		                   ifgrad = 1 -> do compute 
c     ifhess        : flag for computing Hessian
c	                 	   ifhess = 0 -> don't compute 
c		                   ifhess = 1 -> do compute 
c     sources(2,*)  : location of the sources
c     dipstr        : dipole strengths
c     ns            : number of sources
c     target        : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot   (complex *16)      : calculated potential
c     grad  (complex *16)      : calculated gradient
c     hess  (complex *16)      : calculated Hessian
c
c---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 sources(2,ns),target(2)
      complex *16 pot,grad(2),hess(3),potloc,gradloc(2),hessloc(3)
      complex *16 h0,h1,cd,eye,z
      complex *16 dipstr(ns)
c
      data eye/(0.0d0,1.0d0)/
c
      pot = 0.0d0
      if (ifgrad.eq.1) then
         grad(1) = 0.0d0
         grad(2) = 0.0d0
      endif
      if (ifhess.eq.1) then
         hess(1) = 0.0d0
         hess(2) = 0.0d0
         hess(3) = 0.0d0
      endif
c
      do i = 1,ns
         call lpotgrad2d_dp(ifgrad,ifhess,sources(1,i),dipstr(i),target,
     1        potloc,gradloc,hessloc)
         pot = pot + potloc
         if (ifgrad.eq.1) then
         grad(1) = grad(1) + gradloc(1)
         grad(2) = grad(2) + gradloc(2)
         endif
         if (ifhess.eq.1) then
         hess(1) = hess(1) + hessloc(1)
         hess(2) = hess(2) + hessloc(2)
         hess(3) = hess(3) + hessloc(3)
         endif
      enddo
      return
      end
c
c
c
c
c**********************************************************************
      subroutine lpotgrad2d_dp(ifgrad,ifhess,source,dipstr,target,
     1                pot,grad,hess)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine calculates the potential POT, gradient GRAD
c     and Hessian HESS at the target point TARGET, due to a dipole at 
c     SOURCE. The scaling is that required of the delta function
c     response: i.e.,
c     
c              	pot = 1/z
c		grad = gradient = (d/dx, d/dy)
c		hess = Hessian = (d^2/dx^2, d^2/dxdy, d^2/dy^2)
c
c---------------------------------------------------------------------
c     INPUT:
c
c     ifgrad     : flag for computing gradient
c	                 	ifgrad = 0 -> don't compute 
c		                ifgrad = 1 -> do compute 
c     ifhess     : flag for computing hessian
c	                 	ifhess = 0 -> don't compute 
c		                ifhess = 1 -> do compute 
c     source    : location of the source 
c     dipstr    : dipole strength
c     target    : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot       : calculated potential
c     grad      : calculated gradient
c     hess      : calculated hessian
c
c---------------------------------------------------------------------
      real *8 source(2),target(2)
      complex *16 pot,grad(2),hess(3)
      complex *16 dipstr
      complex *16 z, cd, zk, ima
c
      data ima/(0.0d0,1.0d0)/
c
c ... Calculate offsets and distance
c
      xdiff=target(1)-source(1)
      ydiff=target(2)-source(2)
ccc      rr=xdiff*xdiff+ydiff*ydiff
ccc      r=sqrt(rr)
c
      z=dcmplx(xdiff,ydiff)
c
      pot=dipstr/z
c
      if (ifgrad.eq.1) then
         cd = -dipstr/z**2
         grad(1) = cd
         grad(2) = cd*ima
      endif
c
      if (ifhess.eq.1) then
         cd = 2*dipstr/z**3
         hess(1) = cd
         hess(2) = cd*ima
         hess(3) = -cd
      endif
c
      return
      end
c
c
c
c**********************************************************************
      subroutine lpotgrad2dall_sdp(sources,ns,
     $     ifcharge,charge,ifdipole,dipstr,
     1     target,ifpot,pot,ifgrad,grad,ifhess,hess)
c**********************************************************************
c
c     This subroutine calculates the potential POT, gradient GRAD, and
c     Hessian at the target point TARGET, due to a collection of charges
c     and dipoles at SOURCE(2,ns).  The scaling is that required of the
c     delta function response: i.e.,
c     
c              	pot = charge log(abs(z)) + dipstr 1/z
c		grad = gradient = (d/dx, d/dy)
c		hess = Hessian = (d^2/dx^2, d^2/dxdy, d^2/dy^2)
c
c---------------------------------------------------------------------
c     INPUT:
c
c     ifpot         : flag for computing potential
c	                 	   ifpot = 0 -> don't compute 
c		                   ifpot = 1 -> do compute 
c     ifgrad        : flag for computing gradient
c	                 	   ifgrad = 0 -> don't compute 
c		                   ifgrad = 1 -> do compute 
c     ifhess        : flag for computing Hessian
c	                 	   ifhess = 0 -> don't compute 
c		                   ifhess = 1 -> do compute 
c     sources(2,*)  : location of the sources
c     charge        : charge strengths
c     ns            : number of sources
c     target        : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot   (complex *16)      : calculated potential
c     grad  (complex *16)      : calculated gradient
c     hess  (complex *16)      : calculated Hessian
c
c---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 sources(2,ns),target(2)
      complex *16 pot,grad(2),hess(3),potloc,gradloc(2),hessloc(3)
      complex *16 h0,h1,cd,eye,z
      complex *16 charge(ns),dipstr(ns)
c
      data eye/(0.0d0,1.0d0)/
c
      if (ifpot.eq.1) pot = 0.0d0
      if (ifgrad.eq.1) then
         grad(1) = 0.0d0
         grad(2) = 0.0d0
      endif
      if (ifhess.eq.1) then
         hess(1) = 0.0d0
         hess(2) = 0.0d0
         hess(3) = 0.0d0
      endif
c
      do i = 1,ns
c        call lpotgrad2d_sdp(sources(1,i),
c     $     ifcharge,charge(i),ifdipole,dipstr(i),
c     1     target,ifpot,potloc,ifgrad,gradloc,ifhess,hessloc)
c         if (ifpot.eq.1) pot = pot + potloc
c         if (ifgrad.eq.1) then
c         grad(1) = grad(1) + gradloc(1)
c         grad(2) = grad(2) + gradloc(2)
c         endif
c         if (ifhess.eq.1) then
c         hess(1) = hess(1) + hessloc(1)
c         hess(2) = hess(2) + hessloc(2)
c         hess(3) = hess(3) + hessloc(3)
c         endif
        call lpotgrad2d_sdp_add(sources(1,i),
     $     ifcharge,charge(i),ifdipole,dipstr(i),
     1     target,ifpot,pot,ifgrad,grad,ifhess,hess)
      enddo
      return
      end
c
c
c
c
c**********************************************************************
      subroutine lpotgrad2dall_sdp_add(sources,ns,
     $     ifcharge,charge,ifdipole,dipstr,
     1     target,ifpot,pot,ifgrad,grad,ifhess,hess)
c**********************************************************************
c
c     This subroutine calculates the potential POT, gradient GRAD, and
c     Hessian at the target point TARGET, due to a collection of charges
c     and dipoles at SOURCE(2,ns).  The scaling is that required of the
c     delta function response: i.e.,
c     
c              	pot = charge log(abs(z)) + dipstr 1/z
c		grad = gradient = (d/dx, d/dy)
c		hess = Hessian = (d^2/dx^2, d^2/dxdy, d^2/dy^2)
c
c---------------------------------------------------------------------
c     INPUT:
c
c     ifpot         : flag for computing potential
c	                 	   ifpot = 0 -> don't compute 
c		                   ifpot = 1 -> do compute 
c     ifgrad        : flag for computing gradient
c	                 	   ifgrad = 0 -> don't compute 
c		                   ifgrad = 1 -> do compute 
c     ifhess        : flag for computing Hessian
c	                 	   ifhess = 0 -> don't compute 
c		                   ifhess = 1 -> do compute 
c     sources(2,*)  : location of the sources
c     charge        : charge strengths
c     ns            : number of sources
c     target        : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot   (complex *16)      : calculated potential
c     grad  (complex *16)      : calculated gradient
c     hess  (complex *16)      : calculated Hessian
c
c---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 sources(2,ns),target(2)
      complex *16 pot,grad(2),hess(3),potloc,gradloc(2),hessloc(3)
      complex *16 h0,h1,cd,eye,z
      complex *16 charge(ns),dipstr(ns)
c
      data eye/(0.0d0,1.0d0)/
c
c      if (ifpot.eq.1) pot = 0.0d0
c      if (ifgrad.eq.1) then
c         grad(1) = 0.0d0
c         grad(2) = 0.0d0
c      endif
c      if (ifhess.eq.1) then
c         hess(1) = 0.0d0
c         hess(2) = 0.0d0
c         hess(3) = 0.0d0
c      endif
c
      do i = 1,ns
c        call lpotgrad2d_sdp(sources(1,i),
c     $     ifcharge,charge(i),ifdipole,dipstr(i),
c     1     target,ifpot,potloc,ifgrad,gradloc,ifhess,hessloc)
c         if (ifpot.eq.1) pot = pot + potloc
c         if (ifgrad.eq.1) then
c         grad(1) = grad(1) + gradloc(1)
c         grad(2) = grad(2) + gradloc(2)
c         endif
c         if (ifhess.eq.1) then
c         hess(1) = hess(1) + hessloc(1)
c         hess(2) = hess(2) + hessloc(2)
c         hess(3) = hess(3) + hessloc(3)
c         endif
        call lpotgrad2d_sdp_add(sources(1,i),
     $     ifcharge,charge(i),ifdipole,dipstr(i),
     1     target,ifpot,pot,ifgrad,grad,ifhess,hess)
      enddo
      return
      end
c
c
c
c
c**********************************************************************
      subroutine lpotgrad2d_sdp(source,
     $     ifcharge,charge,ifdipole,dipstr,
     1     target,ifpot,pot,ifgrad,grad,ifhess,hess)
c**********************************************************************
c
c     This subroutine calculates the potential POT, gradient GRAD and
c     Hessian HESS at the target point TARGET, due to a charge and a
c     dipole at SOURCE. The scaling is that required of the delta
c     function response: i.e.,
c     
c              	pot = charge log(abs(z)) + dipstr 1/z
c		grad = gradient = (d/dx, d/dy)
c		hess = Hessian = (d^2/dx^2, d^2/dxdy, d^2/dy^2)
c
c---------------------------------------------------------------------
c     INPUT:
c
c     ifpot      : flag for computing potential
c	                 	ifpot = 0 -> don't compute 
c		                ifpot = 1 -> do compute 
c     ifgrad     : flag for computing gradient
c	                 	ifgrad = 0 -> don't compute 
c		                ifgrad = 1 -> do compute 
c     ifhess     : flag for computing hessian
c	                 	ifhess = 0 -> don't compute 
c		                ifhess = 1 -> do compute 
c     source    : location of the source 
c     charge    : charge strength
c     target    : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot       : calculated potential
c     grad      : calculated gradient
c     hess      : calculated hessian
c
c---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 source(2),target(2)
      complex *16 pot,grad(2),hess(3)
      complex *16 charge,dipstr
      complex *16 z, cd, zk, ima, zinv, zinv2
c
      data ima/(0.0d0,1.0d0)/
c
c ... Calculate offsets and distance
c
      xdiff=target(1)-source(1)
      ydiff=target(2)-source(2)
ccc      rr=xdiff*xdiff+ydiff*ydiff
ccc      r=sqrt(rr)
c
      z=dcmplx(xdiff,ydiff)
c
      if (ifpot.eq.1) then
         pot = 0
      endif
c
      if (ifgrad.eq.1) then
         grad(1) = 0
         grad(2) = 0
      endif
c
      if (ifhess.eq.1) then
         hess(1) = 0
         hess(2) = 0
         hess(3) = 0
      endif

        zinv=1/z
        zinv2=zinv*zinv

      if( ifcharge .eq. 1 ) then
c
      if (ifpot.eq.1) pot=charge*log(abs(z))
c
      if (ifgrad.eq.1) then
         cd = charge*zinv
         grad(1) = cd
         grad(2) = cd*ima
      endif
c
      if (ifhess.eq.1) then
         cd = -charge*zinv2
         hess(1) = cd
         hess(2) = cd*ima
         hess(3) = -cd
      endif
c
      endif
c
c
c
      if( ifdipole .eq. 1 ) then

      if (ifpot.eq.1) pot=pot+dipstr*zinv
c
      if (ifgrad.eq.1) then
         cd = -dipstr*zinv2
         grad(1) = grad(1)+cd
         grad(2) = grad(2)+cd*ima
      endif
c
      if (ifhess.eq.1) then
         cd = 2*dipstr/z**3
         hess(1) = hess(1)+cd
         hess(2) = hess(2)+cd*ima
         hess(3) = hess(3)-cd
      endif
      
      endif
c
c
      return
      end
c
c
c
c**********************************************************************
      subroutine lpotgrad2d_sdp_add(source,
     $     ifcharge,charge,ifdipole,dipstr,
     1     target,ifpot,pot,ifgrad,grad,ifhess,hess)
c**********************************************************************
c
c     This subroutine calculates the potential POT, gradient GRAD and
c     Hessian HESS at the target point TARGET, due to a charge and a
c     dipole at SOURCE. The scaling is that required of the delta
c     function response: i.e.,
c     
c              	pot = charge log(abs(z)) + dipstr 1/z
c		grad = gradient = (d/dx, d/dy)
c		hess = Hessian = (d^2/dx^2, d^2/dxdy, d^2/dy^2)
c
c---------------------------------------------------------------------
c     INPUT:
c
c     ifpot      : flag for computing potential
c	                 	ifpot = 0 -> don't compute 
c		                ifpot = 1 -> do compute 
c     ifgrad     : flag for computing gradient
c	                 	ifgrad = 0 -> don't compute 
c		                ifgrad = 1 -> do compute 
c     ifhess     : flag for computing hessian
c	                 	ifhess = 0 -> don't compute 
c		                ifhess = 1 -> do compute 
c     source    : location of the source 
c     charge    : charge strength
c     target    : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot       : calculated potential
c     grad      : calculated gradient
c     hess      : calculated hessian
c
c---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 source(2),target(2)
      complex *16 pot,grad(2),hess(3)
      complex *16 charge,dipstr
      complex *16 z, cd, zk, ima, zinv, zinv2
c
      data ima/(0.0d0,1.0d0)/
c
c ... Calculate offsets and distance
c
      xdiff=target(1)-source(1)
      ydiff=target(2)-source(2)
ccc      rr=xdiff*xdiff+ydiff*ydiff
ccc      r=sqrt(rr)
c
      z=dcmplx(xdiff,ydiff)
c
        zinv=1/z
        zinv2=zinv*zinv

      if( ifcharge .eq. 1 ) then
c
      if (ifpot.eq.1) pot=pot+charge*log(abs(z))
c
      if (ifgrad.eq.1) then
         cd = charge*zinv
         grad(1) = grad(1)+cd
         grad(2) = grad(2)+cd*ima
      endif
c
      if (ifhess.eq.1) then
         cd = -charge*zinv2
         hess(1) = hess(1)+cd
         hess(2) = hess(2)+cd*ima
         hess(3) = hess(3)-cd
      endif
c
      endif
c
c
c
      if( ifdipole .eq. 1 ) then

      if (ifpot.eq.1) pot=pot+dipstr*zinv
c
      if (ifgrad.eq.1) then
         cd = -dipstr*zinv2
         grad(1) = grad(1)+cd
         grad(2) = grad(2)+cd*ima
      endif
c
      if (ifhess.eq.1) then
         cd = 2*dipstr*zinv*zinv2
         hess(1) = hess(1)+cd
         hess(2) = hess(2)+cd*ima
         hess(3) = hess(3)-cd
      endif
      
      endif
c
c
      return
      end
c
c
c
c
c
c**********************************************************************
c
c     Multipole and local expansion evaluation routines for Cauchy sums
c
c**********************************************************************
      subroutine c2dmpeval(rscale,center,mpole,nterms,ztarg,
     1      ifpot,pot,ifgrad,grad,ifhess,hess)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient of the 
c     potential due to an outgoing partial wave expansion.
c               +nterms
c     pot  =      sum   mpole_n / z^n  + mpole_0 log(abs(z)), if ifpot = 1.
c                n=1  
c     grad  = gradient(pot) if ifgrad = 1.
c     hess = hessian if ifhess = 1.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale :    scaling parameter 
c     center :    expansion center
c     mpole  :    multipole expansion 
c     nterms :    order of the multipole expansion
c     ztarg  :    target location
c     ifpot  :   flag controlling evaluation of potential:
c                   ifpot = 0, do not compute potential.
c                   ifpot = 1, compute potential.
c     ifgrad :   flag controlling evaluation of gradient:
c                   ifgrad = 0, do not compute gradient.
c                   ifgrad = 1, compute gradient.
c     ifhess :   flag for computing Hessian:
c	            ifhess = 0 -> don't compute 
c		    ifhess = 1 -> do compute 
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    potential at ztarg (if requested)
c     grad   :    gradient at ztarg (if requested)
c     hess   :    hessian at ztarg (if requested)
c
c-----------------------------------------------------------------------
c
        complex *16 pot,grad,hess,mpole(0:nterms)
        real *8 center(2),ztarg(2),zdiff(2)
        complex *16 ztemp1, ztemp2
ccc        complex *16, allocatable :: z0pow(:), z0powm(:)
        complex *16 z0pow(0:1000)
c
        complex *16 ima,z,cd
        data ima/(0.0d0,1.0d0)/
c
c
        zdiff(1)=ztarg(1)-center(1)
        zdiff(2)=ztarg(2)-center(2)
ccc        call h2cart2polar(zdiff,r,theta)
        z=dcmplx(zdiff(1),zdiff(2))
c
c
        nmax = nterms + 3
ccc        allocate( z0pow(0:nmax) )
        ztemp1=rscale/z
        ztemp2=ztemp1
        z0pow(0)=1
        do i=1,nmax
        z0pow(i)=ztemp1
        ztemp1=ztemp1*ztemp2
        enddo
c
c
        pot=mpole(0)*log(abs(z))
        do n=1,nterms
        pot=pot+mpole(n)*z0pow(n)
        enddo
c
c
        if( ifgrad .eq. 1 ) then

        rinv=1/rscale
        grad=mpole(0)*z0pow(1)
        do n=1,nterms
        grad=grad-mpole(n)*z0pow(n+1)*n
        enddo
        grad=grad*rinv

        endif
c
c
        if( ifhess .eq. 1 ) then

        rinv2=1/rscale**2
        hess=-mpole(0)*z0pow(2)
        do n=1,nterms
        hess=hess+mpole(n)*z0pow(n+2)*n*(n+1)
        enddo
        hess=hess*rinv2

        endif
c
c
        return
        end
c
c
c
c
c
c**********************************************************************
      subroutine c2dtaeval(rscale,center,mpole,nterms,ztarg,
     1      ifpot,pot,ifgrad,grad,ifhess,hess)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient of the 
c     potential due to an incoming partial wave expansion.
c               +nterms
c     pot  =      sum   mpole_n * z^n  + mpole_0, if ifpot = 1.
c                n=1  
c     grad  = gradient(pot) if ifgrad = 1.
c     hess = hessian if ifhess = 1.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale :    scaling parameter 
c     center :    expansion center
c     mpole  :    multipole expansion 
c     nterms :    order of the local expansion
c     ztarg  :    target location
c     ifpot  :   flag controlling evaluation of potential:
c                   ifpot = 0, do not compute potential.
c                   ifpot = 1, compute potential.
c     ifgrad :   flag controlling evaluation of gradient:
c                   ifgrad = 0, do not compute gradient.
c                   ifgrad = 1, compute gradient.
c     ifhess :   flag for computing Hessian:
c	            ifhess = 0 -> don't compute 
c		    ifhess = 1 -> do compute 
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    potential at ztarg (if requested)
c     grad   :    gradient at ztarg (if requested)
c     hess   :    hessian at ztarg (if requested)
c
c-----------------------------------------------------------------------
c
        complex *16 pot,grad,hess,mpole(0:nterms)
        real *8 center(2),ztarg(2),zdiff(2)
        complex *16 ztemp1, ztemp2
ccc        complex *16, allocatable :: z0pow(:), z0powm(:)
        complex *16 z0pow(0:1000)
c
        complex *16 ima,z,cd
        data ima/(0.0d0,1.0d0)/
c
c
        zdiff(1)=ztarg(1)-center(1)
        zdiff(2)=ztarg(2)-center(2)
ccc        call h2cart2polar(zdiff,r,theta)
        z=dcmplx(zdiff(1),zdiff(2))
c
c
c
        nmax = nterms
ccc        allocate( z0pow(0:nmax) )
        ztemp1=z/rscale
        ztemp2=ztemp1
        z0pow(0)=1
        do i=1,nmax
        z0pow(i)=ztemp1
        ztemp1=ztemp1*ztemp2
        enddo
c
c
c
        pot=mpole(0)
        do n=1,nterms
        pot=pot+mpole(n)*z0pow(n)
        enddo
c
c
        if( ifgrad .eq. 1 ) then

        grad=0

        rinv=1/rscale
        do n=1,nterms
        grad=grad+mpole(n)*z0pow(n-1)*n
        enddo
        grad=grad*rinv

        endif
c
c
        if( ifhess .eq. 1 ) then

        hess=0

        rinv2=1/rscale**2
        do n=2,nterms
        hess=hess+mpole(n)*z0pow(n-2)*n*(n-1)
        enddo
        hess=hess*rinv2
        endif
c
c
        return
        end
c
c
c
c
c
c**********************************************************************
      subroutine c2dmpevalall(rscale,center,mpole,nterms,ztarg,ntarg,
     1      ifpot,pot,ifgrad,grad,ifhess,hess)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient of the 
c     potential due to an outgoing partial wave expansion.
c               +nterms
c     pot  =      sum   mpole_n / z^n  + mpole_0 log(abs(z)) if ifgrad = 1.
c                n=1  
c     grad  = gradient(pot) if ifgrad = 1.
c     hess = hessian if ifhess = 1.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale :    scaling parameter 
c     center :    expansion center
c     mpole  :    multipole expansion 
c     nterms :    order of the multipole expansion
c     ztarg  :    target location
c     ifpot  :   flag controlling evaluation of potential:
c                   ifpot = 0, do not compute potential.
c                   ifpot = 1, compute potential.
c     ifgrad :   flag controlling evaluation of gradient:
c                   ifgrad = 0, do not compute gradient.
c                   ifgrad = 1, compute gradient.
c     ifhess :   flag for computing Hessian:
c	            ifhess = 0 -> don't compute 
c		    ifhess = 1 -> do compute 
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    potential at ztarg (if requested)
c     grad   :    gradient at ztarg (if requested)
c     hess   :    hessian at ztarg (if requested)
c
c-----------------------------------------------------------------------
c
        complex *16 pot(1),grad(1),hess(1),mpole(0:nterms)
        complex *16 potloc,gradloc,hessloc
        real *8 center(2),ztarg(2,ntarg),zdiff(2)
c
        complex *16 ima,z
        data ima/(0.0d0,1.0d0)/
c
        do i = 1,ntarg
        call c2dmpeval(rscale,center,mpole,nterms,ztarg(1,i),
     1     ifpot,potloc,ifgrad,gradloc,ifhess,hessloc)
        if (ifpot.eq.1) pot(i) = pot(i) + potloc
        if (ifgrad.eq.1) then
        grad(i) = grad(i) + gradloc
        endif
        if (ifhess.eq.1) then
        hess(i) = hess(i) + hessloc
        endif
        enddo
c
        return
        end
c
c
c
c**********************************************************************
      subroutine c2dtaevalall(rscale,center,mpole,nterms,ztarg,ntarg,
     1      ifpot,pot,ifgrad,grad,ifhess,hess)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient of the 
c     potential due to an incoming partial wave expansion.
c               +nterms
c     pot  =      sum   mpole_n * z^n  + mpole_0 if ifpot = 1.
c                n=1  
c     grad  = gradient(pot) if ifgrad = 1.
c     hess = hessian if ifhess = 1.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale :    scaling parameter 
c     center :    expansion center
c     mpole  :    multipole expansion 
c     nterms :    order of the local expansion
c     ztarg  :    target location
c     ifpot  :   flag controlling evaluation of potential:
c                   ifpot = 0, do not compute potential.
c                   ifpot = 1, compute potential.
c     ifgrad :   flag controlling evaluation of gradient:
c                   ifgrad = 0, do not compute gradient.
c                   ifgrad = 1, compute gradient.
c     ifhess :   flag for computing Hessian:
c	            ifhess = 0 -> don't compute 
c		    ifhess = 1 -> do compute 
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    potential at ztarg (if requested)
c     grad   :    gradient at ztarg (if requested)
c     hess   :    hessian at ztarg (if requested)
c
c-----------------------------------------------------------------------
c
        complex *16 pot(1),grad(1),hess(1),mpole(0:nterms)
        complex *16 potloc,gradloc,hessloc
        real *8 center(2),ztarg(2,ntarg),zdiff(2)
c
        complex *16 ima,z
        data ima/(0.0d0,1.0d0)/
c
        do i = 1,ntarg
        call c2dtaeval(rscale,center,mpole,nterms,ztarg(1,i),
     1     ifpot,potloc,ifgrad,gradloc,ifhess,hessloc)
        if (ifpot.eq.1) pot(i) = pot(i) + potloc
        if (ifgrad.eq.1) then
        grad(i) = grad(i) + gradloc
        endif
        if (ifhess.eq.1) then
        hess(i) = hess(i) + hessloc
        endif
        enddo
c
        return
        end
c
c
c
c
c
c**********************************************************************
c
c     Cauchy sums, direct evaluation routines
c
c**********************************************************************
      subroutine cpotgrad2dall_sdp(sources,ns,
     $     ifcharge,charge,ifdipole,dipstr,
     1     target,ifpot,pot,ifgrad,grad,ifhess,hess)
c**********************************************************************
c
c     This subroutine calculates the potential POT, gradient GRAD, and
c     Hessian at the target point TARGET, due to a collection of 
c     charges and dipoles at SOURCE(2,ns). 
c     The scaling is that required of the delta function
c     response: i.e.,
c     
c              	pot = charge log(abs(z)) + dipstr 1/z
c		grad = gradient = d/dz
c		hess = Hessian = d^2/dz^2
c
c---------------------------------------------------------------------
c     INPUT:
c
c     ifpot         : flag for computing potential
c	                 	   ifpot = 0 -> don't compute 
c		                   ifpot = 1 -> do compute 
c     ifgrad        : flag for computing gradient
c	                 	   ifgrad = 0 -> don't compute 
c		                   ifgrad = 1 -> do compute 
c     ifhess        : flag for computing Hessian
c	                 	   ifhess = 0 -> don't compute 
c		                   ifhess = 1 -> do compute 
c     sources(2,*)  : location of the sources
c     charge        : charge strengths
c     ns            : number of sources
c     target        : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot   (complex *16)      : calculated potential
c     grad  (complex *16)      : calculated gradient
c     hess  (complex *16)      : calculated Hessian
c
c---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 sources(2,ns),target(2)
      complex *16 pot,grad,hess,potloc,gradloc,hessloc
      complex *16 h0,h1,cd,eye,z
      complex *16 charge(ns),dipstr(ns)
c
      data eye/(0.0d0,1.0d0)/
c
      if (ifpot.eq.1) pot = 0.0d0
      if (ifgrad.eq.1) then
         grad = 0.0d0
      endif
      if (ifhess.eq.1) then
         hess = 0.0d0
      endif
c
      do i = 1,ns
        call cpotgrad2d_sdp_add(sources(1,i),
     $     ifcharge,charge(i),ifdipole,dipstr(i),
     1     target,ifpot,pot,ifgrad,grad,ifhess,hess)
      enddo
      return
      end
c
c
c
c
c**********************************************************************
      subroutine cpotgrad2dall_sdp_add(sources,ns,
     $     ifcharge,charge,ifdipole,dipstr,
     1     target,ifpot,pot,ifgrad,grad,ifhess,hess)
c**********************************************************************
c
c     This subroutine calculates the potential POT, gradient GRAD, and
c     Hessian at the target point TARGET, due to a collection of 
c     charges and dipoles at SOURCE(2,ns). 
c     The scaling is that required of the delta function
c     response: i.e.,
c     
c              	pot = charge log(abs(z)) + dipstr 1/z
c		grad = gradient = d/dz
c		hess = Hessian = d^2/dz^2
c
c---------------------------------------------------------------------
c     INPUT:
c
c     ifpot         : flag for computing potential
c	                 	   ifpot = 0 -> don't compute 
c		                   ifpot = 1 -> do compute 
c     ifgrad        : flag for computing gradient
c	                 	   ifgrad = 0 -> don't compute 
c		                   ifgrad = 1 -> do compute 
c     ifhess        : flag for computing Hessian
c	                 	   ifhess = 0 -> don't compute 
c		                   ifhess = 1 -> do compute 
c     sources(2,*)  : location of the sources
c     charge        : charge strengths
c     ns            : number of sources
c     target        : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot   (complex *16)      : calculated potential
c     grad  (complex *16)      : calculated gradient
c     hess  (complex *16)      : calculated Hessian
c
c---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 sources(2,ns),target(2)
      complex *16 pot,grad,hess,potloc,gradloc,hessloc
      complex *16 h0,h1,cd,eye,z
      complex *16 charge(ns),dipstr(ns)
c
      data eye/(0.0d0,1.0d0)/
c
      do i = 1,ns
        call cpotgrad2d_sdp_add(sources(1,i),
     $     ifcharge,charge(i),ifdipole,dipstr(i),
     1     target,ifpot,pot,ifgrad,grad,ifhess,hess)
      enddo
      return
      end
c
c
c
c
c**********************************************************************
      subroutine cpotgrad2d_sdp(source,
     $     ifcharge,charge,ifdipole,dipstr,
     1     target,ifpot,pot,ifgrad,grad,ifhess,hess)
c**********************************************************************
c
c     This subroutine calculates the potential POT, gradient GRAD and
c     Hessian HESS at the target point TARGET, due to a charge and a
c     dipole at SOURCE. The scaling is that required of the delta
c     function response: i.e.,
c     
c              	pot = charge log(abs(z)) + dipstr 1/z
c		grad = gradient = d/dz
c		hess = Hessian = d^2/dz^2
c
c---------------------------------------------------------------------
c     INPUT:
c
c     ifpot      : flag for computing potential
c	                	ifpot = 0 -> don't compute 
c		                ifpot = 1 -> do compute 
c     ifgrad     : flag for computing gradient
c	                 	ifgrad = 0 -> don't compute 
c		                ifgrad = 1 -> do compute 
c     ifhess     : flag for computing hessian
c	                 	ifhess = 0 -> don't compute 
c		                ifhess = 1 -> do compute 
c     source    : location of the source 
c     charge    : charge strength
c     target    : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot       : calculated potential
c     grad      : calculated gradient
c     hess      : calculated hessian
c
c---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 source(2),target(2)
      complex *16 pot,grad,hess
      complex *16 charge,dipstr
      complex *16 z, cd, zk, ima, zinv, zinv2
c
      data ima/(0.0d0,1.0d0)/
c
c ... Calculate offsets and distance
c
      xdiff=target(1)-source(1)
      ydiff=target(2)-source(2)
ccc      rr=xdiff*xdiff+ydiff*ydiff
ccc      r=sqrt(rr)
c
      z=dcmplx(xdiff,ydiff)
c
      if (ifpot.eq.1) then
         pot = 0
      endif
c
      if (ifgrad.eq.1) then
         grad = 0
      endif
c
      if (ifhess.eq.1) then
         hess = 0
      endif

        zinv=1/z
        zinv2=zinv*zinv

      if( ifcharge .eq. 1 ) then
c
      if (ifpot.eq.1) pot=charge*log(abs(z))
c
      if (ifgrad.eq.1) then
         grad = charge*zinv
      endif
c
      if (ifhess.eq.1) then
         hess = -charge*zinv2
      endif
c
      endif
c
c
c
      if( ifdipole .eq. 1 ) then

      if (ifpot.eq.1) pot=pot+dipstr*zinv
c
      if (ifgrad.eq.1) then
         grad = grad-dipstr*zinv2
      endif
c
      if (ifhess.eq.1) then
         hess = hess+2*dipstr*zinv2*zinv
      endif
      
      endif
c
c
      return
      end
c
c
c
c**********************************************************************
      subroutine cpotgrad2d_sdp_add(source,
     $     ifcharge,charge,ifdipole,dipstr,
     1     target,ifpot,pot,ifgrad,grad,ifhess,hess)
c**********************************************************************
c
c     This subroutine calculates the potential POT, gradient GRAD and
c     Hessian HESS at the target point TARGET, due to a charge and a
c     dipole at SOURCE. The scaling is that required of the delta
c     function response: i.e.,
c     
c              	pot = charge log(abs(z)) + dipstr 1/z
c		grad = gradient = d/dz
c		hess = Hessian = d^2/dz^2
c
c---------------------------------------------------------------------
c     INPUT:
c
c     ifpot      : flag for computing potential
c	                 	ifpot = 0 -> don't compute 
c		                ifpot = 1 -> do compute 
c     ifgrad     : flag for computing gradient
c	                 	ifgrad = 0 -> don't compute 
c		                ifgrad = 1 -> do compute 
c     ifhess     : flag for computing hessian
c	                 	ifhess = 0 -> don't compute 
c		                ifhess = 1 -> do compute 
c     source    : location of the source 
c     charge    : charge strength
c     target    : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot       : calculated potential
c     grad      : calculated gradient
c     hess      : calculated hessian
c
c---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 source(2),target(2)
      complex *16 pot,grad,hess
      complex *16 charge,dipstr
      complex *16 z, cd, zk, ima, zinv, zinv2
c
      data ima/(0.0d0,1.0d0)/
c
c ... Calculate offsets and distance
c
      xdiff=target(1)-source(1)
      ydiff=target(2)-source(2)
ccc      rr=xdiff*xdiff+ydiff*ydiff
ccc      r=sqrt(rr)
c
      z=dcmplx(xdiff,ydiff)
c
        zinv=1/z
        zinv2=zinv*zinv

      if( ifcharge .eq. 1 ) then
c
      if (ifpot.eq.1) pot=pot+charge*log(abs(z))
c
      if (ifgrad.eq.1) then
         grad = grad+charge*zinv
      endif
c
      if (ifhess.eq.1) then
         hess = hess-charge*zinv2
      endif
c
      endif
c
c
c
      if( ifdipole .eq. 1 ) then

      if (ifpot.eq.1) pot=pot+dipstr*zinv
c
      if (ifgrad.eq.1) then
         grad = grad-dipstr*zinv2
      endif
c
      if (ifhess.eq.1) then
         hess = hess+2*dipstr*zinv2*zinv
      endif
      
      endif
c
c
      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine cpotgrad2d_sdp_sym(source1,source2,
     $     ifcharge,charge1,charge2,ifdipole,dipstr1,dipstr2,
     1     ifpot,pot1,pot2,ifgrad,grad1,grad2,ifhess,hess1,hess2)
c**********************************************************************
c
c     This subroutine calculates the potential POT, gradient GRAD and
c     Hessian HESS at the target point TARGET, due to a charge and a
c     dipole at SOURCE. The scaling is that required of the delta
c     function response: i.e.,
c     
c              	pot = charge log(abs(z)) + dipstr 1/z
c		grad = gradient = d/dz
c		hess = Hessian = d^2/dz^2
c
c---------------------------------------------------------------------
c     INPUT:
c
c     ifpot      : flag for computing potential
c	                	ifpot = 0 -> don't compute 
c		                ifpot = 1 -> do compute 
c     ifgrad     : flag for computing gradient
c	                 	ifgrad = 0 -> don't compute 
c		                ifgrad = 1 -> do compute 
c     ifhess     : flag for computing hessian
c	                 	ifhess = 0 -> don't compute 
c		                ifhess = 1 -> do compute 
c     source    : location of the source 
c     charge    : charge strength
c     target    : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot       : calculated potential
c     grad      : calculated gradient
c     hess      : calculated hessian
c
c---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 source1(2),source2(2)
      complex *16 pot1,grad1,hess1
      complex *16 pot2,grad2,hess2
      complex *16 charge1,dipstr1
      complex *16 charge2,dipstr2
      complex *16 z, cd, zk, ima, zinv, zinv2
c
      data ima/(0.0d0,1.0d0)/
c
c ... Calculate offsets and distance
c
      xdiff=source2(1)-source1(1)
      ydiff=source2(2)-source1(2)
ccc      rr=xdiff*xdiff+ydiff*ydiff
ccc      r=sqrt(rr)
c
      z=dcmplx(xdiff,ydiff)
c
      if (ifpot.eq.1) then
         pot1 = 0
         pot2 = 0
      endif
c
      if (ifgrad.eq.1) then
         grad1 = 0
         grad2 = 0
      endif
c
      if (ifhess.eq.1) then
         hess1 = 0
         hess2 = 0
      endif

        zinv=1/z
        zinv2=zinv*zinv

      if( ifcharge .eq. 1 ) then
c
      if (ifpot.eq.1) then
        cd = log(abs(z))
        pot2 = charge1*cd
        pot1 = charge2*cd
      endif
c
      if (ifgrad.eq.1) then
         cd = zinv
         grad2 = +charge1*cd
         grad1 = -charge2*cd
      endif
c
      if (ifhess.eq.1) then
         cd = zinv2
         hess2 = -charge1*cd
         hess1 = -charge2*cd
      endif
c
      endif
c
c
c
      if( ifdipole .eq. 1 ) then

      if (ifpot.eq.1) then
         cd = zinv
         pot2 = pot2+dipstr1*cd
         pot1 = pot1-dipstr2*cd
      endif
c
      if (ifgrad.eq.1) then
         cd = zinv2
         grad2 = grad2-dipstr1*cd
         grad1 = grad1-dipstr2*cd
      endif
c
      if (ifhess.eq.1) then
         cd = 2*zinv2*zinv
         hess2 = hess2+dipstr1*cd
         hess1 = hess1-dipstr2*cd
      endif
      
      endif
c
c
      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine cpotgrad2d_sdp_sym_add(source1,source2,
     $     ifcharge,charge1,charge2,ifdipole,dipstr1,dipstr2,
     1     ifpot,pot1,pot2,ifgrad,grad1,grad2,ifhess,hess1,hess2)
c**********************************************************************
c
c     This subroutine calculates the potential POT, gradient GRAD and
c     Hessian HESS at the target point TARGET, due to a charge and a
c     dipole at SOURCE. The scaling is that required of the delta
c     function response: i.e.,
c     
c              	pot = charge log(abs(z)) + dipstr 1/z
c		grad = gradient = d/dz
c		hess = Hessian = d^2/dz^2
c
c---------------------------------------------------------------------
c     INPUT:
c
c     ifpot      : flag for computing potential
c	                	ifpot = 0 -> don't compute 
c		                ifpot = 1 -> do compute 
c     ifgrad     : flag for computing gradient
c	                 	ifgrad = 0 -> don't compute 
c		                ifgrad = 1 -> do compute 
c     ifhess     : flag for computing hessian
c	                 	ifhess = 0 -> don't compute 
c		                ifhess = 1 -> do compute 
c     source    : location of the source 
c     charge    : charge strength
c     target    : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot       : calculated potential
c     grad      : calculated gradient
c     hess      : calculated hessian
c
c---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 source1(2),source2(2)
      complex *16 pot1,grad1,hess1
      complex *16 pot2,grad2,hess2
      complex *16 charge1,dipstr1
      complex *16 charge2,dipstr2
      complex *16 z, cd, zk, ima, zinv, zinv2
c
      data ima/(0.0d0,1.0d0)/
c
c ... Calculate offsets and distance
c
      xdiff=source2(1)-source1(1)
      ydiff=source2(2)-source1(2)
ccc      rr=xdiff*xdiff+ydiff*ydiff
ccc      r=sqrt(rr)
c
      z=dcmplx(xdiff,ydiff)
c
        zinv=1/z
        zinv2=zinv*zinv

      if( ifcharge .eq. 1 ) then
c
      if (ifpot.eq.1) then
        cd = log(abs(z))
        pot2 = pot2+charge1*cd
        pot1 = pot1+charge2*cd
      endif
c
      if (ifgrad.eq.1) then
         cd = zinv
         grad2 = grad2+charge1*cd
         grad1 = grad1-charge2*cd
      endif
c
      if (ifhess.eq.1) then
         cd = zinv2
         hess2 = hess2-charge1*cd
         hess1 = hess1-charge2*cd
      endif
c
      endif
c
c
c
      if( ifdipole .eq. 1 ) then

      if (ifpot.eq.1) then
         cd = zinv
         pot2 = pot2+dipstr1*cd
         pot1 = pot1-dipstr2*cd
      endif
c
      if (ifgrad.eq.1) then
         cd = zinv2
         grad2 = grad2-dipstr1*cd
         grad1 = grad1-dipstr2*cd
      endif
c
      if (ifhess.eq.1) then
         cd = 2*zinv2*zinv
         hess2 = hess2+dipstr1*cd
         hess1 = hess1-dipstr2*cd
      endif
      
      endif
c
c
      return
      end
c
c
c
c**********************************************************************
c
c       Complex valued Laplace particle direct evaluation routines
c
c**********************************************************************
      subroutine rcpotgrad2d_sdp(source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     1     target,ifpot,pot,ifgrad,grad,ifhess,hess)
c**********************************************************************
c
c     This subroutine calculates the potential POT, gradient GRAD and
c     Hessian HESS at the target point TARGET, due to a charge and a
c     dipole at SOURCE. The scaling is that required of the delta
c     function response: i.e.,
c     
c     pot = charge log(r) + dipstr (dipvec \dot \grad log(r) )
c     grad = gradient = (d/dx, d/dy) 
c     hess = Hessian = (d^2/dx^2, d^2/dxdy, d^2/dy^2) 
c
c---------------------------------------------------------------------
c     INPUT:
c
c     ifpot      : flag for computing potential
c	                 	ifpot = 0 -> don't compute 
c		                ifpot = 1 -> do compute 
c     ifgrad     : flag for computing gradient
c	                 	ifgrad = 0 -> don't compute 
c		                ifgrad = 1 -> do compute 
c     ifhess     : flag for computing hessian
c	                 	ifhess = 0 -> don't compute 
c		                ifhess = 1 -> do compute 
c     source    : location of the source 
c     charge    : charge strength
c     dipstr    : dipole strength
c     dipvec    : dipole orientation vector
c     target    : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot       : calculated potential
c     grad      : calculated gradient
c     hess      : calculated hessian
c
c---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 source(2),target(2),dipvec(2)
      complex *16 pot,grad(2),hess(3)
      complex *16 charge,dipstr
      complex *16 z, cd, zk, ima, zinv, zinv2
c
      data ima/(0.0d0,1.0d0)/
c
c ... Calculate offsets and distance
c
      xdiff=target(1)-source(1)
      ydiff=target(2)-source(2)
      rr=xdiff*xdiff+ydiff*ydiff
      r=sqrt(rr)
c
      z=dcmplx(xdiff,ydiff)
c
      if (ifpot.eq.1) then
         pot = 0
      endif
c
      if (ifgrad.eq.1) then
         grad(1) = 0
         grad(2) = 0
      endif
c
      if (ifhess.eq.1) then
         hess(1) = 0
         hess(2) = 0
         hess(3) = 0
      endif

c
      if( ifcharge .eq. 1 ) then
c
      if (ifpot.eq.1) pot=charge*log(r)
c
      if (ifgrad.eq.1) then
         cd = charge/rr
         grad(1) = xdiff*cd
         grad(2) = ydiff*cd
      endif
c
      if (ifhess.eq.1) then
         cd = charge/rr**2
         hess(1) = cd*(rr-2*xdiff*xdiff)
         hess(2) = cd*(  -2*xdiff*ydiff)
         hess(3) = cd*(rr-2*ydiff*ydiff)
      endif
c
      endif
c
c
c
      if( ifdipole .eq. 1 ) then

      if (ifpot.eq.1) then
         cd=dipstr/rr
         pot=pot-cd*(xdiff*dipvec(1)+ydiff*dipvec(2))
      endif
c
      if (ifgrad.eq.1) then
         cd = dipstr/rr**2
         derxx=+rr-2*xdiff*xdiff
         derxy=   -2*xdiff*ydiff
         deryy=+rr-2*ydiff*ydiff
         grad(1) = grad(1)-cd*(derxx*dipvec(1)+derxy*dipvec(2))
         grad(2) = grad(2)-cd*(derxy*dipvec(1)+deryy*dipvec(2))
      endif
c
      if (ifhess.eq.1) then
         cd = dipstr/rr**3
         derxxx=-6*xdiff*rr+8*xdiff**3
         derxxy=-2*ydiff*rr+8*xdiff**2*ydiff
         derxyy=-2*xdiff*rr+8*xdiff*ydiff**2
         deryyy=-6*ydiff*rr+8*ydiff**3
         hess(1) = hess(1)-cd*(derxxx*dipvec(1)+derxxy*dipvec(2))
         hess(2) = hess(2)-cd*(derxxy*dipvec(1)+derxyy*dipvec(2))
         hess(3) = hess(3)-cd*(derxyy*dipvec(1)+deryyy*dipvec(2))
      endif
      
      endif
c
c
      return
      end
c
c
c
c**********************************************************************
      subroutine rcpotgrad2d_sdp_add(source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     1     target,ifpot,pot,ifgrad,grad,ifhess,hess)
c**********************************************************************
c
c     This subroutine calculates the potential POT, gradient GRAD and
c     Hessian HESS at the target point TARGET, due to a charge and a
c     dipole at SOURCE. The scaling is that required of the delta
c     function response: i.e.,
c     
c     pot = charge log(r) + dipstr (dipvec \dot \grad log(r) )
c     grad = gradient = (d/dx, d/dy) 
c     hess = Hessian = (d^2/dx^2, d^2/dxdy, d^2/dy^2) 
c
c---------------------------------------------------------------------
c     INPUT:
c
c     ifpot      : flag for computing potential
c	                 	ifpot = 0 -> don't compute 
c		                ifpot = 1 -> do compute 
c     ifgrad     : flag for computing gradient
c	                 	ifgrad = 0 -> don't compute 
c		                ifgrad = 1 -> do compute 
c     ifhess     : flag for computing hessian
c	                 	ifhess = 0 -> don't compute 
c		                ifhess = 1 -> do compute 
c     source    : location of the source 
c     charge    : charge strength
c     dipstr    : dipole strength
c     dipvec    : dipole orientation vector
c     target    : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot       : calculated potential
c     grad      : calculated gradient
c     hess      : calculated hessian
c
c---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 source(2),target(2),dipvec(2)
      complex *16 pot,grad(2),hess(3)
      complex *16 charge,dipstr
      complex *16 z, cd, zk, ima, zinv, zinv2
c
      data ima/(0.0d0,1.0d0)/
c
c ... Calculate offsets and distance
c
      xdiff=target(1)-source(1)
      ydiff=target(2)-source(2)
      rr=xdiff*xdiff+ydiff*ydiff
      r=sqrt(rr)
c
c
      if( ifcharge .eq. 1 ) then
c
      if (ifpot.eq.1) pot=pot+charge*log(r)
c
      if (ifgrad.eq.1) then
         cd = charge/rr
         grad(1) = grad(1) + xdiff*cd
         grad(2) = grad(2) + ydiff*cd
      endif
c
      if (ifhess.eq.1) then
         cd = charge/rr**2
         hess(1) = hess(1) + cd*(rr-2*xdiff*xdiff)
         hess(2) = hess(2) + cd*(  -2*xdiff*ydiff)
         hess(3) = hess(3) + cd*(rr-2*ydiff*ydiff)
      endif
c
      endif
c
c
c
      if( ifdipole .eq. 1 ) then

      if (ifpot.eq.1) then
         cd=dipstr/rr
         pot=pot-cd*(xdiff*dipvec(1)+ydiff*dipvec(2))
      endif
c
      if (ifgrad.eq.1) then
         cd = dipstr/rr**2
         derxx=+rr-2*xdiff*xdiff
         derxy=   -2*xdiff*ydiff
         deryy=+rr-2*ydiff*ydiff
         grad(1) = grad(1)-cd*(derxx*dipvec(1)+derxy*dipvec(2))
         grad(2) = grad(2)-cd*(derxy*dipvec(1)+deryy*dipvec(2))
      endif
c
      if (ifhess.eq.1) then
         cd = dipstr/rr**3
         derxxx=-6*xdiff*rr+8*xdiff**3
         derxxy=-2*ydiff*rr+8*xdiff**2*ydiff
         derxyy=-2*xdiff*rr+8*xdiff*ydiff**2
         deryyy=-6*ydiff*rr+8*ydiff**3
         hess(1) = hess(1)-cd*(derxxx*dipvec(1)+derxxy*dipvec(2))
         hess(2) = hess(2)-cd*(derxxy*dipvec(1)+derxyy*dipvec(2))
         hess(3) = hess(3)-cd*(derxyy*dipvec(1)+deryyy*dipvec(2))
      endif
      
      endif
c
c
      return
      end
c
c
c
c**********************************************************************
c
c       Real valued Laplace particle direct evaluation routines
c
c**********************************************************************
      subroutine rpotgrad2d_sdp(source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     1     target,ifpot,pot,ifgrad,grad,ifhess,hess)
c**********************************************************************
c
c     This subroutine calculates the potential POT, gradient GRAD and
c     Hessian HESS at the target point TARGET, due to a charge and a
c     dipole at SOURCE. The scaling is that required of the delta
c     function response: i.e.,
c     
c     pot = charge log(r) + dipstr (dipvec \dot \grad log(r) )
c     grad = gradient = (d/dx, d/dy) 
c     hess = Hessian = (d^2/dx^2, d^2/dxdy, d^2/dy^2) 
c
c---------------------------------------------------------------------
c     INPUT:
c
c     ifpot      : flag for computing potential
c	                 	ifpot = 0 -> don't compute 
c		                ifpot = 1 -> do compute 
c     ifgrad     : flag for computing gradient
c	                 	ifgrad = 0 -> don't compute 
c		                ifgrad = 1 -> do compute 
c     ifhess     : flag for computing hessian
c	                 	ifhess = 0 -> don't compute 
c		                ifhess = 1 -> do compute 
c     source    : location of the source 
c     charge    : charge strength
c     dipstr    : dipole strength
c     dipvec    : dipole orientation vector
c     target    : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot       : calculated potential
c     grad      : calculated gradient
c     hess      : calculated hessian
c
c---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 source(2),target(2),dipvec(2)
      real *8 pot,grad(2),hess(3)
      real *8 charge,dipstr
      real *8 z, cd, zk, ima, zinv, zinv2
c
      data ima/(0.0d0,1.0d0)/
c
c ... Calculate offsets and distance
c
      xdiff=target(1)-source(1)
      ydiff=target(2)-source(2)
      rr=xdiff*xdiff+ydiff*ydiff
      r=sqrt(rr)
c
      if (ifpot.eq.1) then
         pot = 0
      endif
c
      if (ifgrad.eq.1) then
         grad(1) = 0
         grad(2) = 0
      endif
c
      if (ifhess.eq.1) then
         hess(1) = 0
         hess(2) = 0
         hess(3) = 0
      endif

c
      if( ifcharge .eq. 1 ) then
c
      if (ifpot.eq.1) pot=charge*log(r)
c
      if (ifgrad.eq.1) then
         cd = charge/rr
         grad(1) = xdiff*cd
         grad(2) = ydiff*cd
      endif
c
      if (ifhess.eq.1) then
         cd = charge/rr**2
         hess(1) = cd*(rr-2*xdiff*xdiff)
         hess(2) = cd*(  -2*xdiff*ydiff)
         hess(3) = cd*(rr-2*ydiff*ydiff)
      endif
c
      endif
c
c
c
      if( ifdipole .eq. 1 ) then

      if (ifpot.eq.1) then
         cd=dipstr/rr
         pot=pot-cd*(xdiff*dipvec(1)+ydiff*dipvec(2))
      endif
c
      if (ifgrad.eq.1) then
         cd = dipstr/rr**2
         derxx=+rr-2*xdiff*xdiff
         derxy=   -2*xdiff*ydiff
         deryy=+rr-2*ydiff*ydiff
         grad(1) = grad(1)-cd*(derxx*dipvec(1)+derxy*dipvec(2))
         grad(2) = grad(2)-cd*(derxy*dipvec(1)+deryy*dipvec(2))
      endif
c
      if (ifhess.eq.1) then
         cd = dipstr/rr**3
         derxxx=-6*xdiff*rr+8*xdiff**3
         derxxy=-2*ydiff*rr+8*xdiff**2*ydiff
         derxyy=-2*xdiff*rr+8*xdiff*ydiff**2
         deryyy=-6*ydiff*rr+8*ydiff**3
         hess(1) = hess(1)-cd*(derxxx*dipvec(1)+derxxy*dipvec(2))
         hess(2) = hess(2)-cd*(derxxy*dipvec(1)+derxyy*dipvec(2))
         hess(3) = hess(3)-cd*(derxyy*dipvec(1)+deryyy*dipvec(2))
      endif
      
      endif
c
c
      return
      end
c
c
c
c**********************************************************************
      subroutine rpotgrad2d_sdp_add(source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     1     target,ifpot,pot,ifgrad,grad,ifhess,hess)
c**********************************************************************
c
c     This subroutine calculates the potential POT, gradient GRAD and
c     Hessian HESS at the target point TARGET, due to a charge and a
c     dipole at SOURCE. The scaling is that required of the delta
c     function response: i.e.,
c     
c     pot = charge log(r) + dipstr (dipvec \dot \grad log(r) )
c     grad = gradient = (d/dx, d/dy) 
c     hess = Hessian = (d^2/dx^2, d^2/dxdy, d^2/dy^2) 
c
c---------------------------------------------------------------------
c     INPUT:
c
c     ifpot      : flag for computing potential
c	                 	ifpot = 0 -> don't compute 
c		                ifpot = 1 -> do compute 
c     ifgrad     : flag for computing gradient
c	                 	ifgrad = 0 -> don't compute 
c		                ifgrad = 1 -> do compute 
c     ifhess     : flag for computing hessian
c	                 	ifhess = 0 -> don't compute 
c		                ifhess = 1 -> do compute 
c     source    : location of the source 
c     charge    : charge strength
c     dipstr    : dipole strength
c     dipvec    : dipole orientation vector
c     target    : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot       : calculated potential
c     grad      : calculated gradient
c     hess      : calculated hessian
c
c---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 source(2),target(2),dipvec(2)
      real *8 pot,grad(2),hess(3)
      real *8 charge,dipstr
      real *8 z, cd, zk, ima, zinv, zinv2
c
      data ima/(0.0d0,1.0d0)/
c
c ... Calculate offsets and distance
c
      xdiff=target(1)-source(1)
      ydiff=target(2)-source(2)
      rr=xdiff*xdiff+ydiff*ydiff
      r=sqrt(rr)
c
c
      if( ifcharge .eq. 1 ) then
c
      if (ifpot.eq.1) pot=pot+charge*log(r)
c
      if (ifgrad.eq.1) then
         cd = charge/rr
         grad(1) = grad(1) + xdiff*cd
         grad(2) = grad(2) + ydiff*cd
      endif
c
      if (ifhess.eq.1) then
         cd = charge/rr**2
         hess(1) = hess(1) + cd*(rr-2*xdiff*xdiff)
         hess(2) = hess(2) + cd*(  -2*xdiff*ydiff)
         hess(3) = hess(3) + cd*(rr-2*ydiff*ydiff)
      endif
c
      endif
c
c
c
      if( ifdipole .eq. 1 ) then

      if (ifpot.eq.1) then
         cd=dipstr/rr
         pot=pot-cd*(xdiff*dipvec(1)+ydiff*dipvec(2))
      endif
c
      if (ifgrad.eq.1) then
         cd = dipstr/rr**2
         derxx=+rr-2*xdiff*xdiff
         derxy=   -2*xdiff*ydiff
         deryy=+rr-2*ydiff*ydiff
         grad(1) = grad(1)-cd*(derxx*dipvec(1)+derxy*dipvec(2))
         grad(2) = grad(2)-cd*(derxy*dipvec(1)+deryy*dipvec(2))
      endif
c
      if (ifhess.eq.1) then
         cd = dipstr/rr**3
         derxxx=-6*xdiff*rr+8*xdiff**3
         derxxy=-2*ydiff*rr+8*xdiff**2*ydiff
         derxyy=-2*xdiff*rr+8*xdiff*ydiff**2
         deryyy=-6*ydiff*rr+8*ydiff**3
         hess(1) = hess(1)-cd*(derxxx*dipvec(1)+derxxy*dipvec(2))
         hess(2) = hess(2)-cd*(derxxy*dipvec(1)+derxyy*dipvec(2))
         hess(3) = hess(3)-cd*(derxyy*dipvec(1)+deryyy*dipvec(2))
      endif
      
      endif
c
c
      return
      end
c
c
c
