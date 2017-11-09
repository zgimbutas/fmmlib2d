cc Copyright (C) 2010: Vladimir Rokhlin
cc Copyright (C) 2009-2012: Leslie Greengard and Zydrunas Gimbutas
cc Contact: greengard@cims.nyu.edu
cc 
cc This software is being released under a modified FreeBSD license
cc (see COPYING in home directory). 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date: 2011-02-22 17:34:23 -0500 (Tue, 22 Feb 2011) $
c    $Revision: 1670 $
c
c
c-----------------------------------------------------------------------------
c
c      h2dterms - determine number of terms in mpole expansions for box
c           of size "size" with Helmholtz parameter zk.
c
c      h2dterms_list2 - build the number of terms table for all boxes 
c           in list 2
c
c      h2dterms_list2e - build the number of terms table for all boxes 
c           in extended list 2
c
c-----------------------------------------------------------------------------
c
c
c
      subroutine h2dterms(size, zk, eps, nterms, ier)
      implicit real *8 (a-h,o-z)
c
c     Determine number of terms in mpole expansions for box of size
c     "size" with Helmholtz parameter zk.
c
c     The method is based on examining the decay of H_n * J_n.
c
c     Maximum number of terms is 10000, which 
c     works for boxes up to 1000 wavelengths in size     
c
c
c-----------------------------------------------------------------------------
c
      integer iscale(0:20000)
      complex *16  zk, z1, z2, z3, jfun(0:20000), ht0,
     1             ht1, ht2, fjder(0:1), ztmp,
     1             hfun(0:20000), fhder(0:1)
c
      ier = 0
c
      z1 = (zk*size)*1.5d0
c
c       the code will run out memory if frequency is too small 
c       set frequency to something more reasonable, nterms is 
c       approximately the same for all small frequencies
c       
ccc        if( abs(z1) .lt. 1d-6 ) z1=1d-6
c
      ntmax = 10000
      ifder = 0
      rscale = 1.0d0
      if (cdabs(zk*size) .lt. 1.0d0) rscale = cdabs(zk*size)
      call h2dall(ntmax,z1,rscale,hfun,ifder,fhder)
c
      z2 = (zk*size) * dsqrt(2d0)/2.d0
c
      ier1 = 0
c
      call jfuns2d(ier1, ntmax, z2, rscale, jfun, ifder, fjder,
     1                 20000, iscale, ntop)
c
c     set error flag if jfuns runs out of memory
c
      if (ier1.eq.8) then 
        ier = 11 
        return
      endif        
c
      xtemp1 = cdabs(jfun(0)*hfun(0))
      xtemp2 = cdabs(jfun(1)*hfun(1))
      xtemp0 = xtemp1+xtemp2
      nterms = 1
      do j = 2, ntmax
        xtemp1 = cdabs(jfun(j)*hfun(j))
        xtemp2 = cdabs(jfun(j-1)*hfun(j-1))
        xtemp = xtemp1+xtemp2
        if(xtemp .lt. eps*xtemp0)then
          nterms = j + 1
          return
        endif
c
      enddo
c
c       ... computational box is too big, set nterms to 10001
c
        ier = 13
        nterms=10001
c
      return
      end
c
c
c
c
c
      subroutine h2dterms_far(size, zk, eps, nterms, ier)
      implicit real *8 (a-h,o-z)
c
c
c     Determine number of terms in mpole expansions for box of size
c     "size" with Helmholtz parameter zk. 
c
c     The method is based on examining the decay of H_n * J_n.
c
c     This routine aasumes slightly larger separation of boxes: the
c     first unit box is located at the origin and the second box is
c     located at (3,0,0).
c
c     Maximum number of terms is 10000, which 
c     works for boxes up to 1000 wavelengths in size     
c
c
c-----------------------------------------------------------------------------
c
      integer iscale(0:20000)
      complex *16  zk, z1, z2, z3, jfun(0:20000), ht0,
     1             ht1, ht2, fjder(0:1), ztmp,
     1             hfun(0:20000), fhder(0:1)
c
      ier = 0
c
      z1 = (zk*size)*2.5d0
c
c       the code will run out memory if frequency is too small 
c       set frequency to something more reasonable, nterms is 
c       approximately the same for all small frequencies
c       
ccc        if( abs(z1) .lt. 1d-6 ) z1=1d-6
c
      ntmax = 10000
      ifder = 0
      rscale = 1.0d0
      if (cdabs(zk*size) .lt. 1.0d0) rscale = cdabs(zk*size)
      call h2dall(ntmax,z1,rscale,hfun,ifder,fhder)
c
      z2 = (zk*size) * dsqrt(2d0)/2.d0
c
      ier1 = 0
c
      call jfuns2d(ier1, ntmax, z2, rscale, jfun, ifder, fjder,
     1                 20000, iscale, ntop)
ccc      call prin2(' jfun is *',jfun,2*ntmax+2)
c
c     set error flag if jfuns runs out of memory
c
      if (ier1.eq.8) then 
        ier = 11 
        return
      endif        
c
      xtemp1 = cdabs(jfun(0)*hfun(0))
      xtemp2 = cdabs(jfun(1)*hfun(1))
      xtemp0 = xtemp1+xtemp2
      nterms = 1
      do j = 2, ntmax
        xtemp1 = cdabs(jfun(j)*hfun(j))
        xtemp2 = cdabs(jfun(j-1)*hfun(j-1))
        xtemp = xtemp1+xtemp2
        if(xtemp .lt. eps*xtemp0)then
          nterms = j + 1
          return
        endif
      enddo
c
c       ... computational box is too big, set nterms to 10001
c
        ier = 13
        nterms=10001
c
      return
      end
c
c
c
c
c
      subroutine h2dterms_list2(size, zk, eps, itable, ier)
      implicit real *8 (a-h,o-z)
c
c
c     Determine number of terms in mpole expansions for box of size
c     "size" with Helmholtz parameter zk.
c
c     The method is based on examining the decay of h_n * j_n.
c
c     Build nterms table for all boxes in list 2
c
c     Maximum number of terms is 1000, which 
c     works for boxes up to 160 wavelengths in size     
c
c
c-----------------------------------------------------------------------------
c
      integer iscale(0:20000)
      complex *16  zk, z1, z2, z3, jfun(0:20000), ht0,
     1             ht1, ht2, fjder(0:1), ztmp,
     1             hfun(0:20000), fhder(0:1)
c
      integer nterms_table(2:3,0:3)
      integer itable(-3:3,-3:3)
c
      ier = 0
c
c
c       first box in x direction
c
ccc      z1 = (zk*size)*1.5d0
c
c       second box in x direction
c
ccc      z1 = (zk*size)*2.5d0
c
c       on the diagonal
c
c      z1 = (zk*size)*dsqrt(2d0)/2.d0*3
c      z1 = (zk*size)*dsqrt(2d0)/2.d0*5
c
        do 1800 ii=2,3
        do 1800 jj=0,3
c
        dx=ii
        dy=jj
c       
        if( dx .gt. 0 ) dx=dx-.5
        if( dy .gt. 0 ) dy=dy-.5
c
        rr=sqrt(dx*dx+dy*dy)
ccc        call prin2('rr=*',rr,1)
ccc        call prin2('rr=*',sqrt(2.0d0)/2*5,1)
c
        z1 = (zk*size)*rr

c       the code will run out memory if frequency is too small 
c       set frequency to something more reasonable, nterms is 
c       approximately the same for all small frequencies
c  
ccc      if( abs(z1) .lt. 1d-6 ) z1=1d-6
c
      ntmax = 10000
      ifder = 0
      rscale = 1.0d0
      if (cdabs(zk*size) .lt. 1.0d0) rscale = cdabs(zk*size)
      call h2dall(ntmax,z1,rscale,hfun,ifder,fhder)
c
      z2 = (zk*size) * dsqrt(2d0)/2.d0
c
      ier1 = 0
c
      call jfuns2d(ier1, ntmax, z2, rscale, jfun, ifder, fjder,
     1                 20000, iscale, ntop)
ccc      call prin2(' jfun is *',jfun,2*ntmax+2)
c
c     set error flag if jfuns runs out of memory
c
      if (ier1.eq.8) then 
        ier = 11 
        return
      endif        
c
      xtemp1 = cdabs(jfun(0)*hfun(0))
      xtemp2 = cdabs(jfun(1)*hfun(1))
      xtemp0 = xtemp1+xtemp2
      nterms = 1
      do j = 2, ntmax
        xtemp1 = cdabs(jfun(j)*hfun(j))
        xtemp2 = cdabs(jfun(j-1)*hfun(j-1))
        xtemp = xtemp1+xtemp2
        if(xtemp .lt. eps*xtemp0)then
          nterms = j + 1
          goto 1600
        endif
      enddo
c
c       ... computational box is too big, set nterms to 10001
c
        ier = 13
        nterms=10001
c
 1600   continue
c
        nterms_table(ii,jj)=nterms
c
 1800   continue
c
ccc        call prinf('nterms=*',nterms_table,2*4*4)
c
c       build the rank table for all boxes in list 2
c
        do i=-3,3
        do j=-3,3
        itable(i,j)=0
        enddo
        enddo
c
        do 2200 i=-3,3
        do 2200 j=-3,3
c
        if( abs(i) .gt. 1 ) then
        itable(i,j)=nterms_table(abs(i),abs(j))
        else if( abs(j) .gt. 1) then
        itable(i,j)=nterms_table(abs(j),abs(i))
        endif
c
 2200   continue
c
      return
      end
c
c
c
c
c
      subroutine h2dterms_list2e(size, zk, eps, itable, ier)
      implicit real *8 (a-h,o-z)
c
c
c     Determine number of terms in mpole expansions for box of size
c     "size" with Helmholtz parameter zk.
c
c     The method is based on examining the decay of h_n * j_n.
c
c     Build nterms table for all boxes in extended list 2
c
c     Maximum number of terms is 10000, which 
c     works for boxes up to 1000 wavelengths in size     
c
c
c-----------------------------------------------------------------------------
c
      integer iscale(0:20000)
      complex *16  zk, z1, z2, z3, jfun(0:20000), ht0,
     1             ht1, ht2, fjder(0:1), ztmp,
     1             hfun(0:20000), fhder(0:1)
c
      integer nterms_table(2:7,0:7)
      integer itable(-7:7,-7:7)
c
      ier = 0
c
c
c       first box in x direction
c
ccc      z1 = (zk*size)*1.5d0
c
c       second box in x direction
c
ccc      z1 = (zk*size)*2.5d0
c
c       on the diagonal
c
c      z1 = (zk*size)*dsqrt(2d0)/2.d0*3
c      z1 = (zk*size)*dsqrt(2d0)/2.d0*5
c
        do 1800 ii=2,7
        do 1800 jj=0,7
c
        dx=ii
        dy=jj
c       
        if( dx .gt. 0 ) dx=dx-.5
        if( dy .gt. 0 ) dy=dy-.5
c
        rr=sqrt(dx*dx+dy*dy)
ccc        call prin2('rr=*',rr,1)
ccc        call prin2('rr=*',sqrt(2.0d0)/2*5,1)
c
        z1 = (zk*size)*rr

c       the code will run out memory if frequency is too small 
c       set frequency to something more reasonable, nterms is 
c       approximately the same for all small frequencies
c  
ccc      if( abs(z1) .lt. 1d-6 ) z1=1d-6
c
      ntmax = 10000
      ifder = 0
      rscale = 1.0d0
      if (cdabs(zk*size) .lt. 1.0d0) rscale = cdabs(zk*size)
      call h2dall(ntmax,z1,rscale,hfun,ifder,fhder)
c
      z2 = (zk*size) * dsqrt(2d0)/2.d0
c
      ier1 = 0
c
      call jfuns2d(ier1, ntmax, z2, rscale, jfun, ifder, fjder,
     1                 20000, iscale, ntop)
ccc      call prin2(' jfun is *',jfun,2*ntmax+2)
c
c     set error flag if jfuns runs out of memory
c
      if (ier1.eq.8) then 
        ier = 11 
        return
      endif        
c
      xtemp1 = cdabs(jfun(0)*hfun(0))
      xtemp2 = cdabs(jfun(1)*hfun(1))
      xtemp0 = xtemp1+xtemp2
      nterms = 1
      do j = 2, ntmax
        xtemp1 = cdabs(jfun(j)*hfun(j))
        xtemp2 = cdabs(jfun(j-1)*hfun(j-1))
        xtemp = xtemp1+xtemp2
        if(xtemp .lt. eps*xtemp0)then
          nterms = j + 1
          goto 1600
        endif
      enddo
c
c       ... computational box is too big, set nterms to 10001
c
        ier = 13
        nterms=10001
c
 1600   continue
c
        nterms_table(ii,jj)=nterms
c
 1800   continue
c
ccc        call prinf('nterms=*',nterms_table,2*4*4)
c
c       build the rank table for all boxes in extended list 2
c
        do i=-7,7
        do j=-7,7
        itable(i,j)=0
        enddo
        enddo
c
        do 2200 i=-7,7
        do 2200 j=-7,7
c
        if( abs(i) .gt. 2 ) then
        itable(i,j)=nterms_table(abs(i),abs(j))
        else if( abs(j) .gt. 2) then
        itable(i,j)=nterms_table(abs(j),abs(i))
        endif
c
 2200   continue
c
      return
      end
c
c
c
c
c
      subroutine h2dterms_eval(itype, size, zk, eps, nterms, ier)
      implicit real *8 (a-h,o-z)
c
c
c     Determine number of terms in mpole expansions for box of size
c     "size" with Helmholtz parameter zk.
c
c     The method is based on examining the decay of h_n * j_n.
c
c     Maximum number of terms is 10000, which 
c     works for boxes up to 1000 wavelengths in size     
c
c
c-----------------------------------------------------------------------------
c
      integer iscale(0:20000)
      complex *16  zk, z1, z2, z3, jfun(0:20000), ht0,
     1             ht1, ht2, fjder(0:1), ztmp,
     1             hfun(0:20000), fhder(0:1)
c
      ier = 0
c
      z1 = (zk*size)*1.5d0
c
c       the code will run out memory if frequency is too small 
c       set frequency to something more reasonable, nterms is 
c       approximately the same for all small frequencies
c       
ccc      if( abs(z1) .lt. 1d-6 ) z1=1d-6
c
      ntmax = 1000
      ifder = 0
      rscale = 1.0d0
      if (cdabs(zk*size) .lt. 1.0d0) rscale = cdabs(zk*size)
      call h2dall(ntmax,z1,rscale,hfun,ifder,fhder)
c
        z2 = (zk*size) * dsqrt(2d0)/2.d0
c
c       corners included
        if( itype .eq. 1 ) z2 = (zk*size) * dsqrt(2d0)/2.d0
c       edges included, no corners
        if( itype .eq. 2 ) z2 = (zk*size) * dsqrt(1d0)/2.d0
c       center only
        if( itype .eq. 3 ) z2 = (zk*size) * 1.0d0/2.d0
c       center only, small interior sphere
        if( itype .eq. 4 ) z2 = (zk*size) * 0.8d0/2.d0
c
      ier1 = 0
c
      call jfuns2d(ier1, ntmax, z2, rscale, jfun, ifder, fjder,
     1                 2000, iscale, ntop)
ccc      call prin2(' jfun is *',jfun,2*ntmax+2)
c
c     set error flag if jfuns runs out of memory
c
      if (ier1.eq.8) then 
        ier = 11 
        return
      endif        
c
      xtemp1 = cdabs(jfun(0)*hfun(0))
      xtemp2 = cdabs(jfun(1)*hfun(1))
      xtemp0 = xtemp1+xtemp2
      nterms = 1
      do j = 2, ntmax
        xtemp1 = cdabs(jfun(j)*hfun(j))
        xtemp2 = cdabs(jfun(j-1)*hfun(j-1))
        xtemp = xtemp1+xtemp2
        if(xtemp .lt. eps*xtemp0)then
          nterms = j + 1
          return
        endif
c
      enddo
c
c       ... computational box is too big, set nterms to 10001
c
        ier = 13
        nterms=10001
c
      return
      end
c
c
c
c
c

