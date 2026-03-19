cc Copyright (C) 2010-2011: Leslie Greengard and Zydrunas Gimbutas
cc SPDX-License-Identifier: GPL-2.0-or-later
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        this is the end of the debugging code and the beginning
c        of the routines for Helmholtz FMM in R^2
c
c        FAR-FIELD + EXPONENTIAL VERSION with nd-density support
c
c        This file contains hfmm2d_list2_far_exp which performs Steps 3-5
c        of the FMM using far-field signature functions for M2L at
c        coarse levels AND plane wave exponential expansions for M2L
c        at fine levels.
c        Supports vectorized evaluations with nd independent densities.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine hfmm2d_list2_far_exp(nd,
     $     zk,bsize,nlev,laddr,scale,nterms,rmlexp,iaddr,epsfmm,
     $     timeinfo,wlists,mptemp,lmptemp,
     $     ifprune_list2,nboxes,iprec,nterms_far,nterms_tra,
     $     min_depth)
        implicit real *8 (a-h,o-z)
c
        integer iaddr(4,*),laddr(2,*)
        integer nterms(0:*),nterms_far(0:*),nterms_tra(0:*)
        complex *16 rmlexp(*)
        real *8 scale(0:*)
        integer itable(-3:3,-3:3)
c
        integer list(10 000)
c
        integer box(15)
        real *8 bsize(0:200)
        real *8 center0(2),corners0(2,4)
c
        integer box1(15)
        real *8 center1(2),corners1(2,4)
c
        real *8 wlists(*)
        complex *16 mptemp(lmptemp)
        complex *16 zk,zkbsize
c
        complex *16 ptemp,ftemp(2),htemp(3)
c
        real *8 timeinfo(10)
c
        integer cache_init(-3:3,-3:3,0:200)
        integer, allocatable :: it_addr(:,:,:)
        complex *16, allocatable :: trans_cache(:)
c
        integer, allocatable :: nexp(:)
        complex *16, allocatable :: zcc(:)
c
        complex *16, allocatable :: ws(:), xs(:)
        complex *16, allocatable :: cshift(:,:,:),cscale(:)
c
        integer, allocatable :: iq_addr(:,:)
        complex *16, allocatable :: quadr(:)

        integer, allocatable :: ie_addr(:,:)
        complex *16, allocatable :: expio(:)
        complex *16, allocatable :: hexp_tmp(:)
        complex *16, allocatable :: lexp_tmp(:)
c
        done = 1.0d0
        pi=4*atan(done)
c
c
         t1=second()
C$        t1=omp_get_wtime()

c       ... preallocate work and indexing arrays for exp
c
        allocate(nexp(0:200))
        allocate(zcc(0:200))

        allocate(xs(400))
        allocate(ws(400))
c
        ifprint = 1
        rsize_wavelengths = abs(zk*bsize(0))/(2*pi)
        wb_cutoff = 8.0d0

        min_exp_depth = log(rsize_wavelengths/wb_cutoff)/log(2.0d0)+2
        if( min_exp_depth .le. 3 ) min_exp_depth = 3

        if( iprec .lt. 0 ) min_exp_depth = nlev+2
        if( iprec .gt. 5 ) min_exp_depth = nlev+2
        if( imag(zk) .lt. 0 ) min_exp_depth = nlev+2

        if( ifprint .eq. 1 ) then
        call prin2('rsize_wavelengths=*',rsize_wavelengths,1)
        call prin2('wb_cutoff=*',wb_cutoff,1)
        call prinf('min_exp_depth=*',min_exp_depth,1)
        endif

        do i = 2,nlev
        call hwtsexp2(ier,iprec,zk*bsize(i),
     $     xs,ws,zcc(i),nexp(i),err)
        call prinf('after hwtsexp2, nexp=*',nexp(i),1)
        enddo

        allocate(iq_addr(8,0:200))
        allocate(ie_addr(2,nboxes))

        do i = 0,200
        do j = 1,8
        iq_addr(j,i) = 0
        enddo
        enddo
c
        do i = 1,nboxes
        ie_addr(1,i) = 0
        ie_addr(2,i) = 0
        enddo

c       ... allocate exponential nodes, weights, precomputed arrays
c
        iptr = 1
        do i = 2,nlev
        iq_addr(1,i) = iptr
        iptr = iptr + nexp(i)
        iq_addr(2,i) = iptr
        iptr = iptr + nexp(i)
        iq_addr(3,i) = iptr
        iptr = iptr + nexp(i)
        iq_addr(4,i) = iptr
        iptr = iptr + nexp(i)*7*7
        iq_addr(5,i) = iptr
        iptr = iptr + nexp(i)*nterms(i)
        iq_addr(6,i) = iptr
        iptr = iptr + nexp(i)*nterms(i)
        iq_addr(7,i) = iptr
        iptr = iptr + nexp(i)*nterms(i)
        iq_addr(8,i) = iptr
        iptr = iptr + nexp(i)*nterms(i)
        enddo

        allocate( quadr(iptr) )

        do i = 2,nlev
        call hwtsexp2(ier,iprec,zk*bsize(i),
     $     quadr(iq_addr(1,i)),quadr(iq_addr(2,i)),
     $     zcc(i),nexp(i),err)
        call h2expini(zk*bsize(i),scale(i),nterms(i),
     $     quadr(iq_addr(1,i)),quadr(iq_addr(2,i)),
     $     zcc(i),nexp(i),
     $     quadr(iq_addr(3,i)),quadr(iq_addr(4,i)),
     $     quadr(iq_addr(5,i)),quadr(iq_addr(6,i)),
     $     quadr(iq_addr(7,i)),quadr(iq_addr(8,i)))
        enddo

c       ... allocate exponentials (nd-vectorized)
c
        iptr=1
        do ilev=2,nlev
        do ibox=laddr(1,ilev+1),laddr(1,ilev+1)+laddr(2,ilev+1)-1
        ie_addr(1,ibox) = iptr
        iptr = iptr + (nexp(ilev)+1)*4*nd
        ie_addr(2,ibox) = iptr
        iptr = iptr + (nexp(ilev)+1)*4*nd
        enddo
        enddo

        allocate( expio(iptr) )

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i = 1,iptr
        expio(i)=0
        enddo
C$OMP END PARALLEL DO

c       ... initialize far field translation cache
c
        allocate(it_addr(-3:3,-3:3,0:200))
c
        do ilev=0,200
        do j = -3,3
        do i = -3,3
        cache_init(i,j,ilev) = 0
        enddo
        enddo
        enddo
c
c       ... allocate far field translation cache
c
        iptr=1
        do ilev=2,min_depth+1
        do j = -3,3
        do i = -3,3
        it_addr(i,j,ilev) = iptr
        iptr = iptr + nterms_far(ilev)
        enddo
        enddo
        enddo
c
        allocate( trans_cache(iptr) )

         t2=second()
C$        t2=omp_get_wtime()
         call prin2('time in list2_far_exp_ini=*',t2-t1,1)

c
        ifprint=1
c
         if (ifprint .ge. 1)
     $     call prinf('=== STEP 3 (merge mp) ===*',i,0)
         t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 3, merge all multipole expansions
c
         do 2300 ilev=nlev,max(3,min_depth+1),-1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level0,level,npts,nkids,radius)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1)
C$OMP$PRIVATE(mptemp,lused,ier,i,j,ptemp,ftemp,cd)
         do 2200 ibox=laddr(1,ilev),laddr(1,ilev)+laddr(2,ilev)-1
c
         call d2tgetb(ier,ibox,box,center0,corners0,wlists)
         call d2tnkids(box,nkids)
c
         if( box(10) .eq. 0 ) goto 2200
c
         if (nkids .ne. 0) then
         level0=box(1)
         if( level0 .ge. 2 ) then
            radius = (corners0(1,1) - center0(1))**2
            radius = radius + (corners0(2,1) - center0(2))**2
            radius = sqrt(radius)
c
            call h2dzero(nd,rmlexp(iaddr(1,ibox)),nterms(level0))
c
            do 2100 i = 1,4
               jbox = box(4+i)
               if (jbox.eq.0) goto 2100
               call d2tgetb(ier,jbox,box1,center1,corners1,wlists)
               level1=box1(1)
               call h2dmpmp(nd,zk,scale(level1),center1,
     1            rmlexp(iaddr(1,jbox)),nterms(level1),scale(level0),
     1            center0,mptemp,nterms(level0))
               call h2dadd(nd,mptemp,rmlexp(iaddr(1,ibox)),
     1            nterms(level0))
 2100       continue
            endif
         endif
 2200    continue
C$OMP END PARALLEL DO
 2300    continue
c
c
        if (ifprint .ge. 1)
     $     call prinf('=== STEP 3af (mp to far) ===*',i,0)
c
c       ... step 3af, form far field signatures
c
         if( min_depth .ge. 2 ) then
         do ilev=min_depth+1,min_depth+1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level0,ier)
         do ibox=laddr(1,ilev),laddr(1,ilev)+laddr(2,ilev)-1
         call d2tgetb(ier,ibox,box,center0,corners0,wlists)
         if( box(10) .gt. 0 ) then
         level0=box(1)
         if( level0 .ge. 2 ) then
         call h2dzero_far(rmlexp(iaddr(3,ibox)),nterms_far(level0))
         call h2dmpfar(rmlexp(iaddr(1,ibox)),nterms(level0),
     $      rmlexp(iaddr(3,ibox)),nterms_far(level0))
         endif
         endif
         enddo
C$OMP END PARALLEL DO
         enddo
         endif
c
c
         if (ifprint .ge. 1)
     $     call prinf('=== STEP 3bf (merge far) ===*',i,0)
c
         do 2320 ilev=min_depth,3,-1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level0,nkids)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1,ier,i)
         do 2220 ibox=laddr(1,ilev),laddr(1,ilev)+laddr(2,ilev)-1
         call d2tgetb(ier,ibox,box,center0,corners0,wlists)
         call d2tnkids(box,nkids)
         if( box(10) .eq. 0 ) goto 2220
         if (nkids .ne. 0) then
         level0=box(1)
         if( level0 .ge. 2 ) then
            call h2dzero_far(rmlexp(iaddr(3,ibox)),nterms_far(level0))
            do 2120 i = 1,4
               jbox = box(4+i)
               if (jbox.eq.0) goto 2120
               call d2tgetb(ier,jbox,box1,center1,corners1,wlists)
               level1=box1(1)
               call h2dfarfar(zk,
     1            center1,rmlexp(iaddr(3,jbox)),nterms_far(level1),
     1            center0,rmlexp(iaddr(3,ibox)),nterms_far(level0))
 2120       continue
            endif
         endif
 2220    continue
C$OMP END PARALLEL DO
 2320    continue
c
         t2=second()
C$        t2=omp_get_wtime()
         timeinfo(3)=t2-t1
c
c
        if (ifprint .ge. 1)
     $     call prinf('=== STEP 4 (mp to lo) ===*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 4, standard M2L between far-field and exponential regimes
c
        do 4310 ilev=min_depth+2,min(min_exp_depth-1,nlev+1)
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level0,list,nlist,itype)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1,ifdirect2)
C$OMP$PRIVATE(mptemp,ier,ilist,if_use_trunc)
C$OMP$SCHEDULE(DYNAMIC)
        do 4210 ibox=laddr(1,ilev),laddr(1,ilev)+laddr(2,ilev)-1
        call d2tgetb(ier,ibox,box,center0,corners0,wlists)
        level0=box(1)
        if (level0 .ge. 2) then
           itype=2
           call d2tgetl(ier,ibox,itype,list,nlist,wlists)
           do 4151 ilist=1,nlist
              jbox=list(ilist)
              call d2tgetb(ier,jbox,box1,center1,corners1,wlists)
              if( box1(10) .eq. 0 ) goto 4151
              if (jbox.eq.0) goto 4151
              if ((box(12).eq.0).and.(ifprune_list2.eq.1))
     $           goto 4151
              level1=box1(1)
              call h2dzero(nd,mptemp,nterms(level1))
              call h2dmploc(nd,zk,scale(level1),center1,
     1           rmlexp(iaddr(1,jbox)),nterms(level1),
     $           scale(level0),
     1           center0,mptemp,nterms(level0))
              call h2dadd(nd,mptemp,rmlexp(iaddr(2,ibox)),
     1           nterms(level0))
 4151         continue
        endif
 4210   continue
C$OMP END PARALLEL DO
 4310   continue

c
c       ... step 4a, form outgoing exponential expansions
c
        if(ifprint .ge. 1)
     $     call prinf('=== STEP 4a (form expout) ====*',i,0)
        t1a=second()
C$        t1a=omp_get_wtime()
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level,npts,zkbsize,idim)
C$OMP$PRIVATE(hexp_tmp)
C$OMP$SCHEDULE(DYNAMIC)
        do 1200 ibox=1,nboxes
        call d2tgetb(ier,ibox,box,center0,corners0,wlists)
        call d2tnkids(box,nkids)
        level=box(1)
        if( level .lt. min_exp_depth-1 ) goto 1200
        if( box(10) .eq. 0 ) goto 1200
c
        zkbsize = zk*bsize(level)
        allocate(hexp_tmp(-nterms(level):nterms(level)))
        do idim=1,nd
        do j=-nterms(level),nterms(level)
           hexp_tmp(j) = rmlexp(iaddr(1,ibox)+(idim-1)
     $        +(j+nterms(level))*nd)
        enddo
        call h2dmpexp_fast(zkbsize,
     $     scale(level),center0,hexp_tmp,nterms(level),
     $     expio(ie_addr(1,ibox)+(idim-1)*4*(nexp(level)+1)),
     $     nexp(level),
     $     quadr(iq_addr(1,level)),quadr(iq_addr(2,level)),
     $     zcc(level),quadr(iq_addr(3,level)),
     $     quadr(iq_addr(5,level)),quadr(iq_addr(6,level)))
        enddo
        deallocate(hexp_tmp)
c
 1200    continue
C$OMP END PARALLEL DO
        t2a=second()
C$        t2a=omp_get_wtime()
c
c       ... step 4b, exponential-to-exponential translations
c
        if (ifprint .ge. 1)
     $     call prinf('=== STEP 4b (expout to expinc) ===*',i,0)

        do 4300 ilev=max(min_depth+2,min_exp_depth),nlev+1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level,list,nlist,itype)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1)
C$OMP$PRIVATE(ier,ilist,ii,jj,idim)
C$OMP$SCHEDULE(DYNAMIC)
        do 4200 ibox=laddr(1,ilev),laddr(1,ilev)+laddr(2,ilev)-1
        call d2tgetb(ier,ibox,box,center0,corners0,wlists)
        level=box(1)
        if (level .ge. 2) then
           itype=2
           call d2tgetl(ier,ibox,itype,list,nlist,wlists)
           do 4150 ilist=1,nlist
              jbox=list(ilist)
              call d2tgetb(ier,jbox,box1,center1,corners1,wlists)
              if( box1(10) .eq. 0 ) goto 4150
              if (jbox.eq.0) goto 4150
              if ((box(12).eq.0).and.(ifprune_list2.eq.1))
     $           goto 4150
              ii=box(2)-box1(2)
              jj=box(3)-box1(3)
              do idim=1,nd
              if( ii .ge. +2) then
              call h2dexpexp_e(
     $           expio(ie_addr(1,jbox)+(idim-1)*4*(nexp(level)+1)),
     $           expio(ie_addr(2,ibox)+(idim-1)*4*(nexp(level)+1)),
     $           nexp(level),quadr(iq_addr(4,level)),ii,jj)
              else if ( ii .le. -2 ) then
              call h2dexpexp_w(
     $           expio(ie_addr(1,jbox)+(idim-1)*4*(nexp(level)+1)),
     $           expio(ie_addr(2,ibox)+(idim-1)*4*(nexp(level)+1)),
     $           nexp(level),quadr(iq_addr(4,level)),-ii,-jj)
              else if ( jj .ge. +2 ) then
              call h2dexpexp_n(
     $           expio(ie_addr(1,jbox)+(idim-1)*4*(nexp(level)+1)),
     $           expio(ie_addr(2,ibox)+(idim-1)*4*(nexp(level)+1)),
     $           nexp(level),quadr(iq_addr(4,level)),jj,-ii)
              else if ( jj .le. -2 ) then
              call h2dexpexp_s(
     $           expio(ie_addr(1,jbox)+(idim-1)*4*(nexp(level)+1)),
     $           expio(ie_addr(2,ibox)+(idim-1)*4*(nexp(level)+1)),
     $           nexp(level),quadr(iq_addr(4,level)),-jj,ii)
              endif
              enddo
 4150       continue
        endif
 4200   continue
C$OMP END PARALLEL DO
 4300   continue
c
c       ... step 4c, convert incoming exponential expansions to local
c
        if(ifprint .ge. 1)
     $     call prinf('=== STEP 4c (form expinc) ====*',i,0)
        t3a=second()
C$        t3a=omp_get_wtime()
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level,nkids,idim)
C$OMP$PRIVATE(lexp_tmp,zkbsize,ier,j)
C$OMP$SCHEDULE(DYNAMIC)
        do 1400 ibox=1,nboxes
        call d2tgetb(ier,ibox,box,center0,corners0,wlists)
        call d2tnkids(box,nkids)
        level=box(1)
        if( level .lt. min_exp_depth-1 ) goto 1400
c
        zkbsize = zk*bsize(level)
        allocate(lexp_tmp(-nterms(level):nterms(level)))
        do idim=1,nd
        do j=-nterms(level),nterms(level)
           lexp_tmp(j) = 0.0d0
        enddo
        call h2dexpta_fast(zkbsize,
     $     scale(level),center0,lexp_tmp,nterms(level),
     $     expio(ie_addr(2,ibox)+(idim-1)*4*(nexp(level)+1)),
     $     nexp(level),
     $     quadr(iq_addr(1,level)),quadr(iq_addr(2,level)),
     $     quadr(iq_addr(7,level)),quadr(iq_addr(8,level)))
        do j=-nterms(level),nterms(level)
           rmlexp(iaddr(2,ibox)+(idim-1)+(j+nterms(level))*nd) =
     $        rmlexp(iaddr(2,ibox)+(idim-1)+(j+nterms(level))*nd)
     $        + lexp_tmp(j)
        enddo
        enddo
        deallocate(lexp_tmp)
c
 1400    continue
C$OMP END PARALLEL DO
        t4a=second()
C$        t4a=omp_get_wtime()
c
c       ... step 4f_init, initialize far field signature translation ops
c
        if (ifprint .ge. 1)
     $     call prinf('=== STEP 4f_init (far to pwl) ===*',i,0)

        do 8301 ilev=3,min_depth+1
        do 8201 ibox=laddr(1,ilev),laddr(1,ilev)+laddr(2,ilev)-1
        call d2tgetb(ier,ibox,box,center0,corners0,wlists)
        level0=box(1)
        if (level0 .ge. 2) then
           itype=2
           call d2tgetl(ier,ibox,itype,list,nlist,wlists)
           do 8151 ilist=1,nlist
              jbox=list(ilist)
              call d2tgetb(ier,jbox,box1,center1,corners1,wlists)
              if( box1(10) .eq. 0 ) goto 8151
              if (jbox.eq.0) goto 8151
              if ((box(12).eq.0).and.(ifprune_list2.eq.1))
     $           goto 8151
              level1=box1(1)
              ii=box(2)-box1(2)
              jj=box(3)-box1(3)
              if( cache_init(ii,jj,level1) .eq. 0) then
              cache_init(ii,jj,level1)=1
              endif
 8151         continue
        endif
 8201   continue
 8301   continue

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(center0,center1,level1)
C$OMP$SCHEDULE(DYNAMIC) COLLAPSE(3)
        do ilev = 3,min_depth+1
        do ii = -3,3
        do jj = -3,3
        level1 = ilev - 1
        center0(1) = 0
        center0(2) = 0
        center1(1) = -bsize(level1)*ii
        center1(2) = -bsize(level1)*jj
        if( cache_init(ii,jj,level1) .eq. 1 ) then
              call h2dfarpwl_init(zk,center1,center0,
     $           nterms_tra(level1),nterms_far(level1),
     $           trans_cache(it_addr(ii,jj,level1)))
        endif
        enddo
        enddo
        enddo
C$OMP END PARALLEL DO

c       ... step 4f, convert far field signatures to incoming plane ones
c
        if (ifprint .ge. 1)
     $     call prinf('=== STEP 4f (far to pwl) ===*',i,0)

        do 8300 ilev=3,min_depth+1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level0,list,nlist,itype)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1)
C$OMP$PRIVATE(ier,ilist,ii,jj)
C$OMP$SCHEDULE(DYNAMIC)
        do 8200 ibox=laddr(1,ilev),laddr(1,ilev)+laddr(2,ilev)-1
        call d2tgetb(ier,ibox,box,center0,corners0,wlists)
        level0=box(1)
        if (level0 .ge. 2) then
           itype=2
           call d2tgetl(ier,ibox,itype,list,nlist,wlists)
           do 8150 ilist=1,nlist
              jbox=list(ilist)
              call d2tgetb(ier,jbox,box1,center1,corners1,wlists)
              if( box1(10) .eq. 0 ) goto 8150
              if (jbox.eq.0) goto 8150
              if ((box(12).eq.0).and.(ifprune_list2.eq.1))
     $           goto 8150
              level1=box1(1)
              ii=box(2)-box1(2)
              jj=box(3)-box1(3)
              call h2dfarpwl(zk,
     $           center1,rmlexp(iaddr(3,jbox)),nterms_far(level1),
     $           center0,rmlexp(iaddr(4,ibox)),nterms_far(level1),
     $           trans_cache(it_addr(ii,jj,level1)))
 8150       continue
        endif
 8200   continue
C$OMP END PARALLEL DO
 8300   continue
c
        t2=second()
C$        t2=omp_get_wtime()
        timeinfo(4)=t2-t1
c


        if (ifprint .ge. 1)
     $     call prinf('=== STEP 5bf (split pwl) ===*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 5bf, split all incoming plane signatures
c
        do 5320 ilev=3,min_depth
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level0,nkids)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1,ier,i)
        do 5220 ibox=laddr(1,ilev),laddr(1,ilev)+laddr(2,ilev)-1
        call d2tgetb(ier,ibox,box,center0,corners0,wlists)
        call d2tnkids(box,nkids)
        if (nkids .ne. 0) then
            level0=box(1)
            if (level0 .ge. 2) then
               do 5120 i = 1,4
	          jbox = box(4+i)
	          if (jbox.eq.0) goto 5120
                  call d2tgetb(ier,jbox,box1,center1,corners1,wlists)
                  level1=box1(1)
                  call h2dpwlpwl(zk,
     1               center0,rmlexp(iaddr(4,ibox)),nterms_far(level0),
     $               center1,rmlexp(iaddr(4,jbox)),nterms_far(level1))
 5120          continue
            endif
        endif
 5220   continue
C$OMP END PARALLEL DO
 5320   continue
c
c
        if (ifprint .ge. 1)
     $     call prinf('=== STEP 5af (pwl to lo) ===*',i,0)
c
        if( min_depth .ge. 2 ) then
        do ilev=min_depth+1,min_depth+1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level0,ier)
         do ibox=laddr(1,ilev),laddr(1,ilev)+laddr(2,ilev)-1
         call d2tgetb(ier,ibox,box,center0,corners0,wlists)
         level0=box(1)
         if( level0 .ge. 2 ) then
         call h2dpwlta(rmlexp(iaddr(2,ibox)),nterms(level0),
     $      rmlexp(iaddr(4,ibox)),nterms_far(level0))
         endif
         enddo
C$OMP END PARALLEL DO
         enddo
         endif
c
        if (ifprint .ge. 1)
     $     call prinf('=== STEP 5 (split lo) ===*',i,0)
c
c       ... step 5, split all local expansions
c
        do 5300 ilev=max(3,min_depth+1),nlev
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level0,level,npts,nkids,radius)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1)
C$OMP$PRIVATE(mptemp,lused,ier,i,j,ptemp,ftemp,cd)
        do 5200 ibox=laddr(1,ilev),laddr(1,ilev)+laddr(2,ilev)-1
c
        call d2tgetb(ier,ibox,box,center0,corners0,wlists)
        call d2tnkids(box,nkids)
c
        if (nkids .ne. 0) then
            level0=box(1)
            if (level0 .ge. 2) then
               do 5100 i = 1,4
	          jbox = box(4+i)
	          if (jbox.eq.0) goto 5100
                  call d2tgetb(ier,jbox,box1,center1,corners1,wlists)
                  level1=box1(1)
                  call h2dlocloc(nd,zk,scale(level0),center0,
     1               rmlexp(iaddr(2,ibox)),nterms(level0),
     1               scale(level1),center1,mptemp,nterms(level1))
                  call h2dadd(nd,mptemp,rmlexp(iaddr(2,jbox)),
     1   	       nterms(level1))
 5100          continue
            endif
        endif
 5200   continue
C$OMP END PARALLEL DO
 5300   continue
c
        t2=second()
C$        t2=omp_get_wtime()
        timeinfo(5)=t2-t1
c
        return
        end
c
c
c
c
c
