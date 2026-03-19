cc Copyright (C) 2010-2011: Leslie Greengard and Zydrunas Gimbutas
cc SPDX-License-Identifier: GPL-2.0-or-later
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        this is the end of the debugging code and the beginning
c        of the routines for Helmholtz FMM in R^2
c
c        FAR-FIELD VERSION with nd-density support
c
c        This file contains hfmm2d_list2_far which performs Steps 3-5
c        of the FMM using far-field signature functions for M2L
c        translations at coarse levels.
c        Supports vectorized evaluations with nd independent densities.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine hfmm2d_list2_far(nd,
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
        complex *16 zk
c
        complex *16 ptemp,ftemp(2),htemp(3)
c
        real *8 timeinfo(10)
c
        integer cache_init(-3:3,-3:3,0:200)
        integer, allocatable :: it_addr(:,:,:)
        complex *16, allocatable :: trans_cache(:)
c
c       ... scratch arrays for far-field nd loops
c
        complex *16, allocatable :: fartemp1(:)
        complex *16, allocatable :: fartemp2(:)
        complex *16, allocatable :: mpscr(:)
c
c
c       ... initialize far field translation cache
c
        allocate(it_addr(-3:3,-3:3,0:200))
c
        do ilev=0,200
c
c       ... allocate memory for outgoing exponentials
c
        do j = -3,3
        do i = -3,3
c
        cache_init(i,j,ilev) = 0
c
        enddo
        enddo

        enddo

c
c       ... allocate far field translation cache
c
        iptr=1
        do ilev=2,min_depth
c
c       ... allocate memory for outgoing exponentials
c
        do j = -3,3
        do i = -3,3

        it_addr(i,j,ilev) = iptr
        iptr = iptr + nterms_far(ilev)
c
        enddo
        enddo

        enddo
c
ccc        call prinf('iptr=*',iptr,1)
        allocate( trans_cache(iptr) )

c
c       ... find max nterms and max nterms_far for scratch allocation
c
        nmax = 0
        nfar_max = 0
        do ilev = 0, nlev
        if( nterms(ilev) .gt. nmax ) nmax = nterms(ilev)
        if( nterms_far(ilev) .gt. nfar_max )
     $     nfar_max = nterms_far(ilev)
        enddo
c
c       ... allocate scratch arrays for far-field nd loops
c
        allocate( fartemp1(nfar_max) )
        allocate( fartemp2(nfar_max) )
        allocate( mpscr(2*nmax+1) )
c

c
c     ifprint is an internal information printing flag.
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c     Prints timing breakdown, list information, and other things if ifprint=2.
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
ccc         do 2200 ibox=nboxes,1,-1
         do 2300 ilev=nlev,max(3,min_depth+1),-1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level0,level,npts,nkids,radius)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1)
C$OMP$PRIVATE(mptemp,lused,ier,i,j,ptemp,ftemp,cd)
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(4)
         do 2200 ibox=laddr(1,ilev),laddr(1,ilev)+laddr(2,ilev)-1
c
         call d2tgetb(ier,ibox,box,center0,corners0,wlists)
         call d2tnkids(box,nkids)
c
c       ... prune all sourceless boxes
c
         if( box(10) .eq. 0 ) goto 2200
c
         if (nkids .ne. 0) then
c
         level0=box(1)
         if( level0 .ge. 2 ) then
ccc         if (level0 .ge. 0) then
            radius = (corners0(1,1) - center0(1))**2
            radius = radius + (corners0(2,1) - center0(2))**2
            radius = sqrt(radius)
c
            if( ifprint .ge. 2 ) then
               call prin2('radius=*',radius,1)
               call prinf('ibox=*',ibox,1)
               call prinf('box=*',box,15)
               call prinf('nkids=*',nkids,1)
            endif
c
c       ... merge multipole expansions of the kids
c
            call h2dzero(nd,rmlexp(iaddr(1,ibox)),nterms(level0))
            if (ifprint .ge. 2) then
               call prin2('center0=*',center0,2)
            endif
c
            do 2100 i = 1,4
               jbox = box(4+i)
               if (jbox.eq.0) goto 2100
               call d2tgetb(ier,jbox,box1,center1,corners1,wlists)
               if (ifprint .ge. 2) then
               call prinf('jbox=*',jbox,1)
               call prin2('center1=*',center1,2)
               endif
               level1=box1(1)
               call h2dmpmp(nd,zk,scale(level1),center1,
     1            rmlexp(iaddr(1,jbox)),nterms(level1),scale(level0),
     1            center0,mptemp,nterms(level0))
               call h2dadd(nd,mptemp,rmlexp(iaddr(1,ibox)),
     1            nterms(level0))
 2100       continue
            if (ifprint .ge. 2) then
            call prinf('=============*',x,0)
            endif
c       ... mark the local expansion of all kids and the parent
c
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
C$OMP$PRIVATE(ibox,box,center0,corners0,level0,level,npts,nkids,radius)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1)
C$OMP$PRIVATE(mpscr,fartemp1,lused,ier,i,j,k,idim,ptemp,ftemp,cd)
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(4)
         do ibox=laddr(1,ilev),laddr(1,ilev)+laddr(2,ilev)-1
c
         call d2tgetb(ier,ibox,box,center0,corners0,wlists)
         call d2tnkids(box,nkids)
c
c       ... prune all sourceless boxes
c
         if( box(10) .gt. 0 ) then
c
         level0=box(1)
         if( level0 .ge. 2 ) then
         call h2dzero_far_nd(nd,
     $      rmlexp(iaddr(3,ibox)),nterms_far(level0))
c
         do idim=1,nd
c
c       ... extract density idim's multipole to contiguous scratch
c
         do j=-nterms(level0),nterms(level0)
         mpscr(j+nterms(level0)+1) =
     $      rmlexp(iaddr(1,ibox)+(idim-1)+(j+nterms(level0))*nd)
         enddo
c
         call h2dmpfar(mpscr,nterms(level0),
     $      fartemp1,nterms_far(level0))
c
c       ... scatter far-field signature back (nd-interleaved)
c
         do k=1,nterms_far(level0)
         rmlexp(iaddr(3,ibox)+(idim-1)+(k-1)*nd) = fartemp1(k)
         enddo
c
         enddo
c
         endif
c
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
c       ... step 3bf, merge all far field signatures
c
         do 2320 ilev=min_depth,3,-1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level0,level,npts,nkids,radius)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1)
C$OMP$PRIVATE(fartemp1,fartemp2)
C$OMP$PRIVATE(lused,ier,i,j,k,idim,ptemp,ftemp,cd)
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(4)
         do 2220 ibox=laddr(1,ilev),laddr(1,ilev)+laddr(2,ilev)-1
c
         call d2tgetb(ier,ibox,box,center0,corners0,wlists)
         call d2tnkids(box,nkids)
c
c       ... prune all sourceless boxes
c
         if( box(10) .eq. 0 ) goto 2220
c
         if (nkids .ne. 0) then
c
         level0=box(1)
         if( level0 .ge. 2 ) then
c
c       ... merge far field signatures of the kids
c
            call h2dzero_far_nd(nd,
     $         rmlexp(iaddr(3,ibox)),nterms_far(level0))
c
            do 2120 i = 1,4
               jbox = box(4+i)
               if (jbox.eq.0) goto 2120
               call d2tgetb(ier,jbox,box1,center1,corners1,wlists)
               if (ifprint .ge. 2) then
               call prinf('jbox=*',jbox,1)
               call prin2('center1=*',center1,2)
               endif
               level1=box1(1)
c
               do idim=1,nd
c
c       ... extract density idim's far-field from child jbox
c
               do k=1,nterms_far(level1)
               fartemp1(k) =
     $            rmlexp(iaddr(3,jbox)+(idim-1)+(k-1)*nd)
               enddo
c
c       ... extract density idim's far-field from parent ibox
c
               do k=1,nterms_far(level0)
               fartemp2(k) =
     $            rmlexp(iaddr(3,ibox)+(idim-1)+(k-1)*nd)
               enddo
c
               call h2dfarfar(zk,
     1            center1,fartemp1,nterms_far(level1),
     1            center0,fartemp2,nterms_far(level0))
c
c       ... scatter back to parent ibox
c
               do k=1,nterms_far(level0)
               rmlexp(iaddr(3,ibox)+(idim-1)+(k-1)*nd) =
     $            fartemp2(k)
               enddo
c
               enddo
c
 2120       continue

            endif
         endif
 2220    continue
C$OMP END PARALLEL DO
 2320    continue
c
ccc        call prinf('=== UPWARD PASS COMPLETE ===*',i,0)
c
c------------------------------------------------------------
c      DEBUGGING SEGMENT - once all multipole expansions are merged
c      to top level, one can compare it to the direct formation of the
c      expansion at the top level from the source locations.
c
ccc        call prinm(rmlexp(iaddr(1,1)),nterms(0))
c
ccc        call h2dformmp(ier,scale(0),source,charge,n,
ccc     1  	center,nterms(0),mptemp)
c
ccc        call prinm(mptemp,nterms(0))
c
ccc        call h2dmperr(rmlexp(iaddr(1,1)),mptemp,nterms(0),d)
ccc        call prin2('error in upward pass=*',d,1)
c
ccc        pause
ccc        stop
c      END DEBUGGING SEGMENT
c------------------------------------------------------------
c
         t2=second()
C$        t2=omp_get_wtime()
ccc        call prin2('time=*',t2-t1,1)
         timeinfo(3)=t2-t1
c
c
c
        if (ifprint .ge. 1)
     $     call prinf('=== STEP 4 (mp to lo) ===*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 4, convert multipole expansions into the local ones
c
        do 4300 ilev=min_depth+2,nlev+1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level0,list,nlists,nlist,itype)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1,ifdirect2,radius)
C$OMP$PRIVATE(mptemp,lused,ier,i,j,ptemp,ftemp,htemp,cd,ilist)
C$OMP$PRIVATE(if_use_trunc,nterms_trunc,ii,jj)
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1)
        do 4200 ibox=laddr(1,ilev),laddr(1,ilev)+laddr(2,ilev)-1
        call d2tgetb(ier,ibox,box,center0,corners0,wlists)
        if (ifprint .ge. 2) then
           call prinf('ibox=*',ibox,1)
           call prinf('box=*',box,15)
        endif
        level0=box(1)
        if (level0 .ge. 2) then
c
c       ... retrieve list #2
c
           itype=2
           call d2tgetl(ier,ibox,itype,list,nlist,wlists)
           if (ifprint .ge. 2) then
              call prinf('list2=*',list,nlist)
           endif
c
c       ... prune all sourceless boxes
c
ccc           if( box(10) .eq. 0 ) nlist=0

c       ... for all pairs in list #2, apply the translation operator
c
           do 4150 ilist=1,nlist
              jbox=list(ilist)
              call d2tgetb(ier,jbox,box1,center1,corners1,wlists)
              if( box1(10) .eq. 0 ) goto 4150
              if (jbox.eq.0) goto 4150
              if ((box(12).eq.0).and.(ifprune_list2.eq.1))
     $           goto 4150
c
c       ... convert multipole expansions for all boxes in list 2 in local exp
c       ... if source is childless, evaluate directly (if cheaper)
c

              level1=box1(1)
              ifdirect2 = 0
c
              if_use_trunc = 0
c
              if (ifdirect2 .eq. 0) then
              if( if_use_trunc .eq. 0) then
                 call h2dzero(nd,mptemp,nterms(level1))
                 call h2dmploc(nd,zk,scale(level1),center1,
     1              rmlexp(iaddr(1,jbox)),nterms(level1),
     $              scale(level0),
     1              center0,mptemp,nterms(level0))
                 call h2dadd(nd,mptemp,rmlexp(iaddr(2,ibox)),
     1              nterms(level0))
c
              else

              ii=box1(2)-box(2)
              jj=box1(3)-box(3)
              nterms_trunc=itable(ii,jj)
              nterms_trunc=min(nterms(level0),nterms_trunc)
              nterms_trunc=min(nterms(level1),nterms_trunc)

              call h2dmploc_add_trunc(nd,zk,scale(level1),center1,
     $           rmlexp(iaddr(1,jbox)),nterms(level1),nterms_trunc,
     $           scale(level0),center0,rmlexp(iaddr(2,ibox)),
     $           nterms(level0),nterms_trunc)

              endif
              endif

 4150       continue
        endif
 4200   continue
C$OMP END PARALLEL DO
 4300   continue
c

        if (ifprint .ge. 1)
     $     call prinf('=== STEP 4f_init (far to pwl) ===*',i,0)
c
c       ... step 4f_init, initialize far field signatures translation ops
c
cc        call prinf('laddr=*',laddr,2*(nlev+1))
cc        call prin2('bsize=*',bsize,(nlev+1))
        do 8301 ilev=3,min_depth+1
        do 8201 ibox=laddr(1,ilev),laddr(1,ilev)+laddr(2,ilev)-1
        call d2tgetb(ier,ibox,box,center0,corners0,wlists)
        if (ifprint .ge. 2) then
           call prinf('ibox=*',ibox,1)
           call prinf('box=*',box,15)
        endif
        level0=box(1)
        if (level0 .ge. 2) then
c
c       ... retrieve list #2
c
           itype=2
           call d2tgetl(ier,ibox,itype,list,nlist,wlists)
           if (ifprint .ge. 2) then
              call prinf('list2=*',list,nlist)
           endif
c
c       ... prune all sourceless boxes
c
ccc           if( box(10) .eq. 0 ) nlist=0

c       ... for all pairs in list #2, apply the translation operator
c
           do 8151 ilist=1,nlist
              jbox=list(ilist)
              call d2tgetb(ier,jbox,box1,center1,corners1,wlists)
              if( box1(10) .eq. 0 ) goto 8151
              if (jbox.eq.0) goto 8151
              if ((box(12).eq.0).and.(ifprune_list2.eq.1))
     $           goto 8151

              level1=box1(1)
c
              ii=box(2)-box1(2)
              jj=box(3)-box1(3)

              if( cache_init(ii,jj,level1) .eq. 0) then
c              call h2dfarpwl_init(zk,center1,center0,
c     $           nterms_tra(level1),nterms_far(level1),
c     $           trans_cache(it_addr(ii,jj,level1)))
              cache_init(ii,jj,level1)=1
              endif

 8151         continue
        endif
 8201   continue
 8301   continue
c


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


        if (ifprint .ge. 1)
     $     call prinf('=== STEP 4f (far to pwl) ===*',i,0)
c
c       ... step 4f, convert far field signatures to incoming plane ones
c
        do 8300 ilev=3,min_depth+1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level0,list,nlists,nlist,itype)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1,ifdirect2,radius)
C$OMP$PRIVATE(fartemp1,fartemp2)
C$OMP$PRIVATE(lused,ier,i,j,k,idim,ptemp,ftemp,htemp,cd,ilist)
C$OMP$PRIVATE(if_use_trunc,nterms_trunc,ii,jj)
C$OMP$SCHEDULE(DYNAMIC)
ccC$OMP$NUM_THREADS(1)
        do 8200 ibox=laddr(1,ilev),laddr(1,ilev)+laddr(2,ilev)-1
        call d2tgetb(ier,ibox,box,center0,corners0,wlists)
        if (ifprint .ge. 2) then
           call prinf('ibox=*',ibox,1)
           call prinf('box=*',box,15)
        endif
        level0=box(1)
        if (level0 .ge. 2) then
c
c       ... retrieve list #2
c
           itype=2
           call d2tgetl(ier,ibox,itype,list,nlist,wlists)
           if (ifprint .ge. 2) then
              call prinf('list2=*',list,nlist)
           endif
c
c       ... prune all sourceless boxes
c
ccc           if( box(10) .eq. 0 ) nlist=0

c       ... for all pairs in list #2, apply the translation operator
c
           do 8150 ilist=1,nlist
              jbox=list(ilist)
              call d2tgetb(ier,jbox,box1,center1,corners1,wlists)
              if( box1(10) .eq. 0 ) goto 8150
              if (jbox.eq.0) goto 8150
              if ((box(12).eq.0).and.(ifprune_list2.eq.1))
     $           goto 8150

              level1=box1(1)
c
              ii=box(2)-box1(2)
              jj=box(3)-box1(3)
c
              do idim=1,nd
c
c       ... extract density idim's far-field from source box jbox
c
              do k=1,nterms_far(level1)
              fartemp1(k) =
     $           rmlexp(iaddr(3,jbox)+(idim-1)+(k-1)*nd)
              enddo
c
c       ... extract density idim's plane-wave from target box ibox
c
              do k=1,nterms_far(level1)
              fartemp2(k) =
     $           rmlexp(iaddr(4,ibox)+(idim-1)+(k-1)*nd)
              enddo
c
              call h2dfarpwl(zk,
     $           center1,fartemp1,nterms_far(level1),
     $           center0,fartemp2,nterms_far(level1),
     $           trans_cache(it_addr(ii,jj,level1)))
c
c       ... scatter plane-wave back to target box ibox
c
              do k=1,nterms_far(level1)
              rmlexp(iaddr(4,ibox)+(idim-1)+(k-1)*nd) =
     $           fartemp2(k)
              enddo
c
              enddo

 8150       continue
        endif
 8200   continue
C$OMP END PARALLEL DO
c        t4=second()
cC$        t4=omp_get_wtime()
c        write(*,*) 'level ', ilev, ' time in list2:', t4-t3
ccc        write(*,*) 'time in list2:', second()-t1
ccc        write(*,*) 'ntops:', ntops
ccc        write(*,*) 'speed:', ntops/(second()-t1)
 8300   continue
c
        t2=second()
C$        t2=omp_get_wtime()
ccc     call prin2('time=*',t2-t1,1)
        timeinfo(4)=t2-t1
c



        if (ifprint .ge. 1)
     $     call prinf('=== STEP 5bf (split pwl) ===*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 5, split all incoming plane signatures
c
        do 5320 ilev=3,min_depth
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level0,level,npts,nkids,radius)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1)
C$OMP$PRIVATE(fartemp1,fartemp2)
C$OMP$PRIVATE(lused,ier,i,j,k,idim,ptemp,ftemp,cd)
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(4)
        do 5220 ibox=laddr(1,ilev),laddr(1,ilev)+laddr(2,ilev)-1
c
        call d2tgetb(ier,ibox,box,center0,corners0,wlists)
        call d2tnkids(box,nkids)
c
        if (nkids .ne. 0) then
            level0=box(1)
            if (level0 .ge. 2) then
c
c       ... split plane-wave signature of the parent box
c
               do 5120 i = 1,4
	          jbox = box(4+i)
	          if (jbox.eq.0) goto 5120
                  call d2tgetb(ier,jbox,box1,center1,corners1,wlists)
                  level1=box1(1)
c
                  do idim=1,nd
c
c       ... extract density idim's plane-wave from parent ibox
c
                  do k=1,nterms_far(level0)
                  fartemp1(k) =
     $               rmlexp(iaddr(4,ibox)+(idim-1)+(k-1)*nd)
                  enddo
c
c       ... extract density idim's plane-wave from child jbox
c
                  do k=1,nterms_far(level1)
                  fartemp2(k) =
     $               rmlexp(iaddr(4,jbox)+(idim-1)+(k-1)*nd)
                  enddo
c
                  call h2dpwlpwl(zk,
     1               center0,fartemp1,nterms_far(level0),
     $               center1,fartemp2,nterms_far(level1))
c
c       ... scatter plane-wave back to child jbox
c
                  do k=1,nterms_far(level1)
                  rmlexp(iaddr(4,jbox)+(idim-1)+(k-1)*nd) =
     $               fartemp2(k)
                  enddo
c
                  enddo
c
 5120          continue
            endif
        endif
c
 5220   continue
C$OMP END PARALLEL DO
 5320   continue
c
c
c
        if (ifprint .ge. 1)
     $     call prinf('=== STEP 5af (pwl to lo) ===*',i,0)
c
c       ... step 5f, convert incoming plane signatures
c
        if( min_depth .ge. 2 ) then
        do ilev=min_depth+1,min_depth+1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level0,level,npts,nkids,radius)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1)
C$OMP$PRIVATE(mpscr,fartemp1)
C$OMP$PRIVATE(lused,ier,i,j,k,idim,ptemp,ftemp,cd)
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(4)
         do ibox=laddr(1,ilev),laddr(1,ilev)+laddr(2,ilev)-1
c
         call d2tgetb(ier,ibox,box,center0,corners0,wlists)
         call d2tnkids(box,nkids)
c
         level0=box(1)
         if( level0 .ge. 2 ) then
c
         do idim=1,nd
c
c       ... extract density idim's local expansion to contiguous scratch
c
         do j=-nterms(level0),nterms(level0)
         mpscr(j+nterms(level0)+1) =
     $      rmlexp(iaddr(2,ibox)+(idim-1)+(j+nterms(level0))*nd)
         enddo
c
c       ... extract density idim's plane-wave to contiguous scratch
c
         do k=1,nterms_far(level0)
         fartemp1(k) =
     $      rmlexp(iaddr(4,ibox)+(idim-1)+(k-1)*nd)
         enddo
c
         call h2dpwlta(mpscr,nterms(level0),
     $      fartemp1,nterms_far(level0))
c
c       ... scatter local expansion back (nd-interleaved)
c
         do j=-nterms(level0),nterms(level0)
         rmlexp(iaddr(2,ibox)+(idim-1)+(j+nterms(level0))*nd) =
     $      mpscr(j+nterms(level0)+1)
         enddo
c
         enddo
c
c         call prin2('rmlexp(iaddr(4,ibox))=*',rmlexp(iaddr(4,ibox)),
c     $      2*nterms_far(level0))
c         call prinf('nterms_far(level0)=*',nterms_far(level0),1)
         endif

         enddo
C$OMP END PARALLEL DO
         enddo
         endif

        if (ifprint .ge. 1)
     $     call prinf('=== STEP 5 (split lo) ===*',i,0)
c
c       ... step 5, split all local expansions
c
ccc        do 5200 ibox=1,nboxes
        do 5300 ilev=max(3,min_depth+1),nlev
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level0,level,npts,nkids,radius)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1)
C$OMP$PRIVATE(mptemp,lused,ier,i,j,ptemp,ftemp,cd)
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(4)
        do 5200 ibox=laddr(1,ilev),laddr(1,ilev)+laddr(2,ilev)-1
c
        call d2tgetb(ier,ibox,box,center0,corners0,wlists)
        call d2tnkids(box,nkids)
c
        if (nkids .ne. 0) then
            level0=box(1)
            if (level0 .ge. 2) then
               if (ifprint .ge. 2) then
                  call prinf('ibox=*',ibox,1)
                  call prinf('box=*',box,15)
                  call prinf('nkids=*',nkids,1)
                  call prin2('center0=*',center0,2)
               endif
c
c       ... split local expansion of the parent box
c
               do 5100 i = 1,4
	          jbox = box(4+i)
	          if (jbox.eq.0) goto 5100
                  call d2tgetb(ier,jbox,box1,center1,corners1,wlists)
                  radius = (corners1(1,1) - center1(1))**2
                  radius = radius + (corners1(2,1) - center1(2))**2
                  radius = sqrt(radius)
                  if (ifprint .ge. 2) then
                     call prinf('jbox=*',jbox,1)
                     call prin2('radius=*',radius,1)
                     call prin2('center1=*',center1,2)
                  endif
                  level1=box1(1)
                  call h2dlocloc(nd,zk,scale(level0),center0,
     1               rmlexp(iaddr(2,ibox)),nterms(level0),
     1               scale(level1),center1,mptemp,nterms(level1))
                  call h2dadd(nd,mptemp,rmlexp(iaddr(2,jbox)),
     1   	       nterms(level1))
 5100          continue
               if (ifprint .ge. 2) call prinf('=============*',x,0)
            endif
        endif
c
        if (nkids .ne. 0) then
            level=box(1)
            if (level .ge. 2) then
               if( ifprint .ge. 2 ) then
                  call prinf('ibox=*',ibox,1)
                  call prinf('box=*',box,15)
                  call prinf('nkids=*',nkids,1)
               endif
            endif
        endif
 5200   continue
C$OMP END PARALLEL DO
 5300   continue
c
        t2=second()
C$        t2=omp_get_wtime()
ccc     call prin2('time=*',t2-t1,1)
        timeinfo(5)=t2-t1
c
        return
        end
c
c
c
c
c
        subroutine h2dmpalloc_far_nd(nd,wlists,iaddr,nboxes,lmptot,
     $     nterms,nterms_far)
        implicit real *8 (a-h,o-z)
        integer box(20)
        integer nterms(0:*),nterms_far(0:*)
        integer iaddr(4,nboxes)
        real *8 center0(2),corners0(2,4)
        real *8 wlists(*)
c
c       ... construct pointer array iaddr for addressing multipole and
c       local expansion
c
        iptr=1
        do ibox=1,nboxes
        call d2tgetb(ier,ibox,box,center0,corners0,wlists)
        level=box(1)
c
c       ... first, allocate memory for the multipole expansion
c
        iaddr(1,ibox)=iptr
        iptr=iptr+(2*nterms(level)+1)*nd
c
c       ... then, allocate memory for the local expansion
c
        iaddr(2,ibox)=iptr
        iptr=iptr+(2*nterms(level)+1)*nd
c
c       ... allocate memory for the far field signature
c
        iaddr(3,ibox)=iptr
        iptr=iptr+(nterms_far(level))*nd
c
c       ... allocate memory for the incoming plane signature
c
        iaddr(4,ibox)=iptr
        iptr=iptr+(nterms_far(level))*nd
c
        enddo
        lmptot = iptr
        return
        end
c
c
c
c
c
        subroutine h2dzero_far(mpfar,nfar)
        implicit real *8 (a-h,o-z)
c
c       ... set far field signature to zero
c
        complex *16 mpfar(nfar)
c
        do n=1,nfar
        mpfar(n)=0
        enddo
c
        return
        end
c
c
c
c
c
        subroutine h2dzero_far_nd(nd,mpfar,nfar)
        implicit real *8 (a-h,o-z)
c
c       ... set nd-interleaved far field signature to zero
c
        complex *16 mpfar(nd,nfar)
c
        do n=1,nfar
        do idim=1,nd
        mpfar(idim,n)=0
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
        subroutine hfmm2dparttree_far(ier,iprec,zk,minlevel,
     $     nsource,source,ntarget,target,
     $     nbox,epsfmm,iisource,iitarget,iwlists,lwlists,
     $     nboxes,laddr,nlev,center,size,
     $     w,lw,lused7)
        implicit real *8 (a-h,o-z)
c
c       Helmholtz FMM in R^2: build the quad-tree
c
c     INPUT PARAMETERS:
c
c       zk: complex *16: Helmholtz parameter
c       nsource: integer:  number of sources
c       source: real *8 (2,n):  source locations
c       w: real *8 (lw): workspace
c       lw:  length of workspace
c
c     OUTPUT PARAMETERS:
c
c       ier   =  error return code
c       lused7 = the amount of workspace used
c
c
        real *8 source(2,*),target(2,*)
c
        real *8 center(2)
c
        integer laddr(2,200)
        integer box(15)
        real *8 center0(2),corners0(2,4)
c
        integer box1(15)
        real *8 center1(2),corners1(2,4)
c
        complex *16 zk
        real *8 w(*)
c
        ier=0
c
        done=1
        pi=4*atan(done)
c
        lused7=0
        ifprint=0
c
        iisource=1
        lused7=lused7+nsource
        if (ifprint.eq.1) call prinf('lused7=*',lused7,1)
        if (lused7 .ge. lw) ier=128
        if( ier .ne. 0 ) return
c
        iitarget=iisource+nsource
        lused7=lused7+ntarget
        if (ifprint.eq.1) call prinf('lused7=*',lused7,1)
        if (lused7 .ge. lw) ier=128
        if( ier .ne. 0 ) return
c
        iwlists=iisource+lused7+10
c
c       ... construct the adaptive FMM quad-tree structure
c
        call d2tstrcr(ier,source,nsource,nbox,
     $     nboxes,w(iisource),laddr,nlev,center,size,
     $     target,ntarget,w(iitarget),w(iwlists),lw-lused7,lused)
c
        if( ier .ne. 0 ) return
c
        minlevel = min(minlevel,nlev)
        ifempty=0
ccc        minlevel=0
        maxlevel=30
        call d2tstrcrem(ier,source,nsource,nbox,
     $     nboxes,w(iisource),laddr,nlev,center,size,
     $     target,ntarget,w(iitarget),w(iwlists),lw-lused7,lused,
     $     ifempty,minlevel,maxlevel)
c
        if( ier .ne. 0 ) return
c
        lwlists=lused
        lused7=lused7+lwlists
        if (lused7 .ge. lw) ier=128
        if( ier .ne. 0 ) return
c
        if (ifprint.eq.1)
     $       call prin2('after d2tstrcr, center=*',center,2)
        if (ifprint.eq.1)
     $       call prin2('after d2tstrcr, size=*',size,1)
        if (ifprint.eq.1)
     $       call prinf('after d2tstrcr, nlev=*',nlev,1)
        if (ifprint.eq.1)
     $       call prinf('after d2tstrcr, nbox=*',nbox,1)
        if (ifprint.eq.1)
     $       call prinf('after d2tstrcr, laddr=*',laddr,2*(nlev+1))
c
ccc        call prinf('after d2tstrcr, isource=*',w(iisource),nsource)
ccc        call prinf('after d2tstrcr, itarget=*',w(iitarget),ntarget)
c
ccc        call d2tprint(w(iwlists),lwlists)
c
c       ... optional, plot the oct-tree in gnuplot compatible format
c
        ifplot = 0
        if (ifplot .eq. 1 .and. nsource .lt. 10000 ) then
c
c       ... plot the boxes
c
        iw=51
        call plot_box2d(iw,center,size)
c
        itag=1
        iw=52
        call plot_label2d(iw,center,size,itag,itag)
c
        iw=60
        call plot_points2d(iw,source,nsource)
c
        iw=63
        call plot_points2d(iw,target,ntarget)
c
        do ibox=1,nboxes
           call d2tgetb(ier,ibox,box,center0,corners0,w(iwlists))
           level=box(1)
           size0=size/2**level
c
           iw=61
           call plot_box2d(iw,center0,size0)
c
           itag=ibox
           iw=62
           call plot_label2d(iw,center0,size0,itag,itag)
        enddo
c
        endif
c
        return
        end
c
c
c
c
c
