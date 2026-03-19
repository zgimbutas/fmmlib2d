ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        this is the end of the debugging code and the beginning
c        of the routines for Laplace FMM in R^2
c
c        EXPONENTIAL VERSION with nd-density support
c
c        This file contains lfmm2d_list2_exp which performs Steps 3-5
c        of the FMM using plane wave expansions for M2L translations.
c        Supports vectorized evaluations with nd independent densities.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine lfmm2d_list2_exp(nd,
     $     bsize,nlev,laddr,scale,nterms,rmlexp,iaddr,epsfmm,
     $     timeinfo,wlists,mptemp,lmptemp,
     $     ifprune_list2,nboxes,iprec)
        implicit real *8 (a-h,o-z)
c
        integer iaddr(2,*),laddr(2,*),nterms(0:*)
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
        complex *16 ptemp,ftemp(2),htemp(3)
c
        real *8 timeinfo(10)
c
        real *8, allocatable :: carray(:,:)

        real *8, allocatable :: ws(:), xs(:), comp(:,:)
        complex *16, allocatable :: cshift(:,:,:)
        complex *16, allocatable :: expout(:,:,:,:), expinc(:,:,:,:)
        complex *16, allocatable :: hexp_tmp(:), lexp_tmp(:)

         t1=second()
C$        t1=omp_get_wtime()

c       
c       ... initialize mpmp, mpta, tata operators
c
        ldc = 100
        allocate( carray(0:ldc,0:ldc) )
        call l2d_init_carray(carray,ldc)

c
c       ... initialize exponential translations
c
        allocate(xs(200))
        allocate(ws(200))
        call lwtsexp2(iprec,xs,ws,nexp)
        call prinf('after lwtsexp2, nexp=*',nexp,1)
ccc        nexp = (nexp/4+1)*4
ccc        call lwtsexp2b(nexp,xs,ws,err)
ccc        call prinf('after lwtsexp2, nexp=*',nexp,1)
ccc        call prin2('after lwtsexp2, xs=*',xs,nexp)
ccc        call prin2('after lwtsexp2, ws=*',ws,nexp)

        nterms0 = nterms(1)
        allocate(comp(0:nterms0,nexp))
        allocate(cshift(nexp,-3:3,-3:3))
        call l2expini(nterms0,xs,ws,nexp,xsum,comp,cshift)

c
c       ... allocate and set all outgoing and incoming exponential
c       expansions to zero
c       ... NOTE: added nd dimension for multi-density support
c
        allocate(expout(0:nexp,1:4,1:nd,1:nboxes))
        allocate(expinc(0:nexp,1:4,1:nd,1:nboxes))

        call zero_fill_cmpl(expout,nboxes*nd*4*(nexp+1))
        call zero_fill_cmpl(expinc,nboxes*nd*4*(nexp+1))

cC$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k)
c        do i = 1,nboxes
c        do k = 1,4
c        do j = 0,nexp
c        expout(j,k,i) = 0
c        enddo
c        enddo
c        enddo
cC$OMP END PARALLEL DO
c
cC$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k)        
c        do i = 1,nboxes
c        do k = 1,4
c        do j = 0,nexp
c        expinc(j,k,i) = 0
c        enddo
c        enddo
c        enddo
cC$OMP END PARALLEL DO
        
         t2=second()
C$        t2=omp_get_wtime()

         call prin2('time in list2_exp_ini=*',t2-t1,1)        

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
         do 2300 ilev=nlev,3,-1
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
            call l2dzero(nd,rmlexp(iaddr(1,ibox)),nterms(level0))
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
               if( nterms(level0)+nterms(level1) .gt. 95 ) then
               call l2dmpmp(nd,scale(level1),center1,
     1            rmlexp(iaddr(1,jbox)),nterms(level1),scale(level0),
     1            center0,mptemp,nterms(level0))
               else
               call l2dmpmp_carray(nd,scale(level1),center1,
     1            rmlexp(iaddr(1,jbox)),nterms(level1),scale(level0),
     1            center0,mptemp,nterms(level0),carray,ldc)
               endif
               call l2dadd(nd,mptemp,rmlexp(iaddr(1,ibox)),
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
        if (ifprint .ge. 1) 
     $     call prinf('=== STEP 4 (mp to lo) ===*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 4, convert multipole expansions into the local ones
c

        if(ifprint .ge. 1) 
     $     call prinf('=== STEP 4a (form expout) ====*',i,0)
        t1a=second()
C$        t1a=omp_get_wtime()
c
c       ... step 4a, form outgoing exponential expansions
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level,npts,i,j,idim)
C$OMP$PRIVATE(hexp_tmp)
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1)
        do 1200 ibox=1,nboxes
c
        call d2tgetb(ier,ibox,box,center0,corners0,wlists)
        call d2tnkids(box,nkids)
c
        level=box(1)
        if( level .lt. 2 ) goto 1200
        npts=box(10)
c
c       ... prune all sourceless boxes
c
        if( box(10) .eq. 0 ) goto 1200
c
c       ... loop over densities and convert each multipole to exponential
c       ... NOTE: now expout has nd dimension, so each density is separate
c
        allocate(hexp_tmp(0:nterms(level)))
        do idim=1,nd
c       ... copy density idim's multipole expansion to contiguous array
        do j=0,nterms(level)
           hexp_tmp(j) = rmlexp(iaddr(1,ibox)+(idim-1)+j*nd)
        enddo
        call l2dmpexp(
     $     scale(level),center0,hexp_tmp,nterms(level),
     $     expout(0,1,idim,ibox),nexp,xs,ws,xsum,comp)
        enddo
        deallocate(hexp_tmp)
c
 1200    continue
C$OMP END PARALLEL DO
c
        if (ifprint .ge. 1) 
     $     call prinf('=== STEP 4b (expout to expinc) ===*',i,0)
        t2a=second()
C$        t2a=omp_get_wtime()

        do 4300 ilev=3,nlev+1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level0,list,nlists,nlist,itype)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1,ifdirect2,radius)
C$OMP$PRIVATE(mptemp,lused,ier,i,j,ptemp,ftemp,htemp,cd,ilist)
C$OMP$PRIVATE(if_use_trunc,nterms_trunc,ii,jj,idim) 
cccC$OMP$SCHEDULE(DYNAMIC)
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
c       ... convert outgoing exp expansions for all boxes in list 2
c       in local exp
c
              ii=box(2)-box1(2)
              jj=box(3)-box1(3)
c
c       ... loop over densities for E2E translations
c
              do idim=1,nd
              if( ii .ge. +2) then
              call l2dexpexp_e(expout(0,1,idim,jbox),
     $           expinc(0,1,idim,ibox),nexp,cshift,ii,jj)
              else if ( ii .le. -2 ) then
              call l2dexpexp_w(expout(0,1,idim,jbox),
     $           expinc(0,1,idim,ibox),nexp,cshift,-ii,-jj)
              else if ( jj .ge. +2 ) then
              call l2dexpexp_n(expout(0,1,idim,jbox),
     $           expinc(0,1,idim,ibox),nexp,cshift,jj,-ii)
              else if ( jj .le. -2 ) then
              call l2dexpexp_s(expout(0,1,idim,jbox),
     $           expinc(0,1,idim,ibox),nexp,cshift,-jj,ii)
              endif
              enddo
c
 4150       continue
        endif
 4200   continue
C$OMP END PARALLEL DO
 4300   continue
c
        if(ifprint .ge. 1) 
     $     call prinf('=== STEP 4c (form expinc) ====*',i,0)
        t3a=second()
C$        t3a=omp_get_wtime()
c
c       ... step 4c, form incoming exponential expansions
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level,npts,nkids,radius,idim)
C$OMP$PRIVATE(mptemp,lused,ier,i,j,ptemp,gtemp,htemp,cd,lexp_tmp)
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1)
        do 1400 ibox=1,nboxes
c
        call d2tgetb(ier,ibox,box,center0,corners0,wlists)
        call d2tnkids(box,nkids)
c
        level=box(1)
        if( level .lt. 2 ) goto 1400
c
c       ... loop over densities and convert each exponential to local
c       ... NOTE: now expinc has nd dimension, so each density is separate
c
        allocate(lexp_tmp(0:nterms(level)))
        do idim=1,nd
c       ... convert exponential to local expansion for this density
        lexp_tmp = 0.0d0
        call l2dexpta(
     $     scale(level),center0,lexp_tmp,nterms(level),
     $     expinc(0,1,idim,ibox),nexp,comp)
c       ... copy back to interleaved array
        do j=0,nterms(level)
           rmlexp(iaddr(2,ibox)+(idim-1)+j*nd) =
     $        rmlexp(iaddr(2,ibox)+(idim-1)+j*nd) + lexp_tmp(j)
        enddo
        enddo
        deallocate(lexp_tmp)
c
 1400    continue
C$OMP END PARALLEL DO
        t4a=second()
C$        t4a=omp_get_wtime()
c
ccc        write(*,*) t2a-t1a, t3a-t2a, t4a-t3a, t4a-t1a
ccc        write(*,*) t2a-t1a, t3a-t2a, t4a-t3a

        t2=second()
C$        t2=omp_get_wtime()
ccc     call prin2('time=*',t2-t1,1)
        timeinfo(4)=t2-t1
c       
        if (ifprint .ge. 1) 
     $     call prinf('=== STEP 5 (split lo) ===*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 5, split all local expansions
c
ccc        do 5200 ibox=1,nboxes
        do 5300 ilev=3,nlev
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
                  if( nterms(level0)+nterms(level1) .gt. 95 ) then
                  call l2dlocloc(nd,scale(level0),center0,
     1               rmlexp(iaddr(2,ibox)),nterms(level0),
     1               scale(level1),center1,mptemp,nterms(level1))
                  else
                  call l2dlocloc_carray(nd,scale(level0),center0,
     1               rmlexp(iaddr(2,ibox)),nterms(level0),
     1               scale(level1),center1,mptemp,nterms(level1),
     1               carray,ldc)
                  endif
                  call l2dadd(nd,mptemp,rmlexp(iaddr(2,jbox)),
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
        subroutine zero_fill_real(a,n)
        implicit real *8 (a-h,o-z)
        real *8 a(*)
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i = 1,n
        a(i) = 0
        enddo
C$OMP END PARALLEL DO
        return
        end
c
c
c
c
c
        subroutine zero_fill_cmpl(a,n)
        implicit real *8 (a-h,o-z)
        complex *16 a(*)
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i = 1,n
        a(i) = 0
        enddo
C$OMP END PARALLEL DO
        return
        end
