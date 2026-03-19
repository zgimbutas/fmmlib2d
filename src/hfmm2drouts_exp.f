cc Copyright (C) 2010: Leslie Greengard and Zydrunas Gimbutas
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
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        this is the end of the debugging code and the beginning
c        of the routines for Helmholtz FMM in R^2
c
c        EXPONENTIAL VERSION with nd-density support
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine hfmm2d_list2_exp(nd,
     $     zk,bsize,nlev,laddr,scale,nterms,rmlexp,iaddr,epsfmm,
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
        complex *16 zk,zkbsize
c
        complex *16 ptemp,ftemp(2),htemp(3)
c
        real *8 timeinfo(10)
c
c
        integer, allocatable :: nexp(:)
        complex *16, allocatable :: zcc(:)
c
        complex *16, allocatable :: ws(:), xs(:)
        complex *16, allocatable :: cshift(:,:,:),cscale(:)
        complex *16, allocatable :: expout(:,:,:), expinc(:,:,:)
c
        integer, allocatable :: iq_addr(:,:)
        complex *16, allocatable :: quadr(:)

        integer, allocatable :: ie_addr(:,:)
        complex *16, allocatable :: expio(:)
c
        done = 1.0d0
        pi=4*atan(done)
c
         t1=second()
C$        t1=omp_get_wtime()

c
c       ... preallocate work and indexing arrays
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
ccc        min_exp_depth = nlev+2

        if( ifprint .eq. 1 ) then
        call prin2('rsize_wavelengths=*',rsize_wavelengths,1)
        call prin2('wb_cufoff=*',wb_cutoff,1)
        call prinf('min_exp_depth=*',min_exp_depth,1)
        endif

        do i = 2,nlev
        call hwtsexp2(ier,iprec,zk*bsize(i),
     $     xs,ws,zcc(i),nexp(i),err)
        call prin2('boxsize (wavelengths)=*',zk*bsize(i)/(2*pi),1)
        call prinf('after hwtsexp2, ier=*',ier,1)
        call prinf('after hwtsexp2, nexp=*',nexp(i),1)
ccc        call prin2('after hwtsexp2, xs=*',xs,2*nexp(i))
ccc        call prin2('after hwtsexp2, ws=*',ws,2*nexp(i))
ccc        call prin2('after hwtsexp2, zcc=*',zcc(i),2)
        enddo

        allocate(iq_addr(8,0:200))
        allocate(ie_addr(2,nboxes))

        do i = 0,200
        iq_addr(1,i) = 0
        iq_addr(2,i) = 0
        iq_addr(3,i) = 0
        iq_addr(4,i) = 0
        iq_addr(5,i) = 0
        iq_addr(6,i) = 0
        iq_addr(7,i) = 0
        iq_addr(8,i) = 0
        enddo
c
        do i = 1,nboxes
        ie_addr(1,i) = 0
        ie_addr(2,i) = 0
        enddo
c

c
c       ... allocate exponential nodes, weights, and precomputed arrays
c
        iptr = 1

        do i = 2,nlev
c
c       ... allocate memory for nodes
c
        iq_addr(1,i) = iptr
        iptr = iptr + nexp(i)
c
c       ... allocate memory for weigths
c
        iq_addr(2,i) = iptr
        iptr = iptr + nexp(i)
c
c       ... allocate memory for cscale
c
        iq_addr(3,i) = iptr
        iptr = iptr + nexp(i)
c
c       ... allocate memory for cshifts
c
        iq_addr(4,i) = iptr
        iptr = iptr + nexp(i)*7*7
c
c       ... allocate memory for zmul12fb
c
        iq_addr(5,i) = iptr
        iptr = iptr + nexp(i)*nterms(i)
c
        iq_addr(6,i) = iptr
        iptr = iptr + nexp(i)*nterms(i)
c
        iq_addr(7,i) = iptr
        iptr = iptr + nexp(i)*nterms(i)
c
        iq_addr(8,i) = iptr
        iptr = iptr + nexp(i)*nterms(i)
c
        enddo

ccc        call prinf('iptr=*',iptr,1)
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

c
c       ... allocate exponentials (nd-vectorized)
c
        iptr=1
        do ilev=2,nlev
        do ibox=laddr(1,ilev+1),laddr(1,ilev+1)+laddr(2,ilev+1)-1
c
c       ... allocate memory for outgoing exponentials
c
        ie_addr(1,ibox) = iptr
        iptr = iptr + (nexp(ilev)+1)*4*nd
c
c       ... allocate memory for incoming exponentials
c
        ie_addr(2,ibox) = iptr
        iptr = iptr + (nexp(ilev)+1)*4*nd
        enddo
        enddo
c
ccc        call prinf('iptr=*',iptr,1)
        allocate( expio(iptr) )

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i = 1,iptr
        expio(i)=0
        enddo
C$OMP END PARALLEL DO

         t2=second()
C$        t2=omp_get_wtime()

         call prin2('time in list2_exp_ini=*',t2-t1,1)

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

c       ... step 4, convert multipole expansions into the local ones
c
c       ... for levels where exponential is not available, use standard M2L
c
        do 4310 ilev=3,min(min_exp_depth-1,nlev+1)
ccc        call h2dterms_list2(bsize(ilev-1),zk,epsfmm, itable, ier)
ccc        call prinf('itable=*',itable,7*7)
c        t3=second()
cC$        t3=omp_get_wtime()
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level0,list,nlists,nlist,itype)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1,ifdirect2,radius)
C$OMP$PRIVATE(mptemp,lused,ier,i,j,ptemp,ftemp,htemp,cd,ilist)
C$OMP$PRIVATE(if_use_trunc,nterms_trunc,ii,jj)
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1)
        do 4210 ibox=laddr(1,ilev),laddr(1,ilev)+laddr(2,ilev)-1
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
           do 4151 ilist=1,nlist
              jbox=list(ilist)
              call d2tgetb(ier,jbox,box1,center1,corners1,wlists)
              if( box1(10) .eq. 0 ) goto 4151
              if (jbox.eq.0) goto 4151
              if ((box(12).eq.0).and.(ifprune_list2.eq.1))
     $           goto 4151
c              radius = (corners1(1,1) - center1(1))**2
c              radius = radius + (corners1(2,1) - center1(2))**2
c              radius = sqrt(radius)
c
c       ... convert multipole expansions for all boxes in list 2 in local exp
c       ... if source is childless, evaluate directly (if cheaper)
c

              level1=box1(1)
c              ifdirect2 = 0
c              if( box1(10) .lt. (nterms(level1)+1)/2 .and.
c     $            box(10) .lt. (nterms(level1)+1)/2  )
c     $               ifdirect2 = 1
c
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
c                 call h2dmploc_add(zk,scale(level1),center1,
c     1              rmlexp(iaddr(1,jbox)),nterms(level1),
c     $              scale(level0),
c     1              center0,rmlexp(iaddr(2,ibox)),nterms(level0))
c
c              call h2dmploc_add_trunc(zk,scale(level1),center1,
c     $           rmlexp(iaddr(1,jbox)),nterms(level1),nterms(level1),
c     $           scale(level0), center0,rmlexp(iaddr(2,ibox)),
c     $           nterms(level0),nterms(level0))
c
c               call prinf('ier=*',ier,1)
c               call prin2('in mploc, mptemp=*',
c     $            mptemp,2*(2*nterms(level0)+1))
c               call prin2('in mploc, rmlexp=*',
c     $            rmlexp(iaddr(1,ibox)),2*(2*nterms(level0)+1))
c               call prin2('in mploc, rmlexp=*',
c     $            rmlexp(iaddr(2,jbox)),2*(2*nterms(level1)+1))
c
              else

              ii=box1(2)-box(2)
              jj=box1(3)-box(3)
              nterms_trunc=itable(ii,jj)
              nterms_trunc=min(nterms(level0),nterms_trunc)
              nterms_trunc=min(nterms(level1),nterms_trunc)

              call h2dmploc_add_trunc(zk,scale(level1),center1,
     $           rmlexp(iaddr(1,jbox)),nterms(level1),nterms_trunc,
     $           scale(level0),center0,rmlexp(iaddr(2,ibox)),
     $           nterms(level0),nterms_trunc)

              endif
              endif

 4151         continue
        endif
 4210   continue
C$OMP END PARALLEL DO
 4310   continue


        if(ifprint .ge. 1)
     $     call prinf('=== STEP 4a (form expout) ====*',i,0)
        t1a=second()
C$        t1a=omp_get_wtime()
c
c       ... step 4a, form outgoing exponential expansions
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level,npts,i,j,zkbsize)
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1)
        do 1200 ibox=1,nboxes
c
        call d2tgetb(ier,ibox,box,center0,corners0,wlists)
        call d2tnkids(box,nkids)
c
        level=box(1)
        if( level .lt. min_exp_depth-1 ) goto 1200
        npts=box(10)
c
c       ... prune all sourceless boxes
c
        if( box(10) .eq. 0 ) goto 1200
c
        zkbsize = zk*bsize(level)
c
c       ... convert multipole to exponential (nd-vectorized)
c
        call h2dmpexp_fast(nd,zkbsize,
     $     scale(level),center0,rmlexp(iaddr(1,ibox)),nterms(level),
     $     expio(ie_addr(1,ibox)),nexp(level),
     $     quadr(iq_addr(1,level)),quadr(iq_addr(2,level)),
     $     zcc(level),quadr(iq_addr(3,level)),
     $     quadr(iq_addr(5,level)),quadr(iq_addr(6,level)))
c
 1200    continue
C$OMP END PARALLEL DO
        t2a=second()
C$        t2a=omp_get_wtime()
c
        if (ifprint .ge. 1)
     $     call prinf('=== STEP 4b (expout to expinc) ===*',i,0)

        do 4300 ilev=max(3,min_exp_depth),nlev+1
c        t3=second()
cC$        t3=omp_get_wtime()
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level0,list,nlists,nlist,itype)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1,ifdirect2,radius)
C$OMP$PRIVATE(mptemp,lused,ier,i,j,ptemp,ftemp,htemp,cd,ilist,level)
C$OMP$PRIVATE(if_use_trunc,nterms_trunc,ii,jj)
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1)
        do 4200 ibox=laddr(1,ilev),laddr(1,ilev)+laddr(2,ilev)-1
        call d2tgetb(ier,ibox,box,center0,corners0,wlists)
        if (ifprint .ge. 2) then
           call prinf('ibox=*',ibox,1)
           call prinf('box=*',box,15)
        endif
        level=box(1)
        if (level .ge. 2) then
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
c       in local exp (nd-vectorized)
c
              ii=box(2)-box1(2)
              jj=box(3)-box1(3)
c
              if( ii .ge. +2) then
              call h2dexpexp_e(nd,
     $           expio(ie_addr(1,jbox)),expio(ie_addr(2,ibox)),
     $           nexp(level),quadr(iq_addr(4,level)),ii,jj)
              else if ( ii .le. -2 ) then
              call h2dexpexp_w(nd,
     $           expio(ie_addr(1,jbox)),expio(ie_addr(2,ibox)),
     $           nexp(level),quadr(iq_addr(4,level)),-ii,-jj)
              else if ( jj .ge. +2 ) then
              call h2dexpexp_n(nd,
     $           expio(ie_addr(1,jbox)),expio(ie_addr(2,ibox)),
     $           nexp(level),quadr(iq_addr(4,level)),jj,-ii)
              else if ( jj .le. -2 ) then
              call h2dexpexp_s(nd,
     $           expio(ie_addr(1,jbox)),expio(ie_addr(2,ibox)),
     $           nexp(level),quadr(iq_addr(4,level)),-jj,ii)
              endif
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
C$OMP$PRIVATE(ibox,box,center0,corners0,level,npts,nkids,radius)
C$OMP$PRIVATE(lused,ier,i,j,ptemp,gtemp,htemp,cd,zkbsize)
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1)
        do 1400 ibox=1,nboxes
c
        call d2tgetb(ier,ibox,box,center0,corners0,wlists)
        call d2tnkids(box,nkids)
c
        level=box(1)
        if( level .lt. min_exp_depth-1 ) goto 1400
c
        zkbsize = zk*bsize(level)
c
c       ... convert exponential to local (nd-vectorized)
c
        call h2dexpta_fast(nd,zkbsize,
     $     scale(level),center0,rmlexp(iaddr(2,ibox)),nterms(level),
     $     expio(ie_addr(2,ibox)),nexp(level),
     $     quadr(iq_addr(1,level)),quadr(iq_addr(2,level)),
     $     quadr(iq_addr(7,level)),quadr(iq_addr(8,level)))
c
 1400    continue
C$OMP END PARALLEL DO
        t4a=second()
C$        t4a=omp_get_wtime()

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
        subroutine zero_fill(a,n)
        implicit real *8 (a-h,o-z)
        complex *16 a(*)
        do i = 1,n
        a(i) = dcmplx(0.0d0,0.0d0)
        enddo
        return
        end
