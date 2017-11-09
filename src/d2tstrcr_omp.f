cc Copyright (C) 2010: Vladimir Rokhlin
cc Copyright (C) 2009-2012: Leslie Greengard and Zydrunas Gimbutas
cc Contact: greengard@cims.nyu.edu
cc 
cc This software is being released under a modified FreeBSD license
cc (see COPYING in home directory). 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date: 2011-04-15 22:14:18 -0400 (Fri, 15 Apr 2011) $
c    $Revision: 1834 $
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       this is the end of the debugging code and the beginning
c       of the actual logic subroutines for the FMM in R^2
c
c       Fortran 95 version
c       
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine d2tstrcr(ier,z,n,nbox,
     1    nboxes,iz,laddr,nlev,center,size,
     $    ztarg,ntarg,iztarg,w,lw,lused777)
        implicit real *8 (a-h,o-z)
        integer iz(*),iztarg(*),w(*),laddr(2,*)
        real *8 z(2,*),ztarg(2,*),center(2),corners(2,4)
        real *8, allocatable :: z0(:,:),z0targ(:,:)
ccc        save
c
c        this subroutine constructs the logical structure for the 
c        fully adaptive FMM in two dimensions and stores it in the
c        array w in the form of a link-list. after that, the user 
c        can obtain the information about various boxes and lists 
c        in it by calling the entries d2tgetb, d2tgetl, d2tlinfo
c        of this subroutine (see).
c
c              note on the list conventions. 
c
c    list 1 of the box ibox - the list of all boxes with which the
c           box ibox interacts directly, including the boxes on the 
c           same level as ibox, boxes on the finer levels, and boxes 
c           on the coarser levels. obviously, list 1 is empty for any
c           box that is not childless.
c
c    list 2 of the box ibox - the list of all boxes with which the
c           box ibox interacts in the regular multipole fashion, i.e. 
c           boxes on the same level as ibox that are separated from it
c           but whose daddies are not separated from the daddy of ibox.
c
c    list 3 of the box ibox - for a childless ibox, the list of all 
c           boxes on the levels finer than that of ibox, which are 
c           separated from ibox, and whose daddys are not separated 
c           from ibox. for a box with children, list 3 is empty.
c           
c    list 4 is dual to list 3, i.e. jbox is on the list 4 of ibox if 
c           and only if ibox is on the list 3 of jbox. 
c
c    list 5 of the box ibox - the list of all boxes at the same level 
c           that are adjacent to the box ibox - the list of colleagues
c
c                            input parameters:
c
c  z - the user-specified points in the space
c  n - the number of elements in array z
c  nbox - the maximum number of points in a box on the finest level
c  lw - the amount of memory in the array w (in integer elements)
c
c                            output parameters:
c
c  ier - error return code
c    ier=0   means successful execution
c    ier=32  means that the amount lw of space in array w
c                 is insufficient
c    ier=16 means that the subroutine attempted to construct more 
c        than 199 levels of subdivision; indicates bad trouble.
c    ier=64  means that the amount lw of space in array w
c                 is severely insufficient
c  nboxes - the total number of boxes created
c  iz - the integer array addressing the particles in 
c         all boxes. 
c       explanation: for a box ibox, the particles living in
c         it are:
c
c         (z(1,j),z(2,j)),(z(1,j+1),z(2,j+1)),
c         (z(1,j+2),z(2,j+2)), . . . 
c         (z(1,j+nj-1),z(2,j+nj-1)),
c         (z(1,j+nj),z(2,j+nj)),
c
c         with j=boxes(9,ibox), and nj=boxes(10,ibox)
c        
c  iztarg - the integer array addressing the targets in 
c         all boxes. 
c
c         (ztarg(1,j),ztag(2,j)),
c         (ztarg(1,j+1),ztarg(2,j+1)),
c         (ztarg(1,j+2),ztarg(2,j+2)), . . . 
c         (ztarg(1,j+nj-1),ztarg(2,j+nj-1)),
c         (ztarg(1,j+nj),ztarg(2,j+nj)),
c
c         with j=boxes(11,ibox), and nj=boxes(12,ibox)
c
c  laddr - an integer array dimensioned (2,nlev), describing the
c         numbers of boxes on various levels of sybdivision, so that
c         the first box on level (i-1) has sequence number laddr(1,i),
c         and there are laddr(2,i) boxes on level i-1
c  nlev - the maximum level number on which any boxes have 
c         been created. the maximum number possible is 200. 
c         it is recommended that the array laddr above be 
c         dimensioned at least (2,200), in case the user underestimates
c         the number of levels required.
c  center - the center of the box on the level 0, containing
c         the whole simulation
c  size - the side of the box on the level 0
c  w - the array containing all tables describing boxes, lists, etc. 
c         it is a link-list (for the most part), and can only be accessed
c         via the entries d2tgetb, d2tgetl, d2tlinfo, of this  subroutine 
c         (see below). the first lused 777 integer locations of 
c         this array should not be altered between the call to this
c         entry and subsequent calls to the entries d2tgetb, d2tgetl,
c         d2tlinfo, of this  subroutine 
c        
c  lused777 - the amount of space in the array w (in integer words)
c        that is occupied by various tables on exit from this 
c        subroutine. this space should not be altered between the
c        call to this entry and subsequent calls to entries d2tgetb,
c        d2tgetl, d2tlinfo, of this  subroutine (see below).
c  
c        . . . construct the quad-tree structure for the user-specified 
c              set of points
c
c       size of real *8 must not exceed the size of two integers
c
        if( n .lt. 1 ) then
c       number of particles less than one, abort
        ier=128
        return
        endif
c
        ier=0
c
        ninire=2
c
        iptr=1
        lptr=500
c
        iiwork=iptr+lptr
        liwork=n+ntarg+4
c     
        iboxes=iiwork+liwork
        lboxes=lw-n-5
        maxboxes=lboxes/15-1
c
c
c        if the memory is insufficient - bomb
c
        if( lw .lt. 12*(n+ntarg) ) then
        ier=64
ccc        call prinf('in d2tstrcr before d2tallb, ier=*',ier,1)
        return
        endif
c
c	 initialize the sorting index 
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)        
        do i=1,n
	iz(i)=i
	enddo
C$OMP END PARALLEL DO
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)        
        do i=1,ntarg
	iztarg(i)=i
	enddo
C$OMP END PARALLEL DO
c
        allocate( z0(2,n) )
        if( ntarg .gt. 0 ) then
        allocate( z0targ(2,ntarg) )
        else
        allocate( z0targ(2,1) )
        endif
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)        
        do i=1,n
	z0(1,i)=z(1,i)
	z0(2,i)=z(2,i)
	enddo
C$OMP END PARALLEL DO
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)        
        do i=1,ntarg
	z0targ(1,i)=ztarg(1,i)
	z0targ(2,i)=ztarg(2,i)
	enddo
C$OMP END PARALLEL DO
c
        t1=second()
C$        t1=omp_get_wtime()
        ifempty=0
        minlevel=0
        maxlevel=100
        call d2tallbem(ier,z0,n,nbox,w(iboxes),maxboxes,
     1    nboxes,iz,laddr,nlev,center,size,w(iiwork),
     $     ifempty,minlevel,maxlevel,z0targ,ntarg,iztarg)
        t2=second()
C$        t2=omp_get_wtime()
        call prin2('time in d2tallbem=*',t2-t1,1)
c
ccc        call d2tallb(ier,z,n,nbox,w(iboxes),maxboxes,
ccc     1    nboxes,iz,laddr,nlev,center,size,w(iiwork) )
c
c        if the memory is insufficient - bomb
c
        if(ier .eq. 0) goto 1100
ccc           call prinf('in d2tstrcr after d2tallb, ier=*',ier,1)
        if(ier .eq. 4) ier=32
        return
 1100 continue
c
c       compress the array w
c   
        nn=nboxes*15
        do 1200 i=1,nn
        w(iiwork+i-1)=w(iboxes+i-1)
 1200 continue
        iboxes=iiwork
 1300 continue
        lboxes=nboxes*15+16
c
c       ... align array for real *16 storage
        lboxes=lboxes+4-mod(lboxes,4)
c
c       construct the centers and the corners for all boxes
c       in the quad-tree
c
        icenters=iboxes+lboxes
        lcenters=(nboxes*2+2)*ninire
c
        icorners=icenters+lcenters
        lcorners=(nboxes*8+2)*ninire
c
        iwlists=icorners+lcorners
        lwlists=lw-iwlists-6
c
        call prinf('lused: ccenters(k)=*', (lcenters+lcorners)/1000,1)
        call prinf('lused(k)=*', (iwlists)/1000,1)
        call d2tcentc(center,size,w(iboxes),nboxes,
     1      w(icenters),w(icorners) )
c
c       now, construct all lists for all boxes
c
        
ccc           call prinf('before d2tlsts, lwlists=*',lwlists,1)
        t1=second()
C$        t1=omp_get_wtime()
        call d2tlsts(ier,w(iboxes),nboxes,w(icorners),
     1        w(iwlists),lwlists,lused)
        t2=second()
C$        t2=omp_get_wtime()
        call prin2('time in d2tlsts=*',t2-t1,1)
c
        lused777=lused+iwlists
ccc        call prinf('after d2tlsts, ier=*',ier,1)
        call prinf('lused(k)=*', (lused777)/1000,1)
c
c       store all pointers
c
        w(1)=nboxes
        w(2)=iboxes
        w(3)=icorners
        w(4)=icenters
        w(5)=iwlists
        w(6)=lused777
c
        w(7)=n
        w(8)=nbox
        w(9)=nlev
        w(10)=ier
        w(11)=0
        w(12)=ifempty
        w(13)=minlevel
        w(14)=maxlevel
        do i=1,200
        w(100+2*i-2)=laddr(1,i)
        w(100+2*i-1)=laddr(2,i)
        enddo
c
        return
        end
c
c
c
c
        subroutine d2tstrcrem(ier,z,n,nbox,
     1    nboxes,iz,laddr,nlev,center,size,
     $    ztarg,ntarg,iztarg,w,lw,lused777,
     $    ifempty,minlevel,maxlevel)
        implicit real *8 (a-h,o-z)
        integer iz(*),iztarg(*),w(*),laddr(2,*)
        real *8 z(2,*),ztarg(2,*),center(2),corners(2,4)
        real *8, allocatable :: z0(:,:),z0targ(:,:)
ccc        save
c
c        this subroutine constructs the logical structure for the 
c        fully adaptive FMM in two dimensions and stores it in the
c        array w in the form of a link-list. after that, the user 
c        can obtain the information about various boxes and lists 
c        in it by calling the entries d2tgetb, d2tgetl, d2tlinfo
c        of this subroutine (see).
c
c              note on the list conventions. 
c
c    list 1 of the box ibox - the list of all boxes with which the
c           box ibox interacts directly, including the boxes on the 
c           same level as ibox, boxes on the finer levels, and boxes 
c           on the coarser levels. obviously, list 1 is empty for any
c           box that is not childless.
c
c    list 2 of the box ibox - the list of all boxes with which the
c           box ibox interacts in the regular multipole fashion, i.e. 
c           boxes on the same level as ibox that are separated from it
c           but whose daddies are not separated from the daddy of ibox.
c
c    list 3 of the box ibox - for a childless ibox, the list of all 
c           boxes on the levels finer than that of ibox, which are 
c           separated from ibox, and whose daddys are not separated 
c           from ibox. for a box with children, list 3 is empty.
c           
c    list 4 is dual to list 3, i.e. jbox is on the list 4 of ibox if 
c           and only if ibox is on the list 3 of jbox. 
c
c    list 5 of the box ibox - the list of all boxes at the same level 
c           that are adjacent to the box ibox - the list of colleagues
c
c                            input parameters:
c
c  z - the user-specified points in the space
c  n - the number of elements in array z
c  nbox - the maximum number of points in a box on the finest level
c  lw - the amount of memory in the array w (in integer elements)
c
c  ifempty - ifempty=0 - remove empty boxes, ifempty=1 - keep empty boxes 
c  minlevel - minimum level of refinement
c  maxlevel - maximum level of refinement
c
c                            output parameters:
c
c  ier - error return code
c    ier=0   means successful execution
c    ier=32  means that the amount lw of space in array w
c                 is insufficient
c    ier=16 means that the subroutine attempted to construct more 
c        than 197 levels of subdivision; indicates bad trouble.
c    ier=64  means that the amount lw of space in array w
c                 is severely insufficient
c  nboxes - the total number of boxes created
c  iz - the integer array addressing the particles in 
c         all boxes. 
c       explanation: for a box ibox, the particles living in
c         it are:
c
c         (z(1,j),z(2,j)),(z(1,j+1),z(2,j+1)),
c         (z(1,j+2),z(2,j+2)), . . . 
c         (z(1,j+nj-1),z(2,j+nj-1)),
c         (z(1,j+nj),z(2,j+nj)),
c
c         with j=boxes(9,ibox), and nj=boxes(10,ibox)
c        
c  iztarg - the integer array addressing the targets in 
c         all boxes. 
c
c         (ztarg(1,j),ztag(2,j)),
c         (ztarg(1,j+1),ztarg(2,j+1)),
c         (ztarg(1,j+2),ztarg(2,j+2)), . . . 
c         (ztarg(1,j+nj-1),ztarg(2,j+nj-1)),
c         (ztarg(1,j+nj),ztarg(2,j+nj)),
c
c         with j=boxes(11,ibox), and nj=boxes(12,ibox)
c
c  laddr - an integer array dimensioned (2,nlev), describing the
c         numbers of boxes on various levels of sybdivision, so that
c         the first box on level (i-1) has sequence number laddr(1,i),
c         and there are laddr(2,i) boxes on level i-1
c  nlev - the maximum level number on which any boxes have 
c         been created. the maximim number possible is 200. 
c         it is recommended that the array laddr above be 
c         dimensioned at least (2,200), in case the user underestimates
c         the number of levels required.
c  center - the center of the box on the level 0, containing
c         the whole simulation
c  size - the side of the box on the level 0
c  w - the array containing all tables describing boxes, lists, etc. 
c         it is a link-list (for the most part), and can only be accessed
c         via the entries d2tgetb, d2tgetl, d2tlinfo, of this  subroutine 
c         (see below). the first lused 777 integer locations of 
c         this array should not be altered between the call to this
c         entry and subsequent calls to the entries d2tgetb, d2tgetl,
c         d2tlinfo, of this  subroutine 
c        
c  lused777 - the amount of space in the array w (in integer words)
c        that is occupied by various tables on exit from this 
c        subroutine. this space should not be altered between the
c        call to this entry and subsequent calls to entries d2tgetb,
c        d2tgetl, d2tlinfo, of this  subroutine (see below).
c  
c        . . . construct the quad-tree structure for the user-specified 
c              set of points
c
c       size of real *8 must not exceed the size of two integers
c
        if( n .lt. 1 ) then
c       number of particles less than one, abort
        ier=128
        return
        endif
c
        ier=0
c
        ninire=2
c
        iptr=1
        lptr=500
c
        iiwork=iptr+lptr
        liwork=n+ntarg+4
c     
        iboxes=iiwork+liwork
        lboxes=lw-n-5
        maxboxes=lboxes/15-1
c
c
c        if the memory is insufficient - bomb
c
        if( lw .lt. 12*(n+ntarg) ) then
        ier=64
ccc        call prinf('in d2tstrcr before d2tallb, ier=*',ier,1)
        return
        endif
c
c	 initialize the sorting index 
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)        
        do i=1,n
	iz(i)=i
	enddo
C$OMP END PARALLEL DO
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)        
        do i=1,ntarg
	iztarg(i)=i
	enddo
C$OMP END PARALLEL DO
c
        allocate( z0(2,n) )
        if( ntarg .gt. 0 ) then
        allocate( z0targ(2,ntarg) )
        else
        allocate( z0targ(2,1) )
        endif
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)        
        do i=1,n
	z0(1,i)=z(1,i)
	z0(2,i)=z(2,i)
	enddo
C$OMP END PARALLEL DO
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)        
        do i=1,ntarg
	z0targ(1,i)=ztarg(1,i)
	z0targ(2,i)=ztarg(2,i)
	enddo
C$OMP END PARALLEL DO
c
        t1=second()
C$        t1=omp_get_wtime()
c        ifempty=0
c        minlevel=0
c        maxlevel=100
        call d2tallbem(ier,z0,n,nbox,w(iboxes),maxboxes,
     1    nboxes,iz,laddr,nlev,center,size,w(iiwork),
     $     ifempty,minlevel,maxlevel,z0targ,ntarg,iztarg)
        t2=second()
C$        t2=omp_get_wtime()
        call prin2('time in d2tallbem=*',t2-t1,1)
c
ccc        call d2tallb(ier,z,n,nbox,w(iboxes),maxboxes,
ccc     1    nboxes,iz,laddr,nlev,center,size,w(iiwork) )
c
c        if the memory is insufficient - bomb
c
        if(ier .eq. 0) goto 1100
ccc           call prinf('in d2tstrcr after d2tallb, ier=*',ier,1)
        if(ier .eq. 4) ier=32
        return
 1100 continue
c
c       compress the array w
c   
        nn=nboxes*15
        do 1200 i=1,nn
        w(iiwork+i-1)=w(iboxes+i-1)
 1200 continue
        iboxes=iiwork
 1300 continue
        lboxes=nboxes*15+16
c
c       ... align array for real *16 storage
        lboxes=lboxes+4-mod(lboxes,4)
c
c       construct the centers and the corners for all boxes
c       in the quad-tree
c
        icenters=iboxes+lboxes
        lcenters=(nboxes*2+2)*ninire
c
        icorners=icenters+lcenters
        lcorners=(nboxes*8+2)*ninire
c
        iwlists=icorners+lcorners
        lwlists=lw-iwlists-6
c
        call prinf('lused: ccenters(k)=*', (lcenters+lcorners)/1000,1)
        call prinf('lused(k)=*', (iwlists)/1000,1)
        call d2tcentc(center,size,w(iboxes),nboxes,
     1      w(icenters),w(icorners) )
c
c       now, construct all lists for all boxes
c
        
ccc           call prinf('before d2tlsts, lwlists=*',lwlists,1)
        t1=second()
C$        t1=omp_get_wtime()
        call d2tlsts(ier,w(iboxes),nboxes,w(icorners),
     1        w(iwlists),lwlists,lused)
        t2=second()
C$        t2=omp_get_wtime()
        call prin2('time in d2tlsts=*',t2-t1,1)
c
        lused777=lused+iwlists
ccc        call prinf('after d2tlsts, ier=*',ier,1)
        call prinf('lused(k)=*', (lused777)/1000,1)
c
c       store all pointers
c
        w(1)=nboxes
        w(2)=iboxes
        w(3)=icorners
        w(4)=icenters
        w(5)=iwlists
        w(6)=lused777
c
        w(7)=n
        w(8)=nbox
        w(9)=nlev
        w(10)=ier
        w(11)=0
        w(12)=ifempty
        w(13)=minlevel
        w(14)=maxlevel
        do i=1,200
        w(100+2*i-2)=laddr(1,i)
        w(100+2*i-1)=laddr(2,i)
        enddo
c
        return
        end
c
c
c
c
        subroutine d2tnkids(box,nkids)
        implicit real *8 (a-h,o-z)
        integer box(15)
c       
        nkids=0
        do ikid=1,4
        if( box(4+ikid) .ne. 0 ) nkids=nkids+1
        enddo
        return
        end
c
c
c
c
c
        subroutine d2tprint(w,lw)
        implicit real *8 (a-h,o-z)
        integer w(*),box(15)
        real *8 center0(2),corners0(2,4)
c
c       ... retrieve stored quad-tree parameters 
c
        ibox=1
        call d2tgetb(ier,ibox,box,center0,corners0,w)
c
        call prin2('after d2tstrcr, center0=*',center0,3)
c
        size0=corners0(1,3)-corners0(1,1)
        call prin2('after d2tstrcr, size0=*',size0,1)
c
        nlev0=w(9)
        call prinf('after d2tstrcr, nlev0=*',nlev0,1)
c
        nbox0=w(8)
        call prinf('after d2tstrcr, nbox0=*',nbox0,1)        
c
        call prinf('after d2tstrcr, laddr0=*',w(100),2*(nlev0+1))
c
        return
        end
c
c
c
        subroutine d2trestore(nboxes,laddr,nlev,center,size,w,lw)
        implicit real *8 (a-h,o-z)
        integer w(*),box(15)
        real *8 center(2),center0(2),corners0(2,4)
        integer laddr(2,200)
c
c       ... retrieve stored quad-tree parameters 
c
        nboxes=w(1)
c
        do i=1,nlev+1
        laddr(1,i)=w(100+2*i-2)
        laddr(2,i)=w(100+2*i-1)
        enddo
c
        nlev=w(9)
c
        ibox=1
        call d2tgetb(ier,ibox,box,center0,corners0,w)
c
        center(1)=center0(1)
        center(2)=center0(2)
c
        size=corners0(1,3)-corners0(1,1)
c
        return
        end
c
c
c
c
        subroutine d2tgetb(ier,ibox,box,center,corners,w)
        implicit real *8 (a-h,o-z)
        integer w(*),box(15),nums(*)
        real *8 center(2),corners(2,4)
c
c        this entry returns to the user the characteristics of
c        user-specified box  ibox.  
c
c                     input parameters:
c
c  ibox - the box number for which the information is desired
c  w - storage area as created by the entry d2tstrcr (see above)
c
c                     output parameters:
c
c  ier - the error return code.
c    ier=0 means successful execution
c    ier=4 means that ibox is either greater than the number of boxes
c           in the structure or less than 1.
c  box - an integer array dimensioned box(15). its elements describe 
c        the box number ibox, as follows:
c
c       1. level - the level of subdivision on which this box 
c             was constructed; 
c       2, 3  - the coordinates of this box among  all
c             boxes on this level
c       4 - the daddy of this box, identified by it address
c             in array boxes
c       5,6,7,8 - the  list of children of this box 
c             (eight of them, and the child is identified by its address
c             in the array boxes; if a box has only one child, only the
c             first of the four child entries is non-zero, etc.)
c       9 - the location in the array iz of the particles 
c             living in this box
c       10 - the number of particles living in this box
c       11 - the location in the array iztarg of the targets
c             living in this box
c       12 - the number of targets living in this box
c       13 - source box type: 0 - empty, 1 - leaf node, 2 - sub-divided
c       14 - target box type: 0 - empty, 1 - leaf node, 2 - sub-divided
c       15 - reserved for future use
c  center - the center of the box number ibox 
c  corners - the corners of the box number ibox 
c
c       . . . return to the user all information about the box ibox
c 
        nboxes=w(1)
        iboxes=w(2)
        icorners=w(3)
        icenters=w(4)
        iwlists=w(5)
c 
        ier=0
        if( (ibox .ge. 1)  .and. (ibox .le. nboxes) ) goto 2100
        ier=4
        return
 2100 continue
c
        ibox0=iboxes+(ibox-1)*15-1
        do 2200 i=1,15
        box(i)=w(ibox0+i)
 2200 continue
c
c      return to the user the center and the corners of the box ibox
c
        call d2tcpcc(w(icenters),w(icorners),ibox,center,corners) 
c
        return
c
c
c
c
         entry d2tgetl(ier,ibox,itype,list,nlist,w)
c
c  ibox - the box number for which the information is desired
c  itype - the type of the desired list for the box ibox
c  w - storage area as created by the entry d2tstrcr (see above)
c
c                     output parameters:
c
c  ier - the error return code.
c    ier=0 means successful execution
c    ier=4 means that the list  itype  for the box  ibox  is empty
c  list - the list  itype  for the box  ibox 
c  nlist - the number of elements in array  list
c
c       return to the user the list number itype for the box ibox
c
        iwlists=w(5)
c
        call d2tlinkretr(ier,itype,ibox,list,nlist,w(iwlists),lused)
        return
c
c
c
c
        entry d2tlinfo(w,lused77,nums)
c
c       this entry returns to the user some of the information
c       about the storage area. namely, it returns the total 
c       amount lused777 of memory utilized in the array w (in integer 
c       locations), and the integer array nums containing the numbers
c       of elements in each of the lists
c
        iwlists=w(5)
c       
        call d2tlinkinfo(w(iwlists),lused77,nums)
        return
        end
c
c
c
c
c
        subroutine d2tlsts(ier,boxes,nboxes,corners,w,lw,lused)
        implicit real *8 (a-h,o-z)
        integer boxes(15,*),collkids(50000),w(*),
     1      dadcolls(2000),list5(20000),stack(60000)
        real *8 corners(2,4,*)
ccc        save
c
c        this subroutine constructs all lists for all boxes 
c        and stores them in the storage area w in the form
c        of a link list. the resulting data can be accessed 
c        by calls to various entries of the subroutine d2tlinkinit (see).
c
c                          input parameters:
c
c  boxes - an integer array dimensioned (15,nboxes), as created by 
c        the subroutine d2tallb (see).  each 15-element column
c         describes one box, as follows:
c
c
c       1. level - the level of subdivision on which this box 
c             was constructed; 
c       2, 3  - the coordinates of this box among  all
c             boxes on this level
c       4 - the daddy of this box, identified by it address
c             in array boxes
c       5,6,7,8 - the  list of children of this box 
c             (eight of them, and the child is identified by its address
c             in the array boxes; if a box has only one child, only the
c             first of the four child entries is non-zero, etc.)
c       9 - the location in the array iz of the particles 
c             living in this box
c       10 - the number of particles living in this box
c       11 - the location in the array iztarg of the targets
c             living in this box
c       12 - the number of targets living in this box
c       13 - source box type: 0 - empty, 1 - leaf node, 2 - sub-divided
c       14 - target box type: 0 - empty, 1 - leaf node, 2 - sub-divided
c       15 - reserved for future use
c
c  nboxes - the total number of boxes created
c  corners - the array of corners of all the boxes to in array boxes
c  lw - the total amount of storage in array w (in integer words)
c 
c              output parameters:
c
c  ier - the error return code.
c    ier=0 means successful execution
c    ier=32 means that the amount lw of space in array w is
c           insufficient. it is a fatal error.
c  w - storage area containing all lists for all boxes in 
c        the form of link-lists, accessible by the subroutine 
c        d2tlinkretr (see).
c  lused - the amount of space in the array w (in integer words)
c        that is occupied by various tables on exit from this 
c        subroutine. this space should not be altered between the
c        call to this subroutine and subsequent calls to the 
c        entries d2tlinkretr, etc. of the subroutine d2tlinkinit (see).
c
c       . . . initialize the storage-retrieval routine for all 
c             boxes
c
ccc         call prinf('in d2tlsts, nboxes=*',nboxes,1)
c       
         lused=0
         lused2=0
c
        ier=0
        ntypes=5
        call d2tlinkinit(jer,nboxes,ntypes,w,lw)
cccc         call prinf('in d2tlsts after d2tlinkinit, ier=*',ier,1)
c
c        construct lists 5,2 for all boxes
c
        do 2000 ibox=2,nboxes
c
c       find this guy's daddy
c
        idad=boxes(4,ibox)
c
c       find daddy's collegues, including daddy himself
c
        dadcolls(1)=idad
        itype5=5
        itype2=2
        call d2tlinkretr(jer,itype5,idad,dadcolls(2),ncolls,w,lused)
        ncolls=ncolls+1
c
c        find the children of the daddy's collegues
c
        nkids=0
        do 1600 i=1,ncolls
        icoll=dadcolls(i)
        do 1400 j=1,4
        kid=boxes(4+j,icoll)
        if(kid .le. 0) goto 1600
        if(kid .eq. ibox) goto 1400
        nkids=nkids+1
        collkids(nkids)=kid
 1400 continue
 1600 continue
c
c       sort the kids of the daddy's collegues into the 
c       lists 2, 5 of the box ibox
c
        nlist1=1
        do 1800 i=1,nkids
c
c       check if this kid is touching the box ibox
c
        kid=collkids(i)
ccc        call d2tifint(corners(1,1,kid),corners(1,1,ibox),ifinter)
        call d2tifint2(boxes(1,kid),boxes(1,ibox),ifinter)
c
        if(ifinter .eq. 1)
     1    call d2tlinkstor(ier,itype5,ibox,kid,nlist1,w,lused)
c
c        if storage capacity has been exceeed - bomb
c
        if(ier .eq. 32) return        
        if(ifinter .eq. 0)
     1    call d2tlinkstor(ier,itype2,ibox,kid,nlist1,w,lused)
 1800 continue
c
c        if storage capacity has been exceeed - bomb
c
        if(ier .eq. 32) return        
 2000 continue
c
ccc        call prinf('constructed lists 5,2; lused=*',lused,1)
c
c       now, construct lists 1, 3
c
        do 3000 i=1,nboxes
c
c       if this box has kids - its lists 1, 3 are empty;
c       do not construct them
c
        if(boxes(5,i) .gt. 0) goto 3000
c
c       do not construct lists 1, 3 for the main box
c
        if(boxes(1,i) .eq. 0) goto 3000
c
        call d2tlinkretr(jer,itype5,i,list5,nlist,w,lused)  
c
        if(jer .eq. 4) goto 3000
c
        do 2200 j=1,nlist
        jbox=list5(j)
        call d2tlst31(ier,i,jbox,boxes,nboxes,
     1    corners,w,stack,lused)
c
c        if storage capacity has been exceeded - bomb
c
        if(ier .eq. 32) return        
 2200 continue
 3000 continue
c
ccc        call prinf('constructed lists 1,3, lused=*',lused,1)
c
        if( 1 .eq. 2 ) then
c
c       copy all elements of lists 1, 2, 3, and 5 while skipping list 4
c       this is not needed, d2tlst31 is skipping list 4 anyway
c
        ntypes5=5
        call d2tlinkinit(jer,nboxes,ntypes5,w(lused+1),lw-(lused+5))
        do 3600 ibox=1,nboxes
        do 2400 itype=1,5
c
        call d2tlinkretr(jer,itype,ibox,list5,nlist,w,lused)
        if(jer .eq. 4) goto 2400
        call d2tlinkstor(ier,itype,ibox,list5,nlist,w(lused+1),lused2)
c
c        if storage capacity has been exceeed - bomb
c
        if(ier .eq. 32) return
 2400 continue
 3600 continue
c
c       compress array w
c
        do 4000 i=1,lused2
        w(i)=w(lused+i)
 4000 continue
        lused=lused2
c
        endif
c
c        finally, construct the lists 4 for all boxes 
c        that need them
c
        itype3=3
        itype4=4
        nlist1=1
        do 4400 ibox=1,nboxes
c
        call d2tlinkretr(jer,itype3,ibox,list5,nlist,w,lused)
        if(jer .eq. 4) goto 4400
        do 4200 j=1,nlist
        call d2tlinkstor(ier,itype4,list5(j),ibox,nlist1,w,lused2)
c
c        if storage capacity has been exceeed - bomb
c
        if(ier .eq. 32) return        
 4200 continue
 4400 continue
        lused=lused2

ccc         call prinf('exiting d2tlsts, lused=*',lused,1)
        return
        end        
c
c
c
c
c
        subroutine d2tlst31(ier,ibox,jbox0,boxes,nboxes,
     1    corners,w,stack,lused)
        implicit real *8 (a-h,o-z)
        integer w(*)
        real *8 corners(2,4,*)
        integer boxes(15,*),stack(3,*)
        data itype1/1/,itype2/2/,itype3/3/,itype4/4/,itype5/5/,
     1      nlist1/1/
ccc        save
c
c       this subroutine constructs all elements of lists 1 and 3 
c       resulting from the subdivision of one element of list 5
c       of the box ibox. all these elements of lists 1, 3 are 
c       stored in the link-lists by the subroutine linstro (see)
c
c        input parameters:
c
c  ibox - the box whose lists are being constructed
c  
c  jbox0 - the element of list 5 of the box ibox that is being
c          subdivided
c  boxes - the array boxes as created by the subroutine d2tallb (see)
c  nboxes - the number of boxes in array boxes
c  corners - the array of corners of all the boxes to in array boxes
c  w - the storage area formatted by the subroutine d2tlinkinit (see) 
c          to be used to store the elements of lists 1, 3 constructed
c          by this subroutine. obviously, by this time, it contains
c          plenty of other lists.
c  
c                      output parameters:
c
c  ier - the error return code.
c    ier=0 means successful execution
c    ier=32 means that the amount lw of space in array w is
c           insufficient. it is a fatal error.
c  w - the augmented storage area, containing all the boxes just 
c          created, in addition to whatever had been stored previously
c  lused - the total length of array w (in integer words) used 
c          on exit from this subroutine
c
c                      work arrays:
c
c  stack - must be at least 600 integer locations long
c
c       . . . starting with the initial element of list 5, 
c             subdivide the boxes recursively and store 
c             the pieces where they belong
c
c        . . . initialize the process
c
        jbox=jbox0
        istack=1
        stack(1,1)=1
        stack(2,1)=jbox
c
        nsons=0
        do 1200 j=5,8
c
        if(boxes(j,jbox) .gt. 0) nsons=nsons+1
 1200 continue
c
        stack(3,1)=nsons
c
c       . . . move up and down the stack, generating the elements 
c             of lists 1, 3 for the box jbox, as appropriate
c
        do 5000 ijk=1,1 000 000 000
c
c       if this box is separated from ibox - store it in list 3;
c       enter this fact in the daddy's table; pass control
c       to the daddy
c
ccc        call d2tifint(corners(1,1,ibox),corners(1,1,jbox),ifinter)
        call d2tifint2(boxes(1,ibox),boxes(1,jbox),ifinter)
c
        if(ifinter .eq. 1) goto 2000
        call d2tlinkstor(ier,itype3,ibox,jbox,nlist1,w,lused)
c
c        if storage capacity has been exceeed - bomb
c
        if(ier .eq. 32) return        
        istack=istack-1
        stack(3,istack)=stack(3,istack)-1
        jbox=stack(2,istack)
        goto 5000
 2000 continue
c
c       this box is not separated from ibox. if it is childless 
c       - enter it in list 1; enter this fact in the daddy's table; 
c       pass control to the daddy
c       
        if(boxes(5,jbox) .ne. 0) goto 3000
        call d2tlinkstor(ier,itype1,ibox,jbox,nlist1,w,lused)
c
c        if storage capacity has been exceeed - bomb
c
        if(ier .eq. 32) return        
c
c       . . . entered jbox in the list1 of ibox; if jbox
c             is on the finer level than ibox - enter ibox
c             in the list 1 of jbox
c
        if(boxes(1,jbox) .eq. boxes(1,ibox)) goto 2400
        call d2tlinkstor(ier,itype1,jbox,ibox,nlist1,w,lused)
c
c        if storage capacity has been exceeed - bomb
c
        if(ier .eq. 32) return        
 2400 continue
c
c       if we have processed the whole box jbox0, get out
c       of the subroutine
c
        if(jbox .eq. jbox0) return
c
        istack=istack-1
        stack(3,istack)=stack(3,istack)-1
        jbox=stack(2,istack)
        goto 5000
 3000 continue
c
c       this box is not separated from ibox, and has children. if 
c       the number of unprocessed sons of this box is zero 
c       - pass control to his daddy
c
        nsons=stack(3,istack)
        if(nsons .ne. 0) goto 4000
c
        if(jbox .eq. jbox0) return
c
        istack=istack-1
        stack(3,istack)=stack(3,istack)-1
        jbox=stack(2,istack)
        goto 5000
 4000 continue
c
c       this box is not separated from ibox; it has sons, and
c       not all of them have been processed. construct the stack
c       element for the appropriate son, and pass the control
c       to him. 
c
        jbox=boxes(4+nsons,jbox)
        istack=istack+1
c
        nsons=0
        do 4600 j=5,8
c
        if(boxes(j,jbox) .gt. 0) nsons=nsons+1
 4600 continue
c
        stack(1,istack)=istack
        stack(2,istack)=jbox
        stack(3,istack)=nsons
c
 5000 continue
c
        return
        end
c
c
c
c
c
        subroutine d2tifint(c1,c2,ifinter)
        implicit real *8 (a-h,o-z)
        real *8 c1(2,4),c2(2,4)
ccc        save
c
c        this subroutine determines if two boxes in the square 
c        intersect or touch.
c
c                input parameters:
c
c  c1 - the four corners of the first box
c  c2 - the four corners of the second box
c
c                output parametes:
c 
c  ifinter - the indicator.
c      ifinter=1 means that the boxes intersect
c      ifinter=0 means that the boxes do not intersect
c
c       . . . find the maximum and minimum coordinates
c             for both boxes
c
        xmin1=c1(1,1)
        ymin1=c1(2,1)
c
        xmax1=c1(1,1)
        ymax1=c1(2,1)
c
        xmin2=c2(1,1)
        ymin2=c2(2,1)
c
        xmax2=c2(1,1)
        ymax2=c2(2,1)
c
c        xmin1=1.0d50
c        ymin1=1.0d50
c
c        xmax1=-1.0d50
c        ymax1=-1.0d50
c
c        xmin2=1.0d50
c        ymin2=1.0d50
c
c        xmax2=-1.0d50
c        ymax2=-1.0d50
c
        do 1200 i=1,4
c
        if(xmin1 .gt. c1(1,i)) xmin1=c1(1,i)
        if(ymin1 .gt. c1(2,i)) ymin1=c1(2,i)
c
        if(xmax1 .lt. c1(1,i)) xmax1=c1(1,i)
        if(ymax1 .lt. c1(2,i)) ymax1=c1(2,i)
c
c
        if(xmin2 .gt. c2(1,i)) xmin2=c2(1,i)
        if(ymin2 .gt. c2(2,i)) ymin2=c2(2,i)
c
        if(xmax2 .lt. c2(1,i)) xmax2=c2(1,i)
        if(ymax2 .lt. c2(2,i)) ymax2=c2(2,i)
c
 1200 continue
c        
c        decide if the boxes intersect
c
        eps=xmax1-xmin1
        if(eps .gt. xmax2-xmin2) eps=xmax2-xmin2
        if(eps .gt. ymax2-ymin2) eps=ymax2-ymin2
c
        if(eps .gt. ymax1-ymin1) eps=ymax1-ymin1
c
        eps=eps/10000
c
        ifinter=1
        if(xmin1 .gt. xmax2+eps) ifinter=0
        if(xmin2 .gt. xmax1+eps) ifinter=0
c
        if(ymin1 .gt. ymax2+eps) ifinter=0
        if(ymin2 .gt. ymax1+eps) ifinter=0
c
        return
        end
c
c
c
c
c
        subroutine d2tifint2(box1,box2,ifinter)
        implicit real *8 (a-h,o-z)
        integer box1(15),box2(15)
c
c        this subroutine determines if two boxes in the square 
c        intersect or touch.
c
c                input parameters:
c
c  box1 - the integer array describing the first box
c  box1 - the integer array describing the second box
c
c       integer arrays are dimensioned (15), as produced by d2tallb
c
c                output parametes:
c 
c  ifinter - the indicator.
c      ifinter=1 means that the boxes intersect
c      ifinter=0 means that the boxes do not intersect
c
c
        ifinter=1
c
        do i=1,2
        level1=box1(1)
        level2=box2(1)
        ip1=box1(i+1)-1
        ip2=box2(i+1)-1
	if( (ip1+1)*2**(level2-level1) .lt. (ip2  ) ) ifinter=0
	if( (ip1  )*2**(level2-level1) .gt. (ip2+1) ) ifinter=0
        if( ifinter .eq. 0 ) return
        enddo
c
        return
        end
c
c
c
c
c
        subroutine d2tcpcc(centers,corners,ibox,center,corner)
        implicit real *8 (a-h,o-z)
        real *8 centers(2,1),corners(2,4,1),center(2),corner(2,4)
ccc        save
c
        center(1)=centers(1,ibox)        
        center(2)=centers(2,ibox)        
c
        do 1200 i=1,4
        corner(1,i)=corners(1,i,ibox)        
        corner(2,i)=corners(2,i,ibox)        
c
 1200 continue
c
        return
        end

c
c
c
c
c
        subroutine d2tallbem(ier,z,n,nbox,boxes,maxboxes,
     1    nboxes,iz,laddr,nlev,center0,size,iwork,
     $    ifempty,minlevel,maxlevel,ztarg,ntarg,iztarg)
        implicit real *8 (a-h,o-z)
        integer boxes(15,*),iz(*),iztarg(*),laddr(2,*),iwork(*),
     1      iisons(4),jjsons(4)
        integer, allocatable :: is(:,:), ns(:,:)
        integer, allocatable :: istarg(:,:), nstarg(:,:)
        integer, allocatable :: son_idx(:)        
        real *8 z(2,1),ztarg(2,1),center0(2),center(2)
        data jjsons/1,1,2,2/,
     1      iisons/1,2,1,2/
ccc        save
c
c        this subroutine constructs a quad-tree corresponding
c        to the user-specified collection of points in the plane
c
c              input parameters:
c
c  z - the set of points in the plane
c  n - the number of elements in z
c
c  ztarg - the user-specified targets in the space
c  ntarg - the number of targets in array ztarg
c
c  nbox - the maximum number of points permitted in a box on 
c        the finest level. in other words, a box will be further
c        subdivided if it contains more than nbox points.
c  maxboxes - the maximum total number of boxes the subroutine 
c        is permitted to create. if the points z are such that 
c        more boxes are needed, the error return code ier is
c        set to 4, and the execution of the subroutine is
c        terminated.
c  
c              output parameters:
c
c  ier - the error return code.
c    ier=0 means successful execution
c    ier=4 means that the subroutine attempted to create more 
c        than maxboxes boxes
c    ier=16 means that the subroutine attempted to construct more 
c        than 197 levels of subdivision.
c  boxes - an integer array dimensioned (15,nboxes). each 15-element
c        column describes one box, as follows:
c
c       1. level - the level of subdivision on which this box 
c             was constructed; 
c       2, 3  - the coordinates of this box among  all
c             boxes on this level
c       4 - the daddy of this box, identified by it address
c             in array boxes
c       5,6,7,8 - the  list of children of this box 
c             (eight of them, and the child is identified by its address
c             in the array boxes; if a box has only one child, only the
c             first of the four child entries is non-zero, etc.)
c       9 - the location in the array iz of the particles 
c             living in this box
c       10 - the number of particles living in this box
c       11 - the location in the array iztarg of the targets
c             living in this box
c       12 - the number of targets living in this box
c       13 - source box type: 0 - empty, 1 - leaf node, 2 - sub-divided
c       14 - target box type: 0 - empty, 1 - leaf node, 2 - sub-divided
c       15 - reserved for future use
c
c    important warning: the array boxes has to be dimensioned 
c                       at least (15,maxboxes)!! otherwise, 
c                       the subroutine is likely to bomb, since
c                       it assumes that it has that much space!!!!
c  nboxes - the total number of boxes created
c  iz - the integer array addressing the particles in 
c         all boxes. 
c       explanation: for a box ibox, the particles living in
c         it are:
c
c         (z(1,j),z(2,j)),(z(1,j+1),z(2,j+1)),
c         (z(1,j+2),z(2,j+2)), . . . 
c         (z(1,j+nj-1),z(2,j+nj-1)),
c         (z(1,j+nj),z(2,j+nj)),
c
c         with j=boxes(9,ibox), and nj=boxes(10,ibox)
c        
c  iztarg - the integer array addressing the targets in 
c         all boxes. 
c
c         (ztarg(1,j),ztag(2,j)),
c         (ztarg(1,j+1),ztarg(2,j+1)),
c         (ztarg(1,j+2),ztarg(2,j+2)), . . . 
c         (ztarg(1,j+nj-1),ztarg(2,j+nj-1)),
c         (ztarg(1,j+nj),ztarg(2,j+nj)),
c
c         with j=boxes(11,ibox), and nj=boxes(12,ibox)
c
c  laddr - an integer array dimensioned (2,numlev), containing
c         the map of array boxes, as follows:
c       laddr(1,i) is the location in array boxes of the information
c         pertaining to level=i-1
c       laddr(2,i) is the number of boxes created on the level i-1
c
c  nlev - the maximum level number on which any boxes have 
c         been created. the maximim number possible is 200. 
c         it is recommended that the array laddr above be 
c         dimensioned at least (2,200), in case the user underestimates
c         the number of levels required.
c  center0 - the center of the box on the level 0, containing
c         the whole simulation
c  size - the side of the box on the level 0
c
c                  work arrays:
c
c  iwork - must be at least n+2 integer*4 elements long.
c
c
c      ... allocate work arrays
c       
        allocate( is(4,maxboxes),ns(4,maxboxes),
     $      istarg(4,maxboxes),nstarg(4,maxboxes),son_idx(maxboxes))

        do i=1,maxboxes
        ns(1,i)=0
        ns(2,i)=0
        ns(3,i)=0
        ns(4,i)=0
        nstarg(1,i)=0
        nstarg(2,i)=0
        nstarg(3,i)=0
        nstarg(4,i)=0
        enddo

c        . . . find the main box containing the whole picture
c
cccc          call prinf('in d2tallb, maxboxes=*',maxboxes,1)
c
        ier=0
        xmin=z(1,1)
        xmax=z(1,1)
        ymin=z(2,1)
        ymax=z(2,1)
c
c        xmin=1.0d50
c        xmax=-xmin
c        ymin=1.0d50
c        ymax=-ymin
c
        do 1100 i=1,n
        if(z(1,i) .lt. xmin) xmin=z(1,i)
        if(z(1,i) .gt. xmax) xmax=z(1,i)
        if(z(2,i) .lt. ymin) ymin=z(2,i)
        if(z(2,i) .gt. ymax) ymax=z(2,i)
 1100 continue
        do 1150 i=1,ntarg
        if(ztarg(1,i) .lt. xmin) xmin=ztarg(1,i)
        if(ztarg(1,i) .gt. xmax) xmax=ztarg(1,i)
        if(ztarg(2,i) .lt. ymin) ymin=ztarg(2,i)
        if(ztarg(2,i) .gt. ymax) ymax=ztarg(2,i)
 1150 continue
        size=xmax-xmin
        sizey=ymax-ymin
        if(sizey .gt. size) size=sizey
c
        center0(1)=(xmin+xmax)/2
        center0(2)=(ymin+ymax)/2
c
ccc         call prin2('in d2tallb, center0=*',center0,3)
ccc         call prin2('in d2tallb, size=*',size,1)
c
        boxes(1,1)=0
        boxes(2,1)=1
        boxes(3,1)=1
        boxes(4,1)=0
        boxes(5,1)=0     
        boxes(6,1)=0     
        boxes(7,1)=0     
        boxes(8,1)=0     
        boxes(9,1)=1
        boxes(10,1)=n
        boxes(11,1)=1
        boxes(12,1)=ntarg
        if( n .le. 0 ) boxes(13,1)=0
        if( ntarg .le. 0 ) boxes(14,1)=0
        if( n .gt. 0 ) boxes(13,1)=1
        if( ntarg .gt. 0 ) boxes(14,1)=1
        boxes(15,1)=0
c
        laddr(1,1)=1
        laddr(2,1)=1
c
        do 1200 i=1,n 
        iz(i)=i
 1200 continue
c
        do 1250 i=1,ntarg
        iztarg(i)=i
 1250 continue
c
c       recursively (one level after another) subdivide all 
c       boxes till none are left with more than nbox particles
c
        maxson=maxboxes
c
ccc         call prinf('in d2tallb, maxson=*',maxson,1) 
c
        maxlev=198
        if( maxlevel .le. maxlev ) maxlev=maxlevel

c
        ison=1
        nlev=0
cccc         call prinf('in d2tallb, nbox=*',nbox,1) 
cccc         call prinf('in d2tallb, n=*',n,1) 
        do 3000 level=0,maxlev-1
cccc          call prinf('in d2tallb, level=*',level,1) 
        laddr(1,level+2)=laddr(1,level+1)+laddr(2,level+1)
        nlevson=0
        idad0=laddr(1,level+1)
        idad1=idad0+laddr(2,level+1)-1
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(idad,numpdad,numtdad)
C$OMP$PRIVATE(ii,jj,iiz,nz,iiztarg,nztarg,center)
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1)
        do 2000 idad=idad0,idad1
c
c       subdivide the box number idad (if needed)
c
        numpdad=boxes(10,idad)
        numtdad=boxes(12,idad)
c
c       ... refine on sources only
ccc        if(numpdad .le. nbox .and. level .ge. minlevel ) goto 2000
c
c       ... not a leaf node on sources
ccc        if( numpdad .gt. nbox ) then
ccc        if( boxes(13,idad) .eq. 1 ) boxes(13,idad)=2
ccc        if( boxes(14,idad) .eq. 1 ) boxes(14,idad)=2
ccc        endif
c
c       ... refine on targets only
cccc        if(numtdad .le. nbox .and. level .ge. minlevel ) goto 2000
c
cccc        if( numtdad .gt. nbox ) then
cccc        if( boxes(13,idad) .eq. 1 ) boxes(13,idad)=2
cccc        if( boxes(14,idad) .eq. 1 ) boxes(14,idad)=2
cccc        endif
c
c       ... refine on both sources and targets
        if(numpdad .le. nbox .and. numtdad .le. nbox .and.
     $     level .ge. minlevel ) goto 2000
c
c       
c       ... not a leaf node on sources or targets
        if( numpdad .gt. nbox .or. numtdad .gt. nbox ) then
        if( boxes(13,idad) .eq. 1 ) boxes(13,idad)=2
        if( boxes(14,idad) .eq. 1 ) boxes(14,idad)=2
        endif
c
c
        ii=boxes(2,idad)
        jj=boxes(3,idad)
        call d2tcentf(center0,size,level,ii,jj,center)
c
        iiz=boxes(9,idad)
        nz=boxes(10,idad)
c
        iiztarg=boxes(11,idad)
        nztarg=boxes(12,idad)
c
cccc        call prinf('before d2tsepa1, nz=*',nz,1)
c
        call d2tsepa1(center,z(1,iiz),iz(iiz),nz,iwork(iiz),
     1    is(1,idad),ns(1,idad))
c
        call d2tsepa1(center,ztarg(1,iiztarg),iztarg(iiztarg),
     $     nztarg,iwork(n+iiztarg),istarg(1,idad),nstarg(1,idad))
c
cccc        call prinf('after d2tsepa1, is=*',is,4)
cccc        call prinf('after d2tsepa1, ns=*',ns,4)
c
 2000   continue
C$OMP END PARALLEL DO
c
        ison_cnt=ison
        do 2010 idad=idad0,idad1 
c
        son_idx(idad)=0
        do 1610 i=1,4
        if(ns(i,idad) .eq. 0 .and. nstarg(i,idad) .eq. 0 .and. 
     $     ifempty .ne. 1) goto 1610
        if( son_idx(idad) .eq. 0 ) son_idx(idad)=ison_cnt
        nlevson=nlevson+1
        nlev=level+1
        ison_cnt=ison_cnt+1
c
c       . . . if the user-provided array boxes is too
c             short - bomb out
c
        if(ison_cnt .le. maxson) goto 1410
        ier=4
        return
 1410 continue
c
 1610   continue
 2010   continue
c
ccc        write(*,*) level, ison, ison_cnt
c
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(idad)
C$OMP$PRIVATE(ii,jj,iiz,nz,iiztarg,nztarg)
C$OMP$PRIVATE(idadson,ison,i,lll,iison,jjson)
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1)
        do 2020 idad=idad0,idad1 
c
        ii=boxes(2,idad)
        jj=boxes(3,idad)
c
        iiz=boxes(9,idad)
        nz=boxes(10,idad)
c
        iiztarg=boxes(11,idad)
        nztarg=boxes(12,idad)
c
c       store in array boxes the sons obtained by the routine 
c       d2tsepa1
c
        if(son_idx(idad) .eq. 0 ) goto 2020
        ison=son_idx(idad)
         idadson=5
        do 1600 i=1,4
        if(ns(i,idad) .eq. 0 .and. nstarg(i,idad) .eq. 0 .and. 
     $     ifempty .ne. 1) goto 1600
        ison=ison+1
c
c        store in array boxes all information about this son
c
        do 1500 lll=5,8
        boxes(lll,ison)=0
 1500 continue
c
        boxes(1,ison)=level+1
        iison=(ii-1)*2+iisons(i)
        jjson=(jj-1)*2+jjsons(i)
        boxes(2,ison)=iison
        boxes(3,ison)=jjson
        boxes(4,ison)=idad        
c
        boxes(9,ison)=is(i,idad)+iiz-1
        boxes(10,ison)=ns(i,idad)
c
        boxes(11,ison)=istarg(i,idad)+iiztarg-1
        boxes(12,ison)=nstarg(i,idad)
c
        if( ns(i,idad) .le. 0 ) boxes(13,ison)=0
        if( nstarg(i,idad) .le. 0 ) boxes(14,ison)=0
        if( ns(i,idad) .gt. 0 ) boxes(13,ison)=1
        if( nstarg(i,idad) .gt. 0 ) boxes(14,ison)=1
        boxes(15,ison)=0
c
        boxes(idadson,idad)=ison
        idadson=idadson+1
 1600 continue
 2020 continue
C$OMP END PARALLEL DO

        ison=ison_cnt
        nboxes=ison
        laddr(2,level+2)=nlevson
         if(nlevson .eq. 0) goto 4000
         level1=level
 3000 continue
        if( level1 .ge. 197 ) ier=16
 4000 continue
        nboxes=ison
        return
        end
c
c
c
c
c
        subroutine d2tcentc(center0,size,boxes,nboxes,
     1      centers,corners)
        implicit real *8 (a-h,o-z)
        integer boxes(15,*)
        real *8 centers(2,*),corners(2,4,*),center(2),center0(2)
ccc        save
c
c       this subroutine produces arrays of centers and 
c       corners for all boxes in the quad-tree structure.
c
c              input parameters:
c
c  center0 - the center of the box on the level 0, containing
c         the whole simulation
c  size - the side of the box on the level 0
c  boxes - an integer array dimensioned (15,nboxes), as produced 
c        by the subroutine d2tallb (see)
c        column describes one box, as follows:
c  nboxes - the total number of boxes created
c
c  
c              output parameters:
c
c  centers - the centers of all boxes in the array boxes
c  corners - the corners of all boxes in the array boxes
c
c       . . . construct the corners for all boxes
c
        x00=center0(1)-size/2
        y00=center0(2)-size/2
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP$PRIVATE(i,level,side,side2,ii,jj,center)   
        do 1400 i=1,nboxes
        level=boxes(1,i)
        side=size/2**level
        side2=side/2
        ii=boxes(2,i)
        jj=boxes(3,i)
        center(1)=x00+(ii-1)*side+side2
        center(2)=y00+(jj-1)*side+side2        
c
        centers(1,i)=center(1)
        centers(2,i)=center(2)
c
        corners(1,1,i)=center(1)-side/2
        corners(1,2,i)=center(1)+side/2
        corners(1,3,i)=center(1)+side/2
        corners(1,4,i)=center(1)-side/2
c
        corners(2,1,i)=center(2)-side/2
        corners(2,2,i)=center(2)-side/2
        corners(2,3,i)=center(2)+side/2
        corners(2,4,i)=center(2)+side/2
 1400 continue
C$OMP END PARALLEL DO
         return
         end
c
c
c
c
c
        subroutine d2tcentf(center0,size,level,i,j,center) 
        implicit real *8 (a-h,o-z)
        real *8 center(2),center0(2)
ccc        data level0/-1/
ccc        save
c
c       this subroutine finds the center of the box 
c       number (i,j) on the level level. note that the 
c       box on level 0 is assumed to have the center 
c       center0, and the side size
c
ccc        if(level .eq. level0) goto 1200
        side=size/2**level
        side2=side/2
        x0=center0(1)-size/2
        y0=center0(2)-size/2
        level0=level
 1200 continue
        center(1)=x0+(i-1)*side+side2
        center(2)=y0+(j-1)*side+side2
        return
        end

c
c
c
c
c
        subroutine d2tsepa1(cent,z,iz,n,iwork,
     1    is,ns)
        implicit real *8 (a-h,o-z)
        real *8 cent(2),z(2,*)
        integer iz(*),iwork(*),is(*),ns(*)
ccc        save
c
c        this subroutine reorders the particles in a box,
c        so that each of the children occupies a contigious 
c        chunk of array iz
c
c        note that we are using a standard numbering convention 
c        for the children:
c
c
c             3,4   
c                                   
c             1,2
c
c 
c                        input parameters:
c
c  cent - the center of the box to be subdivided
c  z - the list of all points in the box to be subdivided
c  iz - the integer array specifying the transposition already
c       applied to the points z, before the subdivision of 
c       the box into children
c  n - the total number of points in array z
c  
c                        output parameters:
c
c  iz - the integer array specifying the transposition already
c       applied to the points z, after the subdivision of 
c       the box into children
c  is - an integer array of length 8 containing the locations 
c       of the sons in array iz
c  ns - an integer array of length 8 containig the numbers of 
c       elements in the sons
c
c                        work arrays:
c
c  iwork - must be n integer elements long
c
c        . . . separate all particles in this box in x
c
        n1=0
        n2=0
        n3=0
        n4=0
c
        n12=0
        n34=0
c
        itype=2
        thresh=cent(2)
        call d2tsepa0ipt(z,iz,n,itype,thresh,iwork,n12)
        n34=n-n12
c
c       at this point, the contents of sons number 1,2 are in
c       the part of array iz with numbers 1,2,...n12
c       the contents of sons number 3,4  are in
c       the part of array iz with numbers n34+1,n34+2,...,n
c
c        . . . separate the boxes 1, 2 and boxes 3, 4
c
        itype=1
        thresh=cent(1)
        if(n12 .ne. 0) 
     1    call d2tsepa0ipt(z,iz,n12,itype,thresh,iwork,n1)
        n2=n12-n1
c
        if(n34 .ne. 0) 
     1    call d2tsepa0ipt
     $     (z(1,n12+1),iz(n12+1),n34,itype,thresh,iwork(n12+1),n3)
        n4=n34-n3
c
c
c       store the information about the sonnies in appropriate arrays
c
        is(1)=1
        ns(1)=n1
c
        is(2)=is(1)+ns(1)
        ns(2)=n2
c
        is(3)=is(2)+ns(2)
        ns(3)=n3
c
        is(4)=is(3)+ns(3)
        ns(4)=n4
c
c
cccc        call prinf('is as created*',is,4)
cccc        call prinf('ns as created*',ns,4)
        return
        end
c
c
c
c
c
        subroutine d2tsepa0(z,iz,n,itype,thresh,iwork,n1)
        implicit real *8 (a-h,o-z)
        real *8 z(2,*)
        integer iz(*),iwork(*)
ccc        save
c
c       subdivide the points in this box, horizontally
c       or vertically, depending on the parameter itype
c
cccc        call prin2('in d2tsepa0,thresh=*',thresh,1)
c
        i1=0
        i2=0
        do 1400 i=1,n
        j=iz(i)
        if(z(itype,j) .gt. thresh) goto 1200
        i1=i1+1
        iz(i1)=j
        goto 1400
c
 1200 continue
        i2=i2+1
        iwork(i2)=j
 1400 continue
c
        do 1600 i=1,i2
        iz(i1+i)=iwork(i)
 1600 continue
        n1=i1
        return
        end
c
c
c
c
c
        subroutine d2tsepa0ip(z,iz,n,itype,thresh,iwork,n1)
        implicit real *8 (a-h,o-z)
        real *8 z(2,*)
        integer iz(*),iwork(*)
ccc        save
c
c       subdivide the points in this box, horizontally
c       or vertically, depending on the parameter itype
c
c       in place sorting algorithm, slightly slower
c
cccc        call prin2('in d2tsepa0,thresh=*',thresh,1)
c
        i1=1
        i2=n

        do while ( i1 .le. i2 )

        do while ( (z(itype,iz(i1)) .le. thresh) ) 
        i1=i1+1
        if( i1 .gt. i2 ) exit
        enddo

        do while ( (z(itype,iz(i2)) .gt. thresh) ) 
        i2=i2-1
        if( i1 .gt. i2 ) exit
        enddo

        if( i1 .lt. i2 ) then
        i3=iz(i2)
        iz(i2)=iz(i1)
        iz(i1)=i3
        i1=i1+1
        endif

        enddo

        n1=i1-1

c        do i=1,n1
c        if( z(itype,iz(i)) .gt. thresh ) then
c        write(*,*) i1,i2,n1
c        write(*,*) i,n, '+', z(itype,iz(i))
c        pause
c        endif
c        enddo

c        do i=n1+1,n
c        if( z(itype,iz(i)) .le. thresh ) then
c        write(*,*) i1,i2,n1
c        write(*,*) i,n, '-', z(itype,iz(i))
c        pause
c        endif
c        enddo

        return
        end
c
c
c
c
c
        subroutine d2tsepa0ipt(z,iz,n,itype,thresh,iwork,n1)
        implicit real *8 (a-h,o-z)
        real *8 z(2,*)
        integer iz(*),iwork(*)
ccc        save
c
c       subdivide the points in this box, horizontally
c       or vertically, depending on the parameter itype
c
c       in place sorting algorithm, sort particles into boxes directly
c       without indirect addressing, this is much faster due to better
c       data access
c
cccc        call prin2('in d2tsepa0,thresh=*',thresh,1)
c
        i1=1
        i2=n

        do while ( i1 .le. i2 )

        do while ( (z(itype,i1) .le. thresh) ) 
        i1=i1+1
        if( i1 .gt. i2 ) exit
        enddo

        do while ( (z(itype,i2) .gt. thresh) ) 
        i2=i2-1
        if( i1 .gt. i2 ) exit
        enddo

        if( i1 .lt. i2 ) then
        i3=iz(i2)
        iz(i2)=iz(i1)
        iz(i1)=i3
        d=z(1,i2)
        z(1,i2)=z(1,i1)
        z(1,i1)=d
        d=z(2,i2)
        z(2,i2)=z(2,i1)
        z(2,i1)=d
        i1=i1+1
        endif

        enddo

        n1=i1-1

c        do i=1,n1
c        if( z(itype,iz(i)) .gt. thresh ) then
c        write(*,*) i1,i2,n1
c        write(*,*) i,n, '+', z(itype,iz(i))
c        pause
c        endif
c        enddo

c        do i=n1+1,n
c        if( z(itype,iz(i)) .le. thresh ) then
c        write(*,*) i1,i2,n1
c        write(*,*) i,n, '-', z(itype,iz(i))
c        pause
c        endif
c        enddo

        return
        end
c
c
c
c
c
        subroutine d2tlinkinit(ier,nboxes,ntypes,w,lw)
        integer w(*),nums(*),inums(20)
        data iilistad/1/,iilists/2/,inumele/3/,inboxes/4/,
     1      intypes/5/,ilw/6/,iltot/7/,
     2      inums/11,12,13,14,15,16,17,18,19,20,21,22,23,24,
     3            25,26,27,28,29,30/
ccc        save
c
c        this is the initialization entry point for the link-list
c        storage-retrieval facility. it formats the array w, to
c        be used by the entries d2tlinkstor, d2tlinkretr, d2tlinkrem,
c        d2tlinkinfo below. 
c
c                     input parameters:
c
c  nboxes - the number of boxes for which the various types of lists
c        will be stored
c  ntypes - the number of types of lists that will be stored
c  lw - the total amount of space in the area w to be used for storage
c        (in integer locations)
c
c                     output parameters:
c
c  ier - error return code;
c    ier=0 means successful execution
c    ier=1024 means that the amount of space in array w is grossly 
c        insufficient
c  w - the formatted area for storage
c
c
c       . . . allocate memory for the storage facility
c
        ier=0
c
        ilistadd=32
        nlistadd=nboxes*ntypes+10
c
        ilists=ilistadd+nlistadd
        numele=0
        ltot=ilists+numele*3+10
c
         if(ltot+100 .lt. lw) goto 1200
         ier=1024
         return
 1200 continue
c
         do 1400 i=1,20
         w(inums(i))=0
 1400 continue
c
        w(iilistad)=ilistadd
        w(iilists)=ilists
        w(inumele)=numele
        w(inboxes)=nboxes
        w(intypes)=ntypes      
        w(ilw)=lw
        w(iltot)=ltot
c
        call d2tlinkini0(w(ilistadd),nboxes,ntypes)
        return
c
c
c
c
        entry d2tlinkstor(ier,itype,ibox,list,nlist,w,lused)
c
c       this entry stores dynamically a list of positive numbers 
c       in the storage array w.
c
c                      input parameters:
c
c  itype - the type of the elements being stored
c  ibox - the box to which these elements corresponds
c  list - the list of positive integer elements to be stored
c  nlist - the number of elements in the array list
c  w - the storage area used by this subroutine; must be first 
c          formatted by the entry d2tlinkinit of this subroutine
c          (see above)
c
c                      output parameters:
c
c  ier - error return code; 
c         ier=0 means successful execution
c         ier=32 means that the storage area w would be exceeded
c                 by this storage request
c  w - the storage area used by this subroutine
c  lused - the number of integer elements used in array
c          w after this call.
c
c       . . . if this storage request exceeds the available memory - bomb
c
        ier=0
        if(w(iilists)+w(inumele)*3+nlist*3 .lt. w(ilw) ) goto 2200
        ier=32
        return
 2200 continue
c
c       store the user-specified list in array w
c
c       Old versions of gfortran <=4.3 will naively optimize w(inumele)
c       away if w(inumele)=0, at optimization levels greater than 1.
c       We use a temporary variable to fix this nasty compiler bug.
c
        itemp = w(inumele)
        call d2tlinksto0(itype,ibox,list,nlist,w(w(iilistad)),
     1      w(inboxes),w(w(iilists)),itemp,w(inums(itype)) )
        w(inumele)=itemp
c
ccc        call d2tlinksto0(itype,ibox,list,nlist,w(w(iilistad)),
ccc     1      w(inboxes),w(w(iilists)),w(inumele),w(inums(itype)) )
c
c       augment the amount of storage used 
c        
        lused=w(iilists)+w(inumele)*3+10
        return
c
c
c
c
        entry d2tlinkretr(ier,itype,ibox,list,nlist,w,lused)
c
c       this entry retrieves from the storage area  w 
c       a list of positive numbers that has been stored there
c       by the entry d2tlinkstor (see above).
c
c                      input parameters:
c
c  itype - the type of the elements to be retrieved
c  ibox - the box to which these elements correspond
c  w - the storage area from which the information is to be 
c          retrieved
c
c                      output parameters:
c
c  ier - error return code;
c          ier=0 means successful execution, nothing to report
c          ier=4 means that no elements are present of the type
c                  itype and the box ibox
c  list - the list of positive integer elements retrieved 
c  nlist - the number of elements in the array list
c  lused - the number of integer elements used in array
c          w after this call.
c
        call d2tlinkret0(ier,itype,ibox,w(w(iilistad)),
     1      w(w(iilists)),list,w(inboxes),nlist)
c
        lused=w(iilists)+w(inumele)*3+10
c
        return
c
c
c
c
        entry d2tlinkrem(ier,itype,ibox,list,nlist,w,lused)
c
c       this entry deletes elements in array lists corresponding 
c       the user-specified list  list. actually, it does not 
c       delete anything, but rather marks these elements for
c       future destruction by making them negative.
c
c                      input parameters:
c
c  itype - the type of the elements to be destroyed
c  ibox - the box to which these elements correspond
c  list - the list of positive integer elements to be destroyed
c  nlist - the number of elements in the array list
c  w - the storage area from which the information is to be 
c          retrieved
c
c                      output parameters:
c
c  ier - error return code;
c          ier=0 means successful execution, nothing to report
c          ier=4*k means that k of the elements on the user-specified
c                   list  list were not present
c          ier=22 means that  no elements whatsoever were present 
c                   for the type  itype   and the box  ibox
c  w - the storage area from which the information is to be 
c          retrieved
c  lused - the number of integer elements used in array
c          w both before and after this call.
c
c       mark for destruction the user-specified elements 
c
        call d2tlinkrem0(ier,itype,ibox,list,nlist,w(w(iilistad)),
     1      w(inboxes),w(w(iilists)),w(inums(itype)) )
c
        lused=w(iilists)+w(inumele)*3+10
        return
c
c
c
c
        entry d2tlinkinfo(w,lused,nums)
c
c       this entry returns to the user some of the information
c       about the storage area. namely, it returns the total 
c       amount lused of memory utilized in the array w (in integer 
c       locations), and the integer array nums containing the numbers
c       of elements in each of the lists
c
        lused=w(iilists)+w(inumele)*3+10
        ntypes7=w(intypes)
         call prinf('in d2tlinkinfo, lused=*',lused,1)
         call prinf('in d2tlinkinfo, ntypes7=*',ntypes7,1)
         call prinf('in w(inumele)=*',w(inumele),1)
        do 6200 i=1,ntypes7
        nums(i)=w(inums(i))
 6200 continue
        return
        end
c
c
c
c
c
        subroutine d2tlinksto0(itype,ibox,list,nlist,listaddr,
     1      nboxes,lists,numele,numtype)
        integer listaddr(nboxes,*),lists(3,*),list(*)
ccc        save
c
c       this entry stores dynamically a list of positive numbers 
c       in the storage array lists, while entering the information
c       about this event in the array listaddr.
c
c                      input parameters:
c
c  itype - the type of the elements being stored
c  ibox - the box to which these elements corresponds
c  list - the list of positive integer elements to be stored
c  nlist - the number of elements in the array list
c  listaddr - the addressing array for the main storage array lists;
c             it is assumed that it has been formatted by a call 
c             to the entry d2tlinkini0 of this subroutine (see below).
c  nboxes - the total number of boxes indexing the elements 
c           in array lists
c  lists - the main storage area used by this subroutine
c  numele - the number of elements stored in array lists on entry 
c             to this subroutinec

c                      output parameters:
c
c  listaddr - the addressing array for the main storage array lists;
c             it is assumed that it has been formatted by a call 
c             to the entry d2tlinkini0 of this subroutine (see below).
c  lists - the main storage area used by this subroutine
c  numele - the number of elements stored in array lists on exit
c             from this subroutine
c  numtype - the total number of elements in all lists of the 
c           type itype AFTER this call
c
c .........................................................................
c
c        interpretation of the entries in arrays lists, listaddr:
c
c  lists(1,i) - the location in array lists of the preceding
c        element in the list with (box,type) combination as the 
c        user-supplied (ibox,itype)
c     lists(1,i) .leq. 0 means that this is the first element 
c        of its type,
c  lists(2,i) - the box number on the interaction list.
c
c
c  listaddr(ibox,itype) - the address in the array lists of the last element
c        of this type for this box;
c     listaddr(ibox,itype) .leq. 0 means that there are no elements 
c        in the list of this type for this box.
c
c       . . . store the user-supplied list elements in the array lists,
c             and enter information about this change in array listaddr
c
        ilast=listaddr(ibox,itype)
        do 1200 i=1,nlist
        numele=numele+1
        numtype=numtype+1
        lists(1,numele)=ilast
        lists(2,numele)=list(i)
c
        ilast=numele
 1200 continue
        listaddr(ibox,itype)=ilast
c
        return
c
c
c
c
        entry d2tlinkret0(ier,itype,ibox,listaddr,lists,list,
     1      nboxes,nlist)
c
c       this entry retrieves from the main storage array lists
c       a list of positive numbers that has been stored there
c       by the entry d2tlinksto0 (see above).
c
c                      input parameters:
c
c  itype - the type of the elements being to be retrieved
c  ibox - the box to which these element corresponds
c  listaddr - the addressing array for the main storage array lists;
c             it is assumed that it has been formatted by a call 
c             to the entry d2tlinkini0 of this subroutine (see below).
c  nboxes - the total number of boxes indexing the elements 
c           in array lists
c  lists - the main storage area used by this subroutine
c
c                      output parameters:
c
c  ier - error return code;
c          ier=0 means successful execution, nothing to report
c          ier=4 means that no elements are present of the type
c                  itype and the box ibox
c  list - the list of positive integer elements to be stored
c  nlist - the number of elements in the array list
c
c
c       . . . retrieve and store in array list the list of the 
c             type itype for the box ibox
c
        ier=0
        ilast=listaddr(ibox,itype)
        if(ilast .gt. 0) goto 2200
        nlist=0
        ier=4
        return
 2200 continue
c
        nlist=0
        do 2400 i=1,1 000 000 000
c
        if(lists(2,ilast) .le. 0) goto 2300       
        nlist=nlist+1
        list(nlist)=lists(2,ilast)
 2300 continue
        ilast=lists(1,ilast)
        if(ilast .le. 0) goto 2600
 2400 continue
 2600 continue
c
        if(nlist .gt. 0) goto 2650
        ier=4
        return
 2650 continue
c
c        flip the retrieved array
c
        if(nlist .eq. 1) return
        do 2700 i=1,nlist/2
        j=list(i)
        list(i)=list(nlist-i+1)
        list(nlist-i+1)=j
 2700 continue
        return
c
c
c
c
        entry d2tlinkini0(listaddr,nboxes,ntypes)
c
c       this subroutine initializes the array listaddr to be used 
c       later by other entries of this subroutine
c
c                        input parameters:
c
c  nboxes - the number of boxes for which the various types of lists
c        will be stored
c  ntypes - the number of types of lists that will be stored
c
c                        output parameters:
c
c  listaddr - the addressing array for the main storage array lists;
c             it is assumed that it has been formatted by a call 
c             to the entry d2tlinkini0 of this subroutine (see below).
c
c       . . . initialize the array listaddr
c
        do 3000 k=1,ntypes
        do 2800 i=1,nboxes
        listaddr(i,k)=-1
 2800 continue
 3000 continue
        return
c
c
c
c
        entry d2tlinkrem0(ier,itype,ibox,list,nlist,listaddr,
     1      nboxes,lists,numtype)
c
c       this entry deletes elements in array lists corresponding 
c       the user-specified list  list. actually, it does not 
c       delete anything, but rather marks these elements for
c       future destruction by making them negative.
c
c                      input parameters:
c
c  itype - the type of the elements to be destroyed
c  ibox - the box to which these elements correspond
c  list - the list of positive integer elements to be destroyed
c  nlist - the number of elements in the array list
c  listaddr - the addressing array for the main storage array lists;
c             it is assumed that it has been formatted by a call 
c             to the entry d2tlinkini0 of this subroutine (see below).
c  nboxes - the total number of boxes indexing the elements 
c           in array lists
c  lists - the main storage area used by this subroutine
c
c                      output parameters:
c
c  ier - error return code;
c          ier=0 means successful execution, nothing to report
c          ier=4*k means that k of the elements on the user-specified
c                   list  list were not present
c          ier=22 means that  no elements whatsoever were present 
c                   for the type  itype   and the box  ibox
c  listaddr - the addressing array for the main storage array lists;
c             it is assumed that it has been formatted by a call 
c             to the entry d2tlinkini0 of this subroutine (see below).
c  lists - the main storage area used by this subroutine
c
c       . . . mark for destruction the elements of type (itype,ibox)
c             that are in the list  list.
c
        ier=0
        do 4000 i=1,nlist 
        ilast=listaddr(ibox,itype)
        if(ilast. gt. 0) goto 3200
        ier=22
        return
 3200 continue
c
        iffound=0
        do 3600 j=1,1 000 000 000
c
        if(ilast .le. 0) goto 3800
        if(lists(2,ilast) .ne. list(i)) goto 3400
        lists(2,ilast)=-lists(2,ilast)
        numtype=numtype-1
        iffound=1
 3400 continue
        ilast=lists(1,ilast)
 3600 continue
 3800 continue
         if(iffound .eq. 0) ier=ier+4
c
 4000 continue
        return
        end
c
c
c
c
c
        subroutine d2tgetbbox(n,z,center,size,corners)
        implicit real *8 (a-h,o-z)
        real *8 z(2,*),center(2),corners(2,4)
c
c       this subroutine finds the center and 
c       corners for the top level box in the quad-tree structure
c
c
        xmin=z(1,1)
        xmax=z(1,1)
        ymin=z(2,1)
        ymax=z(2,1)
c
c        xmin=1.0d50
c        xmax=-xmin
c        ymin=1.0d50
c        ymax=-ymin
c
        do 1100 i=1,n
        if(z(1,i) .lt. xmin) xmin=z(1,i)
        if(z(1,i) .gt. xmax) xmax=z(1,i)
        if(z(2,i) .lt. ymin) ymin=z(2,i)
        if(z(2,i) .gt. ymax) ymax=z(2,i)
 1100 continue
        size=xmax-xmin
        sizey=ymax-ymin
        if(sizey .gt. size) size=sizey
c
c
        center(1)=(xmin+xmax)/2
        center(2)=(ymin+ymax)/2
c
c
        corners(1,1)=center(1)-size/2
        corners(1,2)=corners(1,1)+size
        corners(1,3)=corners(1,1)+size
        corners(1,4)=corners(1,1)
c
c
        corners(2,1)=center(2)-size/2
        corners(2,2)=corners(2,1)
        corners(2,3)=corners(2,1)+size
        corners(2,4)=corners(2,3)+size
c
        return
        end

