c
c     Test driver for FMM tree reuse functionality
c
c     Demonstrates building tree once and evaluating with different nd values
c
        implicit real *8 (a-h,o-z)
        real *8, allocatable :: source(:,:)
        real *8, allocatable :: target(:,:)
        complex *16, allocatable :: charge1(:,:)
        complex *16, allocatable :: charge4(:,:)
        complex *16, allocatable :: charge8(:,:)
        complex *16, allocatable :: dipstr1(:,:)
        complex *16, allocatable :: dipstr4(:,:)
        complex *16, allocatable :: dipstr8(:,:)
        complex *16, allocatable :: pot1(:,:)
        complex *16, allocatable :: pot4(:,:)
        complex *16, allocatable :: pot8(:,:)
        complex *16, allocatable :: pottarg1(:,:)
        complex *16, allocatable :: pottarg4(:,:)
        complex *16, allocatable :: pottarg8(:,:)
c
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c
        call prini(6,13)
c
c       Set problem parameters
c
        nsource = 1000
        ntarget = 1000
        iprec = 4
c
        print *, 'Tree reuse test: nsource =', nsource
        print *, '                 ntarget =', ntarget
        print *, '                 iprec   =', iprec
c
c       Allocate geometry arrays
c
        allocate(source(2,nsource))
        allocate(target(2,ntarget))
c
c       Generate random source locations on unit circle
c
        done = 1.0d0
        pi = 4*atan(done)
c
        call random_seed()
        do i = 1,nsource
          call random_number(phi)
          phi = phi*2*pi
          source(1,i) = cos(phi)
          source(2,i) = sin(phi)
        enddo
c
c       Generate random target locations offset from sources
c
        do i = 1,ntarget
          call random_number(phi)
          phi = phi*2*pi
          target(1,i) = cos(phi) + 3.0d0
          target(2,i) = sin(phi)
        enddo
c
c       Build tree ONCE (geometry only, independent of nd)
c
        print *
        print *, '=== Building FMM tree (one time cost) ==='
        t1 = second()
        call cfmm2dpart_tree_build(ier,iprec,
     $     nsource,source,ntarget,target,tree_handle)
        t2 = second()
        if (ier .ne. 0) then
          print *, 'ERROR: tree_build failed with ier =', ier
          stop
        endif
        print *, 'Tree built successfully, handle =', tree_handle
        print *, 'Tree build time:', t2-t1, 'sec'
c
c       Test 1: Evaluate with nd=1
c
        print *
        print *, '=== Test 1: Evaluation with nd=1 ==='
        nd = 1
        allocate(charge1(nd,nsource))
        allocate(dipstr1(nd,nsource))
        allocate(pot1(nd,nsource))
        allocate(pottarg1(nd,ntarget))
c
        do i = 1,nsource
          dipstr1(1,i) = dcmplx(dble(i), 0.0d0)
        enddo
c
        ifcharge = 0
        ifdipole = 1
        ifpot = 1
        ifgrad = 0
        ifhess = 0
        ifpottarg = 1
        ifgradtarg = 0
        ifhesstarg = 0
c
        t1 = second()
        call cfmm2dpart_tree_eval(ier,tree_handle,nd,
     $     source,ifcharge,charge1,ifdipole,dipstr1,
     $     ifpot,pot1,ifgrad,pot1,ifhess,pot1,
     $     target,ifpottarg,pottarg1,ifgradtarg,pottarg1,ifhesstarg,
     $     pottarg1)
        t2 = second()
        if (ier .ne. 0) then
          print *, 'ERROR: tree_eval failed with ier =', ier
          stop
        endif
        print *, 'Evaluation time (nd=1):', t2-t1, 'sec'
        print *, 'pot1(1) =', pot1(1,1)
c
c       Test 2: Evaluate with nd=4 (same tree!)
c
        print *
        print *, '=== Test 2: Evaluation with nd=4 (reusing tree) ==='
        nd = 4
        allocate(charge4(nd,nsource))
        allocate(dipstr4(nd,nsource))
        allocate(pot4(nd,nsource))
        allocate(pottarg4(nd,ntarget))
c
        do idim = 1,nd
          do i = 1,nsource
            dipstr4(idim,i) = dcmplx(dble(i+idim), 0.0d0)
          enddo
        enddo
c
        t1 = second()
        call cfmm2dpart_tree_eval(ier,tree_handle,nd,
     $     source,ifcharge,charge4,ifdipole,dipstr4,
     $     ifpot,pot4,ifgrad,pot4,ifhess,pot4,
     $     target,ifpottarg,pottarg4,ifgradtarg,pottarg4,ifhesstarg,
     $     pottarg4)
        t2 = second()
        if (ier .ne. 0) then
          print *, 'ERROR: tree_eval failed with ier =', ier
          stop
        endif
        print *, 'Evaluation time (nd=4):', t2-t1, 'sec'
        print *, 'pot4(1,1) =', pot4(1,1)
        print *, 'pot4(4,1) =', pot4(4,1)
c
c       Test 3: Evaluate with nd=8 (same tree!)
c
        print *
        print *, '=== Test 3: Evaluation with nd=8 (reusing tree) ==='
        nd = 8
        allocate(charge8(nd,nsource))
        allocate(dipstr8(nd,nsource))
        allocate(pot8(nd,nsource))
        allocate(pottarg8(nd,ntarget))
c
        do idim = 1,nd
          do i = 1,nsource
            dipstr8(idim,i) = dcmplx(dble(i+idim), 0.0d0)
          enddo
        enddo
c
        t1 = second()
        call cfmm2dpart_tree_eval(ier,tree_handle,nd,
     $     source,ifcharge,charge8,ifdipole,dipstr8,
     $     ifpot,pot8,ifgrad,pot8,ifhess,pot8,
     $     target,ifpottarg,pottarg8,ifgradtarg,pottarg8,ifhesstarg,
     $     pottarg8)
        t2 = second()
        if (ier .ne. 0) then
          print *, 'ERROR: tree_eval failed with ier =', ier
          stop
        endif
        print *, 'Evaluation time (nd=8):', t2-t1, 'sec'
        print *, 'pot8(1,1) =', pot8(1,1)
        print *, 'pot8(8,1) =', pot8(8,1)
c
c       Free tree
c
        call cfmm2dpart_tree_free(tree_handle)
        print *
        print *, '=== Tree freed successfully ==='
c
        print *
        print *, 'SUCCESS: Tree reuse test completed!'
c
        stop
        end
