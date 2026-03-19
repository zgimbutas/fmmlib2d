c
c     Test driver for Helmholtz FMM tree reuse functionality
c
c     Demonstrates building tree once and evaluating with different nd values
c
        implicit real *8 (a-h,o-z)
        real *8, allocatable :: source(:,:)
        real *8, allocatable :: target(:,:)
        real *8, allocatable :: dipvec(:,:)
        complex *16, allocatable :: charge1(:,:)
        complex *16, allocatable :: charge4(:,:)
        complex *16, allocatable :: dipstr1(:,:)
        complex *16, allocatable :: dipstr4(:,:)
        complex *16, allocatable :: pot1(:,:)
        complex *16, allocatable :: pot4(:,:)
        complex *16, allocatable :: grad1(:,:,:)
        complex *16, allocatable :: grad4(:,:,:)
        complex *16, allocatable :: pottarg1(:,:)
        complex *16, allocatable :: pottarg4(:,:)
c
        complex *16 zk
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c
        integer tree_handle
c
        call prini(6,13)
c
c       Set problem parameters
c
        nsource = 1000
        ntarget = 500
        iprec = 4
c
        done = 1.0d0
        pi = 4*atan(done)
c
        zk = dcmplx(32*2*pi, 0.0d0)
c
        print *, 'Helmholtz tree reuse test'
        print *, '  nsource =', nsource
        print *, '  ntarget =', ntarget
        print *, '  iprec   =', iprec
        call prin2('  zk      =*',zk,2)
c
c       Allocate geometry arrays
c
        allocate(source(2,nsource))
        allocate(target(2,ntarget))
        allocate(dipvec(2,nsource))
c
c       Generate random source locations on unit circle
c
        call random_seed()
        do i = 1,nsource
          call random_number(phi)
          phi = phi*2*pi
          source(1,i) = cos(phi)
          source(2,i) = sin(phi)
          dipvec(1,i) = cos(phi)
          dipvec(2,i) = sin(phi)
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
        call hfmm2dpart_tree_build(ier,iprec,zk,
     $     nsource,source,ntarget,target,tree_handle)
        t2 = second()
        if (ier .ne. 0) then
          print *, 'ERROR: tree_build failed with ier =', ier
          stop
        endif
        print *, 'Tree built successfully, handle =', tree_handle
        call prin2('Tree build time (sec)=*',t2-t1,1)
c
c       Test 1: Evaluate with nd=1
c
        print *
        print *, '=== Test 1: Evaluation with nd=1 ==='
        nd = 1
        allocate(charge1(nd,nsource))
        allocate(dipstr1(nd,nsource))
        allocate(pot1(nd,nsource))
        allocate(grad1(nd,2,nsource))
        allocate(pottarg1(nd,ntarget))
c
        do i = 1,nsource
          charge1(1,i) = dcmplx(dble(i), 0.0d0)
          dipstr1(1,i) = dcmplx(dble(i), 0.0d0)
        enddo
c
        ifcharge = 1
        ifdipole = 1
        ifpot = 1
        ifgrad = 1
        ifhess = 0
        ifpottarg = 1
        ifgradtarg = 0
        ifhesstarg = 0
c
        t1 = second()
        call hfmm2dpart_tree_eval(ier,tree_handle,nd,
     $     source,ifcharge,charge1,ifdipole,dipstr1,dipvec,
     $     ifpot,pot1,ifgrad,grad1,ifhess,grad1,
     $     target,ifpottarg,pottarg1,ifgradtarg,pottarg1,ifhesstarg,
     $     pottarg1)
        t2 = second()
        if (ier .ne. 0) then
          print *, 'ERROR: tree_eval failed with ier =', ier
          stop
        endif
        call prin2('Evaluation time nd=1 (sec)=*',t2-t1,1)
        call prin2('pot1(1,1) =*',pot1(1,1),2)
c
c       Test 2: Evaluate with nd=4 (same tree!)
c
        print *
        print *, '=== Test 2: Evaluation with nd=4 (reusing tree) ==='
        nd = 4
        allocate(charge4(nd,nsource))
        allocate(dipstr4(nd,nsource))
        allocate(pot4(nd,nsource))
        allocate(grad4(nd,2,nsource))
        allocate(pottarg4(nd,ntarget))
c
        do idim = 1,nd
          do i = 1,nsource
            charge4(idim,i) = dcmplx(dble(i+idim), 0.0d0)
            dipstr4(idim,i) = dcmplx(dble(i+idim), 0.0d0)
          enddo
        enddo
c
        t1 = second()
        call hfmm2dpart_tree_eval(ier,tree_handle,nd,
     $     source,ifcharge,charge4,ifdipole,dipstr4,dipvec,
     $     ifpot,pot4,ifgrad,grad4,ifhess,grad4,
     $     target,ifpottarg,pottarg4,ifgradtarg,pottarg4,ifhesstarg,
     $     pottarg4)
        t2 = second()
        if (ier .ne. 0) then
          print *, 'ERROR: tree_eval failed with ier =', ier
          stop
        endif
        call prin2('Evaluation time nd=4 (sec)=*',t2-t1,1)
        call prin2('pot4(1,1) =*',pot4(1,1),2)
        call prin2('pot4(4,1) =*',pot4(4,1),2)
c
c       Verify nd=4 consistency: set 4 identical densities,
c       all 4 outputs should match
c
        print *
        print *, '=== Consistency check: 4 identical densities ==='
        do idim = 1,nd
          do i = 1,nsource
            charge4(idim,i) = dcmplx(dble(i), 0.0d0)
            dipstr4(idim,i) = dcmplx(dble(i), 0.0d0)
          enddo
        enddo
c
        call hfmm2dpart_tree_eval(ier,tree_handle,nd,
     $     source,ifcharge,charge4,ifdipole,dipstr4,dipvec,
     $     ifpot,pot4,ifgrad,grad4,ifhess,grad4,
     $     target,ifpottarg,pottarg4,ifgradtarg,pottarg4,ifhesstarg,
     $     pottarg4)
c
c       Check that pot4(1,:) == pot4(2,:) == pot4(3,:) == pot4(4,:)
c
        d = 0
        a = 0
        do i = 1,nsource
          a = a + abs(pot4(1,i))**2
          do idim = 2,nd
            d = d + abs(pot4(idim,i) - pot4(1,i))**2
          enddo
        enddo
        rerr = sqrt(d) / sqrt(a)
        call prin2('nd-consistency error (should be ~0) =*',rerr,1)
        if (rerr .lt. 1.0d-13) then
          print *, 'PASS: nd-consistency check'
        else
          print *, 'FAIL: nd-consistency check, rerr =', rerr
        endif
c
c       Verify nd=1 vs nd=4 (first density)
c
        d = 0
        a = 0
        do i = 1,nsource
          d = d + abs(pot4(1,i) - pot1(1,i))**2
          a = a + abs(pot1(1,i))**2
        enddo
        rerr = sqrt(d) / sqrt(a)
        call prin2('nd=1 vs nd=4 error (should be ~0) =*',rerr,1)
        if (rerr .lt. 1.0d-13) then
          print *, 'PASS: tree reuse consistency (nd=1 vs nd=4)'
        else
          print *, 'FAIL: tree reuse consistency, rerr =', rerr
        endif
c
c       Free tree
c
        call hfmm2dpart_tree_free(tree_handle)
        print *
        print *, '=== Tree freed successfully ==='
c
        print *
        print *, 'SUCCESS: Helmholtz tree reuse test completed!'
c
        stop
        end
