cc SPDX-License-Identifier: BSD-3-Clause
cc Copyright (C) 2025: Zydrunas Gimbutas
cc See LICENSE file for full license text
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     FMM Tree Handle Manager
c
c     This module manages FMM tree structures that are built once
c     and can be reused for multiple evaluations with different nd
c     values and charge distributions.
c
c     The tree structure contains only geometry-dependent data:
c     - Quad-tree spatial decomposition (wlists)
c     - Particle sorting indices (iisource, itarget)
c     - Level-dependent parameters (scale, bsize, nterms)
c     - Box hierarchy (laddr, nboxes, nlev)
c
c     The tree is INDEPENDENT of nd (number of densities).
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     Tree structure definition
c
      module fmm2d_tree_module
      implicit none

      type fmm2d_tree_data
        logical :: allocated

        ! Geometry parameters (tree-dependent, NO nd)
        integer :: nsource
        integer :: ntarget
        integer :: iprec
        integer :: nboxes
        integer :: nlev

        ! Helmholtz parameter (for Helmholtz FMM only)
        complex*16 :: zk

        ! Tree structure arrays
        real*8, allocatable :: wlists(:)
        real*8 :: center(2)
        real*8 :: size
        integer, allocatable :: laddr(:,:)
        real*8, allocatable :: scale(:)
        real*8, allocatable :: bsize(:)
        integer, allocatable :: nterms(:)

        ! Workspace indices into wlists
        integer :: iisource
        integer :: iitarget
        integer :: iwlists
        integer :: lwlists

        ! Derived parameters
        real*8 :: epsfmm
        integer :: nbox

      end type fmm2d_tree_data

      ! Global tree storage (max 100 concurrent trees)
      integer, parameter :: MAX_TREES = 100
      type(fmm2d_tree_data), save, target :: tree_storage(MAX_TREES)
      logical, save :: tree_initialized = .false.

      contains

c
c     Initialize tree storage on first use
c
      subroutine fmm2d_tree_init()
      implicit none
      integer i

      if (tree_initialized) return

      do i = 1,MAX_TREES
        tree_storage(i)%allocated = .false.
      enddo

      tree_initialized = .true.

      return
      end subroutine fmm2d_tree_init
c
c
c     Allocate a new tree handle
c
      subroutine fmm2d_tree_alloc(tree_handle,ier)
      implicit none
      integer, intent(out) :: tree_handle
      integer, intent(out) :: ier
      integer i

      call fmm2d_tree_init()

      ier = 0
      tree_handle = 0

      do i = 1,MAX_TREES
        if (.not. tree_storage(i)%allocated) then
          tree_storage(i)%allocated = .true.
          tree_handle = i
          return
        endif
      enddo

      ! No free slots
      ier = 1
      return
      end subroutine fmm2d_tree_alloc
c
c
c     Get pointer to tree data by handle
c
      subroutine fmm2d_tree_get(tree_handle,tree,ier)
      implicit none
      integer, intent(in) :: tree_handle
      type(fmm2d_tree_data), pointer, intent(out) :: tree
      integer, intent(out) :: ier

      ier = 0
      nullify(tree)

      if (tree_handle .lt. 1 .or. tree_handle .gt. MAX_TREES) then
        ier = 1
        return
      endif

      if (.not. tree_storage(tree_handle)%allocated) then
        ier = 2
        return
      endif

      tree => tree_storage(tree_handle)

      return
      end subroutine fmm2d_tree_get
c
c
c     Free tree handle and deallocate arrays
c
      subroutine fmm2d_tree_dealloc(tree_handle,ier)
      implicit none
      integer, intent(in) :: tree_handle
      integer, intent(out) :: ier
      type(fmm2d_tree_data), pointer :: tree

      call fmm2d_tree_get(tree_handle,tree,ier)
      if (ier .ne. 0) return

      ! Deallocate arrays
      if (allocated(tree%wlists)) deallocate(tree%wlists)
      if (allocated(tree%laddr)) deallocate(tree%laddr)
      if (allocated(tree%scale)) deallocate(tree%scale)
      if (allocated(tree%bsize)) deallocate(tree%bsize)
      if (allocated(tree%nterms)) deallocate(tree%nterms)

      ! Mark as free
      tree%allocated = .false.

      return
      end subroutine fmm2d_tree_dealloc

      end module fmm2d_tree_module
