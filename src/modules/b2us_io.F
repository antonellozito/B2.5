      module b2us_io

      use b2mod_types
      use b2mod_version
      use b2mod_time
      use b2us_map
      use b2us_geo
      use b2us_plasma

      implicit none

      private

*     b2fstate_initial: open_b2fstate opened b2fstati, b2fstate being
*     empty or absent (b2mn wrote b2fplasma with the plasma of the error
*     that stopped the run, which replaces the plasma of b2fstati)
      logical, save :: b2fstate_initial = .false.
      public :: open_b2fstate, b2fstate_initial

      public :: read_b2fgmtry,  write_b2fgmtry,
     .          read_b2fstate,  write_b2fstate,
     .          read_b2fplasma, write_b2fplasma,
     .          read_b2fgmtry_header

      contains

***********************************************************************
*     Module containing basic routines for reading/writing data files
*     in unstructured format.
*     Note: file headers are read elsewhere
***********************************************************************

      subroutine read_b2fgmtry_header(ninp, version, m, lblgm, isymm)
      implicit none
      integer, intent(in) :: ninp
      character(10), intent(in) :: version
      type(mapping), intent(inout) :: m
      character(120), intent(out) :: lblgm
      integer, intent(out) :: isymm
      integer :: idum(0:9)

      if (version.ge.'03.002.000') then
        m%hasTopologicalData = .false.
        m%topoID = 0
        m%geometryID = 13
        m%nXpt = 0
        m%nOpt = 0
        m%nStr = 0
        m%nTgc = 0
        m%nDiv = 0
        m%nfcDiv = 0

        call cfruin (ninp, 7, idum, 'nCi,nCg,nCv,nFc,nVx,nFs,nFt')
        m%nCi = idum(0)
        m%nCg = idum(1)
        m%nCv = idum(2)
        m%nFc = idum(3)
        m%nVx = idum(4)
        m%nFs = idum(5)
        m%nFt = idum(6)

        call cfruin (ninp, 5, idum,
     &    'nCmxVx,nCmxFc,nVmxCv,nVmxFc,nCmxNv')
        m%nCmxVx = idum(0)
        m%nCmxFc = idum(1)
        m%nVmxCv = idum(2)
        m%nVmxFc = idum(3)
        m%nCmxNv = idum(4)

        call cfruin (ninp, 1, idum, 'isClassicalGrid')
        m%isClassicalGrid = idum(0)
        call cfruin_opt (ninp, 3, idum, 'nx,ny,nncut')
        if (m%isClassicalGrid.eq.1) then
          m%nx = idum(0)
          m%ny = idum(1)
          m%nncut = idum(2)
        else
          m%nx = -1
          m%ny = -1
          m%nncut = -1
        endif

        if (version.ge.'03.002.002') then
          m%hasTopologicalData = .true.
          call cfruin (ninp, 1, idum, 'topoflag')
          m%topoID = idum(0)
          call cfruin (ninp, 1, idum, 'GEOMETRY_ID')
          m%geometryID = idum(0)
          call cfruin (ninp, 6, idum, 'nX,nO,nS,nT,nDiv,nDivFc')
          m%nXpt = idum(0)
          m%nOpt = idum(1)
          m%nStr = idum(2)
          m%nTgc = idum(3)
          m%nDiv = idum(4)
          m%nfcDiv = idum(5)
        endif

        call cfruch (ninp, 120, lblgm, 'label')
        call cfruin (ninp, 1, idum, 'isymm')
        isymm = idum(0)
      else
        call xerrab('read_b2fgmtry_header: old version detected, ' //
     &    'not supported')
      endif

      return
      end subroutine read_b2fgmtry_header

***********************************************************************

      subroutine read_b2fgmtry (ninp, mpg, gm)
      implicit none
      integer, intent(in) :: ninp
      type (mapping), intent (inout) :: mpg
      type (geometry), intent (inout) :: gm

* ..read the mapping and geometry

      call read_mapping (ninp, mpg)
      call read_geometry (ninp, mpg%nCv, mpg%nFc, mpg%nVx, mpg%nFs,
     .  mpg%nFt, mpg%nVmxCv, gm)

      return
      end subroutine read_b2fgmtry

***********************************************************************

      subroutine write_b2fgmtry (nout, mpg, gm)
      implicit none
      integer :: nout
      type (mapping), intent (in) :: mpg
      type (geometry), intent (in) :: gm

* ..write the mapping and geometry

      call write_mapping (nout, mpg)
      call write_geometry (nout, mpg%nCv, mpg%nFc, mpg%nVx, mpg%nFs, gm)

      return
      end subroutine write_b2fgmtry

***********************************************************************

      subroutine read_b2fstate (ninp, nCv, nFc, ns, st)
      implicit none
      integer, intent(in) :: ninp, nCv, nFc, ns
      type(B2State), intent(inout) :: st

* ..read the state

      call read_state (ninp, nCv, nFc, ns, st)

      return
      end subroutine read_b2fstate

***********************************************************************

      subroutine write_b2fstate (nout, nCv, nFc, ns, st)
      implicit none
      integer, intent(in) :: nout, nCv, nFc, ns
      type(B2State), intent(in) :: st

* ..write the state

      call write_state (nout, nCv, nFc, ns, st)

      return
      end subroutine write_b2fstate

***********************************************************************

      subroutine read_b2fplasma (ninp, nCv, nFc, ns, st)
      implicit none
      integer, intent(in) :: ninp, nCv, nFc, ns
      type(B2State), intent(inout) :: st

* ..read the plasma state

      call read_plasma (ninp, nCv, nFc, ns, st)

      return
      end subroutine read_b2fplasma

***********************************************************************

      subroutine write_b2fplasma (nout, nCv, nFc, ns, st)
      implicit none
      integer, intent(in) :: nout, nCv, nFc, ns
      type(B2State), intent(in) :: st

* ..write the plasma state

      call write_plasma (nout, nCv, nFc, ns, st)

      return
      end subroutine write_b2fplasma

***********************************************************************

      subroutine open_b2fstate (nput)
      implicit none
      integer, intent(in) :: nput
      integer nbyte, ios, iline
      logical found
      character*256 fname
      character*80 line
      external cfopen, find_file, xertst

* ..open b2fstate, or b2fstati when b2fstate has no plasma state: empty,
*   absent, or only the header b2mn writes at its start (a run stopped
*   by an error, see b2fstate_initial)

      inquire (file='b2fstate', exist=found, size=nbyte)
      b2fstate_initial = .not.(found .and. nbyte.ne.0)
      if (.not.b2fstate_initial) then
*   ..the header holds a few species blocks, the state starts with na
        b2fstate_initial = .true.
        open (nput, file='b2fstate', status='old', form='formatted',
     &        iostat=ios)
        iline = 0
        do while (ios.eq.0 .and. iline.lt.10000)
          read (nput, '(a)', iostat=ios) line
          iline = iline + 1
          if (ios.eq.0 .and. line(1:4).eq.'*cf:') then
            if (index(line(5:), ' na ').gt.0) then
              b2fstate_initial = .false.
              exit
            end if
          end if
        end do
        close (nput)
      endif
      if (.not.b2fstate_initial) then
        call cfopen (nput,'b2fstate','old','un*formatted')
        return
      endif
      fname = 'b2fstati'
      call find_file (fname, found)
      call xertst (found,
     &  'b2fstate without plasma state, and no b2fstati')
      write (*,'(1x,3a)') 'b2fstate without plasma state (an error ',
     &  'stopped the run): the plasma of b2fplasma, the species of ',
     &  trim(fname)
      call cfopen (nput,trim(fname),'old','un*formatted')

      return
      end subroutine open_b2fstate

      end module b2us_io
