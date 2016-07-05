      MODULE nonlinadj
      !
      ! Module for non-linear adjoint iterations.
      !
      ! Currently includes the 'revolve' algorithm implemented by Oana Marin and
      ! Michel Schanen. Jacopo Canton re-organised everything in this module,
      ! commented and tested.
      !
      ! The 'revolve' algorithm is described in
      !    Griewank, A. & Walther, A. 2000
      !    "Algorithm 799: revolve: an implementation of checkpointing for the
      !    reverse or adjoint mode of computational differentiation"
      !    ACM Transaction on Mathematical Software, Vol. 26, pages 19--45
      !
      ! To works, this module needs to be linked with the C/C++ revolve
      ! library.
      !
      ! jcanton@mech.kth.se

      ! TODO
      ! - move here the namelist for the parameters? (I hate those lists)
      ! - some steps of revolve produce the following error message:
      !   'Problem calculating new DT Discriminant=' with negative discriminant

      implicit none

      ! was used when trying to pass fields between usrdat3 and useric
      !! include 'SIZE_DEF'
      !! include 'SIZE'
      !! real :: vx_nladj(lx1, ly1, lz1, lelv)
      !! real :: vy_nladj(lx1, ly1, lz1, lelv)
      !! real :: vz_nladj(lx1, ly1, lz1, lelv)

      ! Interface to the functions in the revolve library
      !
      interface
         function wrap_revolve(check, capo, fine, snaps, info, nid)
            implicit none
            integer :: wrap_revolve, check, capo, fine
            integer :: snaps, info, nid
         end function wrap_revolve 
         subroutine wrap_revolve_reset()
            implicit none
         end subroutine wrap_revolve_reset
         !FUNCTION adj_stack_empty()
         !   implicit none
         !   integer :: adj_stack_empty
         !END FUNCTION adj_stack_empty 
      end interface

      CONTAINS
!==============================================================================

      subroutine fwd_bkw_rvlv(capo_in, fine_in, snaps_in, info_in)
         !
         ! Integrate the direct problem forward in time and the non-linear
         ! adjoint problem backward using the 'revolve' algorithm.
         !
         implicit none

         include 'SIZE_DEF' ! LDMT1
         include 'SIZE'
         include 'TSTEP_DEF' ! ISTEP
         include 'TSTEP' ! DT
         include 'INPUT_DEF'
         include 'INPUT' ! IFPERT
         include 'SOLN_DEF' ! v[xyz], pr, v[xyz]lag
         include 'SOLN'
         include 'ADJOINT_DEF' ! IFADJ
         include 'ADJOINT'
         include 'USERPAR' ! usr_debug

         integer, intent(in) :: capo_in, fine_in, snaps_in, info_in
         integer :: capo, fine, snaps, info, whatodo, check, oldcapo
         logical :: first
         !
         real :: cfl
         integer :: n
         real, allocatable :: chkp(:,:,:,:)


         ! initialise the revolve algorithm
         !
         ! The first time revolve is called, 'check' must be set to -1 so that
         ! appropriate initialisations can be performed internally.
         ! Later, 'check' represents the number of checkpoints stored on the
         ! stack at any one time
         !
         check = -1
         !
         ! The beginning of the subrange currently being processed is defined
         ! by the state number 'capo'.
         ! The value of the state number 'fine' determines the end of the
         ! subrange currently being processed.
         ! So the pair '(capo, fine)' always represents the initial and final
         ! state of the subsequence of states currently being traversed
         ! backward.
         !
         capo = capo_in ! capo <= fine
         fine = fine_in ! Essentially corresponds to the number of time steps.
         !
         ! In contrast to the first three parameters, 'snaps' is called by
         ! value and remains unaltered during the call.
         ! The variable 'snaps' contains the upper bound on the number of
         ! checkpoints stored at any time and is selected by the user, possibly
         ! with the help of the routine 'expense' not included in this module.
         ! (see the paper).
         ! IT IS ASSUMED THAT SNAPS IS NEVER GREATER THAN 64, OTHERWISE THE
         ! EXECUTION IS STOPPED.
         ! If even more checkpoints are needed, the user must change the source
         ! code.
         !
         snaps = snaps_in
         !
         ! 'info' determines how much information about the actions performed
         ! will be printed.
         ! When info=0, no information is sent to standard output.
         ! When info>0, the first call of revolve produces an output that
         ! contains a prediction of the number of forward steps and the factor
         ! by which the execution will be slower than a regular forward
         ! simulation.
         ! Furthermore, before 'terminate' is returned, the total counts for
         ! the number of forward steps, savings of current states onto the
         ! stack, and actions taken are printed on standard output.
         ! If an error occurs the return value of 'info' contains information
         ! about the reason for the irregular termination of revolve.
         !
         info = info_in
         !
         first = .true.


         ! Allocate checkpointing arrays
         !
         n = nx1*ny1*nz1*nelv
         !
         allocate(chkp(n,3,snaps,nbdinp))


         ! Initialise the algorithm
         !
         whatodo=wrap_revolve(check, capo, fine, snaps, info, nid)
         !
         if (info.gt.0) then
            if (nid.eq.0) write(*,*) 'whatodo = ', whatodo
         endif

         select case(whatodo) 
            case(2)
               if (info.gt.0) then
                  write(*,*) 'Store in checkpoint number ',check, '.'
               endif
               call copy(chkp(:,1,check+1,1),vxlag(1,1,1,1,1),n)
               call copy(chkp(:,2,check+1,1),vylag(1,1,1,1,1),n)
               call copy(chkp(:,3,check+1,1),vzlag(1,1,1,1,1),n)
               call copy(chkp(:,1,check+1,2),vxlag(1,1,1,1,2),n)
               call copy(chkp(:,2,check+1,2),vylag(1,1,1,1,2),n)
               call copy(chkp(:,3,check+1,2),vzlag(1,1,1,1,2),n)
               call copy(chkp(:,1,check+1,3),vx,n)
               call copy(chkp(:,2,check+1,3),vy,n)
               call copy(chkp(:,3,check+1,3),vz,n)
            case default
               if(nid.eq.0) write(*,*) 'Error!'
         end select
         
         ! revolve algorithm
         !
         do while (whatodo .ne. 6) ! whatodo=6=terminate

            oldcapo = capo + 1
            whatodo = wrap_revolve(check, capo, fine, snaps, info, nid)

            select case(whatodo) 
            !
            case(1)
               !
               ! ADVANCE
               ! Integrate the direct problem forward to 'capo'
               if (info.gt.0) then
                  if (nid.eq.0) write(*,*) 'ADVANCE'
                  if (nid.eq.0) write(*,*) 'Advance from', oldcapo,
     $                                     ' to ', capo, '.'
               endif

               time = oldcapo*DT
               ifpert = .false.
               ifadj  = .false.
               !count  = 0 ! TODO not used?
               do istep = oldcapo, capo
                  !
                  call nek_advance
                  !
                  if (usr_debug.gt.0) then
                     call compute_cfl(cfl,vx,vy,vz,dt)
                     if (nid.eq.0) write(*,*) 'cfl = ', cfl 	
                  endif
               enddo

            case(2)
               !
               ! TAKESHOT
               ! Save the current state as a checkpoint
               if (info.gt.0) then
                  if (nid.eq.0) write(*,*) 'TAKESHOT'
                  if (nid.eq.0) write(*,*) 'Store in checkpoint number',
     $                                     check, '.'
               endif
               call copy(chkp(:,1,check+1,1),vxlag(1,1,1,1,1),n)
               call copy(chkp(:,2,check+1,1),vylag(1,1,1,1,1),n)
               call copy(chkp(:,3,check+1,1),vzlag(1,1,1,1,1),n)
               call copy(chkp(:,1,check+1,2),vxlag(1,1,1,1,2),n)
               call copy(chkp(:,2,check+1,2),vylag(1,1,1,1,2),n)
               call copy(chkp(:,3,check+1,2),vzlag(1,1,1,1,2),n)
               call copy(chkp(:,1,check+1,3),vx,n)
               call copy(chkp(:,2,check+1,3),vy,n)
               call copy(chkp(:,3,check+1,3),vz,n)
               !
               ! TODO useless?
               !if(firstvx.eqv..true.) then
               !   firstvx=.false.
               !end if

            case(3) 
               !
               ! FIRSTURN
               ! Advance with recording and perform first reverse step
               if (info.gt.0) then
                  if (nid.eq.0) write(*,*) 'FIRSTURN'
                  if (nid.eq.0) write(*,*) 'Advance with recording and',
     $                                     ' perform first reverse step'
               endif

               if(first.eqv..true.) then 
                  call copy(vxp,vx,n)
                  call copy(vyp,vy,n)
                  call copy(vzp,vz,n)
                  first=.false.
               end if
               if (usr_debug.gt.0) then
                  call outpost(vx,vy,vz,pr,t,'lst')
               endif
               !
               ! reverse
               ifpert = .true.
               npert  = 1
               ifadj  = .true. 
               !
               ! TODO this will be needed for the two-level checkpointing
               !call full_restart_save(istep)
               !t = real(gstep)*real(DT) ! gstep is the current global number of steps
               !tp = real(gstep)*real(DT)
               !gstep=gstep-1
               !
               if (usr_debug.gt.0) then
                  call outpost(vx,vy,vz,pr,t,'fdd')
                  call outpost(vxp,vyp,vzp,prp,tp,'adj')
               endif
               !
               call nek_advance
               !
               if (usr_debug.gt.0) then
                  call compute_cfl(cfl,vx,vy,vz,dt)
                  if (nid.eq.0) write(*,*) 'cfl = ', cfl 	
               endif

            case(4)
               !
               ! YOUTURN
               ! subsequent, combined forward/reverse steps
               ! there actually is no forward step here because we're not doing
               ! AD
               if (info.gt.0) then
                  if (nid.eq.0) write(*,*) 'YOUTURN'
                  if (nid.eq.0) write(*,*) 'Forward and reverse one',
     $                                     ' step'
               endif
               !
               ! forward
               t = capo*DT
               ifpert = .false.
               ! npert = 0
               ifadj  = .false.
               istep=capo
               !
               if (info.gt.0) then
                  if (nid.eq.0) write(*,*) 'istep: ', istep,
     $                                     ' capo: ', capo,
     $                                     ' fine: ', fine,
     $                                     ' time: ', time
               endif
               !
               if (usr_debug.gt.0) then
                  if (mod(istep,100).eq.0) then
                     ! TODO this will be needed for the two-level checkpointing
                     !t = real(gstep)*real(DT)
                     call outpost(vx,vy,vz,pr,t,'fdd')
                  endif
               endif
               !
               ! reverse
               ifpert = .true.
               npert  = 1
               ifadj  = .true. 
               !
               call nek_advance
               !
               if (usr_debug.gt.0) then
                  if (mod(istep,100).eq.0) then
                     !tp = real(gstep)*real(DT)
                     call outpost(vxp,vyp,vzp,prp,tp,'adj')
                  endif
                  call compute_cfl(cfl,vx,vy,vz,dt)
                  if (nid.eq.0) write(*,*) 'cfl = ', cfl 	
               endif
               !
               ! TODO this will be needed for the two-level checkpointing
               !gstep=gstep-1

            case(5)
               !
               ! RESTORE
               if (info.gt.0) then
                  if (nid.eq.0) write(*,*) 'RESTORE'
                  if (nid.eq.0) write(*,*) 'Restore from checkpoint',
     $                                     ' number ', check
               endif
               call copy(vxlag(1,1,1,1,1),chkp(:,1,check+1,1),n)
               call copy(vylag(1,1,1,1,1),chkp(:,2,check+1,1),n)
               call copy(vzlag(1,1,1,1,1),chkp(:,3,check+1,1),n)
               call copy(vxlag(1,1,1,1,2),chkp(:,1,check+1,2),n)
               call copy(vylag(1,1,1,1,2),chkp(:,2,check+1,2),n)
               call copy(vzlag(1,1,1,1,2),chkp(:,3,check+1,2),n)
               call copy(vx,chkp(:,1,check+1,3),n)
               call copy(vy,chkp(:,2,check+1,3),n)
               call copy(vz,chkp(:,3,check+1,3),n)

            case(6)
               !
               ! TERMINATE
               if (info.gt.0) then
                  if (nid.eq.0) write(*,*) 'TERMINATE'
                  if (nid.eq.0) write(*,*) 'Computation done.'
               endif

            case default
               !
               ! Something wrong happened
               if(nid .eq. 0) write(*,*) 'Error!'

            end select

         enddo

         ! TODO this will be needed for the two-level checkpointing
         !tp = real(gstep)*real(DT)

         if (usr_debug.gt.0) then
            call outpost(vxp,vyp,vzp,prp,tp,'rev')
         endif

         deallocate(chkp)
         !
         ! reset revolve so it can be used for the next run
         call wrap_revolve_reset()

         ! TODO useless?
         !ssize=ssize-1

      end subroutine fwd_bkw_rvlv

!==============================================================================
      END MODULE nonlinadj
