      MODULE nonlinadj
      !
      ! Module for non-linear adjoint iterations.
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

      CONTAINS

!==============================================================================
!     MAIN ROUTINE
!==============================================================================

      subroutine compute_nlopt(update_alg, fwd_bkw_alg)
         !
         ! Main routine for the computation of non-linear optimal perturbations
         ! jcanton@mech.kth.se
         !
         ! PREREQUISITES:
         ! - we suppose that an initial velocity field with given kinetic
         !   energy has been stored in v[xyz]
         !
         implicit none

         include 'SIZE_DEF'
         include 'SIZE' ! nid
         include 'INPUT_DEF' ! param()
         include 'INPUT'         
         include 'SOLN_DEF' ! v[xyz]
         include 'SOLN'         
         include 'USERPAR' ! usr_debug, nlopt_res, nlopt_maxit

         integer, intent(in) :: update_alg, fwd_bkw_alg
         integer :: nn
         integer :: it_0, it_end
         real    :: tol
         integer :: ite, maxit

         real, dimension(nlopt_maxit) :: res
         !real :: res

         real :: udx(lx1, ly1, lz1, lelv), ! direct velocity field
     $           udy(lx1, ly1, lz1, lelv),
     $           udz(lx1, ly1, lz1, lelv)
         real :: uax(lx1, ly1, lz1, lelv), ! adjoint velocity field
     $           uay(lx1, ly1, lz1, lelv),
     $           uaz(lx1, ly1, lz1, lelv)


         !---------------------------------------------------------------------
         ! Initialise stuff
         !
         nn    = nx1*ny1*nz1*nelv
         tol   = nlopt_res
         maxit = nlopt_maxit
         !
         it_0   = 0
         it_end = param(11)
         if (usr_debug.gt.0) then
            if (nid.eq.0) write(*,*) 'nsteps = ', it_end
         endif
         !
         ! save initial condition in ud[xyz]
         call copy(udx, vx, nn)
         call copy(udy, vy, nn)
         call copy(udz, vz, nn)

         !---------------------------------------------------------------------
         ! Call forward-backward iterations
         !
         ite = 2
         res(1) = 1e3
         !res = 1e3

         do while (res(ite-1) .ge. tol .and. ite .le. maxit)
         !do while (res .ge. tol .and. ite .le. maxit)
            !
            ! integrate forward and backward
            call nladj_fwd_bkw(it_0, it_end, rvlv_snaps, fwd_bkw_alg)
            !
            ! copy 'initial' adjoint field resulting from fwd_bkw iterations
            call copy(uax, vxp, nn)
            call copy(uay, vyp, nn)
            call copy(uaz, vzp, nn)
            !
            ! update iterate
            call nlopt_update(res(ite),
     $                        udx,udy,udz, uax,uay,uaz,
     $                        update_alg)
            !
            ! copy updated field back
            call copy(vx, udx, nn)
            call copy(vy, udy, nn)
            call copy(vz, udz, nn)
            !
            if (usr_debug.gt.0) then
               call outpost(vx,vy,vz,pr,t,'opt')
            endif
            !
            if (usr_debug.gt.0) then
               if (nid.eq.0) then
                  write(*,*) 'ite ', ite-1, ' res = ', res(ite)
                  !write(*,*) 'ite ', ite, ' res = ', res
               endif
            endif
            !
            ite = ite + 1
            !
         enddo


      end subroutine compute_nlopt

!==============================================================================
!     UPDATE ROUTINES
!==============================================================================

      subroutine nlopt_update(res, udx,udy,udz, uax,uay,uaz, alg)
         !
         ! Update the iterate ud[xyz](n) towards a maximum of the lagrangian
         ! cost function.
         ! This routine is a wrapper for different update algorithms.
         ! jcanton@mech.kth.se
         !
         implicit none

         include 'SIZE_DEF'
         include 'SIZE' ! nid
         include 'USERPAR' ! usr_debug, nlopt_res

         real :: res
         real :: udx(1), udy(1), udz(1)
         real :: uax(1), uay(1), uaz(1)
         integer, intent(in) :: alg

         ! update iterate
         !
         select case(alg)

            case(1)
               !
               ! STEEPEST ASCENT
               if (usr_debug.gt.0 .and. nid.eq.0) then
                  write(*,*) 'calling steepest ascent update' 	
               endif
               call nlopt_stpst_ascnt(res, udx,udy,udz, uax,uay,uaz)

            case default
               if(nid.eq.0) then
                  write(*,*) '*** ERROR ***'
                  write(*,*) 'Unknown "alg" in nlopt_update'
                  write(*,*) '*** ERROR ***'
               endif
               call exitt()

         end select

      end subroutine nlopt_update

!------------------------------------------------------------------------------

      subroutine nlopt_stpst_ascnt(res, udx,udy,udz, uax,uay,uaz)
         !
         ! Update the iterate ud[xyz](n) towards the maximum of kinetic energy
         ! using the steepest ascent algorithm.
         ! jcanton@mech.kth.se
         !
         use misc_stuff ! kinetic_energy(), kinetic_energy_2()

         implicit none

         include 'SIZE_DEF' ! nx1, ny1, nz1, nelv
         include 'SIZE' ! nid
         include 'USERPAR' ! usr_debug, nlopt_e0, stasc_eps

         real    :: res
         real    :: udx(1), udy(1), udz(1) ! direct  velocity field
         real    :: uax(1), uay(1), uaz(1) ! adjoint velocity field
         real    :: eps, lambda_1, lambda_2, lambda
         integer :: nn
         real    :: ek0, ek_d, ek_a, ek_da
         real    :: aa, bb, cc
         real    :: tmpx(lx1, ly1, lz1, lelv),
     $              tmpy(lx1, ly1, lz1, lelv),
     $              tmpz(lx1, ly1, lz1, lelv)

         !---------------------------------------------------------------------
         ! Initialise stuff
         !
         nn  = nx1*ny1*nz1*nelv
         ek0 = nlopt_e0
         eps = stasc_eps

         !---------------------------------------------------------------------
         ! Compute lambda
         !
         ek_d  = kinetic_energy(udx, udy, udz)
         ek_a  = kinetic_energy(uax, uay, uaz)
         ek_da = kinetic_energy_2(udx,udy,udz, uax,uay,uaz)
         !
         aa = eps**2*ek_d
         bb = 2*eps*(ek_d - eps*ek_da)
         cc = eps**2*ek_a - 2*eps*ek_da
         !
         if (usr_debug.gt.0) then
            if (nid.eq.0) write(*,*) 'ek_0  = ', ek0
            if (nid.eq.0) write(*,*) 'ek_d  = ', ek_d
            if (nid.eq.0) write(*,*) 'ek_a  = ', ek_a
            if (nid.eq.0) write(*,*) 'ek_da = ', ek_da
            if (nid.eq.0) write(*,*) 'a = ', aa
            if (nid.eq.0) write(*,*) 'b = ', bb
            if (nid.eq.0) write(*,*) 'c = ', cc
         endif
         !
         lambda_1 = (-bb + sqrt(bb**2 - 4*aa*cc))/(2*aa)
         lambda_2 = (-bb - sqrt(bb**2 - 4*aa*cc))/(2*aa)
         lambda   = lambda_1 ! TODO works with lambda_2 as well (?!?)
         !
         if (usr_debug.gt.0) then
            if (nid.eq.0) write(*,*) 'lambda = ', lambda
         endif

         !---------------------------------------------------------------------
         ! Update velocity field
         ! u(n+1) = u(n) + eps*(lambda*u(n) - v(n))
         ! res    =             lambda*u(n) - v(n)
         !
         call copy(tmpx, udx, nn)
         call copy(tmpy, udy, nn)
         call copy(tmpz, udz, nn)
         !
         call cmult(tmpx, lambda, nn)
         call cmult(tmpy, lambda, nn)
         call cmult(tmpz, lambda, nn)
         !
         call sub2(tmpx, uax, nn)
         call sub2(tmpy, uay, nn)
         call sub2(tmpz, uaz, nn)
         !
         res = 2*kinetic_energy(tmpx,tmpy,tmpz)
         !
         call cmult(tmpx, eps, nn)
         call cmult(tmpy, eps, nn)
         call cmult(tmpz, eps, nn)
         !
         call add2(udx, tmpx, nn)
         call add2(udy, tmpy, nn)
         call add2(udz, tmpz, nn)

         if (usr_debug.gt.0) then
            if (nid.eq.0) write(*,*) ''
         endif

      end subroutine nlopt_stpst_ascnt


!==============================================================================
!     INTEGRATION ROUTINES
!==============================================================================

      subroutine nladj_fwd_bkw(it_0, it_end, snaps, alg)
         !
         ! Integrate the non-linear direct problem forward in time and the
         ! non-linear adjoint problem backward.
         ! The initial condition must be stored in v[xyz] and the resulting
         ! adjoint velocity field at time = 0 will be stored in v[xyz]p.
         ! This routine is a wrapper for different integration algorithms.
         ! jcanton@mech.kth.se
         !
         implicit none

         include 'SIZE_DEF'
         include 'SIZE' ! NID
         include 'USERPAR' ! usr_debug

         integer, intent(in) :: it_0, it_end, snaps, alg

         select case(alg)

            case(1)
               !
               ! REVOLVE ALGORITHM
               if (usr_debug.gt.0 .and. nid.eq.0) then
                  write(*,*) 'calling revolve algorithm' 	
               endif
               call fwd_bkw_rvlv(it_0, it_end, snaps)

            case default
               if(nid.eq.0) then
                  write(*,*) '*** ERROR ***'
                  write(*,*) 'Unknown "alg" in nladj_fwd_bkw'
                  write(*,*) '*** ERROR ***'
               endif
               call exitt()

         end select

      end subroutine nladj_fwd_bkw

!------------------------------------------------------------------------------

      subroutine fwd_bkw_rvlv(capo_in, fine_in, snaps_in)
         !
         ! Integrate the direct problem forward in time and the non-linear
         ! adjoint problem backward using the 'revolve' algorithm.
         ! jcanton@mech.kth.se
         !
         ! Currently includes the 'revolve' algorithm implemented by Oana Marin
         ! and Michel Schanen. Jacopo Canton re-organised everything in this
         ! module, commented and tested.
         !
         ! The 'revolve' algorithm is described in
         !    Griewank, A. & Walther, A. 2000
         !    "Algorithm 799: revolve: an implementation of checkpointing for
         !    the reverse or adjoint mode of computational differentiation"
         !    ACM Transaction on Mathematical Software, Vol. 26, pages 19--45
         !
         ! For this routine to work, this module needs to be linked with the
         ! C/C++ revolve library.
         !
         implicit none

         include 'SIZE_DEF' ! LDMT1
         include 'SIZE'
         include 'TSTEP_DEF' ! ISTEP, NBDINP
         include 'TSTEP' ! DT
         include 'INPUT_DEF'
         include 'INPUT' ! IFPERT
         include 'SOLN_DEF' ! v[xyz], pr, v[xyz]lag
         include 'SOLN'
         include 'ADJOINT_DEF' ! IFADJ
         include 'ADJOINT'
         include 'USERPAR' ! usr_debug

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

         integer, intent(in) :: capo_in, fine_in, snaps_in
         integer :: capo, fine, snaps, info, whatodo, check, oldcapo
         logical :: first
         !
         real :: cfl
         integer :: nn, time_order
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
         info = usr_debug
         !
         first = .true.

         ! Initialise stuff
         !
         time_order = nbdinp
         !
         if (.not.(time_order.eq.1 .or. time_order.eq.3)) then
            if(nid.eq.0) then
               write(*,*) '*** ERROR ***'
               write(*,*) 'Only time order = 1 and 3'
               write(*,*) 'have been implemented in fwd_bkw_rvlv'
               write(*,*) '*** ERROR ***'
            endif
            call exitt()
         endif

         ! Allocate checkpointing arrays
         !
         nn = nx1*ny1*nz1*nelv
         !
         allocate(chkp(nn,3,snaps,time_order))


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
                  if (nid.eq.0) write(*,*) 'Store in checkpoint number',
     $                                     check, '.'
               endif
               !
               select case (time_order)
                  case(1)
                     call copy(chkp(:,1,check+1,1), vx,              nn)
                     call copy(chkp(:,2,check+1,1), vy,              nn)
                     call copy(chkp(:,3,check+1,1), vz,              nn)
                  case(3)
                     call copy(chkp(:,1,check+1,1), vxlag(1,1,1,1,1),nn)
                     call copy(chkp(:,2,check+1,1), vylag(1,1,1,1,1),nn)
                     call copy(chkp(:,3,check+1,1), vzlag(1,1,1,1,1),nn)
                     call copy(chkp(:,1,check+1,2), vxlag(1,1,1,1,2),nn)
                     call copy(chkp(:,2,check+1,2), vylag(1,1,1,1,2),nn)
                     call copy(chkp(:,3,check+1,2), vzlag(1,1,1,1,2),nn)
                     call copy(chkp(:,1,check+1,3), vx,              nn)
                     call copy(chkp(:,2,check+1,3), vy,              nn)
                     call copy(chkp(:,3,check+1,3), vz,              nn)
                  case default
               end select
               !
            case default
               if(nid.eq.0) then
                  write(*,*) '*** ERROR ***'
                  write(*,*) 'Unknown "whatodo" in fwd_bkw_rvlv'
                  write(*,*) '*** ERROR ***'
               endif
               call exitt()
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

               time   = oldcapo*DT
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
               !
               select case (time_order)
                  case(1)
                     call copy(chkp(:,1,check+1,1), vx,              nn)
                     call copy(chkp(:,2,check+1,1), vy,              nn)
                     call copy(chkp(:,3,check+1,1), vz,              nn)
                  case(3)
                     call copy(chkp(:,1,check+1,1), vxlag(1,1,1,1,1),nn)
                     call copy(chkp(:,2,check+1,1), vylag(1,1,1,1,1),nn)
                     call copy(chkp(:,3,check+1,1), vzlag(1,1,1,1,1),nn)
                     call copy(chkp(:,1,check+1,2), vxlag(1,1,1,1,2),nn)
                     call copy(chkp(:,2,check+1,2), vylag(1,1,1,1,2),nn)
                     call copy(chkp(:,3,check+1,2), vzlag(1,1,1,1,2),nn)
                     call copy(chkp(:,1,check+1,3), vx,              nn)
                     call copy(chkp(:,2,check+1,3), vy,              nn)
                     call copy(chkp(:,3,check+1,3), vz,              nn)
                  case default
               end select
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
                  first=.false.
                  ! v[xyz]p = v[xyz]
                  call copy(vxp, vx, nn)
                  call copy(vyp, vy, nn)
                  call copy(vzp, vz, nn)
                  ! TODO change sign? see eq.(33) in Kerswell et al. 2014
                  ! seems not: in their pipe paper they don't
                  !call chsign(vxp, nn)
                  !call chsign(vyp, nn)
                  !call chsign(vzp, nn)
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
               !
               select case (time_order)
                  case(1)
                     call copy(vx,               chkp(:,1,check+1,1),nn)
                     call copy(vy,               chkp(:,2,check+1,1),nn)
                     call copy(vz,               chkp(:,3,check+1,1),nn)
                  case(3)
                     call copy(vxlag(1,1,1,1,1), chkp(:,1,check+1,1),nn)
                     call copy(vylag(1,1,1,1,1), chkp(:,2,check+1,1),nn)
                     call copy(vzlag(1,1,1,1,1), chkp(:,3,check+1,1),nn)
                     call copy(vxlag(1,1,1,1,2), chkp(:,1,check+1,2),nn)
                     call copy(vylag(1,1,1,1,2), chkp(:,2,check+1,2),nn)
                     call copy(vzlag(1,1,1,1,2), chkp(:,3,check+1,2),nn)
                     call copy(vx,               chkp(:,1,check+1,3),nn)
                     call copy(vy,               chkp(:,2,check+1,3),nn)
                     call copy(vz,               chkp(:,3,check+1,3),nn)
                  case default
               end select

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
               if(nid.eq.0) then
                  write(*,*) '*** ERROR ***'
                  write(*,*) 'Unknown "whatodo" in fwd_bkw_rvlv'
                  write(*,*) '*** ERROR ***'
               endif
               call exitt()

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
