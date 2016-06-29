      module postprocess
         !real, external :: dnekclock ! defined in comm_mpi.f
         real, external :: glmin, glmax ! defined in math.f


         real divv
         COMMON /SCRUZ/  DIVV (LX2,LY2,LZ2,LELV)
         real work1(lx1,ly1,lz1,lelv),
     $        work2(lx1,ly1,lz1,lelv),
     $        divm1(lx1,ly1,lz1,lelv)
         integer ntdump
         common /rdump/ ntdump

         real rtime, tic
         save rtime
         data rtime  /0./
         real divmax

         real work3(lx1,ly1,lz1,lelv),
     $        work4(lx1,ly1,lz1,lelv),
     $        v_s(lx1,ly1,lz1,lelv), ! streamwise
     $        omega_s(lx1,ly1,lz1,lelv)
         real vorticity(lx1*ly1*lz1*lelv,3)
         real phi(lx1,ly1,lz1,lelv)
         real h1(lx1,ly1,lz1,lelv)
         real h2(lx1,ly1,lz1,lelv)
         real rhs(lx1,ly1,lz1,lelv)

         integer n, nn, maxit
         real toler, p22 

         n  = nelv*lx1*ly1*lz1
         nn = nelv*nx1*ny1*nz1

         ! Output maximum of divergence
         ! This maps divergence from pressure (Pn-2) to velocity mesh (Pn)
         call mappr(divm1,divv,work1,work2)
         ! and averages the divergence on element boundaries
         ! In contrast to outpost(), which does not average
         ! when writing divergence in pressure ouput 
         ! -> therefore, slight discrepancy < ~ factor 2

         divmax = glmax(divm1, nx1*ny1*nz1*nelv)

         if (nid.eq.0) write(*,*) '||div||_max', Divmax


         if (mod(ISTEP,IOSTEP).eq.0.and. istep.gt.0) then
            !
            ! compute and output some postprocessing quantities

!           ! output divergence (in pressure field)
!           call outpost(vx,vy,vz,divv,t,'div')
!
!           ! Compute lambda_2
!           if (NID.eq.0) then
!              write(6,*) ISTEP,IOSTEP,TIME,' compute lambda2'
!           endif
!           call lambda2(t(1,1,1,1,1))
!
!           ! Output streamfunction phi and velocity v_s of streamwise direction
!           !
!           ! compute vorticity
!           call rzero(vorticity,3*n)
!           call comp_vort3(vorticity,work1,work2,vx,vy,vz)  
!           call copy(vorticity(1,1),work1,n);
!           call copy(vorticity(1,2),work2,n);
!           call copy(vorticity(1,3),work3,n);
!           call chsign(work1,n);
!           call chsign(work2,n);
!           call chsign(work3,n);
!           call rzero(h2,n)
!           call chsign(omega_s,n)
!           !
!           ! multiply with the mass matrix
!           call col3(rhs,omega_s,bm1,n)
!           ! setup Poisson problem: h2=0, h1=1
!           call rone(h1,n)
!           ! set solver tolerance
!           ! and save param22
!           toler = 1e-14
!           p22 = param(22)
!           param(22) = toler
!           maxit = 10000    
!           ! Solve Poisson equation:    
!           call hmholtz('PHI ',phi,rhs,h1,h2,v2mask,vmult,1,
!    $                   toler,maxit,1)
!           ! reset param22
!           param(22) = p22
!
!           call outpost(rhs,phi,v_s,pr,t,'stf')

         endif


!        ! Restart routines
!        ! (uncomment the INCLUDE statement at the beginning)
!        call checkpoint
!        ! This call is already executed in setics()/ic.f
!        ! However, it has the wrong time time-value (=0.0) 
!        ! when using our restarting, which initializes 'time' 
!        if (istep.eq.0.and.timeio.ne.0.0)  then
!           ntdump = int( (time+1e-14)/timeio )
!        endif


      end module postprocess
