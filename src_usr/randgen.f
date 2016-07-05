      MODULE randgen
      !
      ! Contains random number generator stolen from
      ! AUTHORS: Richard Chandler
      !      (richard@stats.ucl.ac.uk)
      !      Paul Northrop
      !      (northrop@stats.ox.ac.uk)
      ! LAST MODIFIED: 26/8/03
      !
      ! jcanton@mech.kth.se

      implicit none

      DOUBLE PRECISION ZBQLIX(43), B, C
      INTEGER I
      DATA (ZBQLIX(I),I=1,43) /8.001441D7,5.5321801D8,
     +      1.69570999D8,2.88589940D8,2.91581871D8,1.03842493D8,
     +      7.9952507D7,3.81202335D8,3.11575334D8,4.02878631D8,
     +      2.49757109D8,1.15192595D8,2.10629619D8,3.99952890D8,
     +      4.12280521D8,1.33873288D8,7.1345525D7,2.23467704D8,
     +      2.82934796D8,9.9756750D7,1.68564303D8,2.86817366D8,
     +      1.14310713D8,3.47045253D8,9.3762426D7 ,1.09670477D8,
     +      3.20029657D8,3.26369301D8,9.441177D6,3.53244738D8,
     +      2.44771580D8,1.59804337D8,2.07319904D8,3.37342907D8,
     +      3.75423178D8,7.0893571D7 ,4.26059785D8,3.95854390D8,
     +      2.0081010D7,5.9250059D7,1.62176640D8,3.20429173D8,
     +      2.63576576D8/
      DATA B / 4.294967291D9 /
      DATA C / 0.0D0 /

      CONTAINS
!==============================================================================


      SUBROUTINE init_random_seed()
      
         implicit none
         integer, allocatable :: seed(:)
         integer :: i, n, un, istat, dt(8), pid, t(2), s
         integer(8) :: count, tms
         
         !write(*,*) '    Generating random seed ...'
      
         call random_seed(size = n)
         allocate(seed(n))
         !  XOR:ing the current time and pid. The PID is useful in case one
         !  launches multiple instances of the same program in parallel.
         !write(*,*) '    using current time and pID'
         call system_clock(count)
         if (count /= 0) then
            t = transfer(count, t)
         else
            call date_and_time(values=dt)
            tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000
     $           + dt(2) * 31_8 * 24 * 60 * 60 * 1000
     $           + dt(3) * 24 * 60 * 60 * 60 * 1000
     $           + dt(5) * 60 * 60 * 1000
     $           + dt(6) * 60 * 1000 + dt(7) * 1000
     $           + dt(8)
            t = transfer(tms, t)
         end if
         s = ieor(t(1), t(2))
         pid = getpid() + 1099279 ! Add a prime
         !write(*,*) '    pid: ', pid
         s = ieor(s, pid)
         if (n >= 3) then
            seed(1) = t(1) + 36269
            seed(2) = t(2) + 72551
            seed(3) = pid
            if (n > 3) then
               seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
            end if
         else
            seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
         end if
         call random_seed(put=seed)
      
         !write(*,*) 'seed: ', seed
         !write(*,*) '    Done.'
         !write(*,*)
      
      END SUBROUTINE init_random_seed

!------------------------------------------------------------------------------

!!      block data ZBQLBD01
!!!
!!!        Initializes seed array etc. for random number generator.
!!!        The values below have themselves been generated using the
!!!        NAG generator.
!!!
!!         implicit none
!!
!!         COMMON /ZBQL0001/ ZBQLIX, B, C
!!         DOUBLE PRECISION ZBQLIX(43), B, C
!!         INTEGER I
!!         DATA (ZBQLIX(I),I=1,43) /8.001441D7,5.5321801D8,
!!     +         1.69570999D8,2.88589940D8,2.91581871D8,1.03842493D8,
!!     +         7.9952507D7,3.81202335D8,3.11575334D8,4.02878631D8,
!!     +         2.49757109D8,1.15192595D8,2.10629619D8,3.99952890D8,
!!     +         4.12280521D8,1.33873288D8,7.1345525D7,2.23467704D8,
!!     +         2.82934796D8,9.9756750D7,1.68564303D8,2.86817366D8,
!!     +         1.14310713D8,3.47045253D8,9.3762426D7 ,1.09670477D8,
!!     +         3.20029657D8,3.26369301D8,9.441177D6,3.53244738D8,
!!     +         2.44771580D8,1.59804337D8,2.07319904D8,3.37342907D8,
!!     +         3.75423178D8,7.0893571D7 ,4.26059785D8,3.95854390D8,
!!     +         2.0081010D7,5.9250059D7,1.62176640D8,3.20429173D8,
!!     +         2.63576576D8/
!!         DATA B / 4.294967291D9 /
!!         DATA C / 0.0D0 /
!!      end

!------------------------------------------------------------------------------

      subroutine zbqlini(seed)
!        To initialize the random number generator - either repeatably or
!        nonrepeatably. Need double precision variables because integer storage
!        can't handle the numbers involved
!
!        ARGUMENTS
!        =========
!        SEED    (integer, input). User-input number which generates elements
!        of the array ZBQLIX, which is subsequently used in the random number
!        generation algorithm. If SEED=0, the array is seeded using the system
!        clock if the FORTRAN implementation allows it.
!
         implicit none
!
!        VARIABLES
!        =========
!        SEED    See above
!        ZBQLIX  Seed array for the random number generator. Defined
!                in ZBQLBD01
!        B,C Used in congruential initialisation of ZBQLIX
!        INIT    Indicates whether generator has already been initialised

         INTEGER, INTENT(IN) :: SEED
         INTEGER :: I, INIT
         !DOUBLE PRECISION :: ZBQLIX(43), B, C
         DOUBLE PRECISION :: TMPVAR1, DSS, DMM, DHH, DDD
         !COMMON /ZBQL0001/ ZBQLIX, B, C
         SAVE INIT
!
!        Ensure we don't call this more than once in a program
!
         IF (INIT.GE.1) THEN
            IF(INIT.EQ.1) THEN
               WRITE(*,*) '************* WARNING *****************'
               WRITE(*,*) 'You have called routine ZBQLINI more'
               WRITE(*,*) 'than once. I''m ignoring any subsequent'
               WRITE(*,*) 'calls.'
               WRITE(*,*) '***************************************'
               INIT = 2
            ENDIF
            RETURN
         ELSE
            INIT = 1
         ENDIF
!
!        If SEED = 0, cat the contents of the clock into a file and transform
!        to obtain ZQBLIX(1), then use a congr.  algorithm to set remaining
!        elements. Otherwise take specified value of SEED.
!
         IF (SEED.EQ.0) THEN
            TMPVAR1 = 1.0
         ELSE
            TMPVAR1 = DMOD(DBLE(SEED), B)
         ENDIF
         ZBQLIX(1) = TMPVAR1

         DO 100 I = 2, 43
            TMPVAR1 = ZBQLIX(I-1)*3.0269D4
            TMPVAR1 = DMOD(TMPVAR1,B)       
            ZBQLIX(I) = TMPVAR1
 100     CONTINUE

      end

!------------------------------------------------------------------------------

      function zbqlu01()

         implicit none

!        Returns a uniform random number between 0 & 1, using a Marsaglia-Zaman
!        type subtract-with-borrow generator.  Uses double precision, rather
!        than integer, arithmetic throughout because MZ's integer constants
!        overflow 32-bit integer storage (which goes from -2^31 to 2^31).
!        Ideally, we would explicitly truncate all integer quantities at each
!        stage to ensure that the double precision representations do not
!        accumulate approximation error; however, on some machines the use of
!        DNINT to accomplish this is *seriously* slow (run-time increased by a
!        factor of about 3). This double precision version has been tested
!        against an integer implementation that uses long integers
!        (non-standard and, again, slow) - the output was identical up to the
!        16th decimal place after 10^10 calls, so we're probably OK ...

         DOUBLE PRECISION :: ZBQLU01
         !DOUBLE PRECISION :: B, C, ZBQLIX(43)
         DOUBLE PRECISION :: X,B2,BINV
         INTEGER CURPOS,ID22,ID43

         !COMMON /ZBQL0001/ ZBQLIX,B,C
         !SAVE /ZBQL0001/
         SAVE CURPOS,ID22,ID43
         DATA CURPOS,ID22,ID43 /1,22,43/

         B2 = B
         BINV = 1.0D0/B
 5       X = ZBQLIX(ID22) - ZBQLIX(ID43) - C
         IF (X.LT.0.0D0) THEN
            X = X + B
            C = 1.0D0
         ELSE
            C = 0.0D0
         ENDIF
         ZBQLIX(ID43) = X

!        Update array pointers. Do explicit check for bounds of each to avoid
!        expense of modular arithmetic. If one of them is 0 the others won't be

         CURPOS = CURPOS - 1
         ID22 = ID22 - 1
         ID43 = ID43 - 1
         IF (CURPOS.EQ.0) THEN
            CURPOS=43
         ELSEIF (ID22.EQ.0) THEN
            ID22 = 43
         ELSEIF (ID43.EQ.0) THEN
            ID43 = 43
         ENDIF

!        The integer arithmetic there can yield X=0, which can cause problems
!        in subsequent routines (e.g. ZBQLEXP). The problem is simply that X is
!        discrete whereas U is supposed to be continuous - hence if X is 0, go
!        back and generate another X and return X/B^2 (etc.), which will be
!        uniform on (0,1/B). 

         IF (X.LT.BINV) THEN
            B2 = B2*B
            GOTO 5
         ENDIF

         zbqlu01 = X/B2
      end

!------------------------------------------------------------------------------

      real function rnd_loc (lower, upper)
         !
         ! returns a normalised random number between lower and upper values

         implicit none

         real, intent(in) :: lower, upper
         real rnd

         rnd = real(zbqlu01())
         rnd_loc = lower + rnd*(upper - lower)

         return
      end function rnd_loc

!==============================================================================
      END MODULE randgen
