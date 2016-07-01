      MODULE misc_stuff
      !
      ! Contains stuff that did not fit anywhere else
      ! jcanton@mech.kth.se

      implicit none


      CONTAINS
!==============================================================================

      real function kinetic_energy(u, v, w)

         implicit none

         include 'SIZE_DEF'
         include 'SIZE'
         include 'MASS_DEF' ! bm1
         include 'MASS'

         real, intent(in) :: u(1), v(1), w(1)
         integer :: n
         real, dimension(:), allocatable :: tmp
         real, external :: glsum ! math.f

         n = nx1*ny1*nz1*nelv

         allocate( tmp(n) )

         call rzero(tmp, n)
         call vdot3(tmp, u,v,w, u,v,w, n)
         call col2(tmp, bm1, n)

         kinetic_energy = 0.5*glsum(tmp, n)

         deallocate( tmp )

         return
      end function kinetic_energy

!==============================================================================
      END MODULE misc_stuff
