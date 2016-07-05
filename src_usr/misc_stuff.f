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
         integer :: nn
         real, dimension(:), allocatable :: tmp
         real, external :: glsum ! math.f

         nn = nx1*ny1*nz1*nelv

         allocate( tmp(nn) )

         call rzero(tmp, nn)
         call vdot3(tmp, u,v,w, u,v,w, nn)
         call col2(tmp, bm1, nn)

         kinetic_energy = 0.5*glsum(tmp, nn)

         deallocate( tmp )

         return
      end function kinetic_energy

!==============================================================================
      END MODULE misc_stuff
