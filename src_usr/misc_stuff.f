      MODULE misc_stuff
      !
      ! Contains stuff that did not fit anywhere else
      ! jcanton@mech.kth.se

      implicit none


      CONTAINS
!==============================================================================

      real function kinetic_energy(vx, vy, vz)
         !
         ! Compute the kinetic energy of the v[xyz] velocity field
         ! jcanton@mech.kth.se
         !
         implicit none

         include 'SIZE_DEF' ! nx1, ny1, nz1, nelv
         include 'SIZE'
         include 'MASS_DEF' ! bm1
         include 'MASS'

         real, intent(in) :: vx(1), vy(1), vz(1)
         integer :: nn
         real, dimension(:), allocatable :: tmp
         real, external :: glsum ! math.f

         nn = nx1*ny1*nz1*nelv

         allocate( tmp(nn) )

         call rzero(tmp, nn)
         call vdot3(tmp, vx,vy,vz, vx,vy,vz, nn)
         call col2(tmp, bm1, nn)

         kinetic_energy = 0.5*glsum(tmp, nn)

         deallocate( tmp )

         return
      end function kinetic_energy

!==============================================================================
      END MODULE misc_stuff
