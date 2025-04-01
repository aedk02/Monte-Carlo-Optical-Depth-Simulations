module packet_dynamics
!--------------------------------------------------------------------------------------
! *** Module for packet emission, scattering and movement ***
!--------------------------------------------------------------------------------------
contains

subroutine EMIT_PACKET(x, y, z, nx, ny, nz)     
    !--------------------------------------------------------------------------------------
    ! ***  Initialises the position and direction of packet from a random seed        ***
    !--------------------------------------------------------------------------------------
      implicit none
      double precision :: x,y,z,nx,ny,nz, phi, theta, ran1, ran2
      double precision :: pi = 4.d0*ATAN(1.d0)    !store a value for pi 
    
    
      !set the initial position of the packet to the orign
      x = 0.d0
      y = 0.d0
      z = 0.d0
    
      call random_number(ran1)
      call random_number(ran2)
    
      !generate random polar and azimuthal angles
      theta = ACOS((2.d0*ran1)-1.d0)
      phi = 2.d0 * pi * ran2
    
      !from these angles, calculate random direction packet moves in
      nx = SIN(theta)*COS(phi)
      ny = SIN(theta) * SIN(phi)
      nz = COS(theta)
      
end subroutine EMIT_PACKET
    
subroutine SCATTER(nx,ny,nz)      
    !--------------------------------------------------------------------------------------
    ! ***         Scatters a packet                  ***
    !--------------------------------------------------------------------------------------
      implicit none
      double precision :: nx,ny,nz, phi, theta, ran1, ran2
      double precision :: pi = 4.d0*ATAN(1.d0)    !store a value for pi 
    
      !call the random number generator 
      call random_number(ran1)
      call random_number(ran2)
    
      !generate random polar and azimuthal angles
      theta = ACOS((2.d0*ran1)-1)
      phi = 2.d0 * pi * ran2
    
      !from these angles, calculate random direction packet moves in
      nx = SIN(theta)*COS(phi)
      ny = SIN(theta) * SIN(phi)
      nz = COS(theta)
      
end subroutine SCATTER
    
    

subroutine MOVE_PACKET(rmax, taumax, j, nx, ny, nz, x, y, z, L)     
    !--------------------------------------------------------------------------------------
    ! *** Moves the packet a random distance determined by optical depth and direction ***
    !--------------------------------------------------------------------------------------
        implicit none
        double precision :: rmax, L, tau, ran
        double precision :: taumax(5)
        integer :: j
        double precision :: nx, ny, nz, x, y, z
    
        call random_number(ran)               !random number generator
        tau = -LOG(ran)                       !generate random optical depth
        L = tau * rmax / taumax(j)            !calculate travel distance from tau
        x = x + L * nx                        !update packet x-position
        y = y + L * ny                        !update packet y-position
        z = z + L * nz                        !update packet z-position
    
end subroutine MOVE_PACKET

end module packet_dynamics