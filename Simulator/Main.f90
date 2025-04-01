program MCRT

    !This program performs a monte carlo simulation for light scattering
    !through a uniform density sphere of radius 1m. The program calculates
    !the time taken for a packet of light to travel through the sphere for
    !different optical depths. 
    
    use Packet_dynamics
    use Histogram
    use Pathlengths

    implicit none 
    double precision :: x,y,z,nx,ny,nz                                  
    double precision :: xintersect, yintersect, zintersect
    double precision :: rmax                               
    double precision :: albedo                             
    double precision :: ran                                
    double precision :: L                             
    integer :: i,j, npackets                           
    double precision,dimension(5) :: taumax     
    character(len=80) :: filename
    double precision :: time_taken, pathlength    
    double precision :: v = 299792458.d0    !speed of light [m/s]
    double precision, allocatable, dimension(:) :: times_taken                      
    integer :: num_bins
    double precision,allocatable, dimension(:) :: freq, bin_edges

    num_bins = 1000                     !number of bins for the histogram plotting
    npackets = 100000                   !number of packets or photons to simulate
    
    rmax = 1.d0                         !radius 1[m]
    albedo = 1.d0                       !albedo of the sphere

    !allocate arrays to correct size
    allocate(times_taken(npackets))
    allocate(freq(num_bins))
    allocate(bin_edges(num_bins))

    !initalise list of tau values to test 
    taumax = [dble(0.1), 1.d0, 10.d0, 30.d0, 100.d0]  !radial optical depths
   
    !main loop over test radial optical depths
    do j = 1,5
        do i = 1, npackets                       !loop over packets
            time_taken = 0.d0                    !counter for time taken to leave sphere
            pathlength = 0.d0                    !counter for pathlength travelled                                              
            call EMIT_PACKET(x,y,z,nx,ny,nz)     !create a packet and direction of travel
            
            call MOVE_PACKET(rmax, taumax, j, nx, ny, nz, x, y, z, L) !move packet in direction of travel

            call UPDATE_PATHLENGTH_FIRST_MOVE(x, y, z, nx, ny, nz, L, rmax, pathlength) !update pathlength based on whether packet is inside or outside sphere
            
            do while ((x**2 + y**2 + z**2) <= (rmax**2)) !whilst packet within sphere
              if (ran < albedo) then                     !condition for scattering
                  call SCATTER(nx, ny, nz)               !generate scattering direction
              else                                       !condition for absorption
                  exit
              end if
              
              call MOVE_PACKET(rmax, taumax, j, nx, ny, nz, x, y, z, L) !scatter packet randomly

              call UPDATE_PATHLENGTH(x, y, z, nx, ny, nz, L, rmax, pathlength, xintersect, yintersect, zintersect) !update pathlength

            end do
            !compute time taken and bin
            time_taken = pathlength/v
            times_taken(i) = time_taken
        end do

         !for this optical depth, bin packets and write to data file
          call BIN(num_bins, npackets, times_taken, bin_edges, freq)

         !write the data to a file
          write(filename, '("tau", I1, ".dat")') j
          open(1, file=filename, status='replace')
          print '("file ", A)', trim(filename)
          do i = 1, num_bins
              write(1, *) bin_edges(i), freq(i)
          end do
          close(1)
        
     end do
         
end program MCRT


