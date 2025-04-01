module pathlengths
    !--------------------------------------------------------------------------------------
    ! *** Module for caluclating pathlength quantities ***
    !--------------------------------------------------------------------------------------
    contains


subroutine FIND_DISTANCE_TO_EDGE(x, y, z, L, nx, ny, nz, xintersect, yintersect, zintersect, distance_travelled)     
    !--------------------------------------------------------------------------------------
    ! *** Calculates distance from current point to where packet exits the sphere     ***
    !--------------------------------------------------------------------------------------
        implicit none
        double precision :: x, y, z, L
        double precision :: nx, ny, nz
        double precision :: x0, y0, z0, a, b, c, t
        double precision :: xintersect, yintersect, zintersect
        double precision :: distance_travelled
    
        !backtrack to previous packet position before crossing sphere boundary
        x0 = x - L * nx
        y0 = y - L * ny
        z0 = z - L * nz
    
        !calculate quadratic components for sphere intersection
        b = 2.d0 * (x0 * nx + y0 * ny + z0 * nz)
        a = nx**2 + ny**2 + nz**2
        c = x0**2 + y0**2 + z0**2 - 1.d0
    
        !solve quadratic to find parameter t to sphere edge
        t = (-b + SQRT(b**2 - 4.d0 * a * c)) / (2.d0 * a)
    
        !calculate exact intersection point at sphere edge
        xintersect = x0 + nx * t
        yintersect = y0 + ny * t
        zintersect = z0 + nz * t
    
        !calculate distance travelled to reach sphere boundary
        distance_travelled = SQRT((nx * t)**2 + (ny * t)**2 + (nz * t)**2)
    
    end subroutine FIND_DISTANCE_TO_EDGE


subroutine UPDATE_PATHLENGTH_FIRST_MOVE(x, y, z, nx, ny, nz, L, rmax, pathlength)
    !--------------------------------------------------------------------------------------
    ! *** Updates the pathlength based on whether the packet is inside or outside the sphere   
    !     on it's first move ***
    !--------------------------------------------------------------------------------------
        implicit none
        double precision, intent(in) :: x, y, z, nx, ny, nz, L, rmax
        double precision, intent(inout) :: pathlength
        double precision :: distance_travelled
    
        if ((x**2 + y**2 + z**2) <= (rmax**2)) then
            distance_travelled = sqrt((L*nx)**2.d0 + (L*ny)**2.d0 + (L*nz)**2.d0)
            pathlength = pathlength + distance_travelled
        else
            pathlength = pathlength + rmax
        end if
    
    end subroutine UPDATE_PATHLENGTH_FIRST_MOVE


subroutine UPDATE_PATHLENGTH(x, y, z, nx, ny, nz, L, rmax, pathlength, xintersect, yintersect, zintersect)
    !--------------------------------------------------------------------------------------
    ! *** Updates pathlength based on whether packet is inside or has exited sphere  ***
    !--------------------------------------------------------------------------------------
        implicit none
        double precision, intent(in) :: x, y, z, nx, ny, nz, L, rmax
        double precision, intent(inout) :: pathlength
        double precision, intent(out) :: xintersect, yintersect, zintersect
        double precision :: distance_travelled
    
        if ((x**2 + y**2 + z**2) <= (rmax**2)) then
            distance_travelled = sqrt((L*nx)**2.d0 + (L*ny)**2.d0 + (L*nz)**2.d0)
            pathlength = pathlength + distance_travelled
        else
            call FIND_DISTANCE_TO_EDGE(x, y, z, L, nx, ny, nz, xintersect, yintersect, zintersect, distance_travelled)
            pathlength = pathlength + distance_travelled
        end if
    
    end subroutine UPDATE_PATHLENGTH
    

end module pathlengths