module histogram
    !--------------------------------------------------------------------------------------
    ! *** Module for histogram binning logic ***
    !--------------------------------------------------------------------------------------
    contains

subroutine BIN(num_bins, num_values, value_array, bin_edges, freq)      
    !--------------------------------------------------------------------------------------
    ! ***  Bins the values in a given array according to how many bins are wanted    ***
    !--------------------------------------------------------------------------------------
        implicit none
        integer, intent(in) :: num_bins, num_values
        integer :: i, bin_index         
        double precision :: bin_width
        double precision,dimension(num_values), intent(in) :: value_array   
        double precision,dimension(num_bins), intent(out) :: freq, bin_edges 
    
        !find the width of the bins
        bin_width = (MAXVAL(value_array))/dble(num_bins)
    
        !set the first bin to starting at 0 and ending at the bin width
        bin_edges(1) = bin_width
        
        !initialse the rest of the bin widths
        do i = 2, num_bins
            bin_edges(i) = bin_edges(i-1)+bin_width 
        end do 
    
        !set the array which stores the number of values in each bin to be all zeros
        do i = 1, num_bins
           freq(i) = 0.d0
        end do
    
       !then count all the frequencies 
        do i = 1, num_values
            bin_index = min(ceiling(value_array(i) / bin_width), num_bins)
            freq(bin_index) = freq(bin_index) + 1
        end do
      
    end subroutine BIN

end module histogram 