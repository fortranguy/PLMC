!> \brief Distribution module

module module_distribution

use data_precisions, only : DP
use data_constants, only : PI
use data_distribution, only : dist_dr

implicit none

private
public sphere_volume

contains

    !> Calculate the volume of the sphere
    
    pure function sphere_volume(iDist)
    
        integer, intent(in) :: iDist    
        real(DP) :: sphere_volume
        
        sphere_volume = 4._DP/3._DP * PI * (real(iDist, DP)*dist_dr)**3
        
    end function sphere_volume
    
end module module_distribution
