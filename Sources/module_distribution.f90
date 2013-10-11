!> \brief Distribution module

module module_distribution

use data_precisions, only : DP
use data_constants, only : PI
use data_distribution, only : dist_dr

implicit none

private
public sphereVol

contains

    !> Calculate the volume of the sphere
    
    pure function sphereVol(iDist)
    
        integer, intent(in) :: iDist    
        real(DP) :: sphereVol
        
        sphereVol = 4._DP/3._DP * PI * (real(iDist, DP)*dist_dr)**3
        
    end function sphereVol
    
end module module_distribution
