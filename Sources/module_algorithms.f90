module module_algorithms

use class_spheres
use class_hardSpheres
use class_interactingSpheres
use class_dipolarSpheres

implicit none

contains

    subroutine polymorph(sph)
    
        class(Spheres), intent(in) :: sph
        
        write(*, *) "sph :", sph%name
    
    end subroutine polymorph

end module module_algorithms
