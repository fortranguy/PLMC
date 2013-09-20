module module_algorithms

use class_spheres
use class_hardSpheres
use class_interactingSpheres
use class_dipolarSpheres

implicit none

contains

    subroutine polymorph(sph)
    
        class(Spheres), intent(in) :: sph
        
        select type (sph)
        
            class is (hardSpheres)
            
                write(*, *) "hard : ", sph%name
        
            class is (interactingSpheres)
            
                write(*, *) "inter : ", sph%name
            
            class is (dipolarSpheres)
            
                write(*, *) "dipol : ", sph%name
        
        end select 

    end subroutine polymorph

end module module_algorithms
