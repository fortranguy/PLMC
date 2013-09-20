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
        
            type is (hardSpheres)
            
                write(*, *) "hard : ", sph%name
        
            type is (interactingSpheres)
            
                write(*, *) "inter : ", sph%name
            
            type is (dipolarSpheres)
            
                write(*, *) "dipol : ", sph%name
        
        end select 

    end subroutine polymorph

end module module_algorithms
