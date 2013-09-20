module spheres_changes

use class_spheres
use class_hardSpheres
use class_dipolarSpheres

implicit none

    type, extends(Spheres), public :: Spheres_changes
    
    contains 

    procedure :: polymorph => Spheres_changes_polymorph
    
    end type Spheres_changes
    
    subroutine Spheres_changes_polymorph(this)

        class(Spheres), intent(in) :: this
        
        select type (this)
        
            !class is (InteractingSpheres)
            
            !   write(*, *) "Hello, I'm inter !"
            
            class default  
            
                write(*, *) "Hello world !"
        
        end select

    end subroutine Spheres_changes_polymorph

end module spheres_changes
