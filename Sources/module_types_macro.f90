!> \brief Definition of derived types (macro)

module module_types_macro

use class_neighbour_cells
use class_small_move
use class_small_rotation
use class_hard_spheres_potential
use class_ewald_summation_real

implicit none

    type :: Hard_Spheres_Macro
        type(Neighbour_Cells) :: same_cells, mix_cells
        type(Small_Move) :: move
        type(Hard_Spheres_Potential) :: hard_potential
    end type Hard_Spheres_Macro
    
    type, public, extends(Hard_Spheres_Macro) :: Dipolar_Hard_Spheres_Macro
        type(Small_Rotation) :: rotation
        type(Ewald_Summation_Real) :: ewald_real      
    end type Dipolar_Hard_Spheres_Macro

end module module_types_macro