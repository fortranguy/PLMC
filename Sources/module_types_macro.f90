!> \brief Definition of derived types (macro)

module module_types_macro

use class_neighbour_cells, only: Neighbour_Cells
use class_small_move, only: Small_Move
use class_small_rotation, only: Small_Rotation
use class_hard_spheres_potential_energy, only: Hard_Spheres_Potential_Energy
use class_ewald_summation_real, only: Ewald_Summation_Real
use class_ewald_summation_reci, only: Ewald_Summation_Reci
use class_ewald_summation_self, only: Ewald_Summation_Self
use class_ewald_summation_bound, only: Ewald_Summation_Bound

implicit none

    type, public :: Hard_Spheres_Macro
        type(Neighbour_Cells) :: same_cells, mix_cells
        type(Small_Move) :: move
        type(Hard_Spheres_Potential_Energy) :: hard_potential_energy
    end type Hard_Spheres_Macro
    
    type, public, extends(Hard_Spheres_Macro) :: Dipolar_Hard_Spheres_Macro
        type(Small_Rotation) :: rotation
        type(Ewald_Summation_Real) :: ewald_real
        type(Ewald_Summation_Reci) :: ewald_reci  
        type(Ewald_Summation_Self) :: ewald_self
        type(Ewald_Summation_Bound) :: ewald_bound
    end type Dipolar_Hard_Spheres_Macro

end module module_types_macro
