module types_markov_chain_explorer_wrapper

use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain
use classes_maximum_box_compression_explorer, only: Abstract_Maximum_Box_Compression_Explorer
use classes_random_coordinates, only: Abstract_Random_Coordinates
use classes_particle_insertion_method, only: Abstract_Particle_Insertion_Method
use classes_volume_change_method, only: Abstract_Volume_Change_Method

implicit none

private

    type, public :: Markov_Chain_Explorer_Wrapper
        class(Abstract_Maximum_Box_Compression_Explorer), allocatable :: &
            maximum_box_compression_explorer
        class(Abstract_Volume_Change_Method), allocatable :: volume_change_method
        class(Abstract_Parallelepiped_Domain), allocatable :: particle_insertion_domain
        class(Abstract_Random_Coordinates), allocatable :: random_position, random_orientation
        class(Abstract_Particle_Insertion_Method), allocatable :: particle_insertion_method
    end type Markov_Chain_Explorer_Wrapper

end module types_markov_chain_explorer_wrapper
