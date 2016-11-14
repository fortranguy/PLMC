module types_markov_chain_explorer_wrapper

use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain
use classes_changed_box_size_ratio, only: Abstract_Changed_Box_Size_Ratio
use classes_maximum_box_compression_explorer, only: Abstract_Maximum_Box_Compression_Explorer
use classes_volume_change_method, only: Abstract_Volume_Change_Method
use classes_particle_insertion_method, only: Abstract_Particle_Insertion_Method
use classes_dipolar_neighbourhoods_visitor, only: Abstract_Dipolar_Neighbourhoods_Visitor

implicit none

private

    type, public :: Markov_Chain_Explorer_Wrapper
        class(Abstract_Maximum_Box_Compression_Explorer), allocatable :: &
            maximum_boxes_compression_explorer(:)
        class(Abstract_Changed_Box_Size_Ratio), allocatable :: changed_boxes_size_ratio(:)
        class(Abstract_Volume_Change_Method), allocatable :: volume_change_method
        class(Abstract_Parallelepiped_Domain), allocatable :: particle_insertion_domains(:)
        class(Abstract_Particle_Insertion_Method), allocatable :: particle_insertion_method
        class(Abstract_Dipolar_Neighbourhoods_Visitor), allocatable :: &
            dipolar_neighbourhoods_visitors(:)
    end type Markov_Chain_Explorer_Wrapper

end module types_markov_chain_explorer_wrapper
