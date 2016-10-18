module types_generating_algorithms_wrapper

use classes_generating_algorithm, only: Abstract_Generating_Algorithm

implicit none

private

    type, public :: Generating_Algorithm_Pointer
        class(Abstract_Generating_Algorithm), pointer :: algorithm => null()
    end type Generating_Algorithm_Pointer

    type, public :: Generating_Algorithms_Wrapper
        class(Abstract_Generating_Algorithm), allocatable :: volume_change
        class(Abstract_Generating_Algorithm), allocatable :: one_particle_translation, &
            one_particle_rotation
        class(Abstract_Generating_Algorithm), allocatable :: one_particle_add, one_particle_remove
        class(Abstract_Generating_Algorithm), allocatable :: two_particles_transmutation, two_particles_switch
    end type Generating_Algorithms_Wrapper

end module types_generating_algorithms_wrapper
