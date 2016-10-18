module procedures_generating_algorithms_factory

use procedures_errors, only: error_exit
use types_physical_model_wrapper, only: Physical_Model_Wrapper
use types_changes_wrapper, only: Changes_Wrapper
use classes_generating_algorithm, only: Abstract_Generating_Algorithm, Generating_Algorithm_Wrapper
use procedures_box_volume_change_factory, only: box_volume_change_create => create
use procedures_box_particle_move_factory, only: box_particle_translation_create => &
    create_translation, box_particle_rotation_create => create_rotation
use procedures_box_particle_exchange_factory, only: box_particle_add_create => create_add, &
    box_particle_remove_create => create_remove
use procedures_box_particles_swap_factory, only: box_particles_transmutation_create => &
    create_transmutation, box_particles_switch_create => create_switch

implicit none

private
public :: create, destroy

contains

    subroutine create(generating_algorithms, physical_model, changes)
        type(Generating_Algorithm_Wrapper), intent(out) :: generating_algorithms(:)
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        type(Changes_Wrapper), intent(in) :: changes

        if (size(generating_algorithms) /= 7) then
            call error_exit("procedures_generating_algorithms_factory: create: "//&
                "size(generating_algorithms) must be 7.")
        end if

        call box_volume_change_create(generating_algorithms(1)%algorithm, physical_model, &
            changes%changed_boxes_size)
        call box_particle_translation_create(generating_algorithms(2)%algorithm, &
            physical_model, changes%gemc_components)
        call box_particle_rotation_create(generating_algorithms(3)%algorithm, &
            physical_model, changes%gemc_components)
        call box_particle_add_create(generating_algorithms(4)%algorithm, physical_model, &
            changes)
        call box_particle_remove_create(generating_algorithms(5)%algorithm, physical_model, &
            changes)
        call box_particles_transmutation_create(generating_algorithms(6)%algorithm, &
            physical_model, changes)
        call box_particles_switch_create(generating_algorithms(7)%algorithm, &
            physical_model, changes)
    end subroutine create

    subroutine destroy(generating_algorithms)
        type(Generating_Algorithm_Wrapper), intent(inout) :: generating_algorithms(:)

        integer :: i_algorithm

        do i_algorithm = size(generating_algorithms), 1, -1
            if (allocated(generating_algorithms(i_algorithm)%algorithm)) then
                call generating_algorithms(i_algorithm)%algorithm%destroy()
                deallocate(generating_algorithms(i_algorithm)%algorithm)
            end if
        end do
    end subroutine destroy

end module procedures_generating_algorithms_factory
