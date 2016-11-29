module procedures_generating_algorithms_factory

use data_strings, only: max_word_length
use json_module, only: json_core, json_value
use procedures_errors, only: error_exit
use types_physical_model_wrapper, only: Physical_Model_Wrapper
use types_changes_wrapper, only: Changes_Wrapper
use classes_generating_algorithm, only: Abstract_Generating_Algorithm, Null_Generating_Algorithm, &
    Generating_Algorithm_Wrapper
use procedures_box_volume_change_factory, only: box_volume_change_create => create
use procedures_boxes_volume_exchange_factory, only: boxes_volume_exchange_create => create
use procedures_boxes_particle_teleportation_factory, only: boxes_particle_teleportation_create => &
    create
use procedures_boxes_particles_swap_factory, only: boxes_particles_swap_create => create
use procedures_box_particle_move_factory, only: box_particle_translation_create => &
    create_translation, box_particle_rotation_create => create_rotation
use procedures_box_particle_exchange_factory, only: box_particle_add_create => create_add, &
    box_particle_remove_create => create_remove
use procedures_box_particles_swap_factory, only: box_particles_transmutation_create => &
    create_transmutation, box_particles_switch_create => create_switch

implicit none

private
public :: create, destroy, add_to_report

contains

    subroutine create(generating_algorithms, physical_model, changes)
        type(Generating_Algorithm_Wrapper), intent(out) :: generating_algorithms(:)
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        type(Changes_Wrapper), intent(in) :: changes

        if (size(generating_algorithms) /= 10) then
            call error_exit("procedures_generating_algorithms_factory: create: "//&
                "the number of generating algorithms is wrong.")
        end if

        call box_volume_change_create(generating_algorithms(1)%algorithm, physical_model, &
            changes%changed_boxes_size)
        call boxes_volume_exchange_create(generating_algorithms(2)%algorithm, physical_model, &
            changes%exchanged_boxes_size)
        call boxes_particle_teleportation_create(generating_algorithms(3)%algorithm, &
            physical_model, changes)
        call boxes_particles_swap_create(generating_algorithms(4)%algorithm, physical_model, &
            changes)
        call box_particle_translation_create(generating_algorithms(5)%algorithm, &
            physical_model, changes%components)
        call box_particle_rotation_create(generating_algorithms(6)%algorithm, &
            physical_model, changes%components)
        call box_particle_add_create(generating_algorithms(7)%algorithm, physical_model, &
            changes)
        call box_particle_remove_create(generating_algorithms(8)%algorithm, physical_model, &
            changes)
        call box_particles_transmutation_create(generating_algorithms(9)%algorithm, &
            physical_model, changes)
        call box_particles_switch_create(generating_algorithms(10)%algorithm, &
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

    !> @note The names must coherent with [[procedures_generating_algorithms_factory:create]].
    !> @todo Find a better alternative?
    subroutine add_to_report(json, algorithms_weight, generating_algorithms)
        type(json_core), intent(inout) :: json
        type(json_value), intent(inout), pointer :: algorithms_weight
        type(Generating_Algorithm_Wrapper), intent(in) :: generating_algorithms(:)

        character(len=max_word_length) :: algorithms_name(size(generating_algorithms))
        integer :: i_algorithm, weight_i

        algorithms_name(1) = "box volume change"
        algorithms_name(2) = "boxes volume exchange"
        algorithms_name(3) = "boxes particle teleportation"
        algorithms_name(4) = "boxes particles swap"
        algorithms_name(5) = "box particle translation"
        algorithms_name(6) = "box particle rotation"
        algorithms_name(7) = "box particle add"
        algorithms_name(8) = "box particle remove"
        algorithms_name(9) = "box particles transmutation"
        algorithms_name(10) = "box particles switch"

        do i_algorithm = 1, size(generating_algorithms)
            weight_i = generating_algorithms(i_algorithm)%algorithm%get_num_choices()
            if (weight_i == 0) cycle
            call json%add(algorithms_weight, algorithms_name(i_algorithm), weight_i)
        end do
    end subroutine add_to_report

end module procedures_generating_algorithms_factory
