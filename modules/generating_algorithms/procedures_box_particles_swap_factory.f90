module procedures_box_particles_swap_factory

use procedures_errors, only: error_exit
use classes_hetero_couples, only: Abstract_Hetero_Couples
use procedures_hetero_couples_factory, only: hetero_couples_create_full => create_full, &
    hetero_couples_create_half => create_half, hetero_couples_destroy => destroy
use classes_tower_sampler, only: Abstract_Tower_Sampler
use procedures_tower_sampler_factory, only:tower_sampler_create => create, tower_sampler_destroy =>&
    destroy
use types_physical_model_wrapper, only: Physical_Model_Wrapper
use types_changes_wrapper, only: Changes_Wrapper
use procedures_changes_factory, only: set_can_translate
use classes_generating_algorithm, only: Abstract_Generating_Algorithm, Null_Generating_Algorithm
use classes_box_particles_swap, only: Abstract_Box_Particles_Swap, Box_Particles_Transmutation, &
    Box_Particles_Switch
use procedures_changes_factory, only: set_can_exchange

implicit none

private
public :: create_transmutation, create_switch

contains

    subroutine create_transmutation(particles_transmutation, physical_model, changes)
        class(Abstract_Generating_Algorithm), allocatable, intent(out) :: particles_transmutation
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        type(Changes_Wrapper), intent(in) :: changes

        class(Abstract_Hetero_Couples), allocatable :: couples(:)
        class(Abstract_Tower_Sampler), allocatable :: selectors(:)
        logical :: can_exchange(size(physical_model%mixture%gemc_components, 1), &
            size(physical_model%mixture%gemc_components, 2))

        call set_can_exchange(can_exchange, physical_model%mixture%gemc_components)
        if (count(can_exchange(:, 1)) > 1) then
            allocate(Box_Particles_Transmutation :: particles_transmutation)
        else
            allocate(Null_Generating_Algorithm :: particles_transmutation)
        end if

        call hetero_couples_create_full(couples, size(can_exchange, 2), size(can_exchange, 1))
        call tower_sampler_create(selectors, size(couples), couples(1)%get_num(), count(can_exchange(:, 1)) > 1)
        select type (particles_transmutation)
            type is (Box_Particles_Transmutation)
                call particles_transmutation%construct(physical_model%environment, physical_model%&
                    mixture, physical_model%short_interactions, physical_model%&
                    gemc_dipolar_interactions_dynamic, physical_model%gemc_dipolar_interactions_static, &
                    changes, can_exchange, couples, selectors)
            type is (Null_Generating_Algorithm)
            class default
                call error_exit("procedures_box_particles_swap_factory: create_transmutation: "//&
                    "particles_transmutation: unknown type.")
        end select
        call tower_sampler_destroy(selectors)
        call hetero_couples_destroy(couples)
    end subroutine create_transmutation

    subroutine create_switch(particles_switch, physical_model, changes)
        class(Abstract_Generating_Algorithm), allocatable, intent(out) :: particles_switch
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        type(Changes_Wrapper), intent(in) :: changes

        class(Abstract_Hetero_Couples), allocatable :: couples(:)
        class(Abstract_Tower_Sampler), allocatable :: selectors(:)
        logical :: can_translate(size(changes%gemc_components, 1), size(changes%gemc_components, 2))

        call set_can_translate(can_translate, changes%gemc_components)
        if (count(can_translate(:, 1)) > 1) then
            allocate(Box_Particles_Switch :: particles_switch)
        else
            allocate(Null_Generating_Algorithm :: particles_switch)
        end if

        call hetero_couples_create_half(couples, size(can_translate, 2), size(can_translate, 1))
        call tower_sampler_create(selectors, size(couples), couples(1)%get_num(), count(can_translate(:, 1)) > 1)
        select type (particles_switch)
            type is (Box_Particles_Switch)
                call particles_switch%construct(physical_model%environment, physical_model%mixture,&
                    physical_model%short_interactions, physical_model%gemc_dipolar_interactions_dynamic,&
                    physical_model%gemc_dipolar_interactions_static, changes, can_translate, couples, &
                    selectors)
            type is (Null_Generating_Algorithm)
            class default
                call error_exit("procedures_box_particles_swap_factory: create_switch: "//&
                    "particles_switch: unknown type.")
        end select
        call tower_sampler_destroy(selectors)
        call hetero_couples_destroy(couples)
    end subroutine create_switch

end module procedures_box_particles_swap_factory
