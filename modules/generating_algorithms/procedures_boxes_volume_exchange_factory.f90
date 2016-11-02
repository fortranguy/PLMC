module procedures_boxes_volume_exchange_factory

use procedures_errors, only: error_exit
use classes_hetero_couples, only: Abstract_Hetero_Couples
use procedures_hetero_couples_factory, only: hetero_couples_create_half => create_half, &
    hetero_couples_destroy => destroy
use classes_tower_sampler, only: Abstract_Tower_Sampler
use procedures_tower_sampler_factory, only:tower_sampler_create => create, tower_sampler_destroy =>&
    destroy
use procedures_mixture_factory, only: set_have_positions
use types_physical_model_wrapper, only: Physical_Model_Wrapper
use classes_exchanged_boxes_size, only: Exchanged_Boxes_Size_Line
use classes_generating_algorithm, only: Abstract_Generating_Algorithm, Null_Generating_Algorithm
use classes_boxes_volume_exchange, only: Boxes_Volume_Exchange

implicit none

private
public :: create

contains

    subroutine create(volume_exchange, physical_model, exchanged_boxes_size)
        class(Abstract_Generating_Algorithm), allocatable, intent(out) :: volume_exchange
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        type(Exchanged_Boxes_Size_Line), intent(in) :: exchanged_boxes_size(:)

        class(Abstract_Hetero_Couples), allocatable :: couples
        class(Abstract_Tower_Sampler), allocatable :: selector
        logical :: have_positions(size(physical_model%mixture%components, 1), &
            size(physical_model%mixture%components, 2))

        call set_have_positions(have_positions, physical_model%mixture%components)
        if (size(have_positions, 2) > 1 .and. any(have_positions)) then
            allocate(Boxes_Volume_Exchange :: volume_exchange)
        else
            allocate(Null_Generating_Algorithm :: volume_exchange)
        end if

        call hetero_couples_create_half(couples, size(have_positions, 2))
        call tower_sampler_create(selector, couples%get_num(), size(have_positions, 2) > 1 .and. &
            any(have_positions))
        select type (volume_exchange)
            type is (Boxes_Volume_Exchange)
                call volume_exchange%construct(physical_model%environment, physical_model%mixture%&
                    components, physical_model%short_interactions, physical_model%&
                    dipolar_interactions_facades, exchanged_boxes_size, have_positions, couples, &
                    selector)
            type is (Null_Generating_Algorithm)
            class default
                call error_exit("procedures_boxes_volume_exchange_factory: create: "//&
                    "volume_exchange: unknown type.")
        end select
        call tower_sampler_destroy(selector)
        call hetero_couples_destroy(couples)
    end subroutine create

end module procedures_boxes_volume_exchange_factory
