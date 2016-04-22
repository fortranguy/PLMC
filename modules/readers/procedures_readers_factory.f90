module procedures_readers_factory

use classes_periodic_box, only: Abstract_Periodic_Box
use types_component_wrapper, only: Component_Wrapper
use procedures_box_size_reader_factory, only: box_size_reader_create => create, &
    box_size_reader_destroy => destroy
use procedures_coordinates_reader_factory, only: coordinates_reader_create => create, &
    coordinates_reader_destroy => destroy
use types_readers_wrapper, only: Readers_Wrapper

implicit none

private
public :: readers_create, readers_destroy

contains

    subroutine readers_create(readers, periodic_box, components)
        type(Readers_Wrapper), intent(out) :: readers
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        type(Component_Wrapper), intent(in) :: components(:)

        call box_size_reader_create(readers%box_size, periodic_box)
        call coordinates_reader_create(readers%components, components)
    end subroutine readers_create

    subroutine readers_destroy(readers)
        type(Readers_Wrapper), intent(inout) :: readers

        call coordinates_reader_destroy(readers%components)
        call box_size_reader_destroy(readers%box_size)
    end subroutine readers_destroy

end module procedures_readers_factory
