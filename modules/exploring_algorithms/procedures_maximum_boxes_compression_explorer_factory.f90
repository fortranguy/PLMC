module procedures_maximum_boxes_compression_explorer_factory

use types_physical_model_wrapper, only: Physical_Model_Wrapper
use classes_maximum_box_compression, only: Abstract_Maximum_Box_Compression
use classes_maximum_box_compression_explorer, only: Abstract_Maximum_Box_Compression_Explorer, &
    Concrete_Maximum_Box_Compression_Explorer, Null_Maximum_Box_Compression_Explorer

implicit none

private
public :: create, destroy

contains

    subroutine create(maximum_boxes_compression_explorer, physical_model, maximum_box_compression, &
        measure_maximum_box_compression)
        class(Abstract_Maximum_Box_Compression_Explorer), allocatable, intent(out) :: &
            maximum_boxes_compression_explorer(:)
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        class(Abstract_Maximum_Box_Compression), intent(in) :: maximum_box_compression
        logical, intent(in) :: measure_maximum_box_compression

        integer :: i_box

        if (measure_maximum_box_compression) then
            allocate(Concrete_Maximum_Box_Compression_Explorer :: &
                maximum_boxes_compression_explorer(size(physical_model%environment%periodic_boxes)))
        else
            allocate(Null_Maximum_Box_Compression_Explorer :: &
                maximum_boxes_compression_explorer(size(physical_model%environment%periodic_boxes)))
        end if

        do i_box = 1, size(maximum_boxes_compression_explorer)
            call maximum_boxes_compression_explorer(i_box)%&
                construct(physical_model%environment%periodic_boxes(i_box), physical_model%mixture%&
                components(:, i_box), physical_model%short_interactions%components_pairs, &
                physical_model%short_interactions%cells(i_box), maximum_box_compression)
        end do
    end subroutine create

    subroutine destroy(maximum_boxes_compression_explorer)
        class(Abstract_Maximum_Box_Compression_Explorer), allocatable, intent(inout) :: &
            maximum_boxes_compression_explorer(:)

        integer :: i_box

        if (allocated(maximum_boxes_compression_explorer)) then
            do i_box = size(maximum_boxes_compression_explorer), 1, -1
                call maximum_boxes_compression_explorer(i_box)%destroy()
            end do
            deallocate(maximum_boxes_compression_explorer)
        end if
    end subroutine destroy

end module procedures_maximum_boxes_compression_explorer_factory
