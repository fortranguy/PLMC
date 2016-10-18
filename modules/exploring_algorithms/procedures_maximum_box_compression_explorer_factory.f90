module procedures_maximum_box_compression_explorer_factory

use types_physical_model_wrapper, only: Physical_Model_Wrapper
use classes_maximum_box_compression, only: Abstract_Maximum_Box_Compression
use classes_maximum_box_compression_explorer, only: Abstract_Maximum_Box_Compression_Explorer, &
    Concrete_Maximum_Box_Compression_Explorer, Null_Maximum_Box_Compression_Explorer

implicit none

private
public :: create, destroy

contains

    subroutine create(maximum_box_compression_explorer, physical_model, maximum_box_compression, &
        measure_maximum_box_compression)
        class(Abstract_Maximum_Box_Compression_Explorer), allocatable, intent(out) :: &
            maximum_box_compression_explorer
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        class(Abstract_Maximum_Box_Compression), intent(in) :: maximum_box_compression
        logical, intent(in) :: measure_maximum_box_compression

        if (measure_maximum_box_compression) then
            allocate(Concrete_Maximum_Box_Compression_Explorer :: maximum_box_compression_explorer)
        else
            allocate(Null_Maximum_Box_Compression_Explorer :: maximum_box_compression_explorer)
        end if
        call maximum_box_compression_explorer%construct(physical_model%environment, physical_model%&
            mixture%gemc_components, physical_model%short_interactions, maximum_box_compression)
    end subroutine create

    subroutine destroy(maximum_box_compression_explorer)
        class(Abstract_Maximum_Box_Compression_Explorer), allocatable, intent(inout) :: &
            maximum_box_compression_explorer

        if (allocated(maximum_box_compression_explorer)) then
            call maximum_box_compression_explorer%destroy()
            deallocate(maximum_box_compression_explorer)
        end if
    end subroutine destroy

end module procedures_maximum_box_compression_explorer_factory
