module procedures_tower_sampler_factory

use classes_tower_sampler, only: Abstract_Tower_Sampler, Concrete_Tower_Sampler, Null_Tower_Sampler

implicit none

private
public :: create, destroy

contains

    subroutine create(tower_sampler, num_candidates, needed)
        class(Abstract_Tower_Sampler), allocatable, intent(out) :: tower_sampler
        integer, intent(in) :: num_candidates
        logical, intent(in) :: needed

        if (needed) then
            allocate(Concrete_Tower_Sampler :: tower_sampler)
        else
            allocate(Null_Tower_Sampler :: tower_sampler)
        end if
        call tower_sampler%construct(num_candidates)
    end subroutine create

    subroutine destroy(tower_sampler)
        class(Abstract_Tower_Sampler), allocatable, intent(inout) :: tower_sampler

        if (allocated(tower_sampler)) then
            call tower_sampler%destroy()
            deallocate(tower_sampler)
        end if
    end subroutine destroy

end module procedures_tower_sampler_factory
