module procedures_tower_sampler_factory

use classes_tower_sampler, only: Abstract_Tower_Sampler, Concrete_Tower_Sampler, Null_Tower_Sampler

implicit none

private
public :: create, destroy

interface create
    module procedure :: create_line
    module procedure :: create_element
end interface create

interface destroy
    module procedure :: destroy_element
    module procedure :: destroy_line
end interface destroy

contains

    subroutine create_line(tower_samplers, num_elements, num_candidates, needed)
        class(Abstract_Tower_Sampler), allocatable, intent(out) :: tower_samplers(:)
        integer, intent(in) :: num_elements, num_candidates
        logical, intent(in) :: needed

        integer :: i_element

        if (needed) then
            allocate(Concrete_Tower_Sampler :: tower_samplers(num_elements))
        else
            allocate(Null_Tower_Sampler :: tower_samplers(num_elements))
        end if

        do i_element = 1, size(tower_samplers)
            call tower_samplers(i_element)%construct(num_candidates)
        end do
    end subroutine create_line

    subroutine destroy_line(tower_samplers)
        class(Abstract_Tower_Sampler), allocatable, intent(inout) :: tower_samplers(:)

        integer :: i_element

        if (allocated(tower_samplers)) then
            do i_element = size(tower_samplers), 1, -1
                call tower_samplers(i_element)%destroy()
            end do
            deallocate(tower_samplers)
        end if
    end subroutine destroy_line

    subroutine create_element(tower_sampler, num_candidates, needed)
        class(Abstract_Tower_Sampler), allocatable, intent(out) :: tower_sampler
        integer, intent(in) :: num_candidates
        logical, intent(in) :: needed

        if (needed) then
            allocate(Concrete_Tower_Sampler :: tower_sampler)
        else
            allocate(Null_Tower_Sampler :: tower_sampler)
        end if
        call tower_sampler%construct(num_candidates)
    end subroutine create_element

    subroutine destroy_element(tower_sampler)
        class(Abstract_Tower_Sampler), allocatable, intent(inout) :: tower_sampler

        if (allocated(tower_sampler)) then
            call tower_sampler%destroy()
            deallocate(tower_sampler)
        end if
    end subroutine destroy_element

end module procedures_tower_sampler_factory
