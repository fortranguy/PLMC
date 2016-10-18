module procedures_real_writer_factory

use types_string_wrapper, only: String_Wrapper
use procedures_environment_inquirers, only: box_size_can_change
use classes_changed_box_size, only: Changed_Box_Size_Line
use classes_real_writer, only: Abstract_Real_Writer, Concrete_Real_Writer, Null_Real_Writer

implicit none

private
public :: create, destroy

interface create
    module procedure :: create_accessible_domains_size
    module procedure :: create_line
    module procedure :: create_element
end interface create

interface destroy
    module procedure :: destroy_element
    module procedure :: destroy_line
end interface destroy

contains

    subroutine create_accessible_domains_size(domains_size, paths, filename, changed_boxes_size)
        class(Abstract_Real_Writer), allocatable, intent(out) :: domains_size(:)
        type(String_Wrapper), intent(in) :: paths(:)
        character(len=*), intent(in) :: filename
        type(Changed_Box_Size_Line), intent(in) :: changed_boxes_size(:)

        logical :: boxes_size_can_change(size(changed_boxes_size))
        integer :: i_box

        do i_box = 1, size(boxes_size_can_change)
            boxes_size_can_change(i_box) = box_size_can_change(changed_boxes_size(i_box)%&
                line(i_box)%changed)
        end do

        if (all(boxes_size_can_change)) then
            allocate(Concrete_Real_Writer :: domains_size(size(paths)))
        else
            allocate(Null_Real_Writer :: domains_size(size(paths)))
        end if

        do i_box = 1, size(domains_size)
            call domains_size(i_box)%construct(paths(i_box)%string//filename)
        end do
    end subroutine create_accessible_domains_size

    subroutine create_line(writers, paths, filename, needed)
        class(Abstract_Real_Writer), allocatable, intent(out) :: writers(:)
        type(String_Wrapper), intent(in) :: paths(:)
        character(len=*), intent(in) :: filename
        logical, intent(in) :: needed

        integer :: i_box

        if (needed) then
            allocate(Concrete_Real_Writer :: writers(size(paths)))
        else
            allocate(Null_Real_Writer :: writers(size(paths)))
        end if

        do i_box = 1, size(writers)
            call writers(i_box)%construct(paths(i_box)%string//filename)
        end do
    end subroutine create_line

    subroutine destroy_line(writer)
        class(Abstract_Real_Writer), allocatable, intent(inout) :: writer(:)

        integer :: i_box

        if (allocated(writer)) then
            do i_box = size(writer), 1, -1
                call writer(i_box)%destroy()
            end do
            deallocate(writer)
        end if
    end subroutine destroy_line

    subroutine create_element(writer, filename, needed)
        class(Abstract_Real_Writer), allocatable, intent(out) :: writer
        character(len=*), intent(in) :: filename
        logical, intent(in) :: needed

        if (needed) then
            allocate(Concrete_Real_Writer :: writer)
        else
            allocate(Null_Real_Writer :: writer)
        end if
        call writer%construct(filename)
    end subroutine create_element

    subroutine destroy_element(writer)
        class(Abstract_Real_Writer), allocatable, intent(inout) :: writer

        if (allocated(writer)) then
            call writer%destroy()
            deallocate(writer)
        end if
    end subroutine destroy_element

end module procedures_real_writer_factory
