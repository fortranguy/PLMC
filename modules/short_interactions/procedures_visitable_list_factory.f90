module procedures_visitable_list_factory

use json_module, only: json_file
use procedures_errors, only: error_exit
use procedures_checks, only: check_data_found
use classes_visitable_list, only: Abstract_Visitable_List, Concrete_Visitable_List, &
    Concrete_Visitable_Array, Null_Visitable_List

implicit none

private
public :: allocate, deallocate

contains

    subroutine allocate(list, interact, generating_data, prefix)
        class(Abstract_Visitable_List), allocatable, intent(out) :: list
        logical, intent(in) :: interact
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field, cells_data_structure
        logical :: data_found

        if (interact) then
            data_field = prefix//"Cells.data structure"
            call generating_data%get(data_field, cells_data_structure, data_found)
            call check_data_found(data_field, data_found)
            select case (cells_data_structure)
                case ("list")
                    allocate(Concrete_Visitable_List :: list)
                case ("array")
                    allocate(Concrete_Visitable_Array :: list)
                case default
                    call error_exit(cells_data_structure//" unknown."&
                        //"Choose between 'list' and 'array'.")
            end select
        else
            allocate(Null_Visitable_List :: list)
        end if
    end subroutine allocate

    subroutine deallocate(list)
        class(Abstract_Visitable_List), allocatable, intent(inout) :: list

        if (allocated(list)) deallocate(list)
    end subroutine deallocate

end module procedures_visitable_list_factory
