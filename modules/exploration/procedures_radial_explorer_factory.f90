module procedures_radial_explorer_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_input_prefixes, only: radial_prefix
use classes_number_to_string, only: Concrete_Number_to_String
use json_module, only: json_file
use procedures_errors, only: error_exit
use procedures_mixture_inquirers, only: property_ij_components => ij_components
use types_component_coordinates_reader_selector, only: Component_Coordinates_Reader_Selector
use classes_radial_explorer, only: Abstract_Radial_Explorer, Intra_Radial_Explorer, &
    Inter_Radial_Explorer
use procedures_checks, only: check_data_found

implicit none

private
public :: create, destroy

contains

    subroutine create(radial_explorer, max_box_size, num_components, exploring_data)
        class(Abstract_Radial_Explorer), allocatable, intent(out) :: radial_explorer
        real(DP), intent(in) :: max_box_size(:)
        integer, intent(in) :: num_components
        type(json_file) :: exploring_data

        type(Component_Coordinates_Reader_Selector) :: selector
        integer :: ij_components(2)
        real(DP) :: delta_distance
        character(len=:), allocatable :: filename, data_field
        logical :: data_found
        type(Concrete_Number_to_String) :: string

        selector%read_positions = .true.
        selector%read_orientations = .false.
        ij_components = property_ij_components(exploring_data, radial_prefix)
        filename = "radial_"//string%get(ij_components(1))//"-"//string%get(ij_components(2))//&
            ".out"
        data_field = radial_prefix//"delta"
        call exploring_data%get(data_field, delta_distance, data_found)
        call check_data_found(data_field, data_found)

        if (ij_components(1) == ij_components(2)) then
            allocate(Intra_Radial_Explorer :: radial_explorer)
        else
            allocate(Inter_Radial_Explorer :: radial_explorer)
        end if

        select type (radial_explorer)
            type is (Intra_Radial_Explorer)
                call radial_explorer%construct(num_components, ij_components(1), selector, &
                    minval(max_box_size) / 2._DP, delta_distance, filename)
            type is (Inter_Radial_Explorer)
                call radial_explorer%construct(num_components, ij_components, selector, &
                    minval(max_box_size) / 2._DP, delta_distance, filename)
            class default
                call error_exit("procedures_radial_explorer_factory: create: type unknown.")
        end select
    end subroutine create

    subroutine destroy(radial_explorer)
        class(Abstract_Radial_Explorer), allocatable, intent(inout) :: radial_explorer

        if (allocated(radial_explorer)) then
            call radial_explorer%destroy()
            deallocate(radial_explorer)
        end if
    end subroutine destroy

end module procedures_radial_explorer_factory
