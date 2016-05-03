module procedures_pair_potential_factory

use json_module, only: json_file
use procedures_checks, only: check_data_found
use classes_min_distance, only: Abstract_Min_Distance
use types_potential_domain, only: Short_Potential_Domain
use classes_potential_expression, only: Abstract_Potential_Expression, Null_Potential_Expression
use classes_pair_potential, only: Abstract_Pair_Potential, Tabulated_Pair_Potential, &
    Raw_Pair_Potential, Null_Pair_Potential

implicit none

private
public :: create, destroy

contains

    subroutine create(pair, min_distance, expression, interact, generating_data, prefix)
        class(Abstract_Pair_Potential), allocatable, intent(out) :: pair
        class(Abstract_Min_Distance), intent(in) :: min_distance
        class(Abstract_Potential_Expression), intent(in) :: expression
        logical, intent(in) :: interact
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        call allocate(pair, interact, generating_data, prefix)
        call construct(pair, min_distance, expression, interact, generating_data, prefix)
    end subroutine create

    subroutine allocate(pair, interact, generating_data, prefix)
        class(Abstract_Pair_Potential), allocatable, intent(out) :: pair
        logical, intent(in) :: interact
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found, tabulated_potential

        if (interact) then
            data_field = prefix//"tabulated"
            call generating_data%get(data_field, tabulated_potential, data_found)
            call check_data_found(data_field, data_found)
            if(tabulated_potential) then
                allocate(Tabulated_Pair_Potential :: pair)
            else
                allocate(Raw_Pair_Potential :: pair)
            end if
        else
            allocate(Null_Pair_Potential :: pair)
        end if
    end subroutine allocate

    subroutine construct(pair, min_distance, expression, interact, generating_data, prefix)
        class(Abstract_Pair_Potential), intent(inout) :: pair
        class(Abstract_Min_Distance), intent(in) :: min_distance
        class(Abstract_Potential_Expression), intent(in) :: expression
        logical, intent(in) :: interact
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        type(Short_Potential_Domain) :: domain

        if (interact) then
            domain%min = min_distance%get()
            select type (expression)
                type is (Null_Potential_Expression)
                    domain%max = domain%min
                class default
                    data_field = prefix//"maximum distance"
                    call generating_data%get(data_field, domain%max, data_found)
                    call check_data_found(data_field, data_found)
            end select
            select type (pair)
                type is (Tabulated_Pair_Potential)
                    data_field = prefix//"delta distance"
                    call generating_data%get(data_field, domain%delta, data_found)
                    call check_data_found(data_field, data_found)
            end select
        end if
        call pair%construct(domain, expression)
    end subroutine construct

    subroutine destroy(pair)
        class(Abstract_Pair_Potential), allocatable, intent(inout) :: pair

        if (allocated(pair)) then
            call pair%destroy()
            deallocate(pair)
        end if
    end subroutine destroy

end module procedures_pair_potential_factory
