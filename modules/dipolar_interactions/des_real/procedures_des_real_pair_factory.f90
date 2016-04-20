module procedures_des_real_pair_factory

use json_module, only: json_file
use procedures_checks, only: check_data_found
use types_potential_domain, only: Dipolar_Potential_Domain
use classes_des_real_pair, only: Abstract_DES_Real_Pair, Tabulated_DES_Real_Pair, &
    Raw_DES_Real_Pair, Null_DES_Real_Pair
use classes_permittivity, only: Abstract_Permittivity
use classes_minimum_distance, only: Abstract_Minimum_Distance
use classes_des_convergence_parameter, only: Abstract_DES_Convergence_Parameter

implicit none

private
public :: des_real_pair_create, des_real_pair_destroy

contains

    subroutine des_real_pair_create(pair, permittivity, min_distance, interact, alpha, input_data, &
        prefix)
        class(Abstract_DES_Real_Pair), allocatable, intent(out) :: pair
        class(Abstract_Permittivity), intent(in) :: permittivity
        class(Abstract_Minimum_Distance), intent(in) :: min_distance
        logical, intent(in) :: interact
        class(Abstract_DES_Convergence_Parameter), intent(in) :: alpha
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        call allocate_pair(pair, interact, input_data, prefix)
        call construct_pair(pair, permittivity, min_distance, interact, alpha, input_data, prefix)
    end subroutine des_real_pair_create

    subroutine allocate_pair(pair, interact, input_data, prefix)
        class(Abstract_DES_Real_Pair), allocatable, intent(out) :: pair
        logical, intent(in) :: interact
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found, tabulated_potential

        if (interact) then
            data_field = prefix//"tabulated"
            call input_data%get(data_field, tabulated_potential, data_found)
            call check_data_found(data_field, data_found)
            if(tabulated_potential) then
                allocate(Tabulated_DES_Real_Pair :: pair)
            else
                allocate(Raw_DES_Real_Pair :: pair)
            end if
        else
            allocate(Null_DES_Real_Pair :: pair)
        end if
    end subroutine allocate_pair

    subroutine construct_pair(pair, permittivity, min_distance, interact, alpha,  input_data, &
        prefix)
        class(Abstract_DES_Real_Pair), intent(inout) :: pair
        class(Abstract_Permittivity), intent(in) :: permittivity
        class(Abstract_Minimum_Distance), intent(in) :: min_distance
        logical, intent(in) :: interact
        class(Abstract_DES_Convergence_Parameter), intent(in) :: alpha
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        type(Dipolar_Potential_Domain) :: domain

        if (interact) then
            domain%min = min_distance%get()
            data_field = prefix//"max distance over box edge"
            call input_data%get(data_field, domain%max_over_box, data_found)
            call check_data_found(data_field, data_found)
            select type (pair)
                type is (Tabulated_DES_Real_Pair)
                    data_field = prefix//"delta distance"
                    call input_data%get(data_field, domain%delta, data_found)
                    call check_data_found(data_field, data_found)
            end select
        end if
        call pair%construct(permittivity, alpha, domain)
    end subroutine construct_pair

    subroutine des_real_pair_destroy(pair)
        class(Abstract_DES_Real_Pair), allocatable, intent(inout) :: pair

        if (allocated(pair)) then
            call pair%destroy()
            deallocate(pair)
        end if
    end subroutine des_real_pair_destroy

end module procedures_des_real_pair_factory
