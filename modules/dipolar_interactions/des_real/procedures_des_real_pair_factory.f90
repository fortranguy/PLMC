module procedures_des_real_pair_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_checks, only: check_data_found
use types_potential_domain, only: Dipolar_Potential_Domain
use classes_box_size_memento, only: Abstract_Box_Size_Memento
use classes_permittivity, only: Abstract_Permittivity
use types_min_distance_wrapper, only: Min_Distances_Line
use classes_des_convergence_parameter, only: Abstract_DES_Convergence_Parameter
use classes_des_real_pair, only: Abstract_DES_Real_Pair, Tabulated_DES_Real_Pair, &
    Raw_DES_Real_Pair, Null_DES_Real_Pair

implicit none

private
public :: create, destroy

contains

    subroutine create(pair, box_size_memento, permittivity, min_distances, dipoles_exist, alpha, &
        generating_data, prefix)
        class(Abstract_DES_Real_Pair), allocatable, intent(out) :: pair
        class(Abstract_Box_Size_Memento), intent(in) :: box_size_memento
        class(Abstract_Permittivity), intent(in) :: permittivity
        type(Min_Distances_Line), intent(in) :: min_distances(:)
        logical, intent(in) :: dipoles_exist
        class(Abstract_DES_Convergence_Parameter), intent(in) :: alpha
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        real(DP) :: min_distance
        integer :: i_component, j_component

        min_distance = min_distances(1)%line(1)%distance%get()
        do j_component = 1, size(min_distances)
            do i_component = 1, size(min_distances(j_component)%line)
                associate(min_distance_ij => min_distances(j_component)%line(i_component)%distance%&
                    get())
                    if (min_distance > min_distance_ij) min_distance = min_distance_ij
                end associate
            end do
        end do

        call allocate(pair, dipoles_exist, generating_data, prefix)
        call construct(pair, box_size_memento, permittivity, min_distance, dipoles_exist, alpha, &
            generating_data, prefix)
    end subroutine create

    subroutine allocate(pair, dipoles_exist, generating_data, prefix)
        class(Abstract_DES_Real_Pair), allocatable, intent(out) :: pair
        logical, intent(in) :: dipoles_exist
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found, tabulated_potential

        if (dipoles_exist) then
            data_field = prefix//"tabulated"
            call generating_data%get(data_field, tabulated_potential, data_found)
            call check_data_found(data_field, data_found)
            if(tabulated_potential) then
                allocate(Tabulated_DES_Real_Pair :: pair)
            else
                allocate(Raw_DES_Real_Pair :: pair)
            end if
        else
            allocate(Null_DES_Real_Pair :: pair)
        end if
    end subroutine allocate

    subroutine construct(pair, box_size_memento, permittivity, min_distance, dipoles_exist, alpha, &
        generating_data, prefix)
        class(Abstract_DES_Real_Pair), intent(inout) :: pair
        class(Abstract_Box_Size_Memento), intent(in) :: box_size_memento
        class(Abstract_Permittivity), intent(in) :: permittivity
        real(DP), intent(in) :: min_distance
        logical, intent(in) :: dipoles_exist
        class(Abstract_DES_Convergence_Parameter), intent(in) :: alpha
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        type(Dipolar_Potential_Domain) :: domain

        if (dipoles_exist) then
            domain%min = min_distance
            data_field = prefix//"max distance / box edge"
            call generating_data%get(data_field, domain%max_over_box, data_found)
            call check_data_found(data_field, data_found)
            select type (pair)
                type is (Tabulated_DES_Real_Pair)
                    data_field = prefix//"delta distance"
                    call generating_data%get(data_field, domain%delta, data_found)
                    call check_data_found(data_field, data_found)
            end select
        end if
        call pair%construct(box_size_memento, permittivity, alpha, domain)
    end subroutine construct

    subroutine destroy(pair)
        class(Abstract_DES_Real_Pair), allocatable, intent(inout) :: pair

        if (allocated(pair)) then
            call pair%destroy()
            deallocate(pair)
        end if
    end subroutine destroy

end module procedures_des_real_pair_factory
