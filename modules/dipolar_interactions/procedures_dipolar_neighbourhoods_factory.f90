module procedures_dipolar_neighbourhoods_factory

use iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_errors, only: warning_continue
use procedures_checks, only: check_data_found
use types_component_wrapper, only: Component_Wrapper
use procedures_mixture_total_moments_factory, only: set_are_dipolar
use classes_dipoles_neighbourhood, only: Concrete_Dipolar_Neighbourhood, &
    Null_Dipolar_Neighbourhood, Dipolar_Neighbourhood_Line

implicit none

private
public :: create, destroy

contains

    subroutine create(neighbourhoods, components, needed, exploring_data, prefix)
        type(Dipolar_Neighbourhood_Line), allocatable, intent(out) :: neighbourhoods(:)
        type(Component_Wrapper), intent(in) :: components(:, :)
        logical, intent(in) :: needed
        type(json_file), intent(inout) :: exploring_data
        character(len=*), intent(in) :: prefix

        logical :: are_dipolar(size(components, 1), size(components, 2))
        real(DP) :: max_distance
        character(len=:), allocatable :: data_field
        logical :: data_found
        integer :: i_box, i_component, j_component

        call set_are_dipolar(are_dipolar, components)
        if (.not.any(are_dipolar) .and. needed) then
            call warning_continue("procedures_dipolar_neighbourhoods_factory: "//&
                "there are no dipoles.")
        end if

        i_box = 1
        allocate(neighbourhoods(size(components, 1)))
        do j_component = 1, size(neighbourhoods)
            allocate(neighbourhoods(j_component)%line(j_component))
            do i_component = 1, size(neighbourhoods(j_component)%line)
                if (needed .and. are_dipolar(i_component, i_box) .and. &
                    are_dipolar(j_component, i_box)) then
                    allocate(Concrete_Dipolar_Neighbourhood :: neighbourhoods(j_component)%&
                        line(i_component)%neighbourhood)
                    data_field = prefix//"maximum distance"
                    call exploring_data%get(data_field, max_distance, data_found)
                    call check_data_found(data_field, data_found)
                else
                    allocate(Null_Dipolar_Neighbourhood :: neighbourhoods(j_component)%&
                        line(i_component)%neighbourhood)
                    max_distance = 0._DP
                end if
                call neighbourhoods(j_component)%line(i_component)%neighbourhood%set(max_distance)
            end do
        end do
    end subroutine create

    subroutine destroy(neighbourhoods)
        type(Dipolar_Neighbourhood_Line), allocatable, intent(inout) :: neighbourhoods(:)

        integer :: i_component, j_component

        if (allocated(neighbourhoods)) then
            do j_component = size(neighbourhoods), 1, -1
                if (allocated(neighbourhoods(j_component)%line)) then
                    do i_component = size(neighbourhoods(j_component)%line), 1, -1
                        deallocate(neighbourhoods(j_component)%line(i_component)%neighbourhood)
                    end do
                end if
                deallocate(neighbourhoods(j_component)%line)
            end do
            deallocate(neighbourhoods)
        end if
    end subroutine destroy

end module procedures_dipolar_neighbourhoods_factory
