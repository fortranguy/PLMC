module procedures_des_convergence_parameter_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_checks, only: check_data_found
use classes_des_convergence_parameter, only: Abstract_DES_Convergence_Parameter, &
    Concrete_DES_Convergence_Parameter, Null_DES_Convergence_Parameter

implicit none

private
public :: create, destroy

contains

    subroutine create(alpha, dipoles_exist, generating_data, prefix)
        class(Abstract_DES_Convergence_Parameter), allocatable, intent(out) :: alpha
        logical, intent(in) :: dipoles_exist
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: alpha_x_box_edge

        if (dipoles_exist) then
            data_field = prefix//"alpha times box edge"
            call generating_data%get(data_field, alpha_x_box_edge, data_found)
            call check_data_found(data_field, data_found)
            allocate(Concrete_DES_Convergence_Parameter :: alpha)
        else
            allocate(Null_DES_Convergence_Parameter :: alpha)
        end if
        call alpha%construct(alpha_x_box_edge)
    end subroutine create

    subroutine destroy(alpha)
        class(Abstract_DES_Convergence_Parameter), allocatable, intent(inout) :: alpha

        if (allocated(alpha)) then
            call alpha%destroy()
            deallocate(alpha)
        end if
    end subroutine destroy

end module procedures_des_convergence_parameter_factory
