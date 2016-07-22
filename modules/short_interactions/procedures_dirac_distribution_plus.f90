module procedures_dirac_distribution_plus_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_errors, only: error_exit
use procedures_checks, only: check_data_found
use classes_dirac_distribution_plus, only: Abstract_Dirac_Distribution_Plus, &
    Rectangular_Dirac_Distribution_Plus, Gaussian_Dirac_Distribution_Plus, &
    Null_Dirac_Distribution_Plus

implicit none

private
public :: create, destroy

contains

    subroutine create(dirac_plus, needed, exploring_data, prefix)
        class(Abstract_Dirac_Distribution_Plus), allocatable, intent(out) :: dirac_plus
        logical, intent(in) :: needed
        type(json_file), intent(inout) :: exploring_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: distribution_name, data_field
        logical :: data_found

        if (needed) then
            data_field = prefix//"Dirac Distribution Plus.name"
            call exploring_data%get(data_field, distribution_name, data_found)
            call check_data_found(data_field, data_found)
            select case (distribution_name)
                case ("rectangular")
                    allocate(Rectangular_Dirac_Distribution_Plus :: dirac_plus)
                case ("gaussian")
                    allocate(Gaussian_Dirac_Distribution_Plus :: dirac_plus)
                case default
                    call error_exit("procedures_dirac_distribution_plus_factory: create: "//&
                        distribution_name//" distribution_name is unknown. Choose: 'rectangular'"//&
                            " or 'gaussian'.")
            end select
        else
            allocate(Null_Dirac_Distribution_Plus :: dirac_plus)
        end if

        select type (dirac_plus)
            type is (Rectangular_Dirac_Distribution_Plus)
                block
                    real(DP) :: delta_distance
                    data_field = prefix//"Dirac Distribution Plus.delta distance"
                    call exploring_data%get(data_field, delta_distance, data_found)
                    call check_data_found(data_field, data_found)
                    call dirac_plus%set(delta_distance)
                end block
            type is (Gaussian_Dirac_Distribution_Plus)
                block
                    real(DP) :: max_distance
                    integer :: max_distance_over_std_dev
                    data_field = prefix//"Dirac Distribution Plus.maximum distance"
                    call exploring_data%get(data_field, max_distance, data_found)
                    call check_data_found(data_field, data_found)
                    data_field = prefix//"Dirac Distribution Plus.maximum distance / "//&
                        "standard deviation"
                    call exploring_data%get(data_field, max_distance_over_std_dev, data_found)
                    call check_data_found(data_field, data_found)
                    call dirac_plus%set(max_distance, max_distance_over_std_dev)
                end block
            type is (Null_Dirac_Distribution_Plus)
            class default
                call error_exit("procedures_dirac_distribution_plus_factory: create: "//&
                        "distribution_name: unknown type.")
        end select
    end subroutine create

    subroutine destroy(dirac_plus)
        class(Abstract_Dirac_Distribution_Plus), allocatable, intent(inout) :: dirac_plus

        if (allocated(dirac_plus)) deallocate(dirac_plus)
    end subroutine destroy

end module procedures_dirac_distribution_plus_factory
