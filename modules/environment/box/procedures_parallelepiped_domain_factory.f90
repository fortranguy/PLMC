module procedures_parallelepiped_domain_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_errors, only: error_exit
use procedures_checks, only: check_data_found
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain, &
    Concrete_Parallelepiped_Domain, Concrete_Box_Domain, Null_Parallelepiped_Domain

implicit none

private
public :: create, destroy

contains

    subroutine create(parallelepiped_domain, periodic_box, needed, input_data, prefix)
        class(Abstract_Parallelepiped_Domain), allocatable, intent(out) :: parallelepiped_domain
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        logical, intent(in) :: needed
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: domain_name, data_field
        real(DP), allocatable :: domain_origin(:), domain_size(:)
        logical :: data_found

        if (needed) then
            data_field = prefix//"Parallelepiped Domain.name"
            call input_data%get(data_field, domain_name, data_found)
            call check_data_found(data_field, data_found)
            select case (domain_name)
                case ("domain")
                    data_field = prefix//"Parallelepiped Domain.origin"
                    call input_data%get(data_field, domain_origin, data_found)
                    call check_data_found(data_field, data_found)
                    data_field = prefix//"Parallelepiped Domain.size"
                    call input_data%get(data_field, domain_size, data_found)
                    call check_data_found(data_field, data_found)
                    allocate(Concrete_Parallelepiped_Domain :: parallelepiped_domain)
                case ("box")
                    allocate(Concrete_Box_Domain :: parallelepiped_domain)
                case default
                    call error_exit(domain_name//" domain_name unknown. "//&
                        "Choose between 'domain' and 'box'.")
            end select
        else
            allocate(Null_Parallelepiped_Domain :: parallelepiped_domain)
        end if
        call parallelepiped_domain%construct(periodic_box, domain_origin, domain_size)
    end subroutine create

    subroutine destroy(parallelepiped_domain)
        class(Abstract_Parallelepiped_Domain), allocatable, intent(inout) :: parallelepiped_domain

        if (allocated(parallelepiped_domain)) then
            call parallelepiped_domain%destroy()
            deallocate(parallelepiped_domain)
        end if
    end subroutine destroy

end module procedures_parallelepiped_domain_factory
