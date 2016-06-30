module procedures_parallelepiped_domain_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_errors, only: error_exit, warning_continue
use procedures_checks, only: check_data_found
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_visitable_walls, only: Abstract_Visitable_Walls
use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain, &
    Concrete_Parallelepiped_Domain, Concrete_Boxed_Domain, Concrete_Walled_Domain, &
    Null_Parallelepiped_Domain
use procedures_property_inquirers, only: use_walls

implicit none

private
public :: create_from_json, create_from_box, create_from_walls, destroy

contains

    subroutine create_from_json(parallelepiped_domain, periodic_box, visitable_walls, needed, &
        input_data, prefix)
        class(Abstract_Parallelepiped_Domain), allocatable, intent(out) :: parallelepiped_domain
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Visitable_Walls), intent(in) :: visitable_walls
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
                case ("boxed")
                    allocate(Concrete_Boxed_Domain :: parallelepiped_domain)
                case ("walled")
                    if (.not.use_walls(visitable_walls)) then
                        call warning_continue("procedures_parallelepiped_domain_factory: "//&
                            "create_from_json: walled: walls are not used.")
                    end if
                    allocate(Concrete_Walled_Domain :: parallelepiped_domain)
                case default
                    call error_exit("procedures_parallelepiped_domain_factory: create_from_json: "&
                        //domain_name//" domain_name unknown. Choose between 'domain', 'boxed' and"&
                        //" 'walled'.")
            end select
        else
            allocate(Null_Parallelepiped_Domain :: parallelepiped_domain)
        end if

        select type (parallelepiped_domain)
            type is (Concrete_Parallelepiped_Domain)
                call parallelepiped_domain%construct(periodic_box, domain_origin, domain_size)
            type is (Concrete_Boxed_Domain)
                call parallelepiped_domain%construct(periodic_box)
            type is (Concrete_Walled_Domain)
                call parallelepiped_domain%construct(periodic_box, visitable_walls)
            type is (Null_Parallelepiped_Domain)
            class default
                call error_exit("procedures_parallelepiped_domain_factory: create_from_json: "&
                    //"parallelepiped_domain: unknown type")
        end select
    end subroutine create_from_json

    subroutine create_from_box(parallelepiped_domain, periodic_box, needed)
        class(Abstract_Parallelepiped_Domain), allocatable, intent(out) :: parallelepiped_domain
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        logical, intent(in) :: needed

        if (needed) then
            allocate(Concrete_Boxed_Domain :: parallelepiped_domain)
        else
            allocate(Null_Parallelepiped_Domain :: parallelepiped_domain)
        end if

        select type (parallelepiped_domain)
            type is (Concrete_Boxed_Domain)
                call parallelepiped_domain%construct(periodic_box)
            type is (Null_Parallelepiped_Domain)
            class default
                call error_exit("procedures_parallelepiped_domain_factory: create_from_box: "&
                    //"parallelepiped_domain: unknown type")
        end select
    end subroutine create_from_box

    subroutine create_from_walls(parallelepiped_domain, periodic_box, visitable_walls, needed)
        class(Abstract_Parallelepiped_Domain), allocatable, intent(out) :: parallelepiped_domain
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Visitable_Walls), intent(in) :: visitable_walls
        logical, intent(in) :: needed

        if (needed) then
            allocate(Concrete_Walled_Domain :: parallelepiped_domain)
        else
            allocate(Null_Parallelepiped_Domain :: parallelepiped_domain)
        end if

        select type (parallelepiped_domain)
            type is (Concrete_Walled_Domain)
                call parallelepiped_domain%construct(periodic_box, visitable_walls)
            type is (Null_Parallelepiped_Domain)
            class default
            call error_exit("procedures_parallelepiped_domain_factory: create_from_walls: "&
                //"parallelepiped_domain: unknown type")
        end select
    end subroutine create_from_walls

    subroutine destroy(parallelepiped_domain)
        class(Abstract_Parallelepiped_Domain), allocatable, intent(inout) :: parallelepiped_domain

        if (allocated(parallelepiped_domain)) then
            call parallelepiped_domain%destroy()
            deallocate(parallelepiped_domain)
        end if
    end subroutine destroy

end module procedures_parallelepiped_domain_factory
