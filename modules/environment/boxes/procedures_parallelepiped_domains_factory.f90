module procedures_parallelepiped_domains_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_errors, only: error_exit, warning_continue
use procedures_checks, only: check_data_found
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_visitable_walls, only: Abstract_Visitable_Walls
use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain, &
    Concrete_Parallelepiped_Domain, Boxed_Parallelepiped_Domain, Walled_Parallelepiped_Domain, &
    Null_Parallelepiped_Domain
use procedures_environment_inquirers, only: use_walls

implicit none

private
public :: create_from_json, create_from_boxes, create_from_walls, destroy_line, destroy_element

contains

    subroutine create_from_json(parallelepiped_domains, periodic_boxes, visitable_walls, needed, &
        input_data, prefix)
        class(Abstract_Parallelepiped_Domain), allocatable, intent(out) :: parallelepiped_domains(:)
        class(Abstract_Periodic_Box), intent(in) :: periodic_boxes(:)
        class(Abstract_Visitable_Walls), intent(in) :: visitable_walls(:)
        logical, intent(in) :: needed
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        integer :: i_box
        character(len=:), allocatable :: domain_name, data_field
        real(DP), allocatable :: domain_origin(:), domain_size(:)
        logical :: data_found

        if (needed) then
            data_field = prefix//"Parallelepiped Domain.name"
            call input_data%get(data_field, domain_name, data_found)
            call check_data_found(data_field, data_found)
            select case (domain_name)
                case ("domain")
                    allocate(Concrete_Parallelepiped_Domain :: &
                        parallelepiped_domains(size(periodic_boxes)))
                case ("boxed")
                    allocate(Boxed_Parallelepiped_Domain :: &
                        parallelepiped_domains(size(periodic_boxes)))
                case ("walled")
                    if (all(.not.use_walls(visitable_walls))) then
                        call warning_continue("procedures_parallelepiped_domains_factory: "//&
                            "create_from_json: walled: walls are not used.")
                    end if
                    allocate(Walled_Parallelepiped_Domain :: &
                        parallelepiped_domains(size(periodic_boxes)))
                case default
                    call error_exit("procedures_parallelepiped_domains_factory: create_from_json: "&
                        //domain_name//" domain_name unknown. Choose between 'domain', 'boxed' and"&
                        //" 'walled'.")
            end select
        else
            allocate(Null_Parallelepiped_Domain :: parallelepiped_domains(size(periodic_boxes)))
        end if

        select type (parallelepiped_domains)
            type is (Concrete_Parallelepiped_Domain)
                data_field = prefix//"Parallelepiped Domain.origin"
                call input_data%get(data_field, domain_origin, data_found)
                call check_data_found(data_field, data_found)
                data_field = prefix//"Parallelepiped Domain.size"
                call input_data%get(data_field, domain_size, data_found)
                call check_data_found(data_field, data_found)
                do i_box = 1, size(parallelepiped_domains)
                    call parallelepiped_domains(i_box)%construct(periodic_boxes(i_box), &
                        domain_origin, domain_size)
                end do
            type is (Boxed_Parallelepiped_Domain)
                do i_box = 1, size(parallelepiped_domains)
                    call parallelepiped_domains(i_box)%construct(periodic_boxes(i_box))
                end do
            type is (Walled_Parallelepiped_Domain)
                do i_box = 1, size(parallelepiped_domains)
                    call parallelepiped_domains(i_box)%construct(periodic_boxes(i_box), &
                        visitable_walls(i_box))
                end do
            type is (Null_Parallelepiped_Domain)
            class default
                call error_exit("procedures_parallelepiped_domains_factory: create_from_json: "&
                    //"parallelepiped_domains: unknown type")
        end select
    end subroutine create_from_json

    subroutine create_from_boxes(parallelepiped_domains, periodic_boxes, needed)
        class(Abstract_Parallelepiped_Domain), allocatable, intent(out) :: parallelepiped_domains(:)
        class(Abstract_Periodic_Box), intent(in) :: periodic_boxes(:)
        logical, intent(in) :: needed

        integer :: i_box

        if (needed) then
            allocate(Boxed_Parallelepiped_Domain :: parallelepiped_domains(size(periodic_boxes)))
        else
            allocate(Null_Parallelepiped_Domain :: parallelepiped_domains(size(periodic_boxes)))
        end if

        select type (parallelepiped_domains)
            type is (Boxed_Parallelepiped_Domain)
                do i_box = 1, size(parallelepiped_domains)
                    call parallelepiped_domains(i_box)%construct(periodic_boxes(i_box))
                end do
            type is (Null_Parallelepiped_Domain)
            class default
                call error_exit("procedures_parallelepiped_domains_factory: create_from_boxes: "&
                    //"parallelepiped_domains: unknown type")
        end select
    end subroutine create_from_boxes

    subroutine create_from_walls(parallelepiped_domains, periodic_boxes, visitable_walls, needed)
        class(Abstract_Parallelepiped_Domain), allocatable, intent(out) :: parallelepiped_domains(:)
        class(Abstract_Periodic_Box), intent(in) :: periodic_boxes(:)
        class(Abstract_Visitable_Walls), intent(in) :: visitable_walls(:)
        logical, intent(in) :: needed

        integer :: i_box

        if (needed) then
            allocate(Walled_Parallelepiped_Domain :: parallelepiped_domains(size(periodic_boxes)))
        else
            allocate(Null_Parallelepiped_Domain :: parallelepiped_domains(size(periodic_boxes)))
        end if

        select type (parallelepiped_domains)
            type is (Walled_Parallelepiped_Domain)
                do i_box = 1, size(parallelepiped_domains)
                    call parallelepiped_domains(i_box)%construct(periodic_boxes(i_box), &
                        visitable_walls(i_box))
                end do
            type is (Null_Parallelepiped_Domain)
            class default
            call error_exit("procedures_parallelepiped_domains_factory: create_from_walls: "&
                //"parallelepiped_domains: unknown type")
        end select
    end subroutine create_from_walls

    subroutine destroy_line(parallelepiped_domains)
        class(Abstract_Parallelepiped_Domain), allocatable, intent(inout) :: &
            parallelepiped_domains(:)

        integer :: i_box

        if (allocated(parallelepiped_domains)) then
            do i_box = size(parallelepiped_domains), 1, -1
                call parallelepiped_domains(i_box)%destroy()
            end do
            deallocate(parallelepiped_domains)
        end if
    end subroutine destroy_line

    subroutine destroy_element(parallelepiped_domain)
        class(Abstract_Parallelepiped_Domain), allocatable, intent(inout) :: parallelepiped_domain

        if (allocated(parallelepiped_domain)) then
            call parallelepiped_domain%destroy()
            deallocate(parallelepiped_domain)
        end if
    end subroutine destroy_element

end module procedures_parallelepiped_domains_factory
