module procedures_parallelepiped_domains_print

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions
use class_parallelepiped_domain, only: Abstract_Parallelepiped_Domain

implicit none

private
public :: print_vertices, parallelepiped_domains_overlap

contains

    pure logical function parallelepiped_domains_overlap(domain_1, domain_2)
        class(Abstract_Parallelepiped_Domain), intent(in) :: domain_1, domain_2

        integer :: i_vertex, j_vertex, k_vertex

        parallelepiped_domains_overlap = .false.
        do k_vertex = 1, 2
            do j_vertex = 1, 2
                do i_vertex = 1, 2
                    parallelepiped_domains_overlap = &
                        domain_2%is_inside(domain_1%get_vertices([i_vertex, j_vertex, k_vertex]))
                    if (parallelepiped_domains_overlap) return
                end do
            end do
        end do
    end function parallelepiped_domains_overlap

    subroutine print_vertices(parallelepiped_domain, domain_name)
        class(Abstract_Parallelepiped_Domain), intent(in) :: parallelepiped_domain
        character(len=*), intent(in) :: domain_name

        integer :: i_vertex, j_vertex, k_vertex
        integer :: vertices_unit

        open(newunit=vertices_unit, recl=4096, file=domain_name//"_vertices.out", action="write")
        do k_vertex = 1, 2
            do j_vertex = 1, 2
                do i_vertex = 1, 2
                    write(vertices_unit, *) &
                        parallelepiped_domain%get_vertices([i_vertex, j_vertex, k_vertex])
                end do
                write(vertices_unit, *)
            end do
            write(vertices_unit, *)
        end do
        close(vertices_unit)
    end subroutine print_vertices

end module procedures_parallelepiped_domains_print

program test_parallelepiped_domains

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists
use class_periodic_box, only: Abstract_Periodic_Box
use class_parallelepiped_domain, only: Abstract_Parallelepiped_Domain
use procedures_box_factory, only: allocate_and_set_periodic_box, &
    allocate_and_construct_parallelepiped_domain
use procedures_parallelepiped_domains_print, only: print_vertices, parallelepiped_domains_overlap

implicit none

    class(Abstract_Periodic_Box), allocatable :: periodic_box
    class(Abstract_Parallelepiped_Domain), allocatable :: parallelepiped_domain_0, &
        parallelepiped_domain_1, parallelepiped_domain_2

    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename

    call json_initialize()
    data_filename = "parallelepiped_domains.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename=data_filename)
    deallocate(data_filename)

    call allocate_and_set_periodic_box(periodic_box, input_data, "Test Parallelepiped Domains")
    call allocate_and_construct_parallelepiped_domain(parallelepiped_domain_0, input_data, &
        "Test Parallelepiped Domains.Domain 0", periodic_box)
    call allocate_and_construct_parallelepiped_domain(parallelepiped_domain_1, input_data, &
        "Test Parallelepiped Domains.Domain 1", periodic_box)
    call allocate_and_construct_parallelepiped_domain(parallelepiped_domain_2, input_data, &
        "Test Parallelepiped Domains.Domain 2", periodic_box)

    call print_vertices(parallelepiped_domain_0, "domain_0")
    call print_vertices(parallelepiped_domain_1, "domain_1")
    call print_vertices(parallelepiped_domain_2, "domain_2")
    write(output_unit, *) "Domains 1-2 overlap:", &
        parallelepiped_domains_overlap(parallelepiped_domain_1, parallelepiped_domain_2)
    write(output_unit, *) "Domains 2-1 overlap:", &
        parallelepiped_domains_overlap(parallelepiped_domain_2, parallelepiped_domain_1)

    call parallelepiped_domain_2%destroy()
    deallocate(parallelepiped_domain_2)
    call parallelepiped_domain_1%destroy()
    deallocate(parallelepiped_domain_1)
    deallocate(periodic_box)
    call input_data%destroy()

end program test_parallelepiped_domains
