module procedures_parallelepiped_domains_print

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions
use class_parallelepiped_domain, only: Abstract_Parallelepiped_Domain

implicit none

private
public :: print_vertices

contains

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
use procedures_box_factory, only: box_factory_create, box_factory_destroy
use procedures_parallelepiped_domains_print, only: print_vertices

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

    call box_factory_create(periodic_box, input_data, "Test Parallelepiped Domains")
    call box_factory_create(parallelepiped_domain_0, input_data, &
        "Test Parallelepiped Domains.Domain 0", periodic_box)
    call box_factory_create(parallelepiped_domain_1, input_data, &
        "Test Parallelepiped Domains.Domain 1", periodic_box)
    call box_factory_create(parallelepiped_domain_2, input_data, &
        "Test Parallelepiped Domains.Domain 2", periodic_box)

    call print_vertices(parallelepiped_domain_0, "domain_0")
    call print_vertices(parallelepiped_domain_1, "domain_1")
    call print_vertices(parallelepiped_domain_2, "domain_2")
    write(output_unit, *) "Domains 1-2 overlap:", &
        parallelepiped_domain_1%overlap(parallelepiped_domain_2)
    write(output_unit, *) "Domains 2-1 overlap:", &
        parallelepiped_domain_2%overlap(parallelepiped_domain_1)

    call box_factory_destroy(parallelepiped_domain_2)
    call box_factory_destroy(parallelepiped_domain_1)
    call box_factory_destroy(parallelepiped_domain_0)
    call box_factory_destroy(periodic_box)
    call input_data%destroy()

end program test_parallelepiped_domains
