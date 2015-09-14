module procedures_parallelepiped_domain_write

use, intrinsic :: iso_fortran_env, only: DP => REAL64

implicit none

private
public write_corners

contains

    subroutine write_corners(box_origin, box_size, unit)
        real(DP), intent(in) :: box_origin(:), box_size(:)
        integer, intent(in) :: unit

        integer i, j, k, sgn_i, sgn_j

        sgn_j = 1
        do k = -1, 1, 2
        sgn_j = -sgn_j
        sgn_i = 1
        do j = -1, 1, 2
        sgn_i = -sgn_i
        do i = -1, 1, 2
            write(unit, *) box_origin + real([sgn_i*i, sgn_j*j, k], DP) * box_size / 2._DP
            write(unit, *)
        end do
        end do
        end do
    end subroutine write_corners

end module procedures_parallelepiped_domain_write

program test_parallelepiped_domain

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use data_geometry, only: num_dimensions
use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists, test_data_found
use class_periodic_box, only: Abstract_Periodic_Box, XYZ_Periodic_Box
use class_parallelepiped_domain, only: Abstract_Parallelepiped_Domain, &
                                       Concrete_Parallelepiped_Domain
use procedures_parallelepiped_domain_write, only: write_corners

implicit none

    class(Abstract_Periodic_Box), allocatable :: periodic_box
    class(Abstract_Parallelepiped_Domain), allocatable :: parallelepiped_domain
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename, data_field
    logical :: found
    real(DP), allocatable :: box_size(:), domain_origin(:), domain_size(:)
    real(DP) :: position(num_dimensions)
    integer :: num_particles, num_particles_inside, i_particle
    real(DP) :: rand(num_dimensions)
    integer :: box_unit, domain_unit, positions_box_unit, positions_domain_unit

    call json_initialize()

    data_filename = "parallelepiped_domain.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)

    data_field = "Periodic Box.size"
    call input_data%get(data_field, box_size, found)
    call test_data_found(data_field, found)
    allocate(XYZ_Periodic_Box :: periodic_box)
    call periodic_box%set_size(box_size)
    deallocate(box_size)
    open(newunit=box_unit, recl=4096, file="box.out", action="write")
    call write_corners([0._DP, 0._DP, 0._DP], periodic_box%get_size(), box_unit)
    close(box_unit)

    data_field = "Parallelepiped Domain.origin"
    call input_data%get(data_field, domain_origin, found)
    call test_data_found(data_field, found)
    data_field = "Parallelepiped Domain.size"
    call input_data%get(data_field, domain_size, found)
    call test_data_found(data_field, found)
    allocate(Concrete_Parallelepiped_Domain :: parallelepiped_domain)
    call parallelepiped_domain%construct(periodic_box, domain_origin, domain_size)
    open(newunit=domain_unit, recl=4096, file="domain.out", action="write")
    call write_corners(domain_origin, domain_size, domain_unit)
    close(domain_unit)

    data_field = "Particles.number"
    call input_data%get(data_field, num_particles, found)
    call test_data_found(data_field, found)
    write(output_unit, *) "Box :"
    open(newunit=positions_box_unit, recl=4096, file="positions_box.out", action="write")
    num_particles_inside = 0
    do i_particle = 1, num_particles
        call random_number(rand)
        position = periodic_box%folded(rand * periodic_box%get_size())
        if (parallelepiped_domain%is_inside(position)) then
            num_particles_inside = num_particles_inside + 1
            write(positions_box_unit, *) position
        end if
    end do
    close(positions_box_unit)
    write(output_unit, *) "    Inside positions ratio =", &
        real(num_particles_inside, DP) /  real(num_particles, DP)
    write(output_unit, *) "    Volumes ratio =", &
        parallelepiped_domain%get_volume() / product(periodic_box%get_size())

    write(output_unit, *) "Domain :"
    open(newunit=positions_domain_unit, recl=4096, file="positions_domain.out", action="write")
    num_particles_inside = 0
    do i_particle = 1, num_particles
        position = parallelepiped_domain%random_position()
        if (parallelepiped_domain%is_inside(position)) then
            num_particles_inside = num_particles_inside + 1
            write(positions_domain_unit, *) position
        end if
    end do
    close(positions_domain_unit)
    write(output_unit, *) "    Inside positions ratio =", &
        real(num_particles_inside, DP) /  real(num_particles, DP)

    call parallelepiped_domain%destroy()
    deallocate(parallelepiped_domain)
    deallocate(periodic_box)
    call input_data%destroy()

end program test_parallelepiped_domain
