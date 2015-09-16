module procedures_visitable_list_sum

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use module_particles, only: Concrete_Particle
use class_positions, only: Abstract_Positions
use class_visitable_list, only: Abstract_Visitable_List

implicit none

private
public sum_energy

contains

    subroutine sum_energy(positions, visitable_list, overlap, energy)
        class(Abstract_Positions), intent(in) :: positions
        class(Abstract_Visitable_List), intent(in) :: visitable_list
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy

        type(Concrete_Particle) :: particle
        real(DP) :: energy_i
        integer :: i_particle

        energy = 0._DP
        do i_particle = 1, positions%get_num()
            particle%i = i_particle
            particle%position = positions%get(particle%i)
            call visitable_list%visit(.true., particle, overlap, energy_i)
            if (overlap) exit
            energy = energy + energy_i
        end do
    end subroutine sum_energy

end module procedures_visitable_list_sum

program test_visitable_list

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use data_geometry, only: num_dimensions
use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists, test_data_found
use procedures_errors, only: error_exit
use class_periodic_box, only: Abstract_Periodic_Box, XYZ_Periodic_Box, XY_Periodic_Box
use class_number, only: Abstract_Number, Concrete_Number
use class_diameter, only: Abstract_Diameter, Concrete_Diameter
use class_positions, only: Abstract_Positions, Concrete_Positions
use types_potential_parameters, only: Abstract_Potential_Parameters, Null_Potential_Parameters, &
    Lennard_Jones_Parameters
use class_potential_expression, only: Abstract_Potential_Expression, Null_Potential_Expression, &
    Lennard_Jones_Expression
use types_potential_domain, only: Concrete_Potential_Domain
use class_pair_potential, only: Abstract_Pair_Potential, Concrete_Pair_Potential, &
    Hard_Pair_Potential
use module_particles, only: Concrete_Particle
use class_visitable_list, only: Abstract_Visitable_List, Concrete_Visitable_List
use procedures_visitable_list_sum, only: sum_energy

implicit none

    class(Abstract_Visitable_List), allocatable :: visitable_list
    class(Abstract_Pair_Potential), allocatable :: pair_potential
    type(Concrete_Potential_Domain) :: potential_domain
    class(Abstract_Potential_Expression), allocatable :: potential_expression
    class(Abstract_Potential_Parameters), allocatable :: potential_parameters
    class(Abstract_Positions), allocatable :: positions
    class(Abstract_Diameter), allocatable :: diameter
    class(Abstract_Number), allocatable :: number
    class(Abstract_Periodic_Box), allocatable :: periodic_box

    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename, data_field, positions_input
    character(len=:), allocatable :: box_name, potential_name
    logical :: data_found
    integer :: positions_unit

    real(DP), allocatable :: box_size(:)
    real(DP) :: position(num_dimensions)
    real(DP) :: energy, diameter_value, min_diameter_factor, rand
    integer :: num_particles, i_particle
    integer :: num_exchanges, i_exchange, num_overwrites, i_overwrite, i_target, i_value
    logical :: overlap

    call json_initialize()
    data_filename = "visitable_list.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    deallocate(data_filename)

    data_field = "Box.name"
    call input_data%get(data_field, box_name, data_found)
    call test_data_found(data_field, data_found)
    select case(box_name)
        case("XYZ")
            allocate(XYZ_Periodic_Box :: periodic_box)
        case("XY")
            allocate(XY_Periodic_Box :: periodic_box)
        case default
            call error_exit(data_field//" unkown.")
    end select
    data_field = "Box.size"
    call input_data%get(data_field, box_size, data_found)
    call test_data_found(data_field, data_found)
    call periodic_box%set(box_size)

    data_field = "Particles.number"
    call input_data%get(data_field, num_particles, data_found)
    call test_data_found(data_field, data_found)
    allocate(Concrete_Number :: number)
    call number%set(num_particles)

    allocate(Concrete_Positions :: positions)
    call positions%construct(periodic_box, number)
    positions_input = "positions.in"
    call test_file_exists(positions_input)
    open(newunit=positions_unit, recl=4096, file=positions_input, action="read")
    do i_particle = 1, positions%get_num()
        read(positions_unit, *) position
        call positions%set(i_particle, position)
    end do
    close(positions_unit)

    data_field = "Particles.diameter"
    call input_data%get(data_field, diameter_value, data_found)
    call test_data_found(data_field, data_found)
    data_field = "Particles.minimum diameter factor"
    call input_data%get(data_field, min_diameter_factor, data_found)
    call test_data_found(data_field, data_found)
    allocate(Concrete_Diameter :: diameter)
    call diameter%set(diameter_value, min_diameter_factor)

    data_field = "Potential.name"
    call input_data%get(data_field, potential_name, data_found)
    call test_data_found(data_field, data_found)

    select case(potential_name)
        case ("HS")
            allocate(Hard_Pair_Potential :: pair_potential)
            allocate(Null_Potential_Parameters :: potential_parameters)
        case ("LJ")
            allocate(Lennard_Jones_Parameters :: potential_parameters)
        case default
            call error_exit(data_field//" unkown.")
    end select

    select type(potential_parameters)
        type is (Null_Potential_Parameters)
            allocate(Null_Potential_Expression :: potential_expression)
        type is (Lennard_Jones_Parameters)
            data_field = "Potential.epsilon"
            call input_data%get(data_field, potential_parameters%epsilon, data_found)
            call test_data_found(data_field, data_found)
            potential_parameters%sigma = diameter%get()
            allocate(Lennard_Jones_Expression :: potential_expression)
    end select
    call potential_expression%set(potential_parameters)

    potential_domain%min = diameter%get_min()
    data_field = "Potential.max distance"
    call input_data%get(data_field, potential_domain%max, data_found)
    call test_data_found(data_field, data_found)
    data_field = "Potential.delta distance"
    call input_data%get(data_field, potential_domain%delta, data_found)
    call test_data_found(data_field, data_found)

    if (.not.allocated(pair_potential)) allocate(Concrete_Pair_Potential :: pair_potential)
    call pair_potential%construct(potential_domain, potential_expression)

    allocate(Concrete_Visitable_List :: visitable_list)
    call visitable_list%construct(periodic_box, positions)
    do i_particle = 1, positions%get_num()
        call visitable_list%allocate(i_particle) !artificial
    end do
    call visitable_list%set(pair_potential)
    call sum_energy(positions, visitable_list, overlap, energy)
    if (overlap) then
        write(output_unit,*) "overlap"
    else
        energy = energy / 2._DP
        write(output_unit, *) "[initial] energy =", energy
    end if

    data_field = "Particles.number of exchanges"
    call input_data%get(data_field, num_exchanges, data_found)
    call test_data_found(data_field, data_found)
    do i_exchange = 1, num_exchanges
        call random_number(rand)
        i_particle = int(real(number%get(), DP) * rand) + 1
        call visitable_list%deallocate(i_particle)
        call visitable_list%allocate(i_particle)
    end do
    call sum_energy(positions, visitable_list, overlap, energy)
    if (overlap) then
        write(output_unit,*) "overlap"
    else
        energy = energy / 2._DP
        write(output_unit, *) "[exchange] energy =", energy
    end if

    data_field = "Memory.number of overwrites"
    call input_data%get(data_field, num_overwrites, data_found)
    call test_data_found(data_field, data_found)
    do i_overwrite = 1, num_overwrites
        call random_number(rand)
        i_target = int(real(number%get(), DP) * rand) + 1
        call random_number(rand)
        i_value = int(real(number%get(), DP) * rand) + 1
        call visitable_list%overwrite(i_target, i_value)
        call visitable_list%overwrite(i_value, i_target)
    end do
    call sum_energy(positions, visitable_list, overlap, energy)
    if (overlap) then
        write(output_unit,*) "overlap"
    else
        energy = energy / 2._DP
        write(output_unit, *) "[overwrite] energy =", energy
    end if

    call visitable_list%destroy()
    deallocate(visitable_list)
    call pair_potential%destroy()
    deallocate(pair_potential)
    deallocate(potential_expression)
    deallocate(potential_parameters)
    deallocate(diameter)
    call positions%destroy()
    deallocate(positions)
    deallocate(number)
    deallocate(periodic_box)
    call input_data%destroy()

end program test_visitable_list
