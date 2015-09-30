module procedures_visitable_list_sum

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use types_particle, only: Concrete_Particle
use class_particles_positions, only: Abstract_Particles_Positions
use class_visitable_list, only: Abstract_Visitable_List
use class_pair_potential, only: Abstract_Pair_Potential

implicit none

private
public sum_energy

contains

    subroutine sum_energy(particles_positions, visitable_list, pair_potential, overlap, energy)
        class(Abstract_Particles_Positions), intent(in) :: particles_positions
        class(Abstract_Visitable_List), intent(in) :: visitable_list
        class(Abstract_Pair_Potential), intent(in) :: pair_potential
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy

        type(Concrete_Particle) :: particle
        real(DP) :: energy_i
        integer :: i_particle

        energy = 0._DP
        do i_particle = 1, particles_positions%get_num()
            particle%same_type = .true.
            particle%i = i_particle
            particle%position = particles_positions%get(particle%i)
            call visitable_list%visit(overlap, energy_i, particle, pair_potential)
            if (overlap) exit
            energy = energy + energy_i
        end do
    end subroutine sum_energy

end module procedures_visitable_list_sum

program test_visitable_list

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists, test_data_found
use procedures_errors, only: error_exit
use class_periodic_box, only: Abstract_Periodic_Box
use procedures_environment_factory, only: environment_factory_create, environment_factory_destroy
use class_particles_number, only: Abstract_Particles_Number
use class_particles_diameter, only: Abstract_Particles_Diameter
use class_particles_positions, only: Abstract_Particles_Positions
use procedures_particles_factory, only: particles_factory_create, particles_factory_set, &
     particles_factory_destroy
use class_potential_expression, only: Abstract_Potential_Expression
use class_pair_potential, only: Abstract_Pair_Potential
use class_particles_potential, only: Abstract_Particles_Potential
use procedures_short_potential_factory, only: short_potential_factory_create, &
    short_potential_factory_destroy
use types_particle, only: Concrete_Particle
use class_visitable_list, only: Abstract_Visitable_List, Concrete_Visitable_List, &
    Concrete_Visitable_Array
use procedures_random, only: random_integer
use procedures_visitable_list_sum, only: sum_energy

implicit none

    class(Abstract_Visitable_List), allocatable :: visitable_list
    class(Abstract_Pair_Potential), allocatable :: pair_potential
    class(Abstract_Potential_Expression), allocatable :: potential_expression
    class(Abstract_Particles_Positions), allocatable :: particles_positions
    class(Abstract_Particles_Diameter), allocatable :: particles_diameter
    class(Abstract_Particles_Number), allocatable :: particles_number
    class(Abstract_Periodic_Box), allocatable :: periodic_box

    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename, data_field
    logical :: data_found

    type(Concrete_Particle) :: particle
    real(DP) :: energy
    integer :: num_exchanges, i_exchange, num_overwrites, i_overwrite
    integer :: i_particle, i_target, i_value
    logical :: overlap

    call json_initialize()
    data_filename = "particles_potential.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    deallocate(data_filename)

    call environment_factory_create(periodic_box, input_data, "Environment.")
    call particles_factory_create(particles_number, input_data, "Particles.")
    call particles_factory_create(particles_diameter, input_data, "Particles.")
    call particles_factory_create(particles_positions, periodic_box, particles_number)
    call particles_factory_set(particles_positions, input_data, "Particles.")
    call short_potential_factory_create(potential_expression, input_data, "Short Potential.", &
        particles_diameter)
    call short_potential_factory_create(pair_potential, input_data, "Short Potential.", &
        particles_diameter, potential_expression)
    call short_potential_factory_create(visitable_list, input_data, "Short Potential.", &
        pair_potential)

    call visitable_list%construct(periodic_box)
    do i_particle = 1, particles_positions%get_num()
        particle%i = i_particle
        particle%position = particles_positions%get(particle%i)
        call visitable_list%add(particle) !artificial
    end do
    call sum_energy(particles_positions, visitable_list, pair_potential, overlap, energy)
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
        particle%i = random_integer(particles_number%get())
        particle%position = particles_positions%get(particle%i)
        call visitable_list%remove(particle%i)
        call visitable_list%add(particle)
    end do
    if (num_exchanges > 0) then
        call sum_energy(particles_positions, visitable_list, pair_potential, overlap, energy)
        if (overlap) then
            write(output_unit,*) "overlap"
        else
            energy = energy / 2._DP
            write(output_unit, *) "[exchange] energy =", energy
        end if
    end if

    data_field = "Memory.number of overwrites"
    call input_data%get(data_field, num_overwrites, data_found)
    call test_data_found(data_field, data_found)
    do i_overwrite = 1, num_overwrites
        i_target = random_integer(particles_number%get())
        i_value = random_integer(particles_number%get())
        particle%i = i_value
        particle%position = particles_positions%get(particle%i)
        call visitable_list%set(i_target, particle)
        i_value = particle%i
        particle%i = i_target
        particle%position = particles_positions%get(particle%i)
        call visitable_list%set(i_value, particle)
    end do
    if (num_overwrites > 0) then
        call sum_energy(particles_positions, visitable_list, pair_potential, overlap, energy)
        if (overlap) then
            write(output_unit,*) "overlap"
        else
            energy = energy / 2._DP
            write(output_unit, *) "[overwrite] energy =", energy
        end if
    end if

    call visitable_list%destroy()
    deallocate(visitable_list)
    call short_potential_factory_destroy(pair_potential)
    call short_potential_factory_destroy(potential_expression)
    call particles_factory_destroy(particles_positions)
    call particles_factory_destroy(particles_diameter)
    call particles_factory_destroy(particles_number)
    call environment_factory_destroy(periodic_box)
    call input_data%destroy()

end program test_visitable_list
