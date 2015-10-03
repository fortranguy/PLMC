module procedures_observable_writers_factory

use json_module, only: json_file
use module_data, only: test_data_found
use class_particles_number, only: Abstract_Particles_Number
use class_particles_diameter, only: Abstract_Particles_Diameter
use class_particles_exchange, only: Abstract_Particles_Exchange
use class_moved_positions, only: Abstract_Moved_Positions
use class_rotated_orientations, only: Abstract_Rotated_Orientations
use class_particles_energy_writer, only: Concrete_Energy_Writer_Selector, &
    Abstract_Particles_Energy_Writer, Concrete_Particles_Energy_Writer, Null_Particles_Energy_Writer
use class_inter_energy_writer, only: Abstract_Inter_Energy_Writer, &
    Concrete_Inter_Energy_Writer, Null_Inter_Energy_Writer
use class_changes_writer, only: Concrete_Changes_Selector, &
    Abstract_Changes_Success_Writer, Concrete_Changes_Success_Writer, Null_Changes_Success_Writer
use procedures_property_inquirers, only: use_walls, particles_exist, mixture_exists, &
    particles_can_move, particles_can_rotate, particles_can_exchange

implicit none

private
public :: observable_writers_factory_create, observable_writers_factory_destroy

interface observable_writers_factory_create
    module procedure :: allocate_and_construct_energy_writer
    module procedure :: allocate_and_construct_inter_energy_writer
    module procedure :: allocate_and_construct_changes_writer
end interface observable_writers_factory_create

interface observable_writers_factory_destroy
    module procedure :: destroy_and_deallocate_changes_writer
    module procedure :: destroy_and_deallocate_inter_energy_writer
    module procedure :: destroy_and_deallocate_particles_energy_writer
end interface observable_writers_factory_destroy

contains

    subroutine allocate_and_construct_energy_writer(particles_energy_writer, particles_number, &
        wall_diameter, filename)
        class(Abstract_Particles_Energy_Writer), allocatable, intent(out) :: particles_energy_writer
        class(Abstract_Particles_Number), intent(in) :: particles_number
        class(Abstract_Particles_Diameter), intent(in) :: wall_diameter
        character(len=*), intent(in) :: filename

        type(Concrete_Energy_Writer_Selector) :: energy_selector

        if (particles_exist(particles_number)) then
            allocate(Concrete_Particles_Energy_Writer :: particles_energy_writer)
        else
            allocate(Null_Particles_Energy_Writer :: particles_energy_writer)
        end if
        energy_selector%write_walls = use_walls(wall_diameter)
        call particles_energy_writer%construct(filename, energy_selector)
    end subroutine allocate_and_construct_energy_writer

    subroutine allocate_and_construct_inter_energy_writer(inter_energy_writer, &
        inter_diameter, filename)
        class(Abstract_Inter_Energy_Writer), allocatable, intent(out) :: inter_energy_writer
        class(Abstract_Particles_Diameter), intent(in) :: inter_diameter
        character(len=*), intent(in) :: filename

        if (mixture_exists(inter_diameter)) then
            allocate(Concrete_Inter_Energy_Writer :: inter_energy_writer)
        else
            allocate(Null_Inter_Energy_Writer :: inter_energy_writer)
        end if
        call inter_energy_writer%construct(filename)
    end subroutine allocate_and_construct_inter_energy_writer

    subroutine allocate_and_construct_changes_writer(changes_success_writer, &
        moved_positions, rotated_orientations, particles_exchange, filename)
        class(Abstract_Changes_Success_Writer), allocatable, intent(out) :: changes_success_writer
        class(Abstract_Moved_Positions), intent(in) :: moved_positions
        class(Abstract_Rotated_Orientations), intent(in) :: rotated_orientations
        class(Abstract_Particles_Exchange), intent(in) :: particles_exchange
        character(len=*), intent(in) :: filename

        type(Concrete_Changes_Selector) :: changes_selector

        if (particles_can_move(moved_positions)) then
            allocate(Concrete_Changes_Success_Writer :: changes_success_writer)
        else
            allocate(Null_Changes_Success_Writer :: changes_success_writer)
        end if
        changes_selector%write_rotations = particles_can_rotate(rotated_orientations)
        changes_selector%write_exchanges = particles_can_exchange(particles_exchange)
        call changes_success_writer%construct(filename, changes_selector)
    end subroutine allocate_and_construct_changes_writer

    subroutine destroy_and_deallocate_particles_energy_writer(particles_energy_writer)
        class(Abstract_Particles_Energy_Writer), allocatable, intent(inout) :: &
            particles_energy_writer

        call particles_energy_writer%destroy()
        if (allocated(particles_energy_writer)) deallocate(particles_energy_writer)
    end subroutine destroy_and_deallocate_particles_energy_writer

    subroutine destroy_and_deallocate_inter_energy_writer(inter_energy_writer)
        class(Abstract_Inter_Energy_Writer), allocatable, intent(inout) :: inter_energy_writer

        call inter_energy_writer%destroy()
        if (allocated(inter_energy_writer)) deallocate(inter_energy_writer)
    end subroutine destroy_and_deallocate_inter_energy_writer

    subroutine destroy_and_deallocate_changes_writer(changes_success_writer)
        class(Abstract_Changes_Success_Writer), allocatable, intent(inout) :: changes_success_writer

        call changes_success_writer%destroy()
        if (allocated(changes_success_writer)) deallocate(changes_success_writer)
    end subroutine destroy_and_deallocate_changes_writer

end module procedures_observable_writers_factory

