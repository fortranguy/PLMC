module procedures_average_num_particles_factory

use data_input_prefixes, only: changes_prefix
use procedures_errors, only: error_exit
use json_module, only: json_file
use procedures_checks, only: check_data_found
use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain
use classes_num_particles, only: Abstract_Num_Particles
use classes_component_chemical_potential, only : Abstract_Component_Chemical_Potential
use classes_average_num_particles, only: Abstract_Average_Num_Particles, Constant_Num_Particles, &
    Constant_Chemical_Potential_Num_Particles, Concrete_Average_Num_Particles, &
    Null_Average_Num_Particles
use procedures_mixture_inquirers, only: component_can_exchange

implicit none

private
public :: create, destroy

contains

    subroutine create(average_num_particles, unique_box, accessible_domain, num_particles, exists, &
        chemical_potential, generating_data)
        class(Abstract_Average_Num_Particles), allocatable, intent(out) :: average_num_particles
        logical, intent(in) :: unique_box
        class(Abstract_Parallelepiped_Domain), intent(in) :: accessible_domain
        class(Abstract_Num_Particles), intent(in) :: num_particles
        logical, intent(in) :: exists
        class(Abstract_Component_Chemical_Potential), intent(in) :: chemical_potential
        type(json_file), intent(inout) :: generating_data

        integer :: accumulation_period
        character(len=:), allocatable :: data_field
        logical :: data_found

        if (exists) then
            if (unique_box) then
                if (component_can_exchange(chemical_potential)) then
                    allocate(Constant_Chemical_Potential_Num_Particles :: average_num_particles)
                else
                    allocate(Constant_Num_Particles :: average_num_particles)
                end if
            else
                allocate(Concrete_Average_Num_Particles :: average_num_particles)
            end if

            data_field = changes_prefix//"Components.accumulation period"
            call generating_data%get(data_field, accumulation_period, data_found)
            call check_data_found(data_field, data_found)
        else
            allocate(Null_Average_Num_Particles :: average_num_particles)
        end if

        select type (average_num_particles)
            type is (Constant_Num_Particles)
                call average_num_particles%construct(num_particles, accumulation_period)
            type is (Constant_Chemical_Potential_Num_Particles)
                call average_num_particles%construct(accessible_domain, chemical_potential, &
                    accumulation_period)
            type is (Concrete_Average_Num_Particles)
                call average_num_particles%construct(num_particles, accumulation_period)
            type is (Null_Average_Num_Particles)
                call average_num_particles%construct()
            class default
                call error_exit("procedures_average_num_particles_factory: create: "//&
                    "average_num_particles: type unknown.")
        end select
    end subroutine create

    subroutine destroy(average_num_particles)
        class(Abstract_Average_Num_Particles), allocatable, intent(inout) :: average_num_particles

        if (allocated(average_num_particles)) then
            call average_num_particles%destroy()
            deallocate(average_num_particles)
        end if
    end subroutine destroy

end module procedures_average_num_particles_factory
