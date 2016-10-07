module procedures_component_average_number_factory

use procedures_errors, only: error_exit
use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain
use classes_num_particles, only: Abstract_Num_Particles
use classes_component_chemical_potential, only : Abstract_Component_Chemical_Potential
use classes_component_average_number, only: Abstract_Component_Average_Number, &
    Constant_Component_Average_Number, Variable_Component_Average_Number, &
    Null_Component_Average_Number
use procedures_mixture_inquirers, only: component_can_exchange

implicit none

private
public :: create, destroy

contains

    subroutine create(average_number, accessible_domain, num_particles, chemical_potential)
        class(Abstract_Component_Average_Number), allocatable, intent(out) :: average_number
        class(Abstract_Parallelepiped_Domain), intent(in) :: accessible_domain
        class(Abstract_Num_Particles), intent(in) :: num_particles
        class(Abstract_Component_Chemical_Potential), intent(in) :: chemical_potential

        if (component_can_exchange(chemical_potential)) then
            allocate(Variable_Component_Average_Number :: average_number)
        else
            allocate(Constant_Component_Average_Number :: average_number)
        end if

        select type (average_number)
            type is (Constant_Component_Average_Number)
                call average_number%construct(num_particles)
            type is (Variable_Component_Average_Number)
                call average_number%construct(accessible_domain, chemical_potential)
            type is (Null_Component_Average_Number)
                call average_number%construct()
            class default
                call error_exit("procedures_component_average_number_factory: create: "//&
                    "average_number: type unknown.")
        end select
    end subroutine create

    subroutine destroy(average_number)
        class(Abstract_Component_Average_Number), allocatable, intent(inout) :: average_number

        if (allocated(average_number)) then
            call average_number%destroy()
            deallocate(average_number)
        end if
    end subroutine destroy

end module procedures_component_average_number_factory
