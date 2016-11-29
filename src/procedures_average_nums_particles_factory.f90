module procedures_average_nums_particles_factory

use procedures_errors, only: error_exit
use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain
use types_component_wrapper, only: Component_Wrapper
use classes_average_num_particles, only: Abstract_Average_Num_Particles, Constant_Num_Particles, &
    Constant_Chemical_Potential_Num_Particles, Equipartition_Num_Particles, &
    Null_Average_Num_Particles

implicit none

private
public :: create, destroy

contains

    subroutine create(average_nums_particles, accessible_domains, components, components_exist, &
        can_exchange)
        class(Abstract_Average_Num_Particles), allocatable, intent(out) :: &
            average_nums_particles(:, :)
        class(Abstract_Parallelepiped_Domain), intent(in) :: accessible_domains(:)
        type(Component_Wrapper), target, intent(in) :: components(:, :)
        logical, intent(in) :: components_exist, can_exchange

        integer :: i_box, i_component

        if (components_exist) then
            if (size(components, 2) == 1) then
                if (can_exchange) then
                    allocate(Constant_Chemical_Potential_Num_Particles :: &
                        average_nums_particles(size(components, 1), size(components, 2)))
                else
                    allocate(Constant_Num_Particles :: &
                        average_nums_particles(size(components, 1), size(components, 2)))
                end if
            else
                allocate(Equipartition_Num_Particles :: &
                    average_nums_particles(size(components, 1), size(components, 2)))
            end if
        else
            allocate(Null_Average_Num_Particles :: &
                average_nums_particles(size(components, 1), size(components, 2)))
        end if

        select type (average_nums_particles)
            type is (Constant_Num_Particles)
                do i_box = 1, size(components, 2)
                    do i_component = 1, size(components, 1)
                        call average_nums_particles(i_component, i_box)%&
                            construct(components(i_component, i_box)%num_particles)
                    end do
                end do
            type is (Constant_Chemical_Potential_Num_Particles)
                do i_box = 1, size(components, 2)
                    do i_component = 1, size(components, 1)
                        call average_nums_particles(i_component, i_box)%&
                            construct(accessible_domains(i_box), components(i_component, i_box)%&
                                chemical_potential)
                    end do
                end do
            type is (Equipartition_Num_Particles)
                do i_box = 1, size(components, 2)
                    do i_component = 1, size(components, 1)
                        call average_nums_particles(i_component, i_box)%&
                            construct(components(:, i_box))
                    end do
                end do
            type is (Null_Average_Num_Particles)
            class default
                call error_exit("procedures_average_nums_particles_factory: create: "//&
                    "average_num_particles: type unknown.")
        end select
    end subroutine create

    subroutine destroy(average_nums_particles)
        class(Abstract_Average_Num_Particles), allocatable, intent(inout) :: &
            average_nums_particles(:, :)

        integer :: i_box, i_component

        if (allocated(average_nums_particles)) then
            do i_box = size(average_nums_particles, 2), 1, -1
                do i_component = size(average_nums_particles, 1), 1, -1
                    call average_nums_particles(i_component, i_box)%destroy()
                end do
            end do
            deallocate(average_nums_particles)
        end if
    end subroutine destroy

end module procedures_average_nums_particles_factory
