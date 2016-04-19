module procedures_dlc_factory

use procedures_errors, only: error_exit
use types_environment_wrapper, only: Environment_Wrapper
use types_component_wrapper, only: Component_Wrapper
use classes_dlc_weight, only: Abstract_DLC_Weight, Concrete_DLC_Weight, Null_DLC_Weight
use classes_dlc_structures, only: Abstract_DLC_Structures, Concrete_DLC_Structures, &
    Null_DLC_Structures
use classes_dlc_visitor, only: Abstract_DLC_Visitor, Concrete_DLC_Visitor, Null_DLC_Visitor
use procedures_property_inquirers, only: periodicity_is_xy

implicit none

private
public :: dlc_create, dlc_destroy

interface dlc_create
    module procedure :: create_visitor
    module procedure :: create_structures
    module procedure :: create_weight
end interface dlc_create

interface dlc_destroy
    module procedure :: destroy_weight
    module procedure :: destroy_structures
    module procedure :: destroy_visitor
end interface dlc_destroy

contains

    subroutine create_visitor(visitor, environment, weight, structures)
        class(Abstract_DLC_Visitor), allocatable, intent(out) :: visitor
        type(Environment_Wrapper), intent(in) :: environment
        class(Abstract_DLC_Weight), intent(in) :: weight
        class(Abstract_DLC_Structures), intent(in) :: structures

        select type (structures)
            type is (Concrete_DLC_Structures)
                allocate(Concrete_DLC_Visitor :: visitor)
            type is (Null_DLC_Structures)
                allocate(Null_DLC_Visitor :: visitor)
            class default
                call error_exit("dlc_create: create_visitor: structures type unknown.")
        end select
        call visitor%construct(environment%periodic_box, environment%reciprocal_lattice, weight, &
            structures)
    end subroutine create_visitor

    subroutine destroy_visitor(visitor)
        class(Abstract_DLC_Visitor), allocatable, intent(inout) :: visitor

        if (allocated(visitor)) then
            call visitor%destroy()
            deallocate(visitor)
        end if
    end subroutine destroy_visitor

    subroutine create_structures(structures, environment, components, are_dipolar)
        class(Abstract_DLC_Structures), allocatable, intent(out) :: structures
        type(Environment_Wrapper), intent(in) :: environment
        type(Component_Wrapper), intent(in) :: components(:)
        logical, intent(in) :: are_dipolar(:)

        if (periodicity_is_xy(environment%periodic_box) .and. any(are_dipolar)) then
            allocate(Concrete_DLC_Structures :: structures)
        else
            allocate(Null_DLC_Structures :: structures)
        end if
        call structures%construct(environment%periodic_box, environment%reciprocal_lattice, &
            components, are_dipolar)
    end subroutine create_structures

    subroutine destroy_structures(structures)
        class(Abstract_DLC_Structures), allocatable, intent(inout) :: structures

        if (allocated(structures)) then
            call structures%destroy()
            deallocate(structures)
        end if
    end subroutine destroy_structures

    subroutine create_weight(weight, environment, dipoles_exist)
        class(Abstract_DLC_Weight), allocatable, intent(out) :: weight
        type(Environment_Wrapper), intent(in) :: environment
        logical, intent(in) :: dipoles_exist

        if (periodicity_is_xy(environment%periodic_box) .and. dipoles_exist) then
            allocate(Concrete_DLC_Weight :: weight)
        else
            allocate(Null_DLC_Weight :: weight)
        end if
        call weight%construct(environment%periodic_box, environment%reciprocal_lattice, &
            environment%permittivity)
    end subroutine create_weight

    subroutine destroy_weight(weight)
        class(Abstract_DLC_Weight), allocatable, intent(inout) :: weight

        if (allocated(weight)) then
            call weight%destroy()
            deallocate(weight)
        end if
    end subroutine destroy_weight

end module procedures_dlc_factory
