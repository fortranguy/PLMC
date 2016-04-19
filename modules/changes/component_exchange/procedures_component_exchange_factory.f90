module procedures_component_exchange_factory

use types_component_wrapper, only: Component_Wrapper
use class_component_exchange, only: Abstract_Component_Exchange, Concrete_Component_Exchange, &
    Null_Component_Exchange
use procedures_property_inquirers, only: component_can_exchange

implicit none

private
public :: component_exchange_create, component_exchange_destroy

contains

    subroutine component_exchange_create(component_exchange, component)
        class(Abstract_Component_Exchange), allocatable, intent(out) :: component_exchange
        type(Component_Wrapper), intent(in) :: component

        if (component_can_exchange(component%chemical_potential)) then
            allocate(Concrete_Component_Exchange :: component_exchange)
        else
            allocate(Null_Component_Exchange :: component_exchange)
        end if
        call component_exchange%construct(component)
    end subroutine component_exchange_create

    subroutine component_exchange_destroy(component_exchange)
        class(Abstract_Component_Exchange), allocatable, intent(inout) :: component_exchange

        if (allocated(component_exchange)) then
            call component_exchange%destroy()
            deallocate(component_exchange)
        end if
    end subroutine component_exchange_destroy

end module procedures_component_exchange_factory
