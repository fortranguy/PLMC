module types_potential_domain_selector

implicit none

private

    type, public :: Concrete_Potential_Domain_Selector
        logical :: check_max_over_box_edge = .false.
        logical :: check_delta = .false.
    end type Concrete_Potential_Domain_Selector

end module types_potential_domain_selector
