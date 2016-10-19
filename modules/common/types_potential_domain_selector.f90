module types_potential_domain_selector

implicit none

private

    type, public :: Concrete_Potential_Domain_Selector
        logical :: check_max = .true.
        logical :: check_max_over_box_edge = .true.
        logical :: check_delta = .true.
    end type Concrete_Potential_Domain_Selector

end module types_potential_domain_selector
