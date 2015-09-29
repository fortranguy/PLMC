module types_change_counter

implicit none

private

    type, public :: Concrete_Change_Counter
        integer :: num_hits = 0
        integer :: num_success = 0
    end type Concrete_Change_Counter

end module types_change_counter
