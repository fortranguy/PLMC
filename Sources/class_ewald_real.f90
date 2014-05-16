module class_ewald_real

implicit none

private

    type, public :: Ewald_Real
    
        real(DP) :: min_distance
        real(DP) :: range_cut
        real(DP) :: delta
        integer :: i_min_distance
        integer :: i_range_cut
    
    contains
    
    end type Ewald_Real

end module class_ewald_real
