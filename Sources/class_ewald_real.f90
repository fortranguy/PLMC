module class_ewald_real

use data_precisions, only: DP
use data_constants, only: PI

implicit none

private

    type, public :: Ewald_Real
    
        real(DP) :: min_distance
        real(DP) :: range_cut
        real(DP) :: delta
        
        integer :: i_min_distance
        integer :: i_range_cut        
        real(DP), dimension(:, :), allocatable :: tabulation
    
    contains
    
        !procedure :: construct => Ewald_Real_construct
        !procedure :: destroy => Ewald_Real_destroy
        
        !procedure, private :: set_parameters => Ewald_Real_set_parameters
        !procedure, private :: set_tabulation => Ewald_Real_set_tabulation
        !procedure, private :: set => Ewald_Real_set
        !procedure :: write => Ewald_Real_write
        !procedure, private :: interpolation => Ewald_Real_interpolation
        !procedure, private :: pair => Ewald_Real_pair
        !procedure :: solo => Ewald_Real_solo
        !procedure, private :: Epot_real => Ewald_Real_Epot_real
    
    end type Ewald_Real
    
contains
    
    

end module class_ewald_real
