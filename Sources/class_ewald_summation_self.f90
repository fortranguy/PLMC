module class_ewald_summation_self

use data_precisions, only: DP
use data_constants, only: PI
use class_hard_spheres, only: Dipolar_Hard_Spheres

implicit none

private

    type, public :: Ewald_Summation_Self
        real(DP) :: alpha
    contains
        procedure :: set_alpha => Ewald_Summation_Self_set_alpha
        procedure :: total => Ewald_Summation_Self_total
        procedure :: solo => Ewald_Summation_Self_solo
    end type Ewald_Summation_Self
    
contains
    
    pure subroutine Ewald_Summation_Self_set_alpha(this, alpha)
    
        class(Ewald_Summation_Self), intent(inout) :: this
        real(DP), intent(in) :: alpha
    
        this%alpha = alpha
    
    end subroutine Ewald_Summation_Self_set_alpha
    
    !> Total self energy
    !> \f[ \frac{2}{3}\frac{\alpha^3}{\sqrt{\pi}} \sum_i \vec{\mu}_i\cdot\vec{\mu}_i \f]
    
    pure function Ewald_Summation_Self_total(this, this_spheres) result(total)
    
        class(Ewald_Summation_Self), intent(in) :: this
        class(Dipolar_Hard_Spheres), intent(in) :: this_spheres
        real(DP) :: total

        integer :: i_particle
        
        total = 0._DP
        do i_particle = 1, this_spheres%get_num_particles()
            total = total + this%solo(this_spheres%get_orientation(i_particle))
        end do
        
    end function Ewald_Summation_Self_total
    
    !> Self energy of 1 dipole
    !> \f[ \frac{2}{3}\frac{\alpha^3}{\sqrt{\pi}} \vec{\mu}_i\cdot\vec{\mu}_i \f]
    
    pure function Ewald_Summation_Self_solo(this, orientation) result(solo)
    
        class(Ewald_Summation_Self), intent(in) :: this
        real(DP), dimension(:), intent(in) :: orientation
        real(DP) :: solo
        
        solo = 2._DP/3._DP * this%alpha**3/sqrt(PI) * dot_product(orientation, orientation)
    
    end function Ewald_Summation_Self_solo

end module class_ewald_summation_self
