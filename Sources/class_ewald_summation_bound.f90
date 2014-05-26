module class_ewald_summation_bound

use data_precisions, only: DP
use data_constants, only: PI
use data_box, only: Ndim
use class_hard_spheres, only: Dipolar_Hard_Spheres

implicit none

private

    type, public :: Ewald_Summation_Bound
        real(DP), dimension(Ndim) :: total_moment
    contains
        procedure :: set_total_moment => Ewald_Summation_Bound_set_total_moment
        procedure :: reset_total_moment => Ewald_Summation_Bound_reset_total_moment
        procedure :: total => Ewald_Summation_Bound_total
        
        procedure :: rotation => Ewald_Summation_Bound_rotation
        procedure :: update_total_moment_rotation => Ewald_Summation_Bound_update_total_moment_rotation
        procedure :: exchange => Ewald_Summation_Bound_exchange
        procedure :: update_total_moment_exchange => Ewald_Summation_Bound_update_total_moment_exchange
    end type Ewald_Summation_Bound

contains

    !> Total dipole moment :
    !> \f[ \vec{M} = \sum_j \vec{\mu}_j \f]
    !> \f[ \vec{M}_\underline{l} = \sum_{j \neq l} \vec{\mu}_j \f]
    
    pure subroutine Ewald_Summation_Bound_set_total_moment(this, this_spheres)
    
        class(Ewald_Summation_Bound), intent(inout) :: this
        class(Dipolar_Hard_Spheres), intent(in) :: this_spheres
        integer :: i_particle
        
        this%total_moment(:) = 0._DP
        do i_particle = 1, this_spheres%get_num_particles()
            this%total_moment(:) = this%total_moment(:) + this_spheres%get_orientation(i_particle)
        end do
        
    end subroutine Ewald_Summation_Bound_set_total_moment

    subroutine Ewald_Summation_Bound_reset_total_moment(this, this_spheres, i_step, modulus_unit)

        class(Ewald_Summation_Bound), intent(inout) :: this
        class(Dipolar_Hard_Spheres), intent(in) :: this_spheres
        integer, intent(in) :: i_step
        integer, intent(in) :: modulus_unit

        real(DP) :: modulus_drifted, modulus_reinit

        modulus_drifted = norm2(this%total_moment(:))
        call this%set_total_moment(this_spheres)
        modulus_reinit = norm2(this%total_moment(:))

        write(modulus_unit, *) i_step, abs(modulus_reinit - modulus_drifted)

    end subroutine Ewald_Summation_Bound_reset_total_moment
    
    !> Total shape dependent term
    !> \f[
    !>      J(\vec{M}, S) = \frac{2\pi}{3V} | \vec{M}|^2
    !> \f]
    
    pure function Ewald_Summation_Bound_total(this, Box_size) result(total)
    
        class(Ewald_Summation_Bound), intent(in) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        real(DP) :: total
        
        total = 2._DP*PI/3._DP / product(Box_size) * dot_product(this%total_moment, this%total_moment)
    
    end function Ewald_Summation_Bound_total

    !> Rotation
    
    !> Difference of Energy
    !> \f[
    !>      \Delta U = \frac{2\pi}{3V} [
    !>                      (\vec{\mu}^\prime_l \cdot \vec{\mu}^\prime_l) -
    !>                      (\vec{\mu}_l \cdot \vec{\mu}_l) +
    !>                      2 (\vec{\mu}^\prime_l - \vec{\mu}_l) \cdot \vec{M}_\underline{l}
    !>                 ]
    !> \f]
    
    pure function Ewald_Summation_Bound_rotation(this, Box_size, old_orientation, new_orientation) &
                  result (rotation)
    
        class(Ewald_Summation_Bound), intent(in) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        real(DP), dimension(:), intent(in) :: old_orientation, new_orientation
        real(DP) :: rotation
        
        rotation = dot_product(new_orientation, new_orientation) - &
                   dot_product(old_orientation, old_orientation) + &
                   2._DP * dot_product(new_orientation - old_orientation, &
                   this%total_moment - old_orientation)
                          
        rotation = 2._DP*PI/3._DP / product(Box_size) * rotation
    
    end function Ewald_Summation_Bound_rotation

    !> Update the total moment
    !> \f[
    !>      \Delta \vec{M} = \vec{\mu}^\prime_l - \vec{\mu}_l
    !> \f]

    pure subroutine Ewald_Summation_Bound_update_total_moment_rotation(this, old_orientation, &
                                                                           new_orientation)
                                                                           
        class(Ewald_Summation_Bound), intent(inout) :: this
        real(DP), dimension(:), intent(in) :: old_orientation, new_orientation

        this%total_moment(:) = this%total_moment(:) + new_orientation(:) - old_orientation(:)

    end subroutine Ewald_Summation_Bound_update_total_moment_rotation
    
    !> Exchange
    
    !> Difference of Energy: add
    !> \f[
    !>      \Delta U_{N \rightarrow N+1} = \frac{2\pi}{3V} [
    !>                                         (\vec{\mu}_{N+1} \cdot \vec{\mu}_{N+1})
    !>                                          +2(\vec{\mu} \cdot \vec{M}_N)
    !>                                     ]
    !> \f]
    
    !> Difference of Energy: remove
    !> \f[
    !>      \Delta U_{N \rightarrow N-1} = \frac{2\pi}{3V} [
    !>                                          (\vec{\mu}_{N+1} \cdot \vec{\mu}_{N+1})
    !>                                          -2(\vec{\mu} \cdot \vec{M}_N)
    !>                                      ]
    !> \f]
    
    pure function Ewald_Summation_Bound_exchange(this, Box_size, orientation) result (exchange)
    
        class(Ewald_Summation_Bound), intent(in) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        real(DP), dimension(:), intent(in) :: orientation
        real(DP) :: exchange
        
        exchange = dot_product(orientation, orientation) + &
                   2._DP * dot_product(orientation, this%total_moment)
                          
        exchange = 2._DP*PI/3._DP / product(Box_size) * exchange
    
    end function Ewald_Summation_Bound_exchange   
    
    !> Exchange

    !> Update the total moment: add
    !> \f[
    !>      \Delta \vec{M} = +\vec{\mu}_l
    !> \f]

    !> Update the total moment: remove
    !> \f[
    !>      \Delta \vec{M} = -\vec{\mu}_l
    !> \f]

    pure subroutine Ewald_Summation_Bound_update_total_moment_exchange(this, orientation)

        class(Ewald_Summation_Bound), intent(inout) :: this
        real(DP), dimension(:), intent(in) :: orientation

        this%total_moment(:) = this%total_moment(:) + orientation(:)

    end subroutine Ewald_Summation_Bound_update_total_moment_exchange 

end module class_ewald_summation_bound
