module class_ewald_summation_bound

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: PI
use data_box, only: num_dimensions
use module_types_micro, only: Particle_Index
use module_geometry, only: geometry
use module_physics_micro, only: exchange_sign
use class_hard_spheres, only: Dipolar_Hard_Spheres

implicit none

private

    type, public :: Ewald_Summation_Bound
        real(DP), dimension(num_dimensions) :: total_moment
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

        real(DP) :: modulus_drifted, modulus_reset

        modulus_drifted = norm2(this%total_moment(:))
        call this%set_total_moment(this_spheres)
        modulus_reset = norm2(this%total_moment(:))

        write(modulus_unit, *) i_step, abs(modulus_reset - modulus_drifted)

    end subroutine Ewald_Summation_Bound_reset_total_moment
    
    !> Total shape dependent term
    !> \f[
    !>      J(\vec{M}, S) = \frac{2\pi}{3V} |\vec{M}|^2
    !> \f]
    !> \f[
    !>      J(\vec{M}, R) = \frac{2\pi}{V} |M_z|^2
    !> \f]
    
    pure function Ewald_Summation_Bound_total(this, Box_size) result(total)
    
        class(Ewald_Summation_Bound), intent(in) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        real(DP) :: total

        if (geometry%bulk) then
            total = 2._DP*PI/3._DP / product(Box_size) * dot_product(this%total_moment, this%total_moment)
        else if(geometry%slab) then
            total = 2._DP*PI / product(Box_size) * this%total_moment(3)**2
        end if
    
    end function Ewald_Summation_Bound_total

    !> Rotation
    
    !> Difference of Energy
    !> Bulk:
    !> \f[
    !>      \Delta U = \frac{2\pi}{3V} [
    !>                      (\vec{\mu}^\prime_l \cdot \vec{\mu}^\prime_l) -
    !>                      (\vec{\mu}_l \cdot \vec{\mu}_l) +
    !>                      2 (\vec{\mu}^\prime_l - \vec{\mu}_l) \cdot \vec{M}_\underline{l}
    !>                 ]
    !> \f]
    !> Slab:
    !> \f[
    !>      \Delta U = \frac{2\pi}{V} [
    !>                      \mu^\prime_{l, z} \mu^\prime_{l, z} - \mu_{l, z} \mu_{l, z} +
    !>                      2 (\mu^\prime_{l, z} - \mu_{l, z}) M_{l, z}
    !>                 ]
    !> \f]
    
    pure function Ewald_Summation_Bound_rotation(this, Box_size, old, new) &
                  result (rotation)
    
        class(Ewald_Summation_Bound), intent(in) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        type(Particle_Index), intent(in) :: old, new
        real(DP) :: rotation

        if (geometry%bulk) then
        
            rotation = dot_product(new%orientation, new%orientation) - &
                       dot_product(old%orientation, old%orientation) + &
                       2._DP * dot_product(new%orientation - old%orientation, &
                                           this%total_moment - old%orientation)

            rotation = 2._DP*PI/3._DP / product(Box_size) * rotation

        else if (geometry%slab) then

            rotation = new%orientation(3)**2 - old%orientation(3)**2 + &
                       2._DP * (new%orientation(3) - old%orientation(3)) * &
                               (this%total_moment(3)-old%orientation(3))

            rotation = 2._DP*PI / product(Box_size) * rotation

        end if
    
    end function Ewald_Summation_Bound_rotation

    !> Update the total moment
    !> \f[
    !>      \Delta \vec{M} = \vec{\mu}^\prime_l - \vec{\mu}_l
    !> \f]

    pure subroutine Ewald_Summation_Bound_update_total_moment_rotation(this, old, new)
                                                                           
        class(Ewald_Summation_Bound), intent(inout) :: this
        type(Particle_Index), intent(in) :: old, new

        this%total_moment(:) = this%total_moment(:) + new%orientation(:) - old%orientation(:)

    end subroutine Ewald_Summation_Bound_update_total_moment_rotation
    
    !> Exchange
    
    !> Difference of Energy: add
    !> Bulk:
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
    !> Slab:
    !> Difference of Energy: add
    !> \f[
    !>      \Delta U_{N \rightarrow N+1} = \frac{2\pi}{V} [
    !>                                          \mu_{N+1, z} \mu_{N+1, z} +
    !>                                          2\mu_{N+1, z} M_z^N
    !>                                     ]
    !> \f]
    !> Difference of Energy: remove
    !> \f[
    !>      \Delta U_{N \rightarrow N-1} = \frac{2\pi}{V} [
    !>                                          \mu_{N, z} \mu_{N, z} -
    !>                                          2\mu_{N, z} M_z^N
    !>                                      ]
    !> \f]
    
    pure function Ewald_Summation_Bound_exchange(this, Box_size, particle) result (exchange)
    
        class(Ewald_Summation_Bound), intent(in) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        type(Particle_Index), intent(in) :: particle
        real(DP) :: exchange

        if (geometry%bulk) then
        
            exchange = dot_product(particle%orientation, particle%orientation) + &
                       2._DP * exchange_sign(particle%add) * &
                       dot_product(particle%orientation, this%total_moment)

            exchange = 2._DP*PI/3._DP / product(Box_size) * exchange

        else if(geometry%slab) then

            exchange = particle%orientation(3)**2 + &
                       2._DP * exchange_sign(particle%add) * &
                       particle%orientation(3) * this%total_moment(3)

            exchange = 2._DP*PI / product(Box_size) * exchange

        end if

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

    pure subroutine Ewald_Summation_Bound_update_total_moment_exchange(this, particle)

        class(Ewald_Summation_Bound), intent(inout) :: this
        type(Particle_Index), intent(in) :: particle

        this%total_moment(:) = this%total_moment(:) + &
                               exchange_sign(particle%add) * particle%orientation(:)

    end subroutine Ewald_Summation_Bound_update_total_moment_exchange

end module class_ewald_summation_bound
