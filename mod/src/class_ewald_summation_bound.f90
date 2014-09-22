module class_ewald_summation_bound

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: PI
use data_box, only: num_dimensions
use module_types_micro, only: Particle_Index
use module_geometry, only: geometry
use module_linear_algebra, only: projector_matrix
use module_physics_micro, only: exchange_sign
use class_hard_spheres, only: Dipolar_Hard_Spheres

implicit none

private

    type, public :: Ewald_Summation_Bound
        real(DP), dimension(num_dimensions) :: total_moment
    contains
        procedure :: set_total_moment => Ewald_Summation_Bound_set_total_moment
        procedure :: get_total_moment => Ewald_Summation_Bound_get_total_moment
        procedure :: reset_total_moment => Ewald_Summation_Bound_reset_total_moment
        procedure :: total_energy => Ewald_Summation_Bound_total_energy
        
        procedure :: solo_field => Ewald_Summation_Bound_solo_field
        
        procedure :: rotation_energy => Ewald_Summation_Bound_rotation_energy
        procedure :: update_total_moment_rotation => Ewald_Summation_Bound_update_total_moment_rotation
        procedure :: exchange_energy => Ewald_Summation_Bound_exchange_energy
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
    
    pure function Ewald_Summation_Bound_get_total_moment(this) result(get_total_moment)
    
        class(Ewald_Summation_Bound), intent(in) :: this
        real(DP), dimension(num_dimensions) :: get_total_moment
        
        get_total_moment(:) = this%total_moment(:)
        
    end function Ewald_Summation_Bound_get_total_moment

    subroutine Ewald_Summation_Bound_reset_total_moment(this, this_spheres, i_step, modulus_unit)

        class(Ewald_Summation_Bound), intent(inout) :: this
        class(Dipolar_Hard_Spheres), intent(in) :: this_spheres
        integer, intent(in) :: i_step
        integer, intent(in), optional :: modulus_unit

        real(DP) :: modulus_drifted, modulus_reset

        modulus_drifted = norm2(this%total_moment(:))
        call this%set_total_moment(this_spheres)
        modulus_reset = norm2(this%total_moment(:))

        if (present(modulus_unit)) then
            write(modulus_unit, *) i_step, abs(modulus_reset - modulus_drifted)
        end if

    end subroutine Ewald_Summation_Bound_reset_total_moment
    
    !> Total shape dependent term
    !> \f[
    !>      J(\vec{M}, S) = \frac{2\pi}{3V} |\vec{M}|^2
    !> \f]
    !> \f[
    !>      J(\vec{M}, R) = \frac{2\pi}{V} |M_z|^2
    !> \f]
    
    pure function Ewald_Summation_Bound_total_energy(this, Box_size) result(total_energy)
    
        class(Ewald_Summation_Bound), intent(in) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        real(DP) :: total_energy

        if (geometry%bulk) then
            total_energy = 2._DP*PI/3._DP / product(Box_size) * &
                           dot_product(this%total_moment, this%total_moment)
        else if(geometry%slab) then
            total_energy = 2._DP*PI / product(Box_size) * this%total_moment(3)**2
        end if
    
    end function Ewald_Summation_Bound_total_energy
    
    !> Field
    !> \f[
    !>      \vec{E}(\vec{M}, S) = - \frac{4\pi}{3V}  \vec{M}
    !> \f]
    !> \f[
    !>      \vec{E}(\vec{M}, R) = - \frac{4\pi}{V}  |\vec{e}_z)(\vec{e}_z| \vec{M}
    !> \f]
    
    pure function Ewald_Summation_Bound_solo_field(this, Box_size) result(solo_field)

        class(Ewald_Summation_Bound), intent(in) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        real(DP), dimension(num_dimensions) :: solo_field
        
        if (geometry%bulk) then
            solo_field(:) = -4._DP/3._DP * PI/product(Box_size) * this%total_moment(:)
        else if(geometry%slab) then
            solo_field(:) = -4._DP * PI/product(Box_size) * &
                            matmul(projector_matrix(num_dimensions, 3), this%total_moment)
        end if

    end function Ewald_Summation_Bound_solo_field

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
    
    pure function Ewald_Summation_Bound_rotation_energy(this, Box_size, old, new) &
                  result (rotation_energy)
    
        class(Ewald_Summation_Bound), intent(in) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        type(Particle_Index), intent(in) :: old, new
        real(DP) :: rotation_energy

        if (geometry%bulk) then
        
            rotation_energy = dot_product(new%orientation, new%orientation) - &
                              dot_product(old%orientation, old%orientation) + &
                              2._DP * dot_product(new%orientation - old%orientation, &
                                                  this%total_moment - old%orientation)

            rotation_energy = 2._DP*PI/3._DP / product(Box_size) * rotation_energy

        else if (geometry%slab) then

            rotation_energy = new%orientation(3)**2 - old%orientation(3)**2 + &
                              2._DP * (new%orientation(3) - old%orientation(3)) * &
                                      (this%total_moment(3) - old%orientation(3))

            rotation_energy = 2._DP*PI / product(Box_size) * rotation_energy

        end if
    
    end function Ewald_Summation_Bound_rotation_energy

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
    
    pure function Ewald_Summation_Bound_exchange_energy(this, Box_size, particle) &
         result (exchange_energy)
    
        class(Ewald_Summation_Bound), intent(in) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        type(Particle_Index), intent(in) :: particle
        real(DP) :: exchange_energy

        if (geometry%bulk) then
        
            exchange_energy = dot_product(particle%orientation, particle%orientation) + &
                              2._DP * exchange_sign(particle%add) * &
                              dot_product(particle%orientation, this%total_moment)

            exchange_energy = 2._DP*PI/3._DP / product(Box_size) * exchange_energy

        else if(geometry%slab) then

            exchange_energy = particle%orientation(3)**2 + &
                              2._DP * exchange_sign(particle%add) * &
                              particle%orientation(3) * this%total_moment(3)

            exchange_energy = 2._DP*PI / product(Box_size) * exchange_energy

        end if

    end function Ewald_Summation_Bound_exchange_energy
    
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
