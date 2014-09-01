module class_external_field

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_box, only: num_dimensions
use module_types_micro, only: Particle_Index
use module_physics_micro, only: exchange_sign

implicit none

private

    type, public :: External_Field
        real(DP), dimension(num_dimensions) :: vector
    contains
        procedure :: set => External_Field_set
        procedure :: total => External_Field_total
        
        procedure :: rotation => External_Field_rotation
        procedure :: exchange => External_Field_exchange
    end type External_Field

contains

    subroutine External_Field_set(this, vector)
    
        class(External_Field), intent(inout) :: this
        real(DP), dimension(num_dimensions) :: vector
        
        this%vector(:) = vector
    
    end subroutine External_Field_set
    
    pure function External_Field_total(this, total_moment) result(total)
    
        class(External_Field), intent(in) :: this
        real(DP), dimension(:), intent(in) :: total_moment
        real(DP) :: total
        
        total = -dot_product(total_moment, this%vector)
    
    end function External_Field_total
    
    !> \f[ \Delta U = -(\vec{\mu}^\prime - \vec{\mu} \cdot \vec{E}) \f]

    pure function External_Field_rotation(this, old_orientation, new_orientation) result(rotation)
    
        class(External_Field), intent(in) :: this
        real(DP), dimension(:), intent(in) :: old_orientation, new_orientation
        real(DP) :: rotation

        rotation = -dot_product(new_orientation - old_orientation, this%vector)

    end function External_Field_rotation
    
    !> \f[ \Delta U_{N\rightarrow{}N+1} = -(\vec{\mu}_{N+1} \cdot \vec{E}) \f]

    pure function External_Field_exchange(this, particle) result(exchange)

        class(External_Field), intent(in) :: this
        type(Particle_Index), intent(in) :: particle
        real(DP) :: exchange

        exchange = -exchange_sign(particle%add) * dot_product(particle%orientation, this%vector)

    end function External_Field_exchange
    
end module class_external_field
