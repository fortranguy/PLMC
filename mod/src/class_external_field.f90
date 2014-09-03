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
        procedure :: total_energy => External_Field_total_energy
        
        procedure :: rotation_energy => External_Field_rotation_energy
        procedure :: exchange_energy => External_Field_exchange_energy
    end type External_Field

contains

    subroutine External_Field_set(this, vector)
    
        class(External_Field), intent(inout) :: this
        real(DP), dimension(num_dimensions) :: vector
        
        this%vector(:) = vector
    
    end subroutine External_Field_set
    
    pure function External_Field_total_energy(this, total_moment) result(total_energy)
    
        class(External_Field), intent(in) :: this
        real(DP), dimension(:), intent(in) :: total_moment
        real(DP) :: total_energy
        
        total_energy = -dot_product(total_moment, this%vector)
    
    end function External_Field_total_energy
    
    !> \f[ \Delta U = -(\vec{\mu}^\prime - \vec{\mu} \cdot \vec{E}) \f]

    pure function External_Field_rotation_energy(this, old, new) result(rotation_energy)
    
        class(External_Field), intent(in) :: this
        type(Particle_Index), intent(in) :: old, new
        real(DP) :: rotation_energy

        rotation_energy = -dot_product(new%orientation - old%orientation, this%vector)

    end function External_Field_rotation_energy
    
    !> \f[ \Delta U_{N\rightarrow{}N+1} = -(\vec{\mu}_{N+1} \cdot \vec{E}) \f]

    pure function External_Field_exchange_energy(this, particle) result(exchange_energy)

        class(External_Field), intent(in) :: this
        type(Particle_Index), intent(in) :: particle
        real(DP) :: exchange_energy

        exchange_energy = -exchange_sign(particle%add) * dot_product(particle%orientation, this%vector)

    end function External_Field_exchange_energy
    
end module class_external_field
