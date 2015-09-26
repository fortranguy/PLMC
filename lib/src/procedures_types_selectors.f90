module procedures_types_selectors

use class_particles_number, only: Abstract_Particles_Number, Null_Particles_Number
use class_particles_diameter, only: Abstract_Particles_Diameter, Null_Particles_Diameter
use class_particles_moment_norm, only: Abstract_Particles_Moment_Norm, Null_Particles_Moment_Norm
use class_particles_positions, only: Abstract_Particles_Positions, Null_Particles_Positions
use class_particles_orientations, only: Abstract_Particles_Orientations, Null_Particles_Orientations
use class_particles_dipolar_moments, only: Abstract_Particles_Dipolar_Moments, &
    Null_Particles_Dipolar_Moments
use class_particles_chemical_potential, only: Abstract_Particles_Chemical_Potential, &
    Null_Particles_Chemical_Potential
use class_particles_exchange, only: Abstract_Particles_Exchange, Null_Particles_Exchange

implicit none

private
public :: particles_exist, particles_have_positions, particles_have_orientations, &
    particles_are_dipolar, particles_have_chemical_potential

interface particles_exist
    module procedure :: particles_exist_from_number
    module procedure :: particles_exist_from_diameter
end interface particles_exist

interface particles_are_dipolar
    module procedure :: particles_are_dipolar_from_norm_and_orientations
    module procedure :: particle_are_dipolar_from_dipolar_moments
end interface particles_are_dipolar

contains

    pure logical function particles_exist_from_number(particles_number) result(particles_exist)
        class(Abstract_Particles_Number), intent(in) :: particles_number

        select type (particles_number)
            type is (Null_Particles_Number)
                particles_exist = .false.
            class default
                particles_exist = .true.
        end select
    end function particles_exist_from_number

    pure logical function particles_exist_from_diameter(particles_diameter) result(particles_exist)
        class(Abstract_Particles_Diameter), intent(in) :: particles_diameter

        select type (particles_diameter)
            type is (Null_Particles_Diameter)
                particles_exist = .false.
            class default
                particles_exist = .true.
        end select
    end function particles_exist_from_diameter

    pure logical function particles_have_moment_norm(particles_moment_norm)
        class(Abstract_Particles_Moment_Norm), intent(in) :: particles_moment_norm

        select type (particles_moment_norm)
            type is (Null_Particles_Moment_Norm)
                particles_have_moment_norm = .false.
            class default
                particles_have_moment_norm = .true.
        end select
    end function particles_have_moment_norm

    pure logical function particles_have_positions(partcles_positions)
        class(Abstract_Particles_Positions), intent(in) :: partcles_positions

        select type (partcles_positions)
            type is (Null_Particles_Positions)
                particles_have_positions = .false.
            class default
                particles_have_positions = .true.
        end select
    end function particles_have_positions

    pure logical function particles_have_orientations(particles_orientations)
        class(Abstract_Particles_Orientations), intent(in) :: particles_orientations

        select type (particles_orientations)
            type is (Null_Particles_Orientations)
                particles_have_orientations = .false.
            class default
                particles_have_orientations = .true.
        end select
    end function particles_have_orientations

    pure logical function particles_are_dipolar_from_norm_and_orientations(particles_moment_norm, &
        particles_orientations) result(particles_are_dipolar)
        class(Abstract_Particles_Moment_Norm), intent(in) :: particles_moment_norm
        class(Abstract_Particles_Orientations), intent(in) :: particles_orientations

        particles_are_dipolar = particles_have_moment_norm(particles_moment_norm) .and. &
            particles_have_orientations(particles_orientations)
    end function particles_are_dipolar_from_norm_and_orientations

    pure logical function particle_are_dipolar_from_dipolar_moments(particles_dipolar_moments) &
        result(particles_are_dipolar)
        class(Abstract_Particles_Dipolar_Moments), intent(in) :: particles_dipolar_moments

        select type (particles_dipolar_moments)
            type is (Null_Particles_Dipolar_Moments)
                particles_are_dipolar = .false.
            class default
                particles_are_dipolar = .true.
        end select
    end function particle_are_dipolar_from_dipolar_moments

    pure logical function particles_have_chemical_potential(particles_chemical_potential)
        class(Abstract_Particles_Chemical_Potential), intent(in) :: particles_chemical_potential

        select type (particles_chemical_potential)
            type is (Null_Particles_Chemical_Potential)
                particles_have_chemical_potential = .false.
            class default
                particles_have_chemical_potential = .true.
        end select
    end function particles_have_chemical_potential

end module procedures_types_selectors
