module class_ewald_summation_real

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_box, only: num_dimensions
use json_module, only: json_file
use module_data, only: test_data_found
use module_types_micro, only: Particle_Index
use module_physics_micro, only: set_discrete_length, PBC_vector, ewald_real_B, ewald_real_C
use class_hard_spheres, only: Dipolar_Hard_Spheres

implicit none

private

    type, public :: Ewald_Summation_Real
    
        real(DP) :: min_distance
        real(DP) :: cutoff
        real(DP) :: delta
        
        integer :: i_min_distance
        integer :: i_cutoff
        real(DP), dimension(:, :), allocatable :: tabulation
    
    contains
    
        procedure :: construct => Ewald_Summation_Real_construct
        procedure, private :: set_parameters => Ewald_Summation_Real_set_parameters
        procedure, private :: set_tabulation => Ewald_Summation_Real_set_tabulation
        procedure :: destroy => Ewald_Summation_Real_destroy
        procedure :: write => Ewald_Summation_Real_write
        
        procedure :: total_energy => Ewald_Summation_Real_total_energy
        procedure :: solo_energy => Ewald_Summation_Real_solo_energy
        procedure, private :: pair_energy => Ewald_Summation_Real_pair_energy
        procedure, private :: interpolation => Ewald_Summation_Real_interpolation
        
        procedure :: solo_field => Ewald_Summation_Real_solo_field
    
    end type Ewald_Summation_Real

contains

    subroutine Ewald_Summation_Real_construct(this, Box_size, alpha, min_distance, data_json)
    
        class(Ewald_Summation_Real), intent(inout) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        real(DP), intent(in) :: alpha
        real(DP), intent(in) :: min_distance
        type(json_file), intent(inout) :: data_json

        call this%set_parameters(Box_size, min_distance, data_json)

        if (allocated(this%tabulation)) deallocate(this%tabulation)
        allocate(this%tabulation(this%i_min_distance:this%i_cutoff, 2))

        call this%set_tabulation(alpha)

    end subroutine Ewald_Summation_Real_construct

    subroutine Ewald_Summation_Real_set_parameters(this, Box_size, min_distance, data_json)
    
        class(Ewald_Summation_Real), intent(inout) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        real(DP), intent(in) :: min_distance
        type(json_file), intent(inout) :: data_json

        character(len=4096) :: data_name
        logical :: found

        real(DP) :: cutoff_factor
        
        this%min_distance = min_distance

        data_name = "Potential Energy.Dipolar Hard Spheres.Ewald summation.real.cut off factor"
        call data_json%get(data_name, cutoff_factor, found)
        call test_data_found(data_name, found)
        this%cutoff = cutoff_factor * Box_size(1)

        data_name = "Potential Energy.Dipolar Hard Spheres.Ewald summation.real.delta"
        call data_json%get(data_name, this%delta, found)
        call test_data_found(data_name, found)
        call set_discrete_length(this%min_distance, this%delta)
        
        this%i_min_distance = int(this%min_distance/this%delta)
        this%i_cutoff = int(this%cutoff/this%delta) + 1

    end subroutine Ewald_Summation_Real_set_parameters

    pure subroutine Ewald_Summation_Real_set_tabulation(this, alpha)

        class(Ewald_Summation_Real), intent(inout) :: this
        real(DP), intent(in) :: alpha

        integer :: i_distance
        real(DP) :: distance_i
        
        ! cut
        do i_distance = this%i_min_distance, this%i_cutoff
            distance_i = real(i_distance, DP)*this%delta
            this%tabulation(i_distance, 1) = ewald_real_B(alpha, distance_i)
            this%tabulation(i_distance, 2) = ewald_real_C(alpha, distance_i)
        end do

        ! shift
        this%tabulation(:, 1) = this%tabulation(:, 1) - this%tabulation(this%i_cutoff, 1)
        this%tabulation(:, 2) = this%tabulation(:, 2) - this%tabulation(this%i_cutoff, 2)

    end subroutine Ewald_Summation_Real_set_tabulation

    subroutine Ewald_Summation_Real_destroy(this)
    
        class(Ewald_Summation_Real), intent(inout) :: this

        if (allocated(this%tabulation)) deallocate(this%tabulation)

    end subroutine Ewald_Summation_Real_destroy

    subroutine Ewald_Summation_Real_write(this, potential_energy_unit)
    
        class(Ewald_Summation_Real), intent(in) :: this
        integer, intent(in) :: potential_energy_unit

        integer :: i_distance
        real(DP) :: distance_i

        do i_distance = this%i_min_distance, this%i_cutoff
            distance_i = real(i_distance, DP)*this%delta
            write(potential_energy_unit, *) distance_i, this%tabulation(i_distance, :)
        end do

    end subroutine Ewald_Summation_Real_write

    pure function Ewald_Summation_Real_total_energy(this, Box_size, this_spheres) &
                  result(total_energy)

        class(Ewald_Summation_Real), intent(in) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        type(Dipolar_Hard_Spheres), intent(in) :: this_spheres
        real(DP) :: total_energy

        integer :: i_particle
        type(Particle_Index) :: particle

        total_energy = 0._DP
        
        do i_particle = 1, this_spheres%get_num_particles()
            particle%number = i_particle
            particle%position(:) = this_spheres%get_position(particle%number)
            particle%orientation(:) = this_spheres%get_orientation(particle%number)
            total_energy = total_energy + this%solo_energy(Box_size, this_spheres, particle)
        end do

        total_energy = total_energy/2._DP

    end function Ewald_Summation_Real_total_energy

    !> Energy of 1 dipole with others

    pure function Ewald_Summation_Real_solo_energy(this, Box_size, this_spheres, particle) &
                  result(solo_energy)

        class(Ewald_Summation_Real), intent(in) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        type(Dipolar_Hard_Spheres), intent(in) :: this_spheres
        type(Particle_Index), intent(in) :: particle
        real(DP) :: solo_energy

        integer :: j_particle
        real(DP), dimension(num_dimensions) :: vector_ij

        solo_energy = 0._DP
        
        do j_particle = 1, this_spheres%get_num_particles()
            if (j_particle /= particle%number) then
            
                vector_ij = PBC_vector(Box_size, &
                                       particle%position, this_spheres%get_position(j_particle))
                solo_energy = solo_energy + this%pair_energy(particle%orientation, &
                                                      this_spheres%get_orientation(j_particle), &
                                                      vector_ij)

            end if
        end do

    end function Ewald_Summation_Real_solo_energy

    !> Between 2 particles
    !> \f[ (\vec{\mu}_i\cdot\vec{\mu}_j) B(r_{ij}) -
    !>     (\vec{\mu}_i\cdot\vec{r}_{ij}) (\vec{\mu}_j\cdot\vec{r}_{ij}) C(r_{ij}) \f]

    pure function Ewald_Summation_Real_pair_energy(this, orientation_i, orientation_j, vector_ij) &
                  result(pair_energy)
                  
        class(Ewald_Summation_Real), intent(in) :: this
        real(DP), dimension(:), intent(in) :: orientation_i, orientation_j
        real(DP), dimension(:), intent(in) :: vector_ij
        real(DP) :: pair_energy

        real(DP), dimension(2) :: coefficient

        coefficient(1) = dot_product(orientation_i, orientation_j)
        coefficient(2) =-dot_product(orientation_i, vector_ij) * dot_product(orientation_j, vector_ij)

        pair_energy = dot_product(coefficient, this%interpolation(norm2(vector_ij)))

    end function Ewald_Summation_Real_pair_energy

    !> Linear interpolation

    pure function Ewald_Summation_Real_interpolation(this, distance) &
                  result(interpolation)
    
        class(Ewald_Summation_Real), intent(in) :: this
        real(DP), intent(in) :: distance
        real(DP), dimension(2) :: interpolation

        integer :: i_distance
        real(DP) :: distance_i

        if (distance < this%cutoff) then
            i_distance = int(distance/this%delta)
            distance_i = real(i_distance, DP)*this%delta
            interpolation(:) = this%tabulation(i_distance, :) + (distance-distance_i)/this%delta * &
                              (this%tabulation(i_distance+1, :) - this%tabulation(i_distance, :))
        else
            interpolation(:) = 0._DP
        end if

    end function Ewald_Summation_Real_interpolation
    
    pure function Ewald_Summation_Real_solo_field(this, Box_size, this_spheres, particle) &
                  result(solo_field)
    
        class(Ewald_Summation_Real), intent(in) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        type(Dipolar_Hard_Spheres), intent(in) :: this_spheres
        type(Particle_Index), intent(in) :: particle
        real(DP), dimension(num_dimensions) :: solo_field
        
        integer :: j_particle
        real(DP), dimension(num_dimensions) :: vector_ij
        
        solo_field(:) = 0._DP
        
        do j_particle = 1, this_spheres%get_num_particles()
            if (j_particle /= particle%number) then
            
                vector_ij = PBC_vector(Box_size, &
                                       particle%position, this_spheres%get_position(j_particle))

            end if
        end do
    
    end function Ewald_Summation_Real_solo_field

end module class_ewald_summation_real
