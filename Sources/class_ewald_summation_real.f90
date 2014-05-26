module class_ewald_summation_real

use data_precisions, only: DP
use data_constants, only: PI
use data_box, only: Ndim
use json_module, only: json_file
use module_types_micro, only: Particle_Index
use module_physics_micro, only: set_discrete_length, PBC_vector, ewald_real_B, ewald_real_C
use module_data, only: test_data_found
use class_hard_spheres, only: Dipolar_Hard_Spheres

implicit none

private

    type, public :: Ewald_Summation_Real
    
        real(DP) :: min_distance
        real(DP) :: range_cut
        real(DP) :: delta
        
        integer :: i_min_distance
        integer :: i_range_cut        
        real(DP), dimension(:, :), allocatable :: tabulation
    
    contains
    
        procedure :: construct => Ewald_Summation_Real_construct
        procedure, private :: set_parameters => Ewald_Summation_Real_set_parameters
        procedure, private :: set_tabulation => Ewald_Summation_Real_set_tabulation
        procedure :: destroy => Ewald_Summation_Real_destroy
        procedure :: write => Ewald_Summation_Real_write
        procedure :: total => Ewald_Summation_Real_total
        procedure :: solo => Ewald_Summation_Real_solo
        procedure, private :: pair => Ewald_Summation_Real_pair
        procedure, private :: interpolation => Ewald_Summation_Real_interpolation
    
    end type Ewald_Summation_Real

contains

    subroutine Ewald_Summation_Real_construct(this, Box_size, alpha, min_distance, json)
        class(Ewald_Summation_Real), intent(inout) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        real(DP), intent(in) :: alpha
        real(DP), intent(in) :: min_distance
        type(json_file), intent(inout) :: json

        call this%set_parameters(Box_size, min_distance, json)

        if (allocated(this%tabulation)) deallocate(this%tabulation)
        allocate(this%tabulation(this%i_min_distance:this%i_range_cut, 2))

        call this%set_tabulation(alpha)

    end subroutine Ewald_Summation_Real_construct

    subroutine Ewald_Summation_Real_set_parameters(this, Box_size, min_distance, json)
        class(Ewald_Summation_Real), intent(inout) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        real(DP), intent(in) :: min_distance
        type(json_file), intent(inout) :: json

        character(len=4096) :: data_name
        logical :: found

        real(DP) :: range_cut_factor
        
        this%min_distance = min_distance

        data_name = "Potential.Dipoles.Ewald summation.real.range cut factor"
        call json%get(data_name, range_cut_factor, found)
        call test_data_found(data_name, found)
        this%range_cut = range_cut_factor * Box_size(1)

        data_name = "Potential.Dipoles.Ewald summation.real.delta"
        call json%get(data_name, this%delta, found)
        call test_data_found(data_name, found)
        call set_discrete_length(this%min_distance, this%delta)
        
        this%i_min_distance = int(this%min_distance/this%delta)
        this%i_range_cut = int(this%range_cut/this%delta) + 1

    end subroutine Ewald_Summation_Real_set_parameters

    pure subroutine Ewald_Summation_Real_set_tabulation(this, alpha)

        class(Ewald_Summation_Real), intent(inout) :: this
        real(DP), intent(in) :: alpha

        integer :: i_distance
        real(DP) :: distance_i
        
        ! cut
        do i_distance = this%i_min_distance, this%i_range_cut
            distance_i = real(i_distance, DP)*this%delta
            this%tabulation(i_distance, 1) = ewald_real_B(alpha, distance_i)
            this%tabulation(i_distance, 2) = ewald_real_C(alpha, distance_i)
        end do

        ! shift
        this%tabulation(:, 1) = this%tabulation(:, 1) - this%tabulation(this%i_range_cut, 1)
        this%tabulation(:, 2) = this%tabulation(:, 2) - this%tabulation(this%i_range_cut, 2)

    end subroutine Ewald_Summation_Real_set_tabulation

    subroutine Ewald_Summation_Real_destroy(this)
        class(Ewald_Summation_Real), intent(inout) :: this

        if (allocated(this%tabulation)) deallocate(this%tabulation)

    end subroutine Ewald_Summation_Real_destroy

    subroutine Ewald_Summation_Real_write(this, potential_unit)
        class(Ewald_Summation_Real), intent(in) :: this
        integer, intent(in) :: potential_unit

        integer :: i_distance
        real(DP) :: distance_i

        do i_distance = this%i_min_distance, this%i_range_cut
            distance_i = real(i_distance, DP)*this%delta
            write(potential_unit, *) distance_i, this%tabulation(i_distance, :)
        end do

    end subroutine Ewald_Summation_Real_write

    pure function Ewald_Summation_Real_total(this, Box_size, this_spheres) result(total)

        class(Ewald_Summation_Real), intent(in) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        type(Dipolar_Hard_Spheres), intent(in) :: this_spheres
        real(DP) :: total

        integer :: i_particle
        type(Particle_Index) :: particle

        total = 0._DP
        
        do i_particle = 1, this_spheres%get_num_particles()        
            particle%number = i_particle
            particle%position(:) = this_spheres%get_position(particle%number)
            particle%orientation(:) = this_spheres%get_orientation(particle%number)
            total = total + this%solo(Box_size, this_spheres, particle)            
        end do

        total = total/2._DP

    end function Ewald_Summation_Real_total

    !> Energy of 1 dipole with others

    pure function Ewald_Summation_Real_solo(this, Box_size, this_spheres, particle) result(solo)

        class(Ewald_Summation_Real), intent(in) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        type(Dipolar_Hard_Spheres), intent(in) :: this_spheres
        type(Particle_Index), intent(in) :: particle
        real(DP) :: solo

        integer :: j_particle
        real(DP), dimension(Ndim) :: vector_ij

        solo = 0._DP
        
        do j_particle = 1, this_spheres%get_num_particles()
            if (j_particle /= particle%number) then
            
                vector_ij = PBC_vector(Box_size, &
                                       particle%position, this_spheres%get_position(j_particle))
                solo = solo + this%pair(particle%orientation, &
                                        this_spheres%get_orientation(j_particle), &
                                        vector_ij)

            end if
        end do

    end function Ewald_Summation_Real_solo

    !> Between 2 particles
    !> \f[ (\vec{\mu}_i\cdot\vec{\mu}_j) B(r_{ij}) -
    !>     (\vec{\mu}_i\cdot\vec{r}_{ij}) (\vec{\mu}_j\cdot\vec{r}_{ij}) C(r_{ij}) \f]

    pure function Ewald_Summation_Real_pair(this, orientation_i, orientation_j, vector_ij) &
                  result(pair)
        class(Ewald_Summation_Real), intent(in) :: this
        real(DP), dimension(:), intent(in) :: orientation_i, orientation_j
        real(DP), dimension(:), intent(in) :: vector_ij
        real(DP) :: pair

        real(DP), dimension(2) :: coefficient

        coefficient(1) = dot_product(orientation_i, orientation_j)
        coefficient(2) =-dot_product(orientation_i, vector_ij) * dot_product(orientation_j, vector_ij)

        pair = dot_product(coefficient, this%interpolation(norm2(vector_ij)))

    end function Ewald_Summation_Real_pair

    !> Linear interpolation

    pure function Ewald_Summation_Real_interpolation(this, distance) result(interpolation)
        class(Ewald_Summation_Real), intent(in) :: this
        real(DP), intent(in) :: distance
        real(DP), dimension(2) :: interpolation

        integer :: i_distance
        real(DP) :: distance_i

        if (distance < this%range_cut) then
            i_distance = int(distance/this%delta)
            distance_i = real(i_distance, DP)*this%delta
            interpolation(:) = this%tabulation(i_distance, :) + (distance-distance_i)/this%delta * &
                              (this%tabulation(i_distance+1, :) - this%tabulation(i_distance, :))
        else
            interpolation(:) = 0._DP
        end if

    end function Ewald_Summation_Real_interpolation

end module class_ewald_summation_real
