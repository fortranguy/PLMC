!> \brief Description of the Dipolar_Hard_Spheres class

module class_dipolar_hard_spheres

use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
use data_precisions, only: DP, real_zero, consist_tiny
use data_constants, only: PI
use data_box, only: Ndim
use json_module, only: json_file
use data_monte_carlo, only: dipol_rotate_delta, dipol_rotate_deltaMax, dipol_rotate_rejectFix
use data_potential, only: dipol_rMin_factor, dipol_real_rCut_factor, dipol_real_dr, dipol_alpha_factor
use data_neighbour_cells, only: NnearCell
use data_distribution, only: snap_ratio
use module_types_micro, only: Box_Dimensions, Particle_Index
use module_physics_micro, only: set_discrete_length, PBC_vector, Box_wave1_sym, Box_wave2_sym, &
                                fourier_i, exchange_sign
use module_data, only: test_data_found
use class_hard_spheres

implicit none

private

    type, extends(Hard_Spheres), public :: Dipolar_Hard_Spheres

        private
        
        ! Particles
        real(DP), dimension(:, :), allocatable, public :: all_orientations
        ! Potential
        real(DP) :: alpha !< coefficient of Ewald summation
        real(DP), dimension(Ndim) :: totalMoment
        
    contains

        !> Construction and destruction of the class
        procedure :: construct => Dipolar_Hard_Spheres_construct
        procedure, private :: set_particles => Dipolar_Hard_Spheres_set_particles
        procedure :: destroy => Dipolar_Hard_Spheres_destroy
        
        !> Write a report of the component in a file
        procedure :: write_report => Dipolar_Hard_Spheres_write_report
        
        !> Accessor & Mutator
        procedure :: get_orientation => Dipolar_Hard_Spheres_get_orientation
        procedure :: set_orientation => Dipolar_Hard_Spheres_set_orientation
        procedure :: Epot_set_alpha => Dipolar_Hard_Spheres_Epot_set_alpha
        
        !> Take a snap shot of the configuration: orientations
        procedure :: write_snap_orientations => Dipolar_Hard_Spheres_write_snap_orientations
        
        !> Potential energy
        !>     Real
        !>     Reciprocal
        !>     Self
        procedure :: Epot_self_solo => Dipolar_Hard_Spheres_Epot_self_solo
        procedure, private :: Epot_self => Dipolar_Hard_Spheres_Epot_self
        !>     Total moment
        procedure, private :: set_totalMoment => Dipolar_Hard_Spheres_set_totalMoment
        procedure :: reset_totalMoment => Dipolar_Hard_Spheres_reset_totalMoment
        procedure :: update_totalMoment_rotate => Dipolar_Hard_Spheres_update_totalMoment_rotate
        !>     Boundary conditions
        procedure :: deltaEpot_bound_rotate => Dipolar_Hard_Spheres_deltaEpot_bound_rotate
        procedure :: deltaEpot_bound_exchange => Dipolar_Hard_Spheres_deltaEpot_bound_exchange
        procedure, private :: Epot_bound => Dipolar_Hard_Spheres_Epot_bound
        !>     Total
        procedure :: set_Epot => Dipolar_Hard_Spheres_set_Epot
        procedure :: Epot_conf => Dipolar_Hard_Spheres_Epot_conf
        
    end type Dipolar_Hard_Spheres
    
contains

    subroutine Dipolar_Hard_Spheres_construct(this, json)
    
        class(Dipolar_Hard_Spheres), intent(out) :: this
        type(json_file), intent(inout) :: json
        
        this%name = "dipol"
        write(output_unit, *) this%name, " class construction"
    
        call this%set_particles(json)        
        call this%Hard_Spheres%set_snap(json)
    
    end subroutine Dipolar_Hard_Spheres_construct

    subroutine Dipolar_Hard_Spheres_set_particles(this, json)
        class(Dipolar_Hard_Spheres), intent(inout) :: this
        type(json_file), intent(inout) :: json
        
        character(len=4096) :: data_name
        logical :: found
        
        this%diameter = 1._DP ! = u_length
        
        data_name = "Particles.Dipoles.number of particles"
        call json%get(data_name, this%num_particles, found)
        call test_data_found(data_name, found)
        allocate(this%all_positions(Ndim, this%num_particles))
        allocate(this%all_orientations(Ndim, this%num_particles))
        
        data_name = "Particles.Dipoles.number of Widom particles"
        call json%get(data_name, this%widom_num_particles, found)
        call test_data_found(data_name, found)
        
    end subroutine Dipolar_Hard_Spheres_set_particles
    
    subroutine Dipolar_Hard_Spheres_destroy(this)    
        class(Dipolar_Hard_Spheres), intent(inout) :: this
        
        call this%Hard_Spheres%destroy()
        if (allocated(this%all_orientations)) deallocate(this%all_orientations)   
    end subroutine Dipolar_Hard_Spheres_destroy
    
    !> Report
    
    subroutine Dipolar_Hard_Spheres_write_report(this, report_unit)
    
        class(Dipolar_Hard_Spheres), intent(in) :: this
        integer, intent(in) :: report_unit
        
        call this%Hard_Spheres%write_report(report_unit)

        write(report_unit, *) "    alpha = ", this%alpha
        
    end subroutine Dipolar_Hard_Spheres_write_report
    
    !> Accessors & Mutators
    
    pure function Dipolar_Hard_Spheres_get_orientation(this, i_particle) result(get_orientation)
        class(Dipolar_Hard_Spheres), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP), dimension(Ndim) :: get_orientation
        
        get_orientation(:) = this%all_orientations(:, i_particle)
    end function Dipolar_Hard_Spheres_get_orientation
    
    subroutine Dipolar_Hard_Spheres_set_orientation(this, i_particle, orientation)
        class(Dipolar_Hard_Spheres), intent(inout) :: this
        integer, intent(in) :: i_particle
        real(DP), dimension(:), intent(in) :: orientation
        
        this%all_orientations(:, i_particle) = orientation(:)
    end subroutine Dipolar_Hard_Spheres_set_orientation
    
    subroutine Dipolar_Hard_Spheres_Epot_set_alpha(this, Box_size)
        class(Dipolar_Hard_Spheres), intent(inout) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        
        this%alpha = dipol_alpha_factor / Box_size(1)
    end subroutine Dipolar_Hard_Spheres_Epot_set_alpha
    
    !> Configuration state: orientations
      
    subroutine Dipolar_Hard_Spheres_write_snap_orientations(this, iStep, snap_unit)
        
        class(Dipolar_Hard_Spheres), intent(in) :: this
        integer, intent(in) :: iStep
        integer, intent(in) :: snap_unit
    
        integer :: i_particle
        
        if (modulo(iStep, this%snap_factor) == 0) then
            do i_particle = 1, this%num_particles
                write(snap_unit, *) this%all_orientations(:, i_particle)
            end do
        end if

    end subroutine Dipolar_Hard_Spheres_write_snap_orientations
    
    ! Self: correction ----------------------------------------------------------------------------
    
    !> Self energy of 1 dipole
    !> \f[ \frac{2}{3}\frac{\alpha^3}{\sqrt{\pi}} \vec{\mu}_i\cdot\vec{\mu}_i \f]
    
    pure function Dipolar_Hard_Spheres_Epot_self_solo(this, mCol) result(Epot_self_solo)
    
        class(Dipolar_Hard_Spheres), intent(in) :: this
        real(DP), dimension(:), intent(in) :: mCol
        real(DP) :: Epot_self_solo
        
        Epot_self_solo = 2._DP/3._DP * this%alpha**3/sqrt(PI) * dot_product(mCol, mCol)
    
    end function Dipolar_Hard_Spheres_Epot_self_solo
    
    !> Total self energy
    !> \f[ \frac{2}{3}\frac{\alpha^3}{\sqrt{\pi}} \sum_i \vec{\mu}_i\cdot\vec{\mu}_i \f]
    
    pure function Dipolar_Hard_Spheres_Epot_self(this) result(Epot_self)
    
        class(Dipolar_Hard_Spheres), intent(in) :: this
        real(DP) :: Epot_self

        integer :: i_particle
        
        Epot_self = 0._DP
        do i_particle = 1, this%num_particles
            Epot_self = Epot_self + this%Epot_self_solo(this%all_orientations(:, i_particle))
        end do
        
    end function Dipolar_Hard_Spheres_Epot_self

    ! Total moment ---------------------------------------------------------------------------------

    !> Total dipole moment :
    !> \f[ \vec{M} = \sum_j \vec{\mu}_j \f]
    !> \f[ \vec{M}_\underline{l} = \sum_{j \neq l} \vec{\mu}_j \f]
    
    pure subroutine Dipolar_Hard_Spheres_set_totalMoment(this)
        class(Dipolar_Hard_Spheres), intent(inout) :: this
        integer :: i_particle
        this%totalMoment(:) = 0._DP
        do i_particle = 1, this%num_particles
            this%totalMoment(:) = this%totalMoment(:) + this%all_orientations(:, i_particle)
        end do
    end subroutine Dipolar_Hard_Spheres_set_totalMoment

    !> Reinitialise the total moment factor and write the drift

    subroutine Dipolar_Hard_Spheres_reset_totalMoment(this, iStep, modulus_unit)

        class(Dipolar_Hard_Spheres), intent(inout) :: this
        integer, intent(in) :: iStep
        integer, intent(in) :: modulus_unit

        real(DP) :: modulus_drifted, modulus_reInit

        modulus_drifted = norm2(this%totalMoment(:))
        call this%set_totalMoment()
        modulus_reInit = norm2(this%totalMoment(:))

        write(modulus_unit, *) iStep, abs(modulus_reInit - modulus_drifted)

    end subroutine Dipolar_Hard_Spheres_reset_totalMoment

    !> Rotation

    !> Update the total moment
    !> \f[
    !>      \Delta \vec{M} = \vec{\mu}^\prime_l - \vec{\mu}_l
    !> \f]

    pure subroutine Dipolar_Hard_Spheres_update_totalMoment_rotate(this, mOld, mNew)

        class(Dipolar_Hard_Spheres), intent(inout) :: this
        real(DP), dimension(:), intent(in) :: mOld, mNew

        this%totalMoment(:) = this%totalMoment(:) + mNew(:) - mOld(:)

    end subroutine Dipolar_Hard_Spheres_update_totalMoment_rotate

    ! Boundary conditions: shape-dependent -------------------------------------------------------
    
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
    
    pure function Dipolar_Hard_Spheres_deltaEpot_bound_exchange(this, Box_size, mCol) &
                  result (deltaEpot_bound_exchange)
    
        class(Dipolar_Hard_Spheres), intent(in) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        real(DP), dimension(:), intent(in) :: mCol
        real(DP) :: deltaEpot_bound_exchange
        
        deltaEpot_bound_exchange = dot_product(mCol, mCol) + &
                                   2._DP * dot_product(mCol, this%totalMoment)
                          
        deltaEpot_bound_exchange = 2._DP*PI/3._DP/product(Box_size) * deltaEpot_bound_exchange
    
    end function Dipolar_Hard_Spheres_deltaEpot_bound_exchange
    
    !> Rotation
    
    !> Difference of Energy
    !> \f[
    !>      \Delta U = \frac{2\pi}{3V} [
    !>                      (\vec{\mu}^\prime_l \cdot \vec{\mu}^\prime_l) -
    !>                      (\vec{\mu}_l \cdot \vec{\mu}_l) +
    !>                      2 (\vec{\mu}^\prime_l - \vec{\mu}_l) \cdot \vec{M}_\underline{l}
    !>                 ]
    !> \f]
    
    pure function Dipolar_Hard_Spheres_deltaEpot_bound_rotate(this, Box_size, mOld, mNew) &
                  result (deltaEpot_bound_rotate)
    
        class(Dipolar_Hard_Spheres), intent(in) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        real(DP), dimension(:), intent(in) :: mOld, mNew
        real(DP) :: deltaEpot_bound_rotate
        
        deltaEpot_bound_rotate = dot_product(mNew, mNew) - dot_product(mOld, mOld) + &
                                 2._DP*dot_product(mNew-mOld, this%totalMoment-mOld)
                          
        deltaEpot_bound_rotate = 2._DP*PI/3._DP/product(Box_size) * deltaEpot_bound_rotate
    
    end function Dipolar_Hard_Spheres_deltaEpot_bound_rotate
    
    !> Total shape dependent term
    !> \f[
    !>      J(\vec{M}, S) = \frac{2\pi}{3V} | \vec{M}|^2
    !> \f]
    
    pure function Dipolar_Hard_Spheres_Epot_bound(this, Box_size) result(Epot_bound)
    
        class(Dipolar_Hard_Spheres), intent(in) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        real(DP) :: Epot_bound
        
        Epot_bound = 2._DP*PI/3._DP/product(Box_size) * dot_product(this%totalMoment, this%totalMoment)
    
    end function Dipolar_Hard_Spheres_Epot_bound
    
    ! Total potential energy ----------------------------------------------------------------------
    
    !> Potential energy initialisation
    
    subroutine Dipolar_Hard_Spheres_set_Epot(this, Box)
    
        class(Dipolar_Hard_Spheres), intent(inout) :: this
        type(Box_Dimensions), intent(in) :: Box
        
        call this%Epot_set_alpha(Box%size)
        call this%set_totalMoment()
        
    end subroutine Dipolar_Hard_Spheres_set_Epot

    !> Total potential energy of a configuration
    
    pure function Dipolar_Hard_Spheres_Epot_conf(this, Box) result(Epot_conf)
    
        class(Dipolar_Hard_Spheres), intent(in) :: this
        type(Box_Dimensions), intent(in) :: Box
        real(DP) :: Epot_conf
        
        Epot_conf = - this%Epot_self() + &
                    this%Epot_bound(Box%size) ! this%Epot_real(Box%size) + this%Epot_reci(Box)
    
    end function Dipolar_Hard_Spheres_Epot_conf

end module class_dipolar_hard_spheres
