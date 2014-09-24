module class_distribution_function

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_precisions, only: real_zero
use json_module, only: json_file
use module_data, only: test_data_found

implicit none

private

    type, public :: Distribution_Function
    
        private
    
        real(DP) :: delta
        
        integer :: num_distribution
        
        real(DP), dimension(:, :), allocatable :: distribution_step
        real(DP), dimension(:, :), allocatable :: distribution_function
        integer :: num_observables
    
    contains
    
        procedure :: construct => Distribution_Function_construct
        procedure :: destroy => Distribution_Function_destroy
        procedure :: write => Distribution_Function_write
        
        procedure :: step_init => Distribution_Function_step_init
        procedure :: particle_set => Distribution_Function_particle_set
        procedure :: step_set => Distribution_Function_step_set
    
    end type Distribution_Function
    
contains

    subroutine Distribution_Function_construct(this, data_post_json, distance_max, num_observables)
    
        class(Distribution_Function), intent(out) :: this
        type(json_file), intent(inout) :: data_post_json
        real(DP), intent(in) :: distance_max
        integer, intent(in), optional :: num_observables
        
        character(len=4096) :: data_name
        logical :: found
        
        data_name = "Distribution.delta"
        call data_post_json%get(data_name, this%delta, found)
        call test_data_found(data_name, found)
        
        this%num_distribution = int(distance_max/this%delta)
        
        this%num_observables = num_observables + 1
        allocate(this%distribution_step(this%num_distribution, this%num_observables))
        allocate(this%distribution_function(this%num_distribution, this%num_observables))
        
        this%distribution_function(:, :) = 0._DP
    
    end subroutine Distribution_Function_construct
    
    subroutine Distribution_Function_destroy(this)
    
        class(Distribution_Function), intent(inout) :: this
        
        if (allocated(this%distribution_step)) deallocate(this%distribution_step)
        if (allocated(this%distribution_function)) deallocate(this%distribution_function)
        
    end subroutine Distribution_Function_destroy
    
    subroutine Distribution_Function_write(this, num_steps, distribution_unit)
    
        class(Distribution_Function), intent(inout) :: this
        integer, intent(in) :: num_steps
        integer, intent(in) :: distribution_unit
        
        real(DP) :: distance_i
        integer :: i_distribution
        
        this%distribution_function(:, :) = this%distribution_function(:, :) / real(num_steps, DP)
        this%distribution_step(:, 1) = this%distribution_step(:, 1) / this%delta
        
        do i_distribution = 1, this%num_distribution
            distance_i = (real(i_distribution, DP) + 0.5_DP) * this%delta
            write(distribution_unit, *) distance_i, this%distribution_step(i_distribution, :)
        end do
    
    end subroutine Distribution_Function_write
    
    subroutine Distribution_Function_step_init(this)
    
        class(Distribution_Function), intent(inout) :: this
        
        this%distribution_step(:, :) = 0._DP
    
    end subroutine Distribution_Function_step_init
    
    subroutine Distribution_Function_particle_set(this, distance, data_vector)
    
        class(Distribution_Function), intent(inout) :: this
        real(DP), intent(in) :: distance
        real(DP), dimension(:), intent(in) :: data_vector
        
        integer :: i_distribution
        
        i_distribution = int(distance/this%delta)
        
        this%distribution_step(i_distribution, 1) = this%distribution_step(i_distribution, 1) + 1._DP
        this%distribution_step(i_distribution, 2:this%num_observables) = &
            this%distribution_step(i_distribution, 2:this%num_observables) + data_vector(:)
        
    end subroutine Distribution_Function_particle_set
    
    subroutine Distribution_Function_step_set(this)
    
        class(Distribution_Function), intent(inout) :: this
        
        integer :: i_observable
        
        do i_observable = 2, this%num_observables
            where(this%distribution_step(:, 1) > real_zero)            
                this%distribution_step(:, i_observable) = this%distribution_step(:, i_observable) / &
                                                          this%distribution_step(:, 1)            
            end where
        end do
        
        this%distribution_function(:, :) = this%distribution_function(:, :) + &
                                           this%distribution_step(:, :)
        
    end subroutine Distribution_Function_step_set

end module class_distribution_function
