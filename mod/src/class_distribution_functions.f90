module class_distribution_function

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use module_data, only: test_data_found

implicit none

private

    type, public :: Distribution_Function
    
        private
    
        real(DP) :: distance_max
        real(DP) :: distance_i
        real(DP) :: delta
        
        integer :: num_distribution
        integer :: i_distribution
        
        real(DP), dimension(:, :), allocatable :: distribution_step
        real(DP), dimension(:, :), allocatable :: distribution_function
        integer :: num_observables
    
    contains
    
        procedure :: construct => Distribution_Function_construct
        procedure :: destroy => Distribution_Function_destroy
        
        procedure :: init_at_step => Distribution_Function_init_at_step
        procedure :: set_at_step => Distribution_Function_set_at_step
    
    end type Distribution_Function
    
contains

    subroutine Distribution_Function_construct(this, data_post_json, distance_max, num_observables)
    
        class(Distribution_Function), intent(out) :: this
        type(json_file), intent(inout) :: data_post_json
        real(DP), intent(in) :: distance_max
        integer, intent(in), optional :: num_observables
        
        character(len=4096) :: data_name
        logical :: found
        
        data_name = "Distribution_Function.delta"
        call data_post_json%get(data_name, this%delta, found)
        call test_data_found(data_name, found)
        
        this%distance_max = distance_max
        this%num_distribution = int(this%distance_max/this%delta)
        
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
    
    subroutine Distribution_Function_init_at_step(this)
    
        class(Distribution_Function), intent(inout) :: this
        
        this%distribution_step(:, :) = 0._DP
    
    end subroutine Distribution_Function_init_at_step
    
    subroutine Distribution_Function_set_at_step(this, distance, vector)
    
        class(Distribution_Function), intent(inout) :: this
        real(DP), intent(in) :: distance
        real(DP), dimension(:), intent(in) :: vector
        
        integer :: i_distribution
        
        i_distribution = int(distance/this%delta)
        
        this%distribution_step(i_distribution, 1) = this%distribution_step(i_distribution, 1) + 1._DP
        this%distribution_step(i_distribution, 2:this%num_observables) = &
            this%distribution_step(i_distribution, 2:this%num_observables) + vector(:)
        
    end subroutine Distribution_Function_set_at_step

end module class_distribution_function
