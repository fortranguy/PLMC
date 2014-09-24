module class_distribution

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file

implicit none

private

    type, public :: Distribution
    
        real(DP) :: distance_max
        real(DP) :: distance_i
        real(DP) :: delta
        
        integer :: num_distribution
        integer :: i_distribution
        
        real(DP), dimension(:, :), allocatable :: distribution_step
        real(DP), dimension(:, :), allocatable :: distribution_function
    
    contains
    
        procedure :: contruct => Distribution_contruct
        procedure :: destroy => Distribution_destroy
    
    end type Distribution
    
contains

    subroutine Distribution_contruct(this, data_post_json, distance_max, num_observables)
    
        class(Distribution), intent(out) :: this
        type(json_file), intent(inout) :: data_post_json
        real(DP), intent(in) :: distance_max
        integer, intent(in), optional :: num_observables
        
        character(len=4096) :: data_name
        logical :: found
        
        data_name = "Distribution.delta"
        call data_post_json%get(data_name, this%delta, found)
        call test_data_found(data_name, found)
        
        this%distance_max = distance_max
        this%num_distribution = int(this%distance_max/this%delta)

        allocate(this%distribution_step(this%num_distribution, num_observables))
        allocate(this%distribution_function(this%num_distribution, num_observables))
    
    end subroutine Distribution_contruct
    
    subroutine Distribution_destroy(this)
    
        class(Distribution), intent(inout) :: this
        
        if (allocated(this%distribution_step)) deallocate(this%distribution_step)
        if (allocated(this%distribution_function)) deallocate(this%distribution_function)
        
    end subroutine Distribution_destroy

end module class_distribution
