module class_identity_matrix

use, intrinsic :: iso_fortran_env, only: DP => REAL64

implicit none
private

    type, public :: Identity_Matrix
        real(DP), dimension(:, :), allocatable :: matrix
    contains
        procedure :: construct => Identity_Matrix_construct
        procedure :: destroy => Identity_Matrix_destroy
        procedure :: get => Identity_Matrix_get
    end type Identity_Matrix
    
contains

    subroutine Identity_Matrix_construct(this, num_dimensions)
    
        class(Identity_Matrix), intent(out) :: this
        integer, intent(in) :: num_dimensions
        
        integer :: i_dim
        
        allocate(this%matrix(num_dimensions, num_dimensions))
        
        this%matrix(:, :) = 0._DP
        forall (i_dim = 1:num_dimensions) this%matrix(i_dim, i_dim) = 1._DP
    
    end subroutine Identity_Matrix_construct
    
    subroutine Identity_Matrix_destroy(this)
    
        class(Identity_Matrix), intent(inout) :: this
        
        if (allocated(this%matrix)) deallocate(this%matrix)
    
    end subroutine Identity_Matrix_destroy
    
    pure function Identity_Matrix_get(this) result(get)
    
        class(Identity_Matrix), intent(in) :: this
        real(DP), dimension(size(this%matrix, 1), size(this%matrix, 2)) :: get
        
        get(:, :) = this%matrix(:, :)
    
    end function Identity_Matrix_get

end module class_identity_matrix
