module module_linear_algebra

use, intrinsic :: iso_fortran_env, only: DP => REAL64

implicit none
private
public identity_matrix, projector_matrix, eigen_symmetric

contains

    pure function identity_matrix(num_dimensions)
    
        integer, intent(in) :: num_dimensions
        real(DP), dimension(num_dimensions, num_dimensions) :: identity_matrix
        
        integer :: i_dim
        
        identity_matrix(:, :) = 0._DP
        forall (i_dim = 1:num_dimensions) identity_matrix(i_dim, i_dim) = 1._DP
    
    end function identity_matrix
    
    pure function projector_matrix(num_dimensions, i_dimension)
    
        integer, intent(in) :: num_dimensions, i_dimension
        real(DP), dimension(num_dimensions, num_dimensions) :: projector_matrix
        
        projector_matrix(:, :) = 0._DP
        projector_matrix(i_dimension, i_dimension) = 1._DP
        
    end function projector_matrix

    subroutine eigen_symmetric(matrix, eigenvalues)
    
        real(DP), dimension(:, :), intent(inout) :: matrix
        real(DP), dimension(:), intent(out) :: eigenvalues
        
        real(DP), dimension(3*size(matrix, 1)-1) :: work_space
        integer :: info
    
        call dsyev('V', 'U', &
                   size(matrix, 1), matrix, size(matrix, 1), eigenvalues, &
                   work_space, 3*size(matrix, 1)-1, info)
    
    end subroutine eigen_symmetric

end module module_linear_algebra
