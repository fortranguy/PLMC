module module_linear_algebra

use, intrinsic :: iso_fortran_env, only: DP => REAL64

implicit none
private
public eigen_symmetric

contains

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
