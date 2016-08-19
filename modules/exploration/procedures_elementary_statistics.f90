module procedures_elementary_statistics

use, intrinsic :: iso_fortran_env, only: DP => REAL64

implicit none

private
public :: average, standard_deviation

contains

    !> \[ \langle x \rangle = \frac{1}{N} \sum_{i = 1}^N x_i \]
    pure real(DP) function average(array)
        real(DP), intent(in) :: array(:)

        average = sum(array) / size(array)
    end function average

    !> \[ \Delta x = \sqrt{\frac{1}{N - 1} \sum_{i = 1}^N (x_i - \langle x \rangle)^2} \]
    pure real(DP) function standard_deviation(array)
        real(DP), intent(in) :: array(:)

        standard_deviation = sqrt(sum((array - average(array))**2) / (size(array) - 1))
    end function standard_deviation

end module procedures_elementary_statistics
