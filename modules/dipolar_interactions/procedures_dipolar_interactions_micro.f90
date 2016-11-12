module procedures_dipolar_interactions_micro

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions, PI

implicit none
private
public :: des_real_B, des_real_C, reci_number_1_sym, reci_number_2_sym, set_fourier, set_exp_kz

contains

    !> \[
    !>      B_\alpha(r) = \frac{\mathrm{erfc}(\alpha r)}{r^3} +
    !>          \frac{2\alpha}{\sqrt{\pi}}\frac{e^{-\alpha^2 r^2}}{r^2}
    !> \]
    pure function des_real_B(alpha, r)
        real(DP) :: des_real_B
        real(DP), intent(in) :: alpha
        real(DP), intent(in) :: r

        des_real_B = erfc(alpha*r)/r**3 + 2._DP*alpha/sqrt(PI) * exp(-alpha**2*r**2) / r**2
    end function des_real_B

    !> \[
    !>      C_\alpha(r) = 3\frac{\mathrm{erfc}(\alpha r)}{r^5} +
    !>          \frac{2\alpha}{\sqrt{\pi}}\left(2\alpha^2 + \frac{3}{r^2}\right)
    !>                                     \frac{e^{-\alpha^2 r^2}}{r^2}
    !> \]
    pure function des_real_C(alpha, r)
        real(DP) :: des_real_C
        real(DP), intent(in) :: alpha
        real(DP), intent(in) :: r

        des_real_C = 3._DP*erfc(alpha*r)/r**5 + &
                       2._DP*alpha/sqrt(PI) * (2._DP*alpha**2+3._DP/r**2) * &
                                              exp(-alpha**2*r**2) / r**2
    end function des_real_C

    !> Symmetry: half wave vectors in do loop: reci_number_1
    pure integer function reci_number_1_sym(reci_numbers, n_3, n_2)
        integer, dimension(:), intent(in) :: reci_numbers
        integer, intent(in) :: n_3, n_2

        if (n_2 == 0 .and. n_3 == 0) then
            reci_number_1_sym = 0
        else
            reci_number_1_sym = reci_numbers(1)
        end if
    end function reci_number_1_sym

    !> Symmetry: half wave vectors in do loop: reci_number_2
    pure integer function reci_number_2_sym(reci_numbers, n_3)
        integer, dimension(:), intent(in) :: reci_numbers
        integer, intent(in) :: n_3

        if (n_3 == 0) then
            reci_number_2_sym = 0
        else
            reci_number_2_sym = reci_numbers(2)
        end if
    end function reci_number_2_sym

    !> Set fourier coefficients \( e^{i k_\mathsf{i} x_\mathsf{i}} \) tabulation
    pure subroutine set_fourier(fourier_position_i, reci_number_i, wave_1_x_position_i)
        integer, intent(in) :: reci_number_i
        complex(DP), dimension(-reci_number_i:reci_number_i), intent(out) :: fourier_position_i
        real(DP), intent(in) :: wave_1_x_position_i

        integer :: n_i
        fourier_position_i(0) = (1._DP, 0._DP)
        fourier_position_i(1) = cmplx(cos(wave_1_x_position_i), sin(wave_1_x_position_i), DP)
        fourier_position_i(-1) = conjg(fourier_position_i(1))
        do n_i = 2, reci_number_i
            fourier_position_i(n_i) = fourier_position_i(n_i-1) * fourier_position_i(1)
            fourier_position_i(-n_i) = conjg(fourier_position_i(n_i))
        end do
    end subroutine set_fourier

    !> Set \( e^{k_{1:2} z} \) tabulation
    pure subroutine set_exp_kz(exp_kz_tab, surface_size, position_3)
        real(DP), dimension(0:, 0:), intent(out) :: exp_kz_tab
        real(DP), intent(in) :: surface_size(:), position_3

        real(DP) :: wave_vector(2)
        integer :: n_1, n_2

        do n_2 = 0, ubound(exp_kz_tab, 2)
            wave_vector(2) = 2._DP*PI * real(n_2, DP) / surface_size(2)
            do n_1 = n_2, ubound(exp_kz_tab, 1)
                wave_vector(1) = 2._DP*PI * real(n_1, DP) / surface_size(1)
                exp_kz_tab(n_1, n_2) = exp(norm2(wave_vector) * position_3)
            end do
            do n_1 = 0, n_2-1
                exp_kz_tab(n_1, n_2) = exp_kz_tab(n_2, n_1)
            end do
        end do
    end subroutine set_exp_kz

end module procedures_dipolar_interactions_micro
