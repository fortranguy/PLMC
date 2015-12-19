module procedures_ewald_micro

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions, PI

implicit none
private
public :: ewald_real_B, ewald_real_C, reci_number_1_sym, reci_number_2_sym, set_fourier, set_exp_n_3

contains

    !> \[
    !>      \frac{(\vec{\mu}_i\cdot\vec{\mu_j})}{r^3} -
    !>     3\frac{(\vec{\mu}_i\cdot\vec{r}_{ij}) (\vec{\mu}_j\cdot\vec{r}_{ij})}{r^5}
    !> \]
    pure function dipolar_pair_energy(orientation_i, orientation_j, vector_ij)
        real(DP) :: dipolar_pair_energy
        real(DP), dimension(:), intent(in) :: orientation_i, orientation_j
        real(DP), dimension(:), intent(in) :: vector_ij

        real(DP) :: distance_ij
        distance_ij = norm2(vector_ij)
        dipolar_pair_energy = dot_product(orientation_i, orientation_j) / distance_ij**3 - &
                              3._DP * dot_product(orientation_i, vector_ij) * &
                                      dot_product(orientation_j, vector_ij) / distance_ij**5
    end function dipolar_pair_energy

    !> \[
    !>      B(r) = \frac{\mathrm{erfc}(\alpha r)}{r^3} +
    !>           2\frac{\alpha}{\sqrt{\pi}}\frac{e^{-\alpha^2 r^2}}{r^2}
    !> \]
    pure function ewald_real_B(alpha, r)
        real(DP) :: ewald_real_B
        real(DP), intent(in) :: alpha
        real(DP), intent(in) :: r

        ewald_real_B = erfc(alpha*r)/r**3 + 2._DP*alpha/sqrt(PI) * exp(-alpha**2*r**2) / r**2
    end function ewald_real_B

    !> \[
    !>      C(r) = 3\frac{\mathrm{erfc}(\alpha r)}{r^5} +
    !>            2\frac{\alpha}{\sqrt{\pi}}\left(2\alpha^2 + \frac{3}{r^2}\right)
    !>                                     \frac{e^{-\alpha^2 r^2}}{r^2}
    !> \]
    pure function ewald_real_C(alpha, r)
        real(DP) :: ewald_real_C
        real(DP), intent(in) :: alpha
        real(DP), intent(in) :: r

        ewald_real_C = 3._DP*erfc(alpha*r)/r**5 + &
                       2._DP*alpha/sqrt(PI) * (2._DP*alpha**2+3._DP/r**2) * &
                                              exp(-alpha**2*r**2) / r**2
    end function ewald_real_C

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

    !> Fourier coefficients (bases) tabulation
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

    !> DLC tabulation
    pure subroutine set_exp_n_3(exp_n_3_tab, reci_number, wave_norm, z)
        integer, dimension(:), intent(in) :: reci_number
        real(DP), dimension(0:reci_number(1), 0:reci_number(2)), intent(out) :: exp_n_3_tab
        real(DP), dimension(-reci_number(1):reci_number(1), &
                            -reci_number(2):reci_number(2)), intent(in) :: wave_norm
        real(DP), intent(in) :: z

        integer :: n_1, n_2
        do n_2 = 0, reci_number(2)
            do n_1 = n_2, reci_number(1)
                exp_n_3_tab(n_1, n_2) = exp(wave_norm(n_1, n_2) * z)
            end do
            do n_1 = 0, n_2-1
                exp_n_3_tab(n_1, n_2) = exp_n_3_tab(n_2, n_1)
            end do
        end do
    end subroutine set_exp_n_3

end module procedures_ewald_micro
