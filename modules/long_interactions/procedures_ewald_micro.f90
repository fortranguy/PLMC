module procedures_ewald_micro

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions, PI

implicit none
private
public :: ewald_real_B, ewald_real_C, &
       reciprocal_size_1_sym, reciprocal_size_2_sym, set_fourier, set_exp_kz

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

    !> Symmetry: half wave vectors in do loop: reciprocal_size_1
    pure function reciprocal_size_1_sym(reciprocal_size, ky, kz)
        integer :: reciprocal_size_1_sym
        integer, dimension(:), intent(in) :: reciprocal_size
        integer, intent(in) :: ky, kz

        if (ky == 0 .and. kz == 0) then
            reciprocal_size_1_sym = 0
        else
            reciprocal_size_1_sym = reciprocal_size(1)
        end if
    end function reciprocal_size_1_sym

    !> Symmetry: half wave vectors in do loop: reciprocal_size_2
    pure function reciprocal_size_2_sym(reciprocal_size, kz)
        integer :: reciprocal_size_2_sym
        integer, dimension(:), intent(in) :: reciprocal_size
        integer, intent(in) :: kz

        if (kz == 0) then
            reciprocal_size_2_sym = 0
        else
            reciprocal_size_2_sym = reciprocal_size(2)
        end if
    end function reciprocal_size_2_sym

    !> Fourier coefficients (bases) tabulation
    pure subroutine set_fourier(fourier_position_i, reci_number_i, wave_1_dot_position_i)
        integer, intent(in) :: reci_number_i
        complex(DP), dimension(-reci_number_i:reci_number_i), intent(out) :: fourier_position_i
        real(DP), intent(in) :: wave_1_dot_position_i

        integer :: n_i
        fourier_position_i(0) = (1._DP, 0._DP)
        fourier_position_i(1) = cmplx(cos(wave_1_dot_position_i), sin(wave_1_dot_position_i), DP)
        fourier_position_i(-1) = conjg(fourier_position_i(1))
        do n_i = 2, reci_number_i
            fourier_position_i(n_i) = fourier_position_i(n_i-1) * fourier_position_i(1)
            fourier_position_i(-n_i) = conjg(fourier_position_i(n_i))
        end do
    end subroutine set_fourier

    !> DLC tabulation
    pure subroutine set_exp_kz(exp_kz_tab, reciprocal_size, wave_norm, z)
        integer, dimension(:), intent(in) :: reciprocal_size
        real(DP), dimension(0:reciprocal_size(1), 0:reciprocal_size(2)), intent(out) :: exp_kz_tab
        real(DP), dimension(-reciprocal_size(1):reciprocal_size(1), &
                            -reciprocal_size(2):reciprocal_size(2)), intent(in) :: wave_norm
        real(DP), intent(in) :: z

        integer :: kx, ky
        do ky = 0, reciprocal_size(2)
            do kx = ky, reciprocal_size(1)
                exp_kz_tab(kx, ky) = exp(wave_norm(kx, ky) * z)
            end do
            do kx = 0, ky-1
                exp_kz_tab(kx, ky) = exp_kz_tab(ky, kx)
            end do
        end do
    end subroutine set_exp_kz

end module procedures_ewald_micro
