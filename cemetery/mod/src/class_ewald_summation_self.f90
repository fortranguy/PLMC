module class_ewald_summation_self

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: PI
use data_box, only: num_dimensions
use class_hard_spheres, only: Dipolar_Hard_Spheres

implicit none

private

    type, public :: Ewald_Summation_Self
        real(DP) :: alpha
    contains
        procedure :: set_alpha => Ewald_Summation_Self_set_alpha
        procedure :: total_energy => Ewald_Summation_Self_total_energy
        procedure, private :: total_energy_solo => Ewald_Summation_Self_total_energy_solo
        procedure, private :: total_energy_field => Ewald_Summation_Self_total_energy_field
        procedure :: solo_energy => Ewald_Summation_Self_solo_energy
        procedure :: solo_field => Ewald_Summation_Self_solo_field
        procedure :: test_field => Ewald_Summation_Self_test_field
    end type Ewald_Summation_Self

contains

    pure subroutine Ewald_Summation_Self_set_alpha(this, alpha)

        class(Ewald_Summation_Self), intent(inout) :: this
        real(DP), intent(in) :: alpha

        this%alpha = alpha

    end subroutine Ewald_Summation_Self_set_alpha

    !> Total self energy
    !> \[ \frac{2}{3}\frac{\alpha^3}{\sqrt{\pi}} \sum_i \vec{\mu}_i\cdot\vec{\mu}_i \]

    pure function Ewald_Summation_Self_total_energy(this, this_spheres, using_field) &
                  result(total_energy)

        class(Ewald_Summation_Self), intent(in) :: this
        class(Dipolar_Hard_Spheres), intent(in) :: this_spheres
        logical, intent(in), optional :: using_field
        real(DP) :: total_energy

        if (present(using_field)) then
            if (using_field) then
                total_energy = this%total_energy_field(this_spheres)
            else
                total_energy = this%total_energy_solo(this_spheres)
            end if
        else
            total_energy = this%total_energy_solo(this_spheres)
        end if

    end function Ewald_Summation_Self_total_energy

    pure function Ewald_Summation_Self_total_energy_solo(this, this_spheres) &
                  result(total_energy_solo)

        class(Ewald_Summation_Self), intent(in) :: this
        class(Dipolar_Hard_Spheres), intent(in) :: this_spheres
        real(DP) :: total_energy_solo

        integer :: i_particle

        total_energy_solo = 0._DP
        do i_particle = 1, this_spheres%get_num_particles()
            total_energy_solo = total_energy_solo + &
                                this%solo_energy(this_spheres%get_orientation(i_particle))
        end do

    end function Ewald_Summation_Self_total_energy_solo

    pure function Ewald_Summation_Self_total_energy_field(this, this_spheres) &
                  result(total_energy_field)

        class(Ewald_Summation_Self), intent(in) :: this
        class(Dipolar_Hard_Spheres), intent(in) :: this_spheres
        real(DP) :: total_energy_field

        integer :: i_particle

        total_energy_field = 0._DP
        do i_particle = 1, this_spheres%get_num_particles()
            total_energy_field = total_energy_field + &
                dot_product(this_spheres%get_orientation(i_particle), &
                            this%solo_field(this_spheres%get_orientation(i_particle)))
        end do

        total_energy_field = total_energy_field / 2._DP

    end function Ewald_Summation_Self_total_energy_field

    !> Self energy of 1 dipole
    !> \[ \frac{2}{3}\frac{\alpha^3}{\sqrt{\pi}} \vec{\mu}_i\cdot\vec{\mu}_i \]

    pure function Ewald_Summation_Self_solo_energy(this, orientation) result(solo_energy)

        class(Ewald_Summation_Self), intent(in) :: this
        real(DP), dimension(:), intent(in) :: orientation
        real(DP) :: solo_energy

        solo_energy = 2._DP/3._DP * this%alpha**3/sqrt(PI) * dot_product(orientation, orientation)

    end function Ewald_Summation_Self_solo_energy

    pure function Ewald_Summation_Self_solo_field(this, orientation) &
                  result(solo_field)

        class(Ewald_Summation_Self), intent(in) :: this
        real(DP), dimension(:), intent(in) :: orientation
        real(DP), dimension(num_dimensions) :: solo_field

        solo_field(:) = 4._DP/3._DP * this%alpha**3/sqrt(PI) * orientation(:)

    end function Ewald_Summation_Self_solo_field

    pure function Ewald_Summation_Self_test_field(this, orientation) &
                  result(test_field)

        class(Ewald_Summation_Self), intent(in) :: this
        real(DP), dimension(:), intent(in) :: orientation
        real(DP), dimension(num_dimensions) :: test_field

        test_field(:) = -2._DP/3._DP * this%alpha**3/sqrt(PI) * orientation(:)

    end function Ewald_Summation_Self_test_field

end module class_ewald_summation_self
