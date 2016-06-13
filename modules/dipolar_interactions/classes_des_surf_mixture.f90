module classes_des_surf_mixture

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_permittivity, only: Abstract_Permittivity
use classes_mixture_total_moment, only: Abstract_Mixture_Total_Moment

implicit none

private

    type, abstract, public :: Abstract_DES_Surf_Mixture
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
        real(DP) :: permittivity = 0._DP
        class(Abstract_Mixture_Total_Moment), pointer :: total_moment => null()
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure(Abstract_visit), deferred :: visit
        procedure(Abstract_visit_transmutation), deferred :: visit_transmutation
        procedure :: visit_rotation => Abstract_visit_rotation
        procedure :: visit_add => Abstract_visit_add
        procedure :: visit_remove => Abstract_visit_remove
        procedure(Abstract_visit_exchange), deferred, private :: visit_exchange
    end type Abstract_DES_Surf_Mixture

    abstract interface

        pure real(DP) function Abstract_visit(this)
        import :: DP, Abstract_DES_Surf_Mixture
            class(Abstract_DES_Surf_Mixture), intent(in) :: this
        end function Abstract_visit

        pure real(DP) function Abstract_visit_transmutation(this, ij_components, dipolar_moment_2, &
            dipolar_moment_1)
        import :: DP, Abstract_DES_Surf_Mixture
            class(Abstract_DES_Surf_Mixture), intent(in) :: this
            integer, intent(in) :: ij_components(:)
            real(DP), intent(in) :: dipolar_moment_2(:), dipolar_moment_1(:)
        end function Abstract_visit_transmutation

        pure real(DP) function Abstract_visit_exchange(this, i_component, dipolar_moment, signed)
        import :: DP, Abstract_DES_Surf_Mixture
            class(Abstract_DES_Surf_Mixture), intent(in) :: this
            integer, intent(in) :: i_component
            real(DP), intent(in) :: dipolar_moment(:), signed
        end function Abstract_visit_exchange

    end interface

    type, extends(Abstract_DES_Surf_Mixture), public :: Spherical_DES_Surf_Mixture
    contains
        procedure :: visit => Spherical_visit
        procedure :: visit_transmutation => Spherical_visit_transmutation
        procedure, private :: visit_exchange => Spherical_visit_exchange
    end type Spherical_DES_Surf_Mixture

    type, extends(Abstract_DES_Surf_Mixture), public :: Rectangular_DES_Surf_Mixture
    contains
        procedure :: visit => Rectangular_visit
        procedure :: visit_transmutation => Rectangular_visit_transmutation
        procedure, private :: visit_exchange => Rectangular_visit_exchange
    end type Rectangular_DES_Surf_Mixture

    type, extends(Abstract_DES_Surf_Mixture), public :: Null_DES_Surf_Mixture
    contains
        procedure :: visit => Null_visit
        procedure :: visit_transmutation => Null_visit_transmutation
        procedure, private :: visit_exchange => Null_visit_exchange
    end type Null_DES_Surf_Mixture

contains

!implementation Abstract_DES_Surf_Mixture

    subroutine Abstract_construct(this, periodic_box, permittivity, total_moment)
        class(Abstract_DES_Surf_Mixture), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Permittivity), intent(in) :: permittivity
        class(Abstract_Mixture_Total_Moment), target, intent(in) :: total_moment

        this%periodic_box => periodic_box
        this%permittivity = permittivity%get()
        this%total_moment => total_moment
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_DES_Surf_Mixture), intent(inout) :: this

        this%total_moment => null()
        this%periodic_box => null()
    end subroutine Abstract_destroy

    pure real(DP) function Abstract_visit_rotation(this, i_component, dipolar_moment_2, &
        dipolar_moment_1) result(delta_energy)
        class(Abstract_DES_Surf_Mixture), intent(in) :: this
        integer, intent(in) :: i_component
        real(DP), intent(in) :: dipolar_moment_2(:), dipolar_moment_1(:)

        delta_energy = this%visit_transmutation([i_component, i_component], dipolar_moment_2, &
            dipolar_moment_1)
    end function Abstract_visit_rotation

    pure real(DP) function Abstract_visit_add(this, i_component, dipolar_moment) result(energy)
        class(Abstract_DES_Surf_Mixture), intent(in) :: this
        integer, intent(in) :: i_component
        real(DP), intent(in) :: dipolar_moment(:)

        energy = this%visit_exchange(i_component, dipolar_moment, +1.0_DP)
    end function Abstract_visit_add

    pure real(DP) function Abstract_visit_remove(this, i_component, dipolar_moment) result(energy)
        class(Abstract_DES_Surf_Mixture), intent(in) :: this
        integer, intent(in) :: i_component
        real(DP), intent(in) :: dipolar_moment(:)

        energy = this%visit_exchange(i_component, dipolar_moment, -1.0_DP)
    end function Abstract_visit_remove

!end implementation Abstract_DES_Surf_Mixture

!implementation Spherical_DES_Surf_Mixture

    !> \[ U = \frac{1}{6\epsilon V} \vec{M}^2 \]
    pure real(DP) function Spherical_visit(this) result(energy)
        class(Spherical_DES_Surf_Mixture), intent(in) :: this

        energy = 1._DP/6._DP/this%permittivity / product(this%periodic_box%get_size()) * &
            dot_product(this%total_moment%get(), this%total_moment%get())
    end function Spherical_visit

    !> \[
    !>      \Delta U = \frac{1}{6\epsilon V} [
    !>                      (\vec{\mu}^\prime_\mathsf{i} \cdot \vec{\mu}^\prime_\mathsf{i}) -
    !>                      (\vec{\mu}_\mathsf{i} \cdot \vec{\mu}_\mathsf{i}) +
    !>                      2 (\vec{\mu}^\prime_\mathsf{i} - \vec{\mu}_\mathsf{i}) \cdot
    !>                          (\vec{M} - \vec{\mu}_\mathsf{i})
    !>                 ]
    !> \]
    pure real(DP) function Spherical_visit_transmutation(this, ij_components, dipolar_moment_2, &
        dipolar_moment_1) result(delta_energy)
        class(Spherical_DES_Surf_Mixture), intent(in) :: this
        integer, intent(in) :: ij_components(:)
        real(DP), intent(in) :: dipolar_moment_2(:), dipolar_moment_1(:)

        delta_energy = 0._DP
        if (.not.(this%total_moment%is_dipolar(ij_components(1)) .or. &
            this%total_moment%is_dipolar(ij_components(2)))) return !shortcut?

        delta_energy = 1._DP/6._DP/this%permittivity / product(this%periodic_box%get_size()) * &
            (dot_product(dipolar_moment_2, dipolar_moment_2) - &
             dot_product(dipolar_moment_1, dipolar_moment_1) + &
             2._DP*dot_product(dipolar_moment_2 - dipolar_moment_1, &
                this%total_moment%get() - dipolar_moment_1))
    end function Spherical_visit_transmutation

    !> \[
    !>      \Delta U = \frac{1}{6\epsilon V} [
    !>          (\vec{\mu} \cdot \vec{\mu}) \pm 2(\vec{\mu} \cdot \vec{M})
    !>      ]
    !> \]
    pure real(DP) function Spherical_visit_exchange(this, i_component, dipolar_moment, signed) &
        result(delta_energy)
        class(Spherical_DES_Surf_Mixture), intent(in) :: this
        integer, intent(in) :: i_component
        real(DP), intent(in) :: dipolar_moment(:), signed

        delta_energy = 0._DP
        if (.not.this%total_moment%is_dipolar(i_component)) return

        delta_energy = 1._DP/6._DP/this%permittivity / product(this%periodic_box%get_size()) * &
            (dot_product(dipolar_moment, dipolar_moment) + &
             signed*2._DP*dot_product(dipolar_moment, this%total_moment%get()))
    end function Spherical_visit_exchange

!end implementation Spherical_DES_Surf_Mixture

!implementation Rectangular_DES_Surf_Mixture

    !> \[ U = \frac{1}{2\epsilon V} M_z^2  \]
    pure real(DP) function Rectangular_visit(this) result(energy)
        class(Rectangular_DES_Surf_Mixture), intent(in) :: this

        real(DP) :: total_moment(num_dimensions)

        total_moment = this%total_moment%get()
        energy = 1._DP/2._DP/this%permittivity / product(this%periodic_box%get_size()) * &
            total_moment(3)**2
    end function Rectangular_visit

    !> \[
    !>      \Delta U =
    !>          \frac{1}{2\epsilon V} [
    !>              \mu^{\prime 2}_z - \mu_z^2 + 2(\mu^\prime_z - \mu_z)(M_z - \mu_z)
    !>          ]
    !> \]
    pure real(DP) function Rectangular_visit_transmutation(this, ij_components, dipolar_moment_2, &
        dipolar_moment_1) result(delta_energy)
        class(Rectangular_DES_Surf_Mixture), intent(in) :: this
        integer, intent(in) :: ij_components(:)
        real(DP), intent(in) :: dipolar_moment_2(:), dipolar_moment_1(:)

        real(DP) :: total_moment(num_dimensions)

        delta_energy = 0._DP
        if (.not.(this%total_moment%is_dipolar(ij_components(1)) .or. &
            this%total_moment%is_dipolar(ij_components(2)))) return !shortcut?

        total_moment = this%total_moment%get()
        delta_energy = 1._DP/2._DP/this%permittivity / product(this%periodic_box%get_size()) * &
            (dipolar_moment_2(3)**2 - dipolar_moment_1(3)**2 + &
             2._DP*(dipolar_moment_2(3) - dipolar_moment_1(3)) * &
             (total_moment(3) - dipolar_moment_1(3)))
    end function Rectangular_visit_transmutation

    !> \[
    !>      \Delta U = \frac{1}{2\epsilon V}[\mu_z^2 \pm 2 \mu_z M_z]
    !> \]
    pure real(DP) function Rectangular_visit_exchange(this, i_component, dipolar_moment, signed) &
        result(delta_energy)
        class(Rectangular_DES_Surf_Mixture), intent(in) :: this
        integer, intent(in) :: i_component
        real(DP), intent(in) :: dipolar_moment(:), signed

        real(DP) :: total_moment(num_dimensions)

        delta_energy = 0._DP
        if (.not.this%total_moment%is_dipolar(i_component)) return

        total_moment = this%total_moment%get()
        delta_energy = 1._DP/2._DP/this%permittivity / product(this%periodic_box%get_size()) * &
            (dipolar_moment(3)**2 + &
             signed*2._DP*dipolar_moment(3) * total_moment(3))
    end function Rectangular_visit_exchange

!end implementation Rectangular_DES_Surf_Mixture

!implementation Null_DES_Surf_Mixture

    pure real(DP) function Null_visit(this) result(energy)
        class(Null_DES_Surf_Mixture), intent(in) :: this
        energy = 0._DP
    end function Null_visit

    pure real(DP) function Null_visit_transmutation(this, ij_components, dipolar_moment_2, &
        dipolar_moment_1) result(delta_energy)
        class(Null_DES_Surf_Mixture), intent(in) :: this
        integer, intent(in) :: ij_components(:)
        real(DP), intent(in) :: dipolar_moment_2(:), dipolar_moment_1(:)
        delta_energy = 0._DP
    end function Null_visit_transmutation

    pure real(DP) function Null_visit_exchange(this, i_component, dipolar_moment, signed) &
        result(delta_energy)
        class(Null_DES_Surf_Mixture), intent(in) :: this
        integer, intent(in) :: i_component
        real(DP), intent(in) :: dipolar_moment(:), signed
        delta_energy = 0._DP
    end function Null_visit_exchange

!end implementation Null_DES_Surf_Mixture

end module classes_des_surf_mixture
