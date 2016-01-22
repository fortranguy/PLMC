module class_ewald_surf_mixture

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use class_periodic_box, only: Abstract_Periodic_Box
use class_permittivity, only: Abstract_Permittivity
use class_mixture_total_moment, only: Abstract_Mixture_Total_Moment
use types_temporary_particle, only: Concrete_Temporary_Particle

implicit none

private

    type, abstract, public :: Abstract_Ewald_Surf_Mixture
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
        real(DP) :: permittivity
        class(Abstract_Mixture_Total_Moment), pointer :: mixture_total_moment => null()
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure(Abstract_visit), deferred :: visit
        procedure(Abstract_visit_rotation), deferred :: visit_rotation
        procedure :: visit_add => Abstract_visit_add
        procedure :: visit_remove => Abstract_visit_remove
        procedure(Abstract_visit_exchange), deferred, private :: visit_exchange
    end type Abstract_Ewald_Surf_Mixture

    abstract interface

        pure real(DP) function Abstract_visit(this)
        import :: DP, Abstract_Ewald_Surf_Mixture
            class(Abstract_Ewald_Surf_Mixture), intent(in) :: this
        end function Abstract_visit

        pure real(DP) function Abstract_visit_rotation(this, i_component, new, old)
        import :: DP, Abstract_Ewald_Surf_Mixture, Concrete_Temporary_Particle
            class(Abstract_Ewald_Surf_Mixture), intent(in) :: this
            integer, intent(in) :: i_component
            type(Concrete_Temporary_Particle), intent(in) :: new, old
        end function Abstract_visit_rotation

        pure real(DP) function Abstract_visit_exchange(this, i_component, particle, signed)
        import :: DP, Abstract_Ewald_Surf_Mixture, Concrete_Temporary_Particle
            class(Abstract_Ewald_Surf_Mixture), intent(in) :: this
            integer, intent(in) :: i_component
            type(Concrete_Temporary_Particle), intent(in) :: particle
            real(DP), intent(in) :: signed
        end function Abstract_visit_exchange

    end interface

    type, extends(Abstract_Ewald_Surf_Mixture), public :: Spheric_Ewald_Surf_Mixture
    contains
        procedure :: visit => Spheric_visit
        procedure :: visit_rotation => Spheric_visit_rotation
        procedure, private :: visit_exchange => Spheric_visit_exchange
    end type Spheric_Ewald_Surf_Mixture

    type, extends(Abstract_Ewald_Surf_Mixture), public :: Rectangular_Ewald_Surf_Mixture
    contains
        procedure :: visit => Rectangular_visit
        procedure :: visit_rotation => Rectangular_visit_rotation
        procedure, private :: visit_exchange => Rectangular_visit_exchange
    end type Rectangular_Ewald_Surf_Mixture

    type, extends(Abstract_Ewald_Surf_Mixture), public :: Null_Ewald_Surf_Mixture
    contains
        procedure :: visit => Null_visit
        procedure :: visit_rotation => Null_visit_rotation
        procedure, private :: visit_exchange => Null_visit_exchange
    end type Null_Ewald_Surf_Mixture

contains

!implementation Abstract_Ewald_Surf_Mixture

    subroutine Abstract_construct(this, periodic_box, permittivity, mixture_total_moment)
        class(Abstract_Ewald_Surf_Mixture), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Permittivity), intent(in) :: permittivity
        class(Abstract_Mixture_Total_Moment), target, intent(in) :: mixture_total_moment

        this%periodic_box => periodic_box
        this%permittivity = permittivity%get()
        this%mixture_total_moment => mixture_total_moment
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Ewald_Surf_Mixture), intent(inout) :: this

        this%mixture_total_moment => null()
        this%periodic_box => null()
    end subroutine Abstract_destroy

    pure real(DP) function Abstract_visit_add(this, i_component, particle) result(energy)
        class(Abstract_Ewald_Surf_Mixture), intent(in) :: this
        integer, intent(in) :: i_component
        type(Concrete_Temporary_Particle), intent(in) :: particle

        energy = this%visit_exchange(i_component, particle, 1.0_DP)
    end function Abstract_visit_add

    pure real(DP) function Abstract_visit_remove(this, i_component, particle) result(energy)
        class(Abstract_Ewald_Surf_Mixture), intent(in) :: this
        integer, intent(in) :: i_component
        type(Concrete_Temporary_Particle), intent(in) :: particle

        energy = this%visit_exchange(i_component, particle, -1.0_DP)
    end function Abstract_visit_remove

!end implementation Abstract_Ewald_Surf_Mixture

!implementation Spheric_Ewald_Surf_Mixture

    !> \[ U = \frac{1}{6\epsilon V} \vec{M}^2 \]
    pure real(DP) function Spheric_visit(this) result(energy)
        class(Spheric_Ewald_Surf_Mixture), intent(in) :: this

        energy = 1._DP/6._DP/this%permittivity / product(this%periodic_box%get_size()) * &
            dot_product(this%mixture_total_moment%get(), this%mixture_total_moment%get())
    end function Spheric_visit

    !> \[
    !>      \Delta U = \frac{1}{6\epsilon V} [
    !>                      (\vec{\mu}^\prime_\mathsf{i} \cdot \vec{\mu}^\prime_\mathsf{i}) -
    !>                      (\vec{\mu}_\mathsf{i} \cdot \vec{\mu}_\mathsf{i}) +
    !>                      2 (\vec{\mu}^\prime_\mathsf{i} - \vec{\mu}_\mathsf{i}) \cdot
    !>                          (\vec{M} - \vec{\mu}_\mathsf{i})
    !>                 ]
    !> \]
    pure real(DP) function Spheric_visit_rotation(this, i_component, new, old) result(delta_energy)
        class(Spheric_Ewald_Surf_Mixture), intent(in) :: this
        integer, intent(in) :: i_component
        type(Concrete_Temporary_Particle), intent(in) :: new, old

        delta_energy = 0._DP
        if (.not.this%mixture_total_moment%is_dipolar(i_component)) return

        delta_energy = 1._DP/6._DP/this%permittivity / product(this%periodic_box%get_size()) * &
            (dot_product(new%dipolar_moment, new%dipolar_moment) - &
             dot_product(old%dipolar_moment, old%dipolar_moment) + &
             2._DP*dot_product(new%dipolar_moment - old%dipolar_moment, &
                this%mixture_total_moment%get() - old%dipolar_moment))
    end function Spheric_visit_rotation

    !> \[
    !>      \Delta U = \frac{1}{6\epsilon V} [
    !>          (\vec{\mu} \cdot \vec{\mu}) \pm 2(\vec{\mu} \cdot \vec{M})
    !>      ]
    !> \]
    pure real(DP) function Spheric_visit_exchange(this, i_component, particle, signed) &
        result(delta_energy)
        class(Spheric_Ewald_Surf_Mixture), intent(in) :: this
        integer, intent(in) :: i_component
        type(Concrete_Temporary_Particle), intent(in) :: particle
        real(DP), intent(in) :: signed

        delta_energy = 0._DP
        if (.not.this%mixture_total_moment%is_dipolar(i_component)) return

        delta_energy = 1._DP/6._DP/this%permittivity / product(this%periodic_box%get_size()) * &
            dot_product(particle%dipolar_moment, particle%dipolar_moment) + &
            signed*2._DP*dot_product(particle%dipolar_moment, this%mixture_total_moment%get())
    end function Spheric_visit_exchange

!end implementation Spheric_Ewald_Surf_Mixture

!implementation Rectangular_Ewald_Surf_Mixture

    !> \[ U = \frac{1}{2\epsilon V} M_z^2  \]
    pure real(DP) function Rectangular_visit(this) result(energy)
        class(Rectangular_Ewald_Surf_Mixture), intent(in) :: this

        real(DP) :: total_moment(num_dimensions)

        total_moment = this%mixture_total_moment%get()
        energy = 1._DP/2._DP/this%permittivity / product(this%periodic_box%get_size()) * &
            total_moment(3)**2
    end function Rectangular_visit

    !> \[
    !>      \Delta U =
    !>          \frac{1}{2\epsilon V} [
    !>              \mu^{\prime 2}_z - \mu_z^2 + 2(\mu^\prime_z - \mu_z)(M_z - \mu_z)
    !>          ]
    !> \]
    pure real(DP) function Rectangular_visit_rotation(this, i_component, new, old) &
        result(delta_energy)
        class(Rectangular_Ewald_Surf_Mixture), intent(in) :: this
        integer, intent(in) :: i_component
        type(Concrete_Temporary_Particle), intent(in) :: new, old

        real(DP) :: total_moment(num_dimensions)

        delta_energy = 0._DP
        if (.not.this%mixture_total_moment%is_dipolar(i_component)) return

        total_moment = this%mixture_total_moment%get()
        delta_energy = 1._DP/2._DP/this%permittivity / product(this%periodic_box%get_size()) * &
            (new%dipolar_moment(3)**2 - old%dipolar_moment(3)**2 + &
             2._DP*(new%dipolar_moment(3) - old%dipolar_moment(3)) * &
             (total_moment(3) - old%dipolar_moment(3)))
    end function Rectangular_visit_rotation

    !> \[
    !>      \Delta U = \frac{1}{2\epsilon V}[\mu_z^2 \pm 2 \mu_z M_z]
    !> \]
    pure real(DP) function Rectangular_visit_exchange(this, i_component, particle, signed) &
        result(delta_energy)
        class(Rectangular_Ewald_Surf_Mixture), intent(in) :: this
        integer, intent(in) :: i_component
        type(Concrete_Temporary_Particle), intent(in) :: particle
        real(DP), intent(in) :: signed

        real(DP) :: total_moment(num_dimensions)

        delta_energy = 0._DP
        if (.not.this%mixture_total_moment%is_dipolar(i_component)) return

        total_moment = this%mixture_total_moment%get()
        delta_energy = 1._DP/2._DP/this%permittivity / product(this%periodic_box%get_size()) * &
            (particle%dipolar_moment(3)**2 + &
             signed*2._DP*particle%dipolar_moment(3) * total_moment(3))
    end function Rectangular_visit_exchange

!end implementation Rectangular_Ewald_Surf_Mixture

!implementation Null_Ewald_Surf_Mixture

    pure real(DP) function Null_visit(this) result(energy)
        class(Null_Ewald_Surf_Mixture), intent(in) :: this
        energy = 0._DP
    end function Null_visit

    pure real(DP) function Null_visit_rotation(this, i_component, new, old) result(delta_energy)
        class(Null_Ewald_Surf_Mixture), intent(in) :: this
        integer, intent(in) :: i_component
        type(Concrete_Temporary_Particle), intent(in) :: new, old
        delta_energy = 0._DP
    end function Null_visit_rotation

    pure real(DP) function Null_visit_exchange(this, i_component, particle, signed) &
        result(delta_energy)
        class(Null_Ewald_Surf_Mixture), intent(in) :: this
        integer, intent(in) :: i_component
        type(Concrete_Temporary_Particle), intent(in) :: particle
        real(DP), intent(in) :: signed
        delta_energy = 0._DP
    end function Null_visit_exchange

!end implementation Null_Ewald_Surf_Mixture

end module class_ewald_surf_mixture
