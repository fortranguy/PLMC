module classes_dipolar_interactions_facade

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use types_component_wrapper, only: Component_Wrapper
use types_dipolar_interactions_dynamic_wrapper, only: Dipolar_Interactions_Dynamic_Wrapper
use types_dipolar_interactions_static_wrapper, only: Dipolar_Interactions_Static_Wrapper
use types_reals_line, only: Reals_Line
use procedures_dipolar_interactions_resetter, only: dipolar_interactions_reset => reset
use procedures_dipolar_interactions_visitor, only: dipolar_interactions_visit => visit

implicit none

private

    type, abstract, public :: Abstract_Dipolar_Interactions_Facade
    private
        type(Dipolar_Interactions_Static_Wrapper), pointer :: dipolar_interactions_static => null()
        logical :: real_pair_must_be_reset = .false.
    contains
        procedure(Abstract_destroy), deferred :: destroy
        procedure :: save => Abstract_save
        procedure :: restore => Abstract_restore
        procedure(Abstract_reset), deferred :: reset
        procedure(Abstract_visit), deferred :: visit
        procedure(Abstract_clone), deferred, private :: clone
    end type Abstract_Dipolar_Interactions_Facade

    abstract interface

        subroutine Abstract_destroy(this)
        import :: Abstract_Dipolar_Interactions_Facade
            class(Abstract_Dipolar_Interactions_Facade), intent(inout) :: this
        end subroutine Abstract_destroy

        subroutine Abstract_clone(this, dipolar_interactions_static_target, &
            dipolar_interactions_static_source)
        import :: Dipolar_Interactions_Static_Wrapper, Abstract_Dipolar_Interactions_Facade
            class(Abstract_Dipolar_Interactions_Facade), intent(in) :: this
            type(Dipolar_Interactions_Static_Wrapper), intent(inout) :: &
                dipolar_interactions_static_target
            type(Dipolar_Interactions_Static_Wrapper), intent(in) :: &
                dipolar_interactions_static_source
        end subroutine Abstract_clone

        subroutine Abstract_reset(this)
        import :: Abstract_Dipolar_Interactions_Facade
            class(Abstract_Dipolar_Interactions_Facade), intent(in) :: this
        end subroutine Abstract_reset

        subroutine Abstract_visit(this, new_energies, new_shared_energy, box_volume_ratio, &
            energies, shared_energy)
        import :: DP, Reals_Line, Abstract_Dipolar_Interactions_Facade
            class(Abstract_Dipolar_Interactions_Facade), intent(in) :: this
            type(Reals_Line), intent(inout) :: new_energies(:)
            real(DP), intent(out) :: new_shared_energy
            real(DP), intent(in) :: box_volume_ratio !! \( \frac{V}{V_0} \)
            type(Reals_Line), intent(in) :: energies(:)
            real(DP), intent(in) :: shared_energy
        end subroutine Abstract_visit

    end interface

    type, extends(Abstract_Dipolar_Interactions_Facade), public :: &
        Scalable_Dipolar_Interactions_Facade
    contains
        procedure :: construct => Scalable_construct
        procedure :: destroy => Scalable_destroy
        procedure :: reset => Scalable_reset
        procedure :: visit => Scalable_visit
        procedure, private :: clone => Scalable_clone
    end type Scalable_Dipolar_Interactions_Facade

    type, extends(Abstract_Dipolar_Interactions_Facade), public :: &
        Unscalable_Dipolar_Interactions_Facade
    private
        type(Component_Wrapper), pointer :: components(:) => null()
        type(Dipolar_Interactions_Dynamic_Wrapper), pointer :: dipolar_interactions_dynamic => &
            null()
    contains
        procedure :: construct => Unscalable_construct
        procedure :: destroy => Unscalable_destroy
        procedure :: reset => Unscalable_reset
        procedure :: visit => Unscalable_visit
        procedure, private :: clone => Unscalable_clone
    end type Unscalable_Dipolar_Interactions_Facade

    type, extends(Abstract_Dipolar_Interactions_Facade), public :: Null_Dipolar_Interactions_Facade
    contains
        procedure :: destroy => Null_destroy
        procedure :: save => Null_save
        procedure :: reset => Null_reset
        procedure :: visit => Null_visit
        procedure, private :: clone => Null_clone
    end type Null_Dipolar_Interactions_Facade

contains

!implementation Abstract_Dipolar_Interactions_Facade

    subroutine Abstract_save(this, dipolar_interactions_static, new_box_volume)
        class(Abstract_Dipolar_Interactions_Facade), intent(inout) :: this
        type(Dipolar_Interactions_Static_Wrapper), intent(inout) :: dipolar_interactions_static
        real(DP), intent(in) :: new_box_volume

        this%real_pair_must_be_reset = new_box_volume > this%dipolar_interactions_static%&
            box_volume_memento_real%get()
        call this%clone(dipolar_interactions_static, this%dipolar_interactions_static)
    end subroutine Abstract_save

    subroutine Abstract_restore(this, dipolar_interactions_static)
        class(Abstract_Dipolar_Interactions_Facade), intent(in) :: this
        type(Dipolar_Interactions_Static_Wrapper), intent(in) :: dipolar_interactions_static

        call this%clone(this%dipolar_interactions_static, dipolar_interactions_static)
    end subroutine Abstract_restore

!end implementation Abstract_Dipolar_Interactions_Facade

!implementation Scalable_Dipolar_Interactions_Facade

    subroutine Scalable_construct(this, dipolar_interactions_static)
        class(Scalable_Dipolar_Interactions_Facade), intent(out) :: this
        type(Dipolar_Interactions_Static_Wrapper), target, intent(in) :: dipolar_interactions_static

        this%dipolar_interactions_static => dipolar_interactions_static
    end subroutine Scalable_construct

    subroutine Scalable_destroy(this)
        class(Scalable_Dipolar_Interactions_Facade), intent(inout) :: this

        this%dipolar_interactions_static => null()
    end subroutine Scalable_destroy

    subroutine Scalable_clone(this, dipolar_interactions_static_target, &
        dipolar_interactions_static_source)
        class(Scalable_Dipolar_Interactions_Facade), intent(in) :: this
        type(Dipolar_Interactions_Static_Wrapper), intent(inout) :: &
            dipolar_interactions_static_target
        type(Dipolar_Interactions_Static_Wrapper), intent(in) :: dipolar_interactions_static_source

        if (this%real_pair_must_be_reset) then
            allocate(dipolar_interactions_static_target%real_pair, &
                source=dipolar_interactions_static_source%real_pair)
        end if
    end subroutine Scalable_clone

    subroutine Scalable_reset(this)
        class(Scalable_Dipolar_Interactions_Facade), intent(in) :: this

        if (this%real_pair_must_be_reset) then
            call this%dipolar_interactions_static%box_volume_memento_real%save()
            call this%dipolar_interactions_static%real_pair%reset()
        end if
    end subroutine Scalable_reset

    subroutine Scalable_visit(this, new_energies, new_shared_energy, box_volume_ratio, energies, &
        shared_energy)
        class(Scalable_Dipolar_Interactions_Facade), intent(in) :: this
        type(Reals_Line), intent(inout) :: new_energies(:)
        real(DP), intent(out) :: new_shared_energy
        real(DP), intent(in) :: box_volume_ratio
        type(Reals_Line), intent(in) :: energies(:)
        real(DP), intent(in) :: shared_energy

        integer :: i_component

        do i_component = 1, size(new_energies)
            new_energies(i_component)%line = energies(i_component)%line / box_volume_ratio
        end do
        new_shared_energy = shared_energy / box_volume_ratio
    end subroutine Scalable_visit

!implementation Scalable_Dipolar_Interactions_Facade

!implementation Unscalable_Dipolar_Interactions_Facade

    subroutine Unscalable_construct(this, components, dipolar_interactions_dynamic, &
        dipolar_interactions_static)
        class(Unscalable_Dipolar_Interactions_Facade), intent(out) :: this
        type(Component_Wrapper), target, intent(in) :: components(:)
        type(Dipolar_Interactions_Dynamic_Wrapper), target, intent(in) :: &
            dipolar_interactions_dynamic
        type(Dipolar_Interactions_Static_Wrapper), target, intent(in) :: dipolar_interactions_static

        this%components => components
        this%dipolar_interactions_dynamic => dipolar_interactions_dynamic
        this%dipolar_interactions_static => dipolar_interactions_static
    end subroutine Unscalable_construct

    subroutine Unscalable_destroy(this)
        class(Unscalable_Dipolar_Interactions_Facade), intent(inout) :: this

        this%dipolar_interactions_static => null()
        this%dipolar_interactions_dynamic => null()
        this%components => null()
    end subroutine Unscalable_destroy

    subroutine Unscalable_clone(this, dipolar_interactions_static_target, &
        dipolar_interactions_static_source)
        class(Unscalable_Dipolar_Interactions_Facade), intent(in) :: this
        type(Dipolar_Interactions_Static_Wrapper), intent(inout) :: &
            dipolar_interactions_static_target
        type(Dipolar_Interactions_Static_Wrapper), intent(in) :: dipolar_interactions_static_source

        allocate(dipolar_interactions_static_target%box_volume_memento_real, &
            source=dipolar_interactions_static_source%box_volume_memento_real)
        if (this%real_pair_must_be_reset) then
            allocate(dipolar_interactions_static_target%real_pair, &
                source=dipolar_interactions_static_source%real_pair)
        end if
        allocate(dipolar_interactions_static_target%box_volume_memento_reci, &
            source=dipolar_interactions_static_source%box_volume_memento_reci)
        allocate(dipolar_interactions_static_target%reci_weight, &
            source=dipolar_interactions_static_source%reci_weight)
        allocate(dipolar_interactions_static_target%reci_structure, &
            source=dipolar_interactions_static_source%reci_structure)
        allocate(dipolar_interactions_static_target%dlc_weight, &
            source=dipolar_interactions_static_source%dlc_weight)
        allocate(dipolar_interactions_static_target%dlc_structures, &
            source=dipolar_interactions_static_source%dlc_structures)
    end subroutine Unscalable_clone

    subroutine Unscalable_reset(this)
        class(Unscalable_Dipolar_Interactions_Facade), intent(in) :: this

        call dipolar_interactions_reset(this%dipolar_interactions_static, this%&
            real_pair_must_be_reset)
    end subroutine Unscalable_reset

    subroutine Unscalable_visit(this, new_energies, new_shared_energy, box_volume_ratio, energies, &
        shared_energy)
        class(Unscalable_Dipolar_Interactions_Facade), intent(in) :: this
        type(Reals_Line), intent(inout) :: new_energies(:)
        real(DP), intent(out) :: new_shared_energy
        real(DP), intent(in) :: box_volume_ratio
        type(Reals_Line), intent(in) :: energies(:)
        real(DP), intent(in) :: shared_energy

        call dipolar_interactions_visit(new_energies, new_shared_energy, this%components, this%&
            dipolar_interactions_dynamic)
    end subroutine Unscalable_visit

!end implementation Unscalable_Dipolar_Interactions_Facade

!implementation Null_Dipolar_Interactions_Facade

    subroutine Null_destroy(this)
        class(Null_Dipolar_Interactions_Facade), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_save(this, dipolar_interactions_static, new_box_volume)
        class(Null_Dipolar_Interactions_Facade), intent(inout) :: this
        type(Dipolar_Interactions_Static_Wrapper), intent(inout) :: dipolar_interactions_static
        real(DP), intent(in) :: new_box_volume
    end subroutine Null_save

    subroutine Null_clone(this, dipolar_interactions_static_target, &
        dipolar_interactions_static_source)
        class(Null_Dipolar_Interactions_Facade), intent(in) :: this
        type(Dipolar_Interactions_Static_Wrapper), intent(inout) :: &
            dipolar_interactions_static_target
        type(Dipolar_Interactions_Static_Wrapper), intent(in) :: dipolar_interactions_static_source
    end subroutine Null_clone

    subroutine Null_reset(this)
        class(Null_Dipolar_Interactions_Facade), intent(in) :: this
    end subroutine Null_reset

    subroutine Null_visit(this, new_energies, new_shared_energy, box_volume_ratio, energies, &
        shared_energy)
        class(Null_Dipolar_Interactions_Facade), intent(in) :: this
        type(Reals_Line), intent(inout) :: new_energies(:)
        real(DP), intent(out) :: new_shared_energy
        real(DP), intent(in) :: box_volume_ratio
        type(Reals_Line), intent(in) :: energies(:)
        real(DP), intent(in) :: shared_energy

        integer :: i_component

        do i_component = 1, size(new_energies)
            new_energies(i_component)%line = 0._DP
        end do
        new_shared_energy = 0._DP
    end subroutine Null_visit

!end implementation Null_Dipolar_Interactions_Facade

end module classes_dipolar_interactions_facade
