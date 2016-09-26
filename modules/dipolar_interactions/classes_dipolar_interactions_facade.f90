module classes_dipolar_interactions_facade

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_box_size_memento, only: Abstract_Box_Size_Memento
use procedures_box_factory, box_destroy => destroy
use classes_des_real_pair, only: Abstract_DES_Real_Pair
use procedures_des_reci_factory, only: des_reci_destroy => destroy
use procedures_dlc_factory, only: dlc_destroy => destroy
use types_component_wrapper, only: Component_Wrapper
use procedures_des_real_factory, only: des_real_destroy => destroy
use types_dipolar_interactions_dynamic_wrapper, only: Dipolar_Interactions_Dynamic_Wrapper
use types_dipolar_interactions_static_wrapper, only: Dipolar_Interactions_Static_Wrapper
use types_reals_line, only: Reals_Line
use procedures_dipolar_interactions_resetter, only: dipolar_interactions_reset => reset, &
    dipolar_interactions_reset_real => reset_real
use procedures_dipolar_interactions_visitor, only: dipolar_interactions_visit => visit

implicit none

private

    type, abstract, public :: Abstract_Dipolar_Interactions_Facade
    private
        type(Component_Wrapper), pointer :: components(:) => null()
        type(Dipolar_Interactions_Dynamic_Wrapper), pointer :: dipolar_interactions_dynamic => &
                null()
        type(Dipolar_Interactions_Static_Wrapper), pointer :: dipolar_interactions_static => null()
        logical :: crop_real_pair = .false.
        logical :: reset_real_pair = .false.
    contains
        procedure(Abstract_destroy), deferred :: destroy
        procedure :: save => Abstract_save
        procedure :: restore => Abstract_restore
        procedure(Abstract_reset), deferred :: reset
        procedure(Abstract_visit), deferred :: visit
        procedure(Abstract_set_real_pair_flags), private, deferred :: set_real_pair_flags
        procedure(Abstract_clone), deferred, private :: clone
        procedure, private :: clone_real => Abstract_clone_real
        procedure(Abstract_target), deferred, private :: target
        procedure, private :: target_real => Abstract_target_real
    end type Abstract_Dipolar_Interactions_Facade

    abstract interface

        subroutine Abstract_destroy(this)
        import :: Abstract_Dipolar_Interactions_Facade
            class(Abstract_Dipolar_Interactions_Facade), intent(inout) :: this
        end subroutine Abstract_destroy

        subroutine Abstract_set_real_pair_flags(this, new_box_size)
        import :: DP, Abstract_Dipolar_Interactions_Facade
            class(Abstract_Dipolar_Interactions_Facade), intent(inout) :: this
            real(DP), intent(in) :: new_box_size(:)
        end subroutine Abstract_set_real_pair_flags

        subroutine Abstract_clone(this, dipolar_interactions_static_target, &
            dipolar_interactions_static_source)
        import :: Dipolar_Interactions_Static_Wrapper, Abstract_Dipolar_Interactions_Facade
            class(Abstract_Dipolar_Interactions_Facade), intent(in) :: this
            type(Dipolar_Interactions_Static_Wrapper), intent(inout) :: &
                dipolar_interactions_static_target
            type(Dipolar_Interactions_Static_Wrapper), intent(in) :: &
                dipolar_interactions_static_source
        end subroutine Abstract_clone

        subroutine Abstract_target(this)
        import :: Abstract_Dipolar_Interactions_Facade
            class(Abstract_Dipolar_Interactions_Facade), intent(in) :: this
        end subroutine Abstract_target

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
        procedure, private :: set_real_pair_flags => Scalable_set_real_pair_flags
        procedure, private :: clone => Scalable_clone
        procedure, private :: target => Scalable_target
    end type Scalable_Dipolar_Interactions_Facade

    type, extends(Abstract_Dipolar_Interactions_Facade), public :: &
        Unscalable_Dipolar_Interactions_Facade
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
    contains
        procedure :: construct => Unscalable_construct
        procedure :: destroy => Unscalable_destroy
        procedure :: reset => Unscalable_reset
        procedure :: visit => Unscalable_visit
        procedure, private :: set_real_pair_flags => Unscalable_set_real_pair_flags
        procedure, private :: clone => Unscalable_clone
        procedure, private :: target => Unscalable_target
    end type Unscalable_Dipolar_Interactions_Facade

    type, extends(Abstract_Dipolar_Interactions_Facade), public :: Null_Dipolar_Interactions_Facade
    contains
        procedure :: destroy => Null_destroy
        procedure :: save => Null_save
        procedure :: reset => Null_reset
        procedure :: visit => Null_visit
        procedure, private :: set_real_pair_flags => Null_set_real_pair_flags
        procedure, private :: clone => Null_clone
        procedure, private :: target => Null_target
    end type Null_Dipolar_Interactions_Facade

contains

!implementation Abstract_Dipolar_Interactions_Facade

    subroutine Abstract_save(this, dipolar_interactions_static, new_box_size)
        class(Abstract_Dipolar_Interactions_Facade), intent(inout) :: this
        type(Dipolar_Interactions_Static_Wrapper), intent(inout) :: dipolar_interactions_static
        real(DP), intent(in) :: new_box_size(:)

        call this%set_real_pair_flags(new_box_size)
        call this%clone(dipolar_interactions_static, this%dipolar_interactions_static)
    end subroutine Abstract_save

    subroutine Abstract_restore(this, dipolar_interactions_static)
        class(Abstract_Dipolar_Interactions_Facade), intent(in) :: this
        type(Dipolar_Interactions_Static_Wrapper), intent(in) :: dipolar_interactions_static

        call this%clone(this%dipolar_interactions_static, dipolar_interactions_static)
        call this%target()
    end subroutine Abstract_restore

    subroutine Abstract_clone_real(this, box_size_memento_target, des_real_pair_target, &
        box_size_memento_source, des_real_pair_source)
        class(Abstract_Dipolar_Interactions_Facade), intent(in) :: this
        class(Abstract_Box_Size_Memento), allocatable, intent(inout) :: box_size_memento_target
        class(Abstract_DES_Real_Pair), allocatable, intent(inout) :: des_real_pair_target
        class(Abstract_Box_Size_Memento), intent(in) :: box_size_memento_source
        class(Abstract_DES_Real_Pair), intent(in) :: des_real_pair_source

        if (this%reset_real_pair) then
            call des_real_destroy(des_real_pair_target)
            call box_destroy(box_size_memento_target)
            allocate(box_size_memento_target, source=box_size_memento_source)
            allocate(des_real_pair_target, source=des_real_pair_source)
        end if
    end subroutine Abstract_clone_real

    subroutine Abstract_target_real(this)
        class(Abstract_Dipolar_Interactions_Facade), intent(in) :: this

        integer :: i_component, j_component

        if (this%reset_real_pair) then
            call this%dipolar_interactions_static%real_pair%target(this%&
                dipolar_interactions_static%box_size_memento_real)
            do j_component = 1, size(this%dipolar_interactions_dynamic%real_components, 2)
                do i_component = 1, size(this%dipolar_interactions_dynamic%real_components, 1)
                    call this%dipolar_interactions_dynamic%&
                        real_components(i_component, j_component)%component%&
                        target(this%dipolar_interactions_static%box_size_memento_real, this%&
                            dipolar_interactions_static%real_pair)
                end do
            end do
        end if
    end subroutine Abstract_target_real

!end implementation Abstract_Dipolar_Interactions_Facade

!implementation Scalable_Dipolar_Interactions_Facade

    subroutine Scalable_construct(this, components, dipolar_interactions_dynamic, &
        dipolar_interactions_static)
        class(Scalable_Dipolar_Interactions_Facade), intent(out) :: this
        type(Component_Wrapper), target, intent(in) :: components(:)
        type(Dipolar_Interactions_Dynamic_Wrapper), target, intent(in) :: &
            dipolar_interactions_dynamic
        type(Dipolar_Interactions_Static_Wrapper), target, intent(in) :: dipolar_interactions_static

        this%components => components
        this%dipolar_interactions_dynamic => dipolar_interactions_dynamic
        this%dipolar_interactions_static => dipolar_interactions_static
    end subroutine Scalable_construct

    subroutine Scalable_destroy(this)
        class(Scalable_Dipolar_Interactions_Facade), intent(inout) :: this

        this%dipolar_interactions_static => null()
        this%dipolar_interactions_dynamic => null()
        this%components => null()
    end subroutine Scalable_destroy

    subroutine Scalable_set_real_pair_flags(this, new_box_size)
        class(Scalable_Dipolar_Interactions_Facade), intent(inout) :: this
        real(DP), intent(in) :: new_box_size(:)

        this%crop_real_pair = .false.
        this%reset_real_pair = product(new_box_size) > product(this%dipolar_interactions_static%&
            box_size_memento_real%get())
    end subroutine Scalable_set_real_pair_flags

    subroutine Scalable_target(this)
        class(Scalable_Dipolar_Interactions_Facade), intent(in) :: this

        call this%target_real()
    end subroutine Scalable_target

    subroutine Scalable_clone(this, dipolar_interactions_static_target, &
        dipolar_interactions_static_source)
        class(Scalable_Dipolar_Interactions_Facade), intent(in) :: this
        type(Dipolar_Interactions_Static_Wrapper), intent(inout) :: &
            dipolar_interactions_static_target
        type(Dipolar_Interactions_Static_Wrapper), intent(in) :: dipolar_interactions_static_source

        call this%clone_real(dipolar_interactions_static_target%box_size_memento_real, &
            dipolar_interactions_static_target%real_pair, dipolar_interactions_static_source%&
            box_size_memento_real, dipolar_interactions_static_source%real_pair)
    end subroutine Scalable_clone

    subroutine Scalable_reset(this)
        class(Scalable_Dipolar_Interactions_Facade), intent(in) :: this

        call dipolar_interactions_reset_real(this%dipolar_interactions_static%&
            box_size_memento_real, this%dipolar_interactions_static%real_pair, this%reset_real_pair)
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

    subroutine Unscalable_construct(this, periodic_box, components, dipolar_interactions_dynamic, &
        dipolar_interactions_static)
        class(Unscalable_Dipolar_Interactions_Facade), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        type(Component_Wrapper), target, intent(in) :: components(:)
        type(Dipolar_Interactions_Dynamic_Wrapper), target, intent(in) :: &
            dipolar_interactions_dynamic
        type(Dipolar_Interactions_Static_Wrapper), target, intent(in) :: dipolar_interactions_static

        this%periodic_box => periodic_box
        this%components => components
        this%dipolar_interactions_dynamic => dipolar_interactions_dynamic
        this%dipolar_interactions_static => dipolar_interactions_static
    end subroutine Unscalable_construct

    subroutine Unscalable_destroy(this)
        class(Unscalable_Dipolar_Interactions_Facade), intent(inout) :: this

        this%dipolar_interactions_static => null()
        this%dipolar_interactions_dynamic => null()
        this%components => null()
        this%periodic_box => null()
    end subroutine Unscalable_destroy

    subroutine Unscalable_set_real_pair_flags(this, new_box_size)
        class(Unscalable_Dipolar_Interactions_Facade), intent(inout) :: this
        real(DP), intent(in) :: new_box_size(:)

        if (new_box_size(1) < this%dipolar_interactions_static%box_size_memento_real%get_edge()) &
            then
            this%crop_real_pair = .true.
            this%reset_real_pair = .false.
        else
            this%crop_real_pair = .false.
            this%reset_real_pair = .true.
        end if
    end subroutine Unscalable_set_real_pair_flags

    subroutine Unscalable_target(this)
        class(Unscalable_Dipolar_Interactions_Facade), intent(in) :: this

        call this%dipolar_interactions_static%box_size_memento_real%target(this%periodic_box)
        call this%target_real()
        call this%dipolar_interactions_static%box_size_memento_reci%target(this%periodic_box)
        call this%dipolar_interactions_static%reci_weight%target(this%dipolar_interactions_static%&
            box_size_memento_reci)
        call this%dipolar_interactions_static%reci_structure%target(this%periodic_box, this%&
            dipolar_interactions_static%box_size_memento_reci, this%components)
        call this%dipolar_interactions_dynamic%reci_visitor%target(this%&
            dipolar_interactions_static%box_size_memento_reci, this%dipolar_interactions_static%&
            reci_weight, this%dipolar_interactions_static%reci_structure)
        call this%dipolar_interactions_static%dlc_weight%target(this%periodic_box)
        call this%dipolar_interactions_static%dlc_structures%target(this%periodic_box, this%&
            components)
        call this%dipolar_interactions_dynamic%dlc_visitor%target(this%dipolar_interactions_static%&
            dlc_weight, this%dipolar_interactions_static%dlc_structures)
    end subroutine Unscalable_target

    subroutine Unscalable_clone(this, dipolar_interactions_static_target, &
        dipolar_interactions_static_source)
        class(Unscalable_Dipolar_Interactions_Facade), intent(in) :: this
        type(Dipolar_Interactions_Static_Wrapper), intent(inout) :: &
            dipolar_interactions_static_target
        type(Dipolar_Interactions_Static_Wrapper), intent(in) :: dipolar_interactions_static_source

        call this%clone_real(dipolar_interactions_static_target%box_size_memento_real, &
            dipolar_interactions_static_target%real_pair, dipolar_interactions_static_source%&
            box_size_memento_real, dipolar_interactions_static_source%real_pair)
        call box_destroy(dipolar_interactions_static_target%box_size_memento_reci)
        allocate(dipolar_interactions_static_target%box_size_memento_reci, &
            source=dipolar_interactions_static_source%box_size_memento_reci)
        call des_reci_destroy(dipolar_interactions_static_target%reci_weight)
        allocate(dipolar_interactions_static_target%reci_weight, &
            source=dipolar_interactions_static_source%reci_weight)
        call des_reci_destroy(dipolar_interactions_static_target%reci_structure)
        allocate(dipolar_interactions_static_target%reci_structure, &
            source=dipolar_interactions_static_source%reci_structure)
        call dlc_destroy(dipolar_interactions_static_target%dlc_weight)
        allocate(dipolar_interactions_static_target%dlc_weight, &
            source=dipolar_interactions_static_source%dlc_weight)
        call dlc_destroy(dipolar_interactions_static_target%dlc_structures)
        allocate(dipolar_interactions_static_target%dlc_structures, &
            source=dipolar_interactions_static_source%dlc_structures)
    end subroutine Unscalable_clone

    subroutine Unscalable_reset(this)
        class(Unscalable_Dipolar_Interactions_Facade), intent(in) :: this

        if (this%crop_real_pair) then
            call this%dipolar_interactions_static%real_pair%crop()
        end if
        call dipolar_interactions_reset(this%dipolar_interactions_static, this%reset_real_pair)
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

    subroutine Null_set_real_pair_flags(this, new_box_size)
        class(Null_Dipolar_Interactions_Facade), intent(inout) :: this
        real(DP), intent(in) :: new_box_size(:)
    end subroutine Null_set_real_pair_flags

    subroutine Null_save(this, dipolar_interactions_static, new_box_size)
        class(Null_Dipolar_Interactions_Facade), intent(inout) :: this
        type(Dipolar_Interactions_Static_Wrapper), intent(inout) :: dipolar_interactions_static
        real(DP), intent(in) :: new_box_size(:)
    end subroutine Null_save

    subroutine Null_target(this)
        class(Null_Dipolar_Interactions_Facade), intent(in) :: this
    end subroutine Null_target

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
