module classes_maximum_box_compression_explorer

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_errors, only: error_exit
use types_environment_wrapper, only: Environment_Wrapper
use types_component_wrapper, only: Component_Wrapper
use types_short_interactions_wrapper, only: Short_Interactions_Wrapper
use procedures_short_interactions_visitor, only: short_interactions_visit_cells => visit_cells
use classes_maximum_box_compression, only: Abstract_Maximum_Box_Compression
use procedures_maximum_box_compression_factory, only: maximum_box_compression_destroy => destroy
use types_exploring_observables_wrapper, only: Exploring_Observables_Wrapper
use classes_exploring_algorithm, only: Abstract_Exploring_Algorithm

implicit none

private

    type, extends(Abstract_Exploring_Algorithm), abstract, public :: &
        Abstract_Maximum_Box_Compression_Explorer
    private
        type(Environment_Wrapper), pointer :: environment => null()
        type(Component_Wrapper), pointer :: components(:, :) => null()
        type(Short_Interactions_Wrapper), pointer :: short_interactions => null()
        class(Abstract_Maximum_Box_Compression), allocatable :: maximum_box_compression
        real(DP) :: min_distance
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: try => Abstract_try
        procedure, private :: set_min_distance => Abstract_set_min_distance
    end type Abstract_Maximum_Box_Compression_Explorer

    type, extends(Abstract_Maximum_Box_Compression_Explorer), public :: &
        Concrete_Maximum_Box_Compression_Explorer
    end type Concrete_Maximum_Box_Compression_Explorer

    type, extends(Abstract_Maximum_Box_Compression_Explorer), public :: &
        Null_Maximum_Box_Compression_Explorer
    contains
        procedure :: Null_construct
        procedure :: Null_destroy
        procedure :: Null_try
    end type Null_Maximum_Box_Compression_Explorer

contains

!implementation Abstract_Maximum_Box_Compression_Explorer

    subroutine Abstract_construct(this, environment, components, short_interactions, &
        maximum_box_compression)
        class(Abstract_Maximum_Box_Compression_Explorer), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        type(Component_Wrapper), target, intent(in) :: components(:, :)
        type(Short_Interactions_Wrapper), target, intent(in) :: short_interactions
        class(Abstract_Maximum_Box_Compression), intent(in) :: maximum_box_compression

        this%environment => environment
        this%components => components
        this%short_interactions => short_interactions
        call this%set_min_distance()
        allocate(this%maximum_box_compression, source=maximum_box_compression)
    end subroutine Abstract_construct

    !> @warning periodic_boxes(1): arbitrary
    subroutine Abstract_set_min_distance(this)
        class(Abstract_Maximum_Box_Compression_Explorer), intent(inout) :: this

        integer :: i_component, j_component
        real(DP) :: min_distance_ij

        this%min_distance = this%environment%periodic_boxes(1)%get_max_distance()
        do j_component = 1, size(this%components, 1)
            do i_component = 1, j_component
                min_distance_ij = this%short_interactions%components_pairs(j_component)%&
                    line(i_component)%potential%get_min_distance()
                if (min_distance_ij < this%min_distance) this%min_distance = min_distance_ij
            end do
        end do
    end subroutine Abstract_set_min_distance

    subroutine Abstract_destroy(this)
        class(Abstract_Maximum_Box_Compression_Explorer), intent(inout) :: this

        call maximum_box_compression_destroy(this%maximum_box_compression)
        this%short_interactions => null()
        this%components => null()
        this%environment => null()
    end subroutine Abstract_destroy

    subroutine Abstract_try(this, observables)
        class(Abstract_Maximum_Box_Compression_Explorer), intent(in) :: this
        type(Exploring_Observables_Wrapper), intent(inout) :: observables

        logical :: overlap
        integer :: i_box
        real(DP) :: min_distance_ratio, max_distance_ratio

        do i_box = 1, size(this%environment%periodic_boxes)
            max_distance_ratio = this%environment%periodic_boxes(i_box)%get_max_distance() / this%min_distance
            call short_interactions_visit_cells(overlap, min_distance_ratio, max_distance_ratio, this%&
                components(:, i_box), this%short_interactions%cells(i_box)%visitable_cells)
            if (overlap) call error_exit("Abstract_Maximum_Box_Compression_Explorer: try: "//&
                "short_interactions_visit_cells: overlap")
            observables%maximum_boxes_compression_delta(i_box) = this%maximum_box_compression%&
                get_delta(min_distance_ratio)
        end do
    end subroutine Abstract_try

!end implementation Abstract_Maximum_Box_Compression_Explorer

!implementation Null_Maximum_Box_Compression_Explorer

    subroutine Null_construct(this, environment, components, short_interactions)
        class(Null_Maximum_Box_Compression_Explorer), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        type(Component_Wrapper), target, intent(in) :: components(:, :)
        type(Short_Interactions_Wrapper), target, intent(in) :: short_interactions
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Maximum_Box_Compression_Explorer), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_try(this, min_distance_ratio)
        class(Null_Maximum_Box_Compression_Explorer), intent(in) :: this
        real(DP), intent(out) :: min_distance_ratio
        min_distance_ratio = 0._DP
    end subroutine Null_try

!implementation Null_Maximum_Box_Compression_Explorer

end module classes_maximum_box_compression_explorer
