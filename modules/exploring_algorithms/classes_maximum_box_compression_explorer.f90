module classes_maximum_box_compression_explorer

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_errors, only: error_exit
use classes_periodic_box, only: Abstract_Periodic_Box
use procedures_box_size, only: box_size_max_distance => max_distance
use types_component_wrapper, only: Component_Wrapper
use classes_pair_potential, only: Pair_Potential_Line
use types_cells_wrapper, only: Cells_Wrapper
use procedures_short_interactions_visitor, only: short_interactions_visit => visit
use classes_maximum_box_compression, only: Abstract_Maximum_Box_Compression
use procedures_maximum_box_compression_factory, only: maximum_box_compression_destroy => destroy

implicit none

private

    type, abstract, public :: Abstract_Maximum_Box_Compression_Explorer
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
        type(Component_Wrapper), pointer :: components(:) => null()
        type(Pair_Potential_Line), pointer :: components_pairs(:) => null()
        type(Cells_Wrapper), pointer :: cells => null()
        class(Abstract_Maximum_Box_Compression), allocatable :: maximum_box_compression
        real(DP) :: min_distance = 0._DP
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: reset => Abstract_set_min_distance
        procedure :: try => Abstract_try
    end type Abstract_Maximum_Box_Compression_Explorer

    type, extends(Abstract_Maximum_Box_Compression_Explorer), public :: &
        Concrete_Maximum_Box_Compression_Explorer
    end type Concrete_Maximum_Box_Compression_Explorer

    type, extends(Abstract_Maximum_Box_Compression_Explorer), public :: &
        Null_Maximum_Box_Compression_Explorer
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: reset => Null_set_min_distance
        procedure :: try => Null_try
    end type Null_Maximum_Box_Compression_Explorer

contains

!implementation Abstract_Maximum_Box_Compression_Explorer

    subroutine Abstract_construct(this, periodic_box, components, components_pairs, cells, &
        maximum_box_compression)
        class(Abstract_Maximum_Box_Compression_Explorer), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        type(Component_Wrapper), target, intent(in) :: components(:)
        type(Pair_Potential_Line), target, intent(in) :: components_pairs(:)
        type(Cells_Wrapper), target, intent(in) :: cells
        class(Abstract_Maximum_Box_Compression), intent(in) :: maximum_box_compression

        this%periodic_box => periodic_box
        this%components => components
        this%components_pairs => components_pairs
        this%cells => cells
        allocate(this%maximum_box_compression, source=maximum_box_compression)
    end subroutine Abstract_construct

    subroutine Abstract_set_min_distance(this)
        class(Abstract_Maximum_Box_Compression_Explorer), intent(inout) :: this

        integer :: i_component, j_component
        real(DP) :: min_distance_ij

        this%min_distance = box_size_max_distance(this%periodic_box%get_size())
        do j_component = 1, size(this%components)
            do i_component = 1, j_component
                min_distance_ij = this%components_pairs(j_component)%&
                    line(i_component)%potential%get_min_distance()
                if (min_distance_ij < this%min_distance) this%min_distance = min_distance_ij
            end do
        end do
    end subroutine Abstract_set_min_distance

    subroutine Abstract_destroy(this)
        class(Abstract_Maximum_Box_Compression_Explorer), intent(inout) :: this

        call maximum_box_compression_destroy(this%maximum_box_compression)
        this%cells => null()
        this%components_pairs => null()
        this%components => null()
        this%periodic_box => null()
    end subroutine Abstract_destroy

    subroutine Abstract_try(this, maximum_box_compression_delta)
        class(Abstract_Maximum_Box_Compression_Explorer), intent(in) :: this
        real(DP), intent(out) :: maximum_box_compression_delta

        logical :: overlap
        real(DP) :: min_distance_ratio, max_distance_ratio

        max_distance_ratio = box_size_max_distance(this%periodic_box%get_size()) / this%min_distance
        call short_interactions_visit(overlap, min_distance_ratio, max_distance_ratio, &
            this%components, this%cells%visitable_cells)
        if (overlap) call error_exit("Abstract_Maximum_Box_Compression_Explorer: try: "//&
            "short_interactions_visit: overlap")
        maximum_box_compression_delta = this%maximum_box_compression%get_delta(min_distance_ratio)
    end subroutine Abstract_try

!end implementation Abstract_Maximum_Box_Compression_Explorer

!implementation Null_Maximum_Box_Compression_Explorer

    subroutine Null_construct(this, periodic_box, components, components_pairs, cells, &
        maximum_box_compression)
        class(Null_Maximum_Box_Compression_Explorer), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        type(Component_Wrapper), target, intent(in) :: components(:)
        type(Pair_Potential_Line), target, intent(in) :: components_pairs(:)
        type(Cells_Wrapper), target, intent(in) :: cells
        class(Abstract_Maximum_Box_Compression), intent(in) :: maximum_box_compression
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Maximum_Box_Compression_Explorer), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_set_min_distance(this)
        class(Null_Maximum_Box_Compression_Explorer), intent(inout) :: this
    end subroutine Null_set_min_distance

    subroutine Null_try(this, maximum_box_compression_delta)
        class(Null_Maximum_Box_Compression_Explorer), intent(in) :: this
        real(DP), intent(out) :: maximum_box_compression_delta
        maximum_box_compression_delta = 0._DP
    end subroutine Null_try

!implementation Null_Maximum_Box_Compression_Explorer

end module classes_maximum_box_compression_explorer
