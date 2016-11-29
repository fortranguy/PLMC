module classes_dipolar_neighbourhoods_visitor

use types_logical_wrapper, only: Logical_Rectangle
use procedures_errors, only: error_exit
use types_component_wrapper, only: Component_Wrapper
use classes_visitable_cells, only: Abstract_Visitable_Cells
use procedures_short_interactions_visitor, only: short_interactions_visit => visit

implicit none

private

    type, abstract, public :: Abstract_Dipolar_Neighbourhoods_Visitor
    private
        type(Component_Wrapper), pointer :: components(:) => null()
        class(Abstract_Visitable_Cells), pointer :: visitable_cells(:, :) => null()
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: try => Abstract_try
    end type Abstract_Dipolar_Neighbourhoods_Visitor

    type, extends(Abstract_Dipolar_Neighbourhoods_Visitor), public :: &
        Concrete_Dipolar_Neighbourhoods_Visitor

    end type Concrete_Dipolar_Neighbourhoods_Visitor

    type, extends(Abstract_Dipolar_Neighbourhoods_Visitor), public :: &
        Null_Dipolar_Neighbourhoods_Visitor
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: try => Null_try
    end type Null_Dipolar_Neighbourhoods_Visitor

contains

!implementation Abstract_Dipolar_Neighbourhoods_Visitor

    subroutine Abstract_construct(this, components, visitable_cells)
        class(Abstract_Dipolar_Neighbourhoods_Visitor), intent(out) :: this
        type(Component_Wrapper), target, intent(in) :: components(:)
        class(Abstract_Visitable_Cells), target, intent(in) :: visitable_cells(:, :)

        this%components => components
        this%visitable_cells => visitable_cells
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Dipolar_Neighbourhoods_Visitor), intent(inout) :: this

        this%visitable_cells => null()
        this%components => null()
    end subroutine Abstract_destroy

    subroutine Abstract_try(this, adjacency_matrices)
        class(Abstract_Dipolar_Neighbourhoods_Visitor), intent(in) :: this
        type(Logical_Rectangle), intent(out) :: adjacency_matrices(:, :)

        logical :: overlap
        integer :: i_component, j_component, nums_dipoles(size(this%components))

        do i_component = 1, size(nums_dipoles)
            nums_dipoles(i_component) = this%components(i_component)%dipole_moments%get_num()
        end do

        do j_component = 1, size(this%visitable_cells, 2)
            do i_component = 1, size(this%visitable_cells, 1)
                if (allocated(adjacency_matrices(i_component, j_component)%rectangle)) &
                    deallocate(adjacency_matrices(i_component, j_component)%rectangle)
                allocate(adjacency_matrices(i_component, j_component)%&
                    rectangle(nums_dipoles(i_component), nums_dipoles(j_component)))
                adjacency_matrices(i_component, j_component)%rectangle = .false.
            end do
        end do

        call short_interactions_visit(overlap, adjacency_matrices, this%components, this%&
            visitable_cells)
        if (overlap) call error_exit("Abstract_Dipolar_Neighbourhoods_Visitor: try: "//&
            "short_interactions_visit: overlap")
    end subroutine Abstract_try

!end implementation Abstract_Dipolar_Neighbourhoods_Visitor

!implementation Null_Dipolar_Neighbourhoods_Visitor

    subroutine Null_construct(this, components, visitable_cells)
        class(Null_Dipolar_Neighbourhoods_Visitor), intent(out) :: this
        type(Component_Wrapper), target, intent(in) :: components(:)
        class(Abstract_Visitable_Cells), target, intent(in) :: visitable_cells(:, :)
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Dipolar_Neighbourhoods_Visitor), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_try(this, adjacency_matrices)
        class(Null_Dipolar_Neighbourhoods_Visitor), intent(in) :: this
        type(Logical_Rectangle), intent(out) :: adjacency_matrices(:, :)
        integer :: i_component, j_component
        do j_component = 1, size(this%visitable_cells, 2)
            do i_component = 1, size(this%visitable_cells, 1)
                if (allocated(adjacency_matrices(i_component, j_component)%rectangle)) &
                    deallocate(adjacency_matrices(i_component, j_component)%rectangle)
                allocate(adjacency_matrices(i_component, j_component)%rectangle(0, 0))
                adjacency_matrices(i_component, j_component)%rectangle = .false.
            end do
        end do
    end subroutine Null_try

!end implementation Null_Dipolar_Neighbourhoods_Visitor

end module classes_dipolar_neighbourhoods_visitor
