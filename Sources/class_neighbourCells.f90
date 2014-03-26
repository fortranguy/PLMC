!> \brief Description of the neighbours class

module class_neighbourCells

use, intrinsic :: iso_fortran_env, only: error_unit
use data_precisions, only: DP, real_zero
use data_box, only: Ndim
use data_neighbourCells, only: NnearCell_dim, NnearCell
use module_types, only: Node, LinkedList
use module_physics_micro, only: index_from_coord, coord_PBC

implicit none
private

    type, public :: NeighbourCells
    
        private
        
        real(DP), dimension(Ndim) :: cell_size
        integer, dimension(Ndim) :: NtotalCell_dim
        integer :: NtotalCell
        integer, dimension(:, :), allocatable, public :: near_among_total
        type(LinkedList), dimension(:), allocatable, public :: beginCells
        type(LinkedList), dimension(:), allocatable :: currentCells
        type(LinkedList), dimension(:), allocatable :: nextCells
        
    contains
    
        procedure :: construct => NeighbourCells_construct
        procedure :: destroy => NeighbourCells_destroy
        
        procedure :: get_cell_size => NeighbourCells_get_cell_size
        procedure :: get_NtotalCell_dim => NeighbourCells_get_NtotalCell_dim
        
        procedure :: alloc_nodes => NeighbourCells_alloc_nodes
        procedure :: alloc_cells => NeighbourCells_alloc_cells
        procedure :: dealloc_nodes => NeighbourCells_dealloc_nodes
        procedure :: dealloc_cells => NeighbourCells_dealloc_cells
        procedure :: check_cellsSize => NeighbourCells_check_cellsSize
        procedure :: index_from_position => NeighbourCells_index_from_position
        procedure :: all_cols_to_cells => NeighbourCells_all_cols_to_cells
        procedure :: remove_col_from_cell => NeighbourCells_remove_col_from_cell
        procedure :: add_col_to_cell => NeighbourCells_add_col_to_cell
        procedure, private :: init_near_among_total => &
                              NeighbourCells_init_near_among_total
        
    end type NeighbourCells
    
contains

    subroutine NeighbourCells_construct(this, Box_size, proposed_cell_size, rCut)
    
        class(NeighbourCells), intent(out) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        real(DP), dimension(:), intent(in) :: proposed_cell_size
        real(DP), intent(in) :: rCut
        
        this%NtotalCell_dim(:) = floor(Box_size(:)/proposed_cell_size(:))
        this%NtotalCell = product(this%NtotalCell_dim)
        this%cell_size(:) = Box_size(:)/real(this%NtotalCell_dim(:), DP)
        
        allocate(this%near_among_total(NnearCell, this%NtotalCell))
            
        call this%check_CellsSize(Box_size, rCut)
        call this%alloc_cells()
        call this%init_near_among_total()
    
    end subroutine NeighbourCells_construct
    
    subroutine NeighbourCells_destroy(this)
    
        class(NeighbourCells), intent(inout) :: this
        
        if (allocated(this%near_among_total)) deallocate(this%near_among_total)
        
        call this%dealloc_cells()
        
    end subroutine NeighbourCells_destroy
    
    !> Accessors:
    
    pure function NeighbourCells_get_cell_size(this) result(get_cell_size)
        class(NeighbourCells), intent(in) :: this
        real(DP), dimension(Ndim) :: get_cell_size
        get_cell_size(:) = this%cell_size(:)
    end function NeighbourCells_get_cell_size
    
    pure function NeighbourCells_get_NtotalCell_dim(this) result(get_NtotalCell_dim)
        class(NeighbourCells), intent(in) :: this
        integer, dimension(Ndim) :: get_NtotalCell_dim
        get_NtotalCell_dim(:) = this%NtotalCell_dim(:)
    end function NeighbourCells_get_NtotalCell_dim

    !> Linked-list allocation
    
    subroutine NeighbourCells_alloc_nodes(this)
    
        class(NeighbourCells), intent(inout) :: this
    
        integer :: iCell
    
        do iCell = 1, this%NtotalCell

            allocate(this%beginCells(iCell)%particle)
            this%currentCells(iCell)%particle => this%beginCells(iCell)%particle
            this%currentCells(iCell)%particle%iCol = 0
            
            allocate(this%nextCells(iCell)%particle)
            this%nextCells(iCell)%particle%iCol = 0
            this%currentCells(iCell)%particle%next => this%nextCells(iCell)%particle
            this%currentCells(iCell)%particle => this%nextCells(iCell)%particle
    
        end do
    
    end subroutine NeighbourCells_alloc_nodes
    
    subroutine NeighbourCells_alloc_cells(this)
    
        class(NeighbourCells), intent(inout) :: this

        allocate(this%beginCells(this%NtotalCell))
        allocate(this%currentCells(this%NtotalCell))
        allocate(this%nextCells(this%NtotalCell))
        
        call this%alloc_nodes()
        
    end subroutine NeighbourCells_alloc_cells
    
    ! Linked-list deallocation
    
    recursive pure subroutine free_link(current)

        type(Node), pointer :: current
        
        if (associated(current%next)) then
            call free_link(current%next)
        end if
        deallocate(current)
        
    end subroutine free_link
    
    pure subroutine NeighbourCells_dealloc_nodes(this)
    
        class(NeighbourCells), intent(inout) :: this
    
        integer :: iCell

        do iCell = 1, this%NtotalCell
            if (associated(this%beginCells(iCell)%particle)) then
                call free_link(this%beginCells(iCell)%particle)
            end if
        end do
    
    end subroutine NeighbourCells_dealloc_nodes
    
    pure subroutine NeighbourCells_dealloc_cells(this)
    
        class(NeighbourCells), intent(inout) :: this
        
        call this%dealloc_nodes()
        
        if (allocated(this%beginCells)) deallocate(this%beginCells)
        if (allocated(this%currentCells)) deallocate(this%currentCells)
        if (allocated(this%nextCells)) deallocate(this%nextCells)
    
    end subroutine NeighbourCells_dealloc_cells
    
    ! NeighbourCells cells size check
    
    subroutine NeighbourCells_check_cellsSize(this, Box_size, rCut)
    
        class(NeighbourCells), intent(in) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        real(DP), intent(in) :: rCut
        
        integer :: jDim
        real(DP) :: Box_size_mod_cell_size
        
        do jDim = 1, Ndim
        
            if (this%cell_size(jDim) < rCut) then
                write(error_unit, *) "    rCut too big in the dimension", jDim, ": "
                write(error_unit, *) "    ", this%cell_size(jDim), "<", rCut
                error stop
            end if
            
            if (this%NtotalCell_dim(jDim) < NnearCell_dim(jDim)) then
                write(error_unit, *) "    Too few cells in the dimension", jDim, ": "
                write(error_unit, *) "    ", this%NtotalCell_dim(jDim), "<", NnearCell_dim(jDim)
                stop
            end if
            
            Box_size_mod_cell_size = modulo(Box_size(jDim), this%cell_size(jDim))
            if (Box_size_mod_cell_size > real_zero .and. &
            abs(Box_size_mod_cell_size - this%cell_size(jDim)) > real_zero) then
                write(error_unit, *) "    Cell size is not a divisor of the system size"
                write(error_unit, *) "    in the dimension", jDim, ": "
                write(error_unit, *) "    Box_size", Box_size(jDim)
                write(error_unit, *) "    cell_size", this%cell_size(jDim)
                write(error_unit, *) "    modulo(Box_size, cell_size) = ", Box_size_mod_cell_size
                stop
            end if
            
        end do
        
    end subroutine NeighbourCells_check_cellsSize
    
    ! Assignment: particle -> cell
    
    pure function NeighbourCells_index_from_position(this, xCol) result(index_from_position)
    
        class(NeighbourCells), intent(in) :: this
        real(DP), dimension(:), intent(in) :: xCol
        integer :: index_from_position
        
        integer, dimension(Ndim) :: cell_coord
    
        cell_coord(:) = int(xCol(:)/this%cell_size(:)) + 1
        index_from_position = cell_coord(1) + this%NtotalCell_dim(1)*(cell_coord(2)-1) + &
                             this%NtotalCell_dim(1)*this%NtotalCell_dim(2)*(cell_coord(3)-1)
    
    end function NeighbourCells_index_from_position
    
    pure subroutine NeighbourCells_all_cols_to_cells(this, Ncol, positions)
    
        class(NeighbourCells), intent(inout) :: this
        integer, intent(in) :: Ncol
        real(DP), dimension(:, :), intent(in) :: positions
    
        integer :: iCol
        integer :: iCell
    
        do iCol = 1, Ncol
    
            iCell = this%index_from_position(positions(:,iCol))
            this%currentCells(iCell)%particle%iCol = iCol
            
            allocate(this%nextCells(iCell)%particle)
            this%nextCells(iCell)%particle%iCol = 0
            this%currentCells(iCell)%particle%next => this%nextCells(iCell)%particle
            this%currentCells(iCell)%particle => this%nextCells(iCell)%particle
            
        end do
        
        do iCell = 1, this%NtotalCell
            this%currentCells(iCell)%particle%next => null()
        end do
        
    end subroutine NeighbourCells_all_cols_to_cells
    
    ! Neighbour cells update
    
    subroutine NeighbourCells_remove_col_from_cell(this, iCol, iCellOld)
    
        class(NeighbourCells), intent(inout) :: this
        integer, intent(in) :: iCol, iCellOld
        
        type(Node), pointer :: current => null()
        type(Node), pointer :: next => null(), previous => null()
    
        previous => this%beginCells(iCellOld)%particle
        current => previous%next
        
        do
        
            next => current%next
        
            if (current%iCol == iCol) then
                previous%next => current%next
                deallocate(current)
                current => next
                exit
            else
                previous => current
            end if
            
            current => next
        
        end do
            
    end subroutine NeighbourCells_remove_col_from_cell
    
    subroutine NeighbourCells_add_col_to_cell(this, iCol, iCellNew)
    
        class(NeighbourCells), intent(inout) :: this
        integer, intent(in) :: iCol, iCellNew
    
        type(Node), pointer :: new => null()
        type(Node), pointer :: next => null(), previous => null()
        
        previous => this%beginCells(iCellNew)%particle
        next => previous%next

        allocate(new)
        new%next => previous%next
        previous%next => new
        new%iCol = iCol

    end subroutine NeighbourCells_add_col_to_cell
    
    ! Neighbour cells initialisation
    
    pure subroutine NeighbourCells_init_near_among_total(this)
    
        class(NeighbourCells), intent(inout) :: this
    
        integer :: iTotalCell, jTotalCell, kTotalCell, totalCell_index
        integer :: iNearCell, jNearCell, kNearCell, nearCell_index
        integer, dimension(Ndim) :: totalCell_coord, nearCell_coord
        
        do kTotalCell = 1, this%NtotalCell_dim(3)
        do jTotalCell = 1, this%NtotalCell_dim(2)
        do iTotalCell = 1, this%NtotalCell_dim(1)
            
            totalCell_index = index_from_coord([iTotalCell, jTotalCell, kTotalCell], &
                                               this%NtotalCell_dim)

            do kNearCell = 1, NnearCell_dim(3)
            do jNearCell = 1, NnearCell_dim(2)
            do iNearCell = 1, NnearCell_dim(1)
            
                nearCell_coord(:) = [iNearCell, jNearCell, kNearCell]
                nearCell_index = index_from_coord(nearCell_coord, NnearCell_dim)
                nearCell_coord(:) = nearCell_coord(:) - NnearCell_dim(:) + 1
                    ! with respect to the center (?) [iTotalCell, jTotalCell, kTotalCell]
                
                totalCell_coord(:) = [iTotalCell, jTotalCell, kTotalCell] + nearCell_coord(:)
                totalCell_coord(:) = coord_PBC(totalCell_coord, this%NtotalCell_dim)
                
                this%near_among_total(nearCell_index, totalCell_index) = &
                    index_from_coord(totalCell_coord, this%NtotalCell_dim)
                    
            end do
            end do
            end do
        
        end do
        end do
        end do
            
    end subroutine NeighbourCells_init_near_among_total

end module class_neighbourCells
