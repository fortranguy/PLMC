!> \brief Description of the neighbours class

module class_neighbourCells

use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
use data_precisions, only : DP, real_zero
use data_box, only : Ndim, Lsize
use data_neighbourCells, only : NnearCell_dim, NnearCell
use module_physics, only : index_from_coord, coord_PBC

implicit none
private
public :: Node

    type Node
        integer :: iCol
        type(Node), pointer :: next => null()
    end type Node
    
    type LinkedList
        type(Node), pointer :: particle => null()
    end type LinkedList

    type, public :: NeighbourCells
    
        private
        
        real(DP), dimension(Ndim) :: cell_size
        integer, dimension(Ndim) :: NtotalCell_dim
        integer :: NtotalCell
        integer, dimension(:, :), allocatable, public :: nearCells_among_totalCells
        type(LinkedList), dimension(:), allocatable, public :: beginCells
        type(LinkedList), dimension(:), allocatable :: currentCells
        type(LinkedList), dimension(:), allocatable :: nextCells        
        
    contains
    
        procedure :: construct => NeighbourCells_construct
        procedure :: destroy => NeighbourCells_destroy
        
        procedure :: get_cell_size => NeighbourCells_get_cell_size
        procedure :: get_NtotalCell_dim => NeighbourCells_get_NtotalCell_dim
        
        procedure :: alloc_cells => NeighbourCells_alloc_cells
        procedure :: dealloc_cells => NeighbourCells_dealloc_cells
        procedure :: check_cellsSize => NeighbourCells_check_cellsSize
        procedure :: index_from_position => NeighbourCells_index_from_position
        procedure :: all_cols_to_cells => NeighbourCells_all_cols_to_cells
        procedure :: remove_col_from_cell => NeighbourCells_remove_col_from_cell
        procedure :: add_col_to_cell => NeighbourCells_add_col_to_cell
        procedure, private :: init_nearCells_among_totalCells => &
                              NeighbourCells_init_nearCells_among_totalCells
        
    end type NeighbourCells
    
contains

    subroutine NeighbourCells_construct(this, proposed_cell_size, rCut)
    
        class(NeighbourCells), intent(out) :: this
        real(DP), dimension(:), intent(in) :: proposed_cell_size
        real(DP), intent(in) :: rCut
        
        integer :: jDim
        
        this%NtotalCell_dim(:) = floor(Lsize(:)/proposed_cell_size(:))
        this%NtotalCell = product(this%NtotalCell_dim)
        this%cell_size(:) = Lsize(:)/real(this%NtotalCell_dim(:), DP)
        
        do jDim=1, Ndim        
            if (proposed_cell_size(jDim) < this%cell_size(jDim)) then            
                write(error_unit, *) "    Warning : cell size in the dimension", jDim, &
                                      "was changed."
                write(error_unit, *) "    ", proposed_cell_size(jDim), "->", this%cell_size(jDim)            
            end if        
        end do
        
        allocate(this%nearCells_among_totalCells(NnearCell, this%NtotalCell))
            
        call this%check_CellsSize(rCut)        
        call this%alloc_cells()
        call this%init_nearCells_among_totalCells()
    
    end subroutine NeighbourCells_construct
    
    pure subroutine NeighbourCells_destroy(this)
    
        class(NeighbourCells), intent(inout) :: this
        
        if (allocated(this%nearCells_among_totalCells)) then
            deallocate(this%nearCells_among_totalCells)
        end if
        
        call this%dealloc_cells()
        
    end subroutine NeighbourCells_destroy
    
    !> Accessor : cell_size
    
    pure function NeighbourCells_get_cell_size(this) result(get_cell_size)
    
        class(NeighbourCells), intent(in) :: this
        real(DP), dimension(Ndim) :: get_cell_size
        
        get_cell_size(:) = this%cell_size(:)
    
    end function NeighbourCells_get_cell_size
    
    !> Accessor : NtotalCell_dim
    
    pure function NeighbourCells_get_NtotalCell_dim(this) result(get_NtotalCell_dim)
    
        class(NeighbourCells), intent(in) :: this    
        integer, dimension(Ndim) :: get_NtotalCell_dim
        
        get_NtotalCell_dim(:) = this%NtotalCell_dim(:)
        
    end function NeighbourCells_get_NtotalCell_dim

    !> Linked-list allocation
    
    pure subroutine NeighbourCells_alloc_cells(this)
    
        class(NeighbourCells), intent(inout) :: this
    
        integer :: iCell

        allocate(this%beginCells(this%NtotalCell))
        allocate(this%currentCells(this%NtotalCell))
        allocate(this%nextCells(this%NtotalCell))
        
        do iCell = 1, this%NtotalCell

            allocate(this%beginCells(iCell)%particle)
            this%currentCells(iCell)%particle => this%beginCells(iCell)%particle
            this%currentCells(iCell)%particle%iCol = 0
            
            allocate(this%nextCells(iCell)%particle)
            this%nextCells(iCell)%particle%iCol = 0            
            this%currentCells(iCell)%particle%next => this%nextCells(iCell)%particle
            this%currentCells(iCell)%particle => this%nextCells(iCell)%particle     
    
        end do
        
    end subroutine NeighbourCells_alloc_cells
    
    ! Linked-list deallocation
    
    recursive pure subroutine free_link(current)

        type(Node), pointer :: current
        
        if (associated(current%next)) then
            call free_link(current%next)
        end if
        deallocate(current)
        
    end subroutine free_link
    
    pure subroutine NeighbourCells_dealloc_cells(this)
    
        class(NeighbourCells), intent(inout) :: this
    
        integer :: iCell

        do iCell = 1, this%NtotalCell
            if (associated(this%beginCells(iCell)%particle)) then
                call free_link(this%beginCells(iCell)%particle)
            end if
        end do
    
    end subroutine NeighbourCells_dealloc_cells
    
    ! NeighbourCells cells size check
    
    subroutine NeighbourCells_check_cellsSize(this, rCut)
    
        class(NeighbourCells), intent(in) :: this
        real(DP), intent(in) :: rCut
        
        integer :: jDim
        
        do jDim = 1, Ndim
        
            if (this%cell_size(jDim) < rCut) then
                write(error_unit, *) "    Warning : big rCut in the dimension", jDim, ":"
                write(error_unit, *) "    ", this%cell_size(jDim), "<", rCut
            end if
            
            if (Lsize(jDim)/2._DP < rCut) then
                write(error_unit, *) "    rCut too large in the dimension", jDim, ":"
                write(error_unit, *) "    ", Lsize(jDim)/2._DP, "<", rCut
                stop
            end if
            
            if (this%NtotalCell_dim(jDim) < NnearCell_dim(jDim)) then
                write(error_unit, *) "    Too few cells in the dimension", jDim, ":"
                write(error_unit, *) "    ", this%NtotalCell_dim(jDim), "<", NnearCell_dim(jDim)
                stop
            end if
            
            if (modulo(Lsize(jDim), this%cell_size(jDim)) > real_zero) then
                if (abs(modulo(Lsize(jDim), this%cell_size(jDim)) - this%cell_size(jDim)) > real_zero) then
                    write(error_unit, *) "    Cell size is not a divisor of the system size"
                    write(error_unit, *) "    in the dimension", jDim, ":"
                    write(error_unit, *) "    Lsize", Lsize(jDim)
                    write(error_unit, *) "    cell_size", this%cell_size(jDim)
                    write(error_unit, *) "    modulo(Lsize, cell_size) = ", &
                                          modulo(Lsize(jDim), this%cell_size(jDim))
                    stop
                end if
            end if
            
        end do
        
    end subroutine NeighbourCells_check_cellsSize
    
    ! Assignment : particle -> cell
    
    pure function NeighbourCells_index_from_position(this, xCol) result(index_from_position)
    
        class(NeighbourCells), intent(in) :: this
        real(DP), dimension(:), intent(in) :: xCol
        integer :: index_from_position
        
        integer, dimension(Ndim) :: cell_coord        
    
        cell_coord(:) = int(xCol(:)/this%cell_size(:)) + 1
        index_from_position = cell_coord(1) + this%NtotalCell_dim(1)*(cell_coord(2)-1) + &
                             this%NtotalCell_dim(1)*this%NtotalCell_dim(2)*(cell_coord(3)-1)
    
    end function NeighbourCells_index_from_position
    
    pure subroutine NeighbourCells_all_cols_to_cells(this, Ncol, X)
    
        class(NeighbourCells), intent(inout) :: this
        integer, intent(in) :: Ncol
        real(DP), dimension(:, :), intent(in) :: X
    
        integer :: iCol
        integer :: iCell
    
        do iCol = 1, Ncol
    
            iCell = this%index_from_position(X(:,iCol))
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
    
    pure subroutine NeighbourCells_init_nearCells_among_totalCells(this)
    
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
                
                this%nearCells_among_totalCells(nearCell_index, totalCell_index) = &
                    index_from_coord(totalCell_coord, this%NtotalCell_dim)
                    
            end do
            end do
            end do
        
        end do
        end do
        end do
            
    end subroutine NeighbourCells_init_nearCells_among_totalCells

end module class_neighbourCells
