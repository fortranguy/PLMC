!> \brief Description of the neighbours class

module class_neighbourCells

use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
use data_precisions, only : DP
use data_cell, only : Ndim, Lsize
use data_neighbourCells, only : NnearCell_dim, NnearCell

implicit none
private
public :: Link

    type Link
        integer :: iCol
        type(Link), pointer :: next => null()
    end type Link
    
    type LinkedList
        type(Link), pointer :: particle => null()
    end type LinkedList

    type, public :: NeighbourCells
    
        private
        
        real(DP), dimension(Ndim), public :: cell_size
        integer, dimension(Ndim), public :: NtotalCell_dim
        integer :: NtotalCell
        integer, dimension(:, :), allocatable, public :: nearCells_from_totalCells
        type(LinkedList), dimension(:), allocatable, public :: beginCells
        type(LinkedList), dimension(:), allocatable :: currentCells
        type(LinkedList), dimension(:), allocatable :: nextCells        
        
    contains
    
        procedure :: construct => NeighbourCells_construct
        procedure :: destroy => NeighbourCells_destroy
        
        procedure :: alloc_cells => NeighbourCells_alloc_cells
        procedure :: dealloc_cells => NeighbourCells_dealloc_cells
        procedure :: check_cellsSize => NeighbourCells_check_cellsSize
        procedure :: position_to_cell => NeighbourCells_position_to_cell
        procedure :: all_cols_to_cells => NeighbourCells_all_cols_to_cells
        procedure :: remove_col_from_cell => NeighbourCells_remove_col_from_cell
        procedure :: add_col_to_cell => NeighbourCells_add_col_to_cell
        procedure, private :: totalCell_coord_to_index => NeighbourCells_cell_coord_to_ind
        procedure, private :: totalCell_PBC => NeighbourCells_totalCell_PBC
        procedure :: nearCells_from_totalCells_init => NeighbourCells_nearCells_from_totalCells_init
        
    end type NeighbourCells
    
contains

    subroutine NeighbourCells_construct(this, cell_size, rCut)
    
        class(NeighbourCells), intent(out) :: this
        real(DP), dimension(:), intent(in) :: cell_size
        real(DP), intent(in) :: rCut
        
        this%cell_size(:) = cell_size(:)
        this%NtotalCell_dim(:) = int(Lsize(:)/this%cell_size(:))
        this%NtotalCell = product(this%NtotalCell_dim)
        allocate(this%nearCells_from_totalCells(NnearCell, this%NtotalCell))
            
        call this%check_CellsSize(rCut)
    
    end subroutine NeighbourCells_construct
    
    subroutine NeighbourCells_destroy(this)
    
        class(NeighbourCells), intent(inout) :: this
        
        if (allocated(this%nearCells_from_totalCells)) then
            deallocate(this%nearCells_from_totalCells)
        end if
        
        call this%dealloc_cells()
        
    end subroutine NeighbourCells_destroy

    ! Linked-list allocation
    
    subroutine NeighbourCells_alloc_cells(this)
    
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
    
    recursive subroutine free_link(current)

        type(Link), pointer :: current
        
        if (associated(current%next)) then
            call free_link(current%next)
        end if
        deallocate(current)
        
    end subroutine free_link
    
    subroutine NeighbourCells_dealloc_cells(this)
    
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
        
        integer :: iDim
        
        do iDim = 1, Ndim
        
            if (this%cell_size(iDim) < rCut .and. this%cell_size(iDim) /= Lsize(iDim)/3._DP) then
                write(error_unit, *) "Warning : big rCut in the dimension", iDim, ":"
                write(error_unit, *) this%cell_size(iDim), "<", rCut
            end if
            
            if (Lsize(iDim)/2._DP*sqrt(3._DP) < rCut) then
                write(error_unit, *) "rCut too large in the dimension", iDim, ":"
                write(error_unit, *) Lsize(iDim)/2._DP*sqrt(3._DP), "<", rCut
                stop
            end if
            
            if (this%NtotalCell_dim(iDim) < NnearCell_dim(iDim)) then
                write(error_unit, *) "Too few cells in the dimension", iDim, ":"
                write(error_unit, *) this%NtotalCell_dim(iDim), "<", NnearCell_dim(iDim)
                stop
            end if
            
            if (modulo(Lsize(iDim), this%cell_size(iDim)) /= 0) then
                write(error_unit, *) "Cell size is not a divisor of the system size"
                write(error_unit, *) "in the dimension", iDim, ":"
                write(error_unit, *) "Lsize", Lsize(iDim)
                write(error_unit, *) "cell_size", this%cell_size(iDim)
                write(error_unit, *) "modulo(Lsize, cell_size) = ", &
                                      modulo(Lsize(iDim), this%cell_size(iDim))
                stop
            end if
            
        end do
        
    end subroutine NeighbourCells_check_cellsSize
    
    ! Assignment : particle -> cell
    
    pure function NeighbourCells_position_to_cell(this, xCol) result(position_to_cell)
    
        class(NeighbourCells), intent(in) :: this
        real(DP), dimension(:), intent(in) :: xCol
        integer :: position_to_cell
        
        integer, dimension(Ndim) :: cell_coord        
    
        cell_coord(:) = int(xCol(:)/this%cell_size(:)) + 1
        position_to_cell = cell_coord(1) + this%NtotalCell_dim(1)*(cell_coord(2)-1) + &
                           this%NtotalCell_dim(1)*this%NtotalCell_dim(2)*(cell_coord(3)-1)
    
    end function NeighbourCells_position_to_cell
    
    subroutine NeighbourCells_all_cols_to_cells(this, Ncol, X)
    
        class(NeighbourCells), intent(inout) :: this
        integer, intent(in) :: Ncol
        real(DP), dimension(:, :), intent(in) :: X
    
        integer :: iCol
        integer :: iCell
    
        do iCol = 1, Ncol
    
            iCell = this%position_to_cell(X(:,iCol))
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
        
        type(Link), pointer :: current => null()
        type(Link), pointer :: next => null(), previous => null()
    
        previous => this%beginCells(iCellOld)%particle
        current => previous%next
        
        do
        
            next => current%next
        
            if ( current%iCol == iCol ) then
            
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
    
        type(Link), pointer :: new => null()
        type(Link), pointer :: next => null(), previous => null()           
          
        
        previous => this%beginCells(iCellNew)%particle
        
        do
        
            next => previous%next
            
            if (.not. associated(next%next)) then
            
                allocate(new)
                new%next => previous%next
                previous%next => new
                new%iCol = iCol
                exit
                
            end if
            
            previous => next
            
        end do
            
    end  subroutine NeighbourCells_add_col_to_cell
    
    ! Neighbour cells ------------------------------------------------------------------------------

    pure function NeighbourCells_cell_coord_to_ind(this, totalCell_coord) &
                  result(totalCell_coord_to_index)
        
        class(NeighbourCells), intent(in) :: this
        integer, dimension(:), intent(in) :: totalCell_coord
        integer :: totalCell_coord_to_index
        
        totalCell_coord_to_index = totalCell_coord(1) + this%NtotalCell_dim(1)*(totalCell_coord(2)-1) &
                                   + this%NtotalCell_dim(1)*this%NtotalCell_dim(2)* &
                                     (totalCell_coord(3)-1)
    
    end function NeighbourCells_cell_coord_to_ind
    
    pure function nearCell_coord_to_index(nearCell_coord)
    
        integer, dimension(:), intent(in) :: nearCell_coord
        integer :: nearCell_coord_to_index
        
        nearCell_coord_to_index = nearCell_coord(1) + NnearCell_dim(1)*(nearCell_coord(2)-1) + &
                                  NnearCell_dim(1)*NnearCell_dim(2)*(nearCell_coord(3)-1)
    
    end function nearCell_coord_to_index
    
    pure function NeighbourCells_totalCell_PBC(this, totalCell_coord) result(totalCell_PBC)
    
        class(NeighbourCells), intent(in) :: this
        integer, dimension(:), intent(in) :: totalCell_coord
        integer, dimension(Ndim) :: totalCell_PBC
        
        totalCell_PBC(:) = totalCell_coord(:)
        
        where (totalCell_PBC(:) < 1)
            totalCell_PBC(:) = totalCell_PBC(:) + this%NtotalCell_dim(:)
        end where
        
        where (totalCell_PBC(:) > this%NtotalCell_dim(:))
            totalCell_PBC(:) = totalCell_PBC(:) - this%NtotalCell_dim(:)
        end where
    
    end function NeighbourCells_totalCell_PBC
    
    subroutine NeighbourCells_nearCells_from_totalCells_init(this)
    
        class(NeighbourCells), intent(inout) :: this
    
        integer :: iTotalCell, jTotalCell, kTotalCell, totalCell_index
        integer :: iNearCell, jNearCell, kNearCell, nearCell_index
        integer, dimension(Ndim) :: totalCell_coord, nearCell_coord
        
        do iTotalCell = 1, this%NtotalCell_dim(1)
        do jTotalCell = 1, this%NtotalCell_dim(2)
        do kTotalCell = 1, this%NtotalCell_dim(3)
            
            totalCell_index = this%totalCell_coord_to_index([iTotalCell, jTotalCell, kTotalCell])

            do iNearCell = 1, NnearCell_dim(1)
            do jNearCell = 1, NnearCell_dim(2)
            do kNearCell = 1, NnearCell_dim(3)
            
                nearCell_coord(:) = [iNearCell, jNearCell, kNearCell]
                nearCell_index = nearCell_coord_to_index(nearCell_coord(:))
                nearCell_coord(:) = nearCell_coord(:) - NnearCell_dim(:) + 1
                    ! with respect to the center (?) [iTotalCell, jTotalCell, kTotalCell]
                
                totalCell_coord(:) = [iTotalCell, jTotalCell, kTotalCell] + nearCell_coord(:)
                
                this%nearCells_from_totalCells(nearCell_index, totalCell_index) = &
                    this%totalCell_coord_to_index(this%totalCell_PBC(totalCell_coord(:)))
                    
            end do
            end do
            end do
        
        end do
        end do
        end do
            
    end subroutine NeighbourCells_nearCells_from_totalCells_init

end module class_neighbourCells
