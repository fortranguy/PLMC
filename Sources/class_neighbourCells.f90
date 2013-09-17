!> \brief Description of the neighbours class

module class_neighbourCells

use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
use data_precisions, only : DP
use data_box, only : Ndim, Lsize
use data_neighbourCells, only : NnearCell_dim, NnearCell
use module_physics, only : index_from_coord, coord_PBC

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
        
        procedure :: getCell_size => NeighbourCells_getCell_size
        procedure :: getNtotalCell_dim => NeighbourCells_getNtotalCell_dim
        
        procedure :: alloc_cells => NeighbourCells_alloc_cells
        procedure :: dealloc_cells => NeighbourCells_dealloc_cells
        procedure :: check_cellsSize => NeighbourCells_check_cellsSize
        procedure :: index_from_position => NeighbourCells_index_from_position
        procedure :: all_cols_to_cells => NeighbourCells_all_cols_to_cells
        procedure :: remove_col_from_cell => NeighbourCells_remove_col_from_cell
        procedure :: add_col_to_cell => NeighbourCells_add_col_to_cell
        procedure :: nearCells_among_totalCells_init => NeighbourCells_nearCells_among_totalCells_init
        
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
            if (proposed_cell_size(jDim) /= this%cell_size(jDim)) then            
                write(error_unit, *) "    Warning : cell size in the dimension", jDim, &
                                      "was changed."
                write(error_unit, *) "    ", proposed_cell_size(jDim), "->", this%cell_size(jDim)            
            end if        
        end do
        
        allocate(this%nearCells_among_totalCells(NnearCell, this%NtotalCell))
            
        call this%check_CellsSize(rCut)        
        call this%alloc_cells()
        call this%nearCells_among_totalCells_init()
    
    end subroutine NeighbourCells_construct
    
    subroutine NeighbourCells_destroy(this)
    
        class(NeighbourCells), intent(inout) :: this
        
        if (allocated(this%nearCells_among_totalCells)) then
            deallocate(this%nearCells_among_totalCells)
        end if
        
        call this%dealloc_cells()
        
    end subroutine NeighbourCells_destroy
    
    !> Accessor : cell_size
    
    pure function NeighbourCells_getCell_size(this) result(getCell_size)
    
        class(NeighbourCells), intent(in) :: this
        real(DP), dimension(Ndim) :: getCell_size
        
        getCell_size(:) = this%cell_size(:)
    
    end function NeighbourCells_getCell_size
    
    !> Accessor : NtotalCell_dim
    
    pure function NeighbourCells_getNtotalCell_dim(this) result(getNtotalCell_dim)
    
        class(NeighbourCells), intent(in) :: this    
        integer, dimension(Ndim) :: getNtotalCell_dim
        
        getNtotalCell_dim(:) = this%NtotalCell_dim(:)
        
    end function NeighbourCells_getNtotalCell_dim

    !> Linked-list allocation
    
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
        
        integer :: jDim
        
        do jDim = 1, Ndim
        
            if (this%cell_size(jDim) < rCut .and. this%cell_size(jDim) /= Lsize(jDim)/3._DP) then
                write(error_unit, *) "  Warning : big rCut in the dimension", jDim, ":"
                write(error_unit, *) "  ", this%cell_size(jDim), "<", rCut
            end if
            
            if (Lsize(jDim)/2._DP < rCut) then
                write(error_unit, *) "  rCut too large in the dimension", jDim, ":"
                write(error_unit, *) "  ", Lsize(jDim)/2._DP, "<", rCut
                stop
            end if
            
            if (this%NtotalCell_dim(jDim) < NnearCell_dim(jDim)) then
                write(error_unit, *) "  Too few cells in the dimension", jDim, ":"
                write(error_unit, *) "  ", this%NtotalCell_dim(jDim), "<", NnearCell_dim(jDim)
                stop
            end if
            
            if (modulo(Lsize(jDim), this%cell_size(jDim)) /= 0) then
                write(error_unit, *) "  Cell size is not a divisor of the system size"
                write(error_unit, *) "  in the dimension", jDim, ":"
                write(error_unit, *) "  Lsize", Lsize(jDim)
                write(error_unit, *) "  cell_size", this%cell_size(jDim)
                write(error_unit, *) "  modulo(Lsize, cell_size) = ", &
                                      modulo(Lsize(jDim), this%cell_size(jDim))
                stop
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
    
    subroutine NeighbourCells_all_cols_to_cells(this, Ncol, X)
    
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
        
        type(Link), pointer :: current => null()
        type(Link), pointer :: next => null(), previous => null()
    
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
    
    ! Neighbour cells initialisation
    
    subroutine NeighbourCells_nearCells_among_totalCells_init(this)
    
        class(NeighbourCells), intent(inout) :: this
    
        integer :: iTotalCell, jTotalCell, kTotalCell, totalCell_index
        integer :: iNearCell, jNearCell, kNearCell, nearCell_index
        integer, dimension(Ndim) :: totalCell_coord, nearCell_coord
        
        do iTotalCell = 1, this%NtotalCell_dim(1)
        do jTotalCell = 1, this%NtotalCell_dim(2)
        do kTotalCell = 1, this%NtotalCell_dim(3)
            
            totalCell_index = index_from_coord([iTotalCell, jTotalCell, kTotalCell], &
                                               this%NtotalCell_dim)

            do iNearCell = 1, NnearCell_dim(1)
            do jNearCell = 1, NnearCell_dim(2)
            do kNearCell = 1, NnearCell_dim(3)
            
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
            
    end subroutine NeighbourCells_nearCells_among_totalCells_init

end module class_neighbourCells
