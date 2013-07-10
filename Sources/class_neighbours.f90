!> \brief Description of the neighbours class

module class_neighbourCells

use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
use data_precisions, only : DP
use data_cell, only : Ndim, Lsize
use data_neighbourCells, only : NnearNeighCell_dim, NnearNeighCell

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

    type, public :: Neighbours
        
        real(DP), dimension(Ndim) :: cell_Lsize
        integer, dimension(Ndim) :: NtotalNeighCell_dim
        integer, dimension(:, :), allocatable :: cell_neighs
        type(LinkedList), dimension(:), allocatable :: cells
        type(LinkedList), dimension(:), allocatable :: cellsNext
        type(LinkedList), dimension(:), allocatable :: cellsBegin
        
    contains
    
        procedure :: construct => Neighbours_construct
        procedure :: destroy => Neighbours_destroy
        
        procedure :: alloc_cells => Neighbours_alloc_cells
        procedure :: dealloc_cells => Neighbours_dealloc_cells
        procedure :: check_cellsSize => Neighbours_check_cellsSize
        procedure :: position_to_cell => Neighbours_position_to_cell
        procedure :: all_cols_to_cells => Neighbours_all_cols_to_cells
        procedure :: remove_col_from_cell => Neighbours_remove_col_from_cell
        procedure :: add_col_to_cell => Neighbours_add_col_to_cell
        procedure, private :: cell_coord_to_ind => Neighbours_cell_coord_to_ind
        procedure, private :: cell_period => Neighbours_cell_period
        procedure :: cell_neighs_init => Neighbours_cell_neighs_init
        
    end type Neighbours
    
contains

    subroutine Neighbours_construct(this, cell_Lsize, rCut)
    
        class(Neighbours), intent(out) :: this
        real(DP), dimension(:), intent(in) :: cell_Lsize
        real(DP), intent(in) :: rCut
        
        this%cell_Lsize(:) = cell_Lsize(:)
        this%NtotalNeighCell_dim(:) = int(Lsize(:)/this%cell_Lsize(:))
        allocate(this%cell_neighs(NnearNeighCell, product(this%NtotalNeighCell_dim)))
            
        call this%check_CellsSize(rCut)
    
    end subroutine Neighbours_construct
    
    subroutine Neighbours_destroy(this)
    
        class(Neighbours), intent(inout) :: this
        
        if (allocated(this%cell_neighs)) then
            deallocate(this%cell_neighs)
        end if
        
        call this%dealloc_cells()
        
    end subroutine Neighbours_destroy

    ! Linked-list allocation
    
    subroutine Neighbours_alloc_cells(this)
    
        class(Neighbours), intent(inout) :: this
    
        integer :: iCell, nCells
        
        nCells = product(this%NtotalNeighCell_dim)

        allocate(this%cellsBegin(nCells))
        allocate(this%cells(nCells))
        allocate(this%cellsNext(nCells))
        
        do iCell = 1, nCells

            allocate(this%cellsBegin(iCell)%particle)
            this%cells(iCell)%particle => this%cellsBegin(iCell)%particle
            this%cells(iCell)%particle%iCol = 0
            
            allocate(this%cellsNext(iCell)%particle)
            this%cellsNext(iCell)%particle%iCol = 0            
            this%cells(iCell)%particle%next => this%cellsNext(iCell)%particle
            this%cells(iCell)%particle => this%cellsNext(iCell)%particle     
    
        end do
        
    end subroutine Neighbours_alloc_cells
    
    ! Linked-list deallocation
    
    recursive subroutine free_link(current)

        type(Link), pointer :: current
        
        if (associated(current%next)) then
            call free_link(current%next)
        end if
        deallocate(current)
        
    end subroutine free_link
    
    subroutine Neighbours_dealloc_cells(this)
    
        class(Neighbours), intent(inout) :: this
    
        integer :: iCell
        integer :: nCells
    
        nCells = product(this%NtotalNeighCell_dim)
        do iCell = 1, nCells
            if (associated(this%cellsBegin(iCell)%particle)) then
                call free_link(this%cellsBegin(iCell)%particle)
            end if
        end do
    
    end subroutine Neighbours_dealloc_cells
    
    ! Neighbours cells size check
    
    subroutine Neighbours_check_cellsSize(this, rCut)
    
        class(Neighbours), intent(in) :: this
        real(DP), intent(in) :: rCut
        
        integer :: iDim
        
        do iDim = 1, Ndim
        
            if (this%cell_Lsize(iDim) < rCut .and. this%cell_Lsize(iDim) /= Lsize(iDim)/3._DP) then
                write(error_unit, *) "Warning : big rCut in the dimension", iDim, ":"
                write(error_unit, *) this%cell_Lsize(iDim), "<", rCut
            end if
            
            if (Lsize(iDim)/2._DP*sqrt(3._DP) < rCut) then
                write(error_unit, *) "rCut too large in the dimension", iDim, ":"
                write(error_unit, *) Lsize(iDim)/2._DP*sqrt(3._DP), "<", rCut
                stop
            end if
            
            if (this%NtotalNeighCell_dim(iDim) < NnearNeighCell_dim(iDim)) then
                write(error_unit, *) "Too few cells in the dimension", iDim, ":"
                write(error_unit, *) this%NtotalNeighCell_dim(iDim), "<", NnearNeighCell_dim(iDim)
                stop
            end if
            
            if (modulo(Lsize(iDim), this%cell_Lsize(iDim)) /= 0) then
                write(error_unit, *) "Cell size is not a divisor of the system size"
                write(error_unit, *) "in the dimension", iDim, ":"
                write(error_unit, *) "Lsize", Lsize(iDim)
                write(error_unit, *) "cell_Lsize", this%cell_Lsize(iDim)
                write(error_unit, *) "modulo(Lsize, cell_Lsize) = ", &
                                      modulo(Lsize(iDim), this%cell_Lsize(iDim))
                stop
            end if
            
        end do
        
    end subroutine Neighbours_check_cellsSize
    
    ! Assignment : particle -> cell
    
    pure function Neighbours_position_to_cell(this, xCol) result(position_to_cell)
    
        class(Neighbours), intent(in) :: this
        real(DP), dimension(:), intent(in) :: xCol
        integer :: position_to_cell
        
        integer, dimension(Ndim) :: cell_coord        
    
        cell_coord(:) = int(xCol(:)/this%cell_Lsize(:)) + 1
        position_to_cell = cell_coord(1) + this%NtotalNeighCell_dim(1)*(cell_coord(2)-1) + &
                           this%NtotalNeighCell_dim(1)*this%NtotalNeighCell_dim(2)*(cell_coord(3)-1)
    
    end function Neighbours_position_to_cell
    
    subroutine Neighbours_all_cols_to_cells(this, Ncol, X)
    
        class(Neighbours), intent(inout) :: this
        integer, intent(in) :: Ncol
        real(DP), dimension(:, :), intent(in) :: X
    
        integer :: iCol
        integer :: iCell, nCells
        
        nCells = product(this%NtotalNeighCell_dim)
    
        do iCol = 1, Ncol
    
            iCell = this%position_to_cell(X(:,iCol))
            this%cells(iCell)%particle%iCol = iCol
            
            allocate(this%cellsNext(iCell)%particle)
            this%cellsNext(iCell)%particle%iCol = 0
            this%cells(iCell)%particle%next => this%cellsNext(iCell)%particle
            this%cells(iCell)%particle => this%cellsNext(iCell)%particle
            
        end do
        
        do iCell = 1, nCells
            
            this%cells(iCell)%particle%next => null()
            
        end do
        
    end subroutine Neighbours_all_cols_to_cells
    
    ! Neighbours cells update
    
    subroutine Neighbours_remove_col_from_cell(this, iCol, iCellOld)
    
        class(Neighbours), intent(inout) :: this
    
        integer, intent(in) :: iCol, iCellOld
        
        type(Link), pointer :: current => null()
        type(Link), pointer :: next => null(), previous => null()
    
        previous => this%cellsBegin(iCellOld)%particle
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
            
    end subroutine Neighbours_remove_col_from_cell   
    
    subroutine Neighbours_add_col_to_cell(this, iCol, iCellNew)
    
        class(Neighbours), intent(inout) :: this
    
        integer, intent(in) :: iCol, iCellNew
    
        type(Link), pointer :: new => null()
        type(Link), pointer :: next => null(), previous => null()           
          
        
        previous => this%cellsBegin(iCellNew)%particle
        
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
            
    end  subroutine Neighbours_add_col_to_cell
    
    ! Neighbour cells ------------------------------------------------------------------------------

    pure function Neighbours_cell_coord_to_ind(this, coord) result(cell_coord_to_ind)
        
        class(Neighbours), intent(in) :: this
        integer, dimension(:), intent(in) :: coord
        integer :: cell_coord_to_ind
        
        cell_coord_to_ind = coord(1) + this%NtotalNeighCell_dim(1)*(coord(2)-1) + &
                            this%NtotalNeighCell_dim(1)*this%NtotalNeighCell_dim(2)*(coord(3)-1)
    
    end function Neighbours_cell_coord_to_ind
    
    pure function cell_neigh_coord_to_ind(neigh_coord)
    
        integer, dimension(:), intent(in) :: neigh_coord        
        integer :: cell_neigh_coord_to_ind
        
        cell_neigh_coord_to_ind = neigh_coord(1) + NnearNeighCell_dim(1)*(neigh_coord(2)-1) + &
                                  NnearNeighCell_dim(1)*NnearNeighCell_dim(2)*(neigh_coord(3)-1)
    
    end function cell_neigh_coord_to_ind
    
    pure function Neighbours_cell_period(this, coord) result(cell_period)
    
        class(Neighbours), intent(in) :: this    
        integer, dimension(:), intent(in) :: coord        
        integer, dimension(Ndim) :: cell_period
        
        cell_period(:) = coord(:)
        
        where (cell_period(:) < 1)
            cell_period(:) = cell_period(:) + this%NtotalNeighCell_dim(:)
        end where
        
        where (cell_period(:) > this%NtotalNeighCell_dim(:))
            cell_period(:) = cell_period(:) - this%NtotalNeighCell_dim(:)
        end where
    
    end function Neighbours_cell_period
    
    subroutine Neighbours_cell_neighs_init(this)
    
        class(Neighbours), intent(inout) :: this 
    
        integer :: i, j, k, ind
        integer :: neigh_i, neigh_j, neigh_k, neigh_ind
        integer, dimension(Ndim) :: coord, neigh_coord
        
        do i = 1, this%NtotalNeighCell_dim(1)
        do j = 1, this%NtotalNeighCell_dim(2)
        do k = 1, this%NtotalNeighCell_dim(3)
            
            ind = this%cell_coord_to_ind([i, j, k])

            do neigh_i = 1, NnearNeighCell_dim(1)
            do neigh_j = 1, NnearNeighCell_dim(2)
            do neigh_k = 1, NnearNeighCell_dim(3)
            
                neigh_coord(:) = [neigh_i, neigh_j, neigh_k]                
                neigh_ind = cell_neigh_coord_to_ind(neigh_coord(:))          
                neigh_coord(:) = neigh_coord(:) - NnearNeighCell_dim(:) + 1
                    ! with respect to the center (?) [i, j, k]
                
                coord(:) = [i, j, k] + neigh_coord(:)
                
                this%cell_neighs(neigh_ind, ind) = this%cell_coord_to_ind(this%cell_period(coord(:)))
                    
            end do
            end do
            end do
        
        end do
        end do
        end do
            
    end subroutine Neighbours_cell_neighs_init

end module class_neighbourCells
