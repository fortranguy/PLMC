module class_neighbours

use data_cell
use data_particles
use data_potentiel
use data_neighbours
use data_mc
use mod_pbc

implicit none

private

    type Link
    integer :: iCol
    type(Link), pointer :: next => null()
    end type Link
    
    type LinkedList
        type(Link), pointer :: particle => null()
    end type LinkedList

    type, public :: Neighbours
        
        real(DP), dimension(Dim), private :: cell_Lsize
        integer, dimension(Dim), private :: cell_coordMax
        integer, dimension(:, :), allocatable :: cell_neighs
        type(LinkedList), allocatable, dimension(:), private :: cells
        type(LinkedList), allocatable, dimension(:), private :: cellsNext
        type(LinkedList), allocatable, dimension(:), private :: cellsBegin
        
    contains
    
        procedure :: alloc_Cells => component_alloc_Cells
        procedure :: dealloc_Cells => component_dealloc_Cells
        procedure :: check_CellsSize => component_check_CellsSize
        procedure :: position_to_cell => component_position_to_cell
        procedure :: all_col_to_cell => component_all_col_to_cell
        procedure :: remove_cell_col => component_remove_cell_col
        procedure :: add_cell_col => component_add_cell_col
        procedure :: cell_coord_to_ind => component_cell_coord_to_ind
        procedure :: cell_period => component_cell_period
        procedure :: ini_cell_neighs => component_ini_cell_neighs
        
    end type Neighbours
    
contains

    ! Linked-list allocation
    
    subroutine component_alloc_Cells(this)
    
        class(Component), intent(inout) :: this
    
        integer :: iCell, nCells
        
        nCells = product(this%cell_coordMax)

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
        
    end subroutine component_alloc_Cells
    
    ! Linked-list deallocation
    
    recursive subroutine libere_chaine(courant)

        type(Link), pointer :: courant
        
        if (associated(courant%next)) then
            call libere_chaine(courant%next)
        end if
        deallocate(courant)
        
    end subroutine libere_chaine
    
    subroutine component_dealloc_Cells(this)
    
        class(Component), intent(inout) :: this
    
        integer :: iCell
        integer :: nCells
    
        nCells = this%cell_coordMax(1) * this%cell_coordMax(2) * &
            this%cell_coordMax(3)
        do iCell = 1, nCells
            if (associated(this%cellsBegin(iCell)%particle)) then
                call libere_chaine(this%cellsBegin(iCell)%particle)
            end if
        end do
    
    end subroutine component_dealloc_Cells
    
    ! Neighbours cells size check
    
    subroutine component_check_CellsSize(this)
    
        class(Component), intent(in) :: this
        
        integer :: iDir
        
        do iDir = 1, Dim
        
            if (this%cell_Lsize(iDir) < this%rcut) then
                write(*, *) "Cellule trop petite dans la direction", iDir, ":"
                write(*, *) this%cell_Lsize(iDir), "<", this%rcut
                stop
            end if
            
            if (this%cell_coordMax(iDir) < cell_neigh_coordMax(iDir)) then
                write(*, *) "Trop peu de cellules dans la direction", iDir, ":"
                write(*, *) this%cell_coordMax(iDir), "<",&
                    cell_neigh_coordMax(iDir)
                stop
            end if
            
        end do
        
    end subroutine component_check_CellsSize
    
    ! Assignment : particle -> cell
    
    function component_position_to_cell(this, xCol) &
        result(position_to_cell)
    
        class(Component), intent(in) :: this
        real(DP), dimension(Dim), intent(in) :: xCol
        
        integer, dimension(Dim) :: cell_coord
        integer :: position_to_cell
    
        cell_coord(:) = int( xCol(:)/this%cell_Lsize(:) ) + 1
        position_to_cell = cell_coord(1) + this%cell_coordMax(1) * &
            (cell_coord(2)-1) + this%cell_coordMax(1) * &
            this%cell_coordMax(2) * (cell_coord(3)-1)
    
    end function component_position_to_cell
    
    subroutine component_all_col_to_cell(this)
    
        class(Component), intent(inout) :: this
    
        integer :: iCol
        integer :: iCell, nCells
        
        nCells = this%cell_coordMax(1) * this%cell_coordMax(2) * &
            this%cell_coordMax(3)
    
        do iCol = 1, this%Ncol
    
            iCell = this%position_to_cell(this%X(:,iCol))
            this%cells(iCell)%particle%iCol = iCol
            
            allocate(this%cellsNext(iCell)%particle)
            this%cellsNext(iCell)%particle%iCol = 0
            this%cells(iCell)%particle%next => &
                this%cellsNext(iCell)%particle
            this%cells(iCell)%particle => this%cellsNext(iCell)%particle
            
        end do
        
        do iCell = 1, nCells
            
            this%cells(iCell)%particle%next => null()
            
        end do
        
    end subroutine component_all_col_to_cell
    
    ! Neighbours cells update
    
    subroutine component_remove_cell_col(this, iCol, iCellBefore)
    
        class(Component), intent(inout) :: this
    
        integer, intent(in) :: iCol, iCellBefore
        
        type(Link), pointer :: courant => null()
        type(Link), pointer :: suivant => null(), precedent => null()
    
        precedent => this%cellsBegin(iCellBefore)%particle
        courant => precedent%next
        
        do
        
            suivant => courant%next
        
            if ( courant%iCol == iCol ) then
            
                precedent%next => courant%next
                deallocate(courant)
                courant => suivant
                exit
                
            else
            
                precedent => courant
                
            end if
            
            courant => suivant
        
        end do
            
    end subroutine component_remove_cell_col   
    
    subroutine component_add_cell_col(this, iCol, iCellAfter)
    
        class(Component), intent(inout) :: this
    
        integer, intent(in) :: iCol, iCellAfter
    
        type(Link), pointer :: nouveau => null()
        type(Link), pointer :: suivant => null(), precedent => null()           
          
        
        precedent => this%cellsBegin(iCellAfter)%particle
        
        do
        
            suivant => precedent%next
            
            if (.not. associated(suivant%next)) then
            
                allocate(nouveau)
                nouveau%next => precedent%next
                precedent%next => nouveau
                nouveau%iCol = iCol
                exit
                
            end if
            
            precedent => suivant
            
        end do
            
    end  subroutine component_add_cell_col
    
! -----------------------------------------------------------------------------
! Neighbours :
! -----------------------------------------------------------------------------
    
    function component_cell_coord_to_ind(this, coord) result(cell_coord_to_ind)
        
        class(Component), intent(in) :: this    
        integer, dimension(Dim), intent(in) :: coord
        
        integer :: cell_coord_to_ind
        
        cell_coord_to_ind = coord(1) + this%cell_coordMax(1)*(coord(2)-1) + &
            this%cell_coordMax(1)*this%cell_coordMax(2)*(coord(3)-1)
    
    end function component_cell_coord_to_ind
    
    function cell_neigh_coord_to_ind(neigh_coord)
    
        integer, dimension(Dim), intent(in) :: neigh_coord
        
        integer :: cell_neigh_coord_to_ind
        
        cell_neigh_coord_to_ind = neigh_coord(1) + &
            cell_neigh_coordMax(1) * (neigh_coord(2)-1) &
            + cell_neigh_coordMax(1) * cell_neigh_coordMax(2) * &
            (neigh_coord(3)-1)
    
    end function cell_neigh_coord_to_ind
    
    function component_cell_period(this, coord) result(cell_period)
    
        class(Component), intent(in) :: this    
        integer, dimension(Dim), intent(in) :: coord
        
        integer, dimension(Dim) :: cell_period
        
        cell_period(:) = coord(:)
        
        where (cell_period(:) < 1)
            cell_period(:) = cell_period(:) + this%cell_coordMax(:)
        end where
        
        where (cell_period(:) > this%cell_coordMax(:))
            cell_period(:) = cell_period(:) - this%cell_coordMax(:)
        end where
    
    end function component_cell_period
    
    subroutine component_ini_cell_neighs(this)
    
        class(Component), intent(inout) :: this 
    
        integer :: i, j, k, ind
        integer :: neigh_i, neigh_j, neigh_k, neigh_ind
        integer, dimension(Dim) :: coord, neigh_coord
        
        do i = 1, this%cell_coordMax(1)
        do j = 1, this%cell_coordMax(2)
        do k = 1, this%cell_coordMax(3)
            
            ind = this%cell_coord_to_ind([i, j, k])

            do neigh_i = 1, cell_neigh_coordMax(1)
            do neigh_j = 1, cell_neigh_coordMax(2)
            do neigh_k = 1, cell_neigh_coordMax(3)
            
                neigh_coord(:) = [neigh_i, neigh_j, neigh_k]                
                neigh_ind = cell_neigh_coord_to_ind(neigh_coord(:))          
                neigh_coord(:) = neigh_coord(:) - cell_neigh_coordMax(:) + 1
                    ! Par rapport au centre [i, j, k]
                
                coord(:) = [i, j, k] + neigh_coord(:)
                
                this%cell_neighs(neigh_ind, ind) = &
                    this%cell_coord_to_ind( this%cell_period(coord(:)) )
                    
            end do
            end do
            end do
        
        end do
        end do
        end do
            
    end subroutine component_ini_cell_neighs

end module class_neighbours
