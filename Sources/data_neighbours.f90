!***********************************************************************
!* MODULE: Neighbours                                                  *
!***********************************************************************
module data_neighbours

use data_constants
use data_potentiel

    ! Type
    type Link
        integer :: iCol
        type(Link), pointer :: next => null()
    end type Link
    
    type LinkedList
        type(Link), pointer :: particle => null()
    end type LinkedList
    

        
    ! LinkedList
    type(LinkedList), allocatable, dimension(:) :: cells, cellsNext
    type(LinkedList), allocatable, dimension(:), protected :: cellsBegin


    
contains
    
    ! Allocation
    
    subroutine alloc_Cells()
    
        integer :: iCell, nCells
        
        nCells = product(cell_coordMax)

        allocate(cellsBegin(nCells))
        allocate(cells(nCells))
        allocate(cellsNext(nCells))
        
        do iCell = 1, nCells

            allocate(cellsBegin(iCell)%particle)
            cells(iCell)%particle => cellsBegin(iCell)%particle
            cells(iCell)%particle%iCol = 0
            
            allocate(cellsNext(iCell)%particle)
            cellsNext(iCell)%particle%iCol = 0            
            cells(iCell)%particle%next => cellsNext(iCell)%particle
            cells(iCell)%particle => cellsNext(iCell)%particle     
    
        end do
        
    end subroutine alloc_Cells
    
    ! Libération
    
    recursive subroutine libere_chaine(courant)
    
        type(Link), pointer :: courant
        
        if (associated(courant%next)) then
            call libere_chaine(courant%next)
        end if
        deallocate(courant)
        
    end subroutine libere_chaine
    
    subroutine dealloc_Cells()
    
        integer :: iCell
        integer :: nCells = cell_coordMax(1)*cell_coordMax(2)*cell_coordMax(3)
    
        do iCell = 1, nCells
            if (associated(cellsBegin(iCell)%particle)) then
                call libere_chaine(cellsBegin(iCell)%particle)
            end if
        end do
    
    end subroutine dealloc_Cells