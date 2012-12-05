!***********************************************************************
!* MODULE: Neighbours                                                  *
!***********************************************************************
module data_neighbours

use data_potentiel

    ! Type
    type Particle
        integer :: iCol
        type(Particle), pointer :: next => null()
    end type Particle
    
    type ContainerParticle
        type(Particle), pointer :: particle => null()
    end type ContainerParticle
    
    ! Cellules
    real(DP), parameter :: cell_Lsize1 = rcut11, cell_Lsize2 = rcut11, &
        cell_Lsize3 = rcut11
    real(DP), dimension(Dim), parameter :: cell_Lsize = &
        [cell_Lsize1, cell_Lsize2, cell_Lsize3]
    integer, parameter :: cell_iMax = int(Lsize1/cell_Lsize1), &
    cell_jMax = int(Lsize2/cell_Lsize2), cell_kMax = int(Lsize3/cell_Lsize3)
    integer, dimension(dim), parameter :: cell_coordMax = &
        [cell_iMax, cell_jMax, cell_kMax]
    type(ContainerParticle), allocatable, dimension(:) :: cells, cellsNext
    type(ContainerParticle), allocatable, dimension(:), protected :: cellsBegin

    ! Voisins
    integer, parameter :: cell_neigh_iMax = 3, &
        cell_neigh_jMax = cell_neigh_iMax, &
        cell_neigh_kMax = cell_neigh_iMax ! évident
    integer, dimension(dim), parameter :: cell_neigh_coordMax = &
        [cell_neigh_iMax, cell_neigh_jMax, cell_neigh_kMax]
    integer, parameter :: cell_neighs_nb = &
        cell_neigh_iMax*cell_neigh_jMax*cell_neigh_kMax ! inclus soi-même
    integer, dimension(cell_neighs_nb, cell_iMax*cell_jMax*cell_kMax) :: &
        cell_neighs
    
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
    
        type(Particle), pointer :: courant
        
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
    
end module data_neighbours
!***********************************************************************
