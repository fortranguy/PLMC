module mod_neighbours

use data_constants
use data_cell
use data_neighbours
use data_particles
use data_potentiel

implicit none

contains

    ! Vérification de la taille des cellules (voisines)
    
    subroutine check_CellsSize()
        
        integer :: iDir
        
        do iDir = 1, Dim
        
            if (cell_Lsize(iDir) < rcut11) then
                write(*, *) "Cellule trop petite dans la direction", iDir, ":"
                write(*, *) cell_Lsize(iDir), "<", rcut11
                stop
            end if
            
            if (cell_coordMax(iDir) < cell_neigh_coordMax(iDir)) then
                write(*, *) "Trop peu de cellules dans la direction", iDir, ":"
                write(*, *) cell_coordMax(iDir), "<", cell_neigh_coordMax(iDir)
                stop
            end if
            
        end do
        
    end subroutine check_CellsSize
    
    ! Assignation : particule -> cellule
    
    function col_to_cell(iCol)
    
        integer, intent(in) :: iCol
        integer, dimension(Dim) :: cell_coord
        integer :: col_to_cell
    
        cell_coord = int( x(:, iCol)/cell_Lsize(:) ) + 1
        col_to_cell = cell_coord(1) + cell_coordMax(1) * (cell_coord(2)-1) + &
            cell_coordMax(1)*cell_coordMax(2)*(cell_coord(3)-1)
    
    end function col_to_cell
    
    function position_to_cell(pos)
    
        real(DP), dimension(Dim), intent(in) :: pos
        integer, dimension(Dim) :: cell_coord
        integer :: position_to_cell
    
        cell_coord = int( pos(:)/cell_Lsize(:) ) + 1
        position_to_cell = cell_coord(1) + cell_coordMax(1) * (cell_coord(2)-1)&
            + cell_coordMax(1)*cell_coordMax(2)*(cell_coord(3)-1)
    
    end function position_to_cell
    
    subroutine all_col_to_cell()
    
        integer :: iCol
        integer :: iCell, nCells = cell_coordMax(1)*cell_coordMax(2)*&
            cell_coordMax(3)
    
        do iCol = 1, Ncol1
    
            iCell = col_to_cell(iCol)            
            cells(iCell)%particle%iCol = iCol
            
            allocate(cellsNext(iCell)%particle)
            cellsNext(iCell)%particle%iCol = 0
            cells(iCell)%particle%next => cellsNext(iCell)%particle
            cells(iCell)%particle => cellsNext(iCell)%particle
            
        end do
        
        do iCell = 1, nCells
            
            cells(iCell)%particle%next => null()
            
        end do
        
    end subroutine all_col_to_cell
    
    ! Mise à jour des TV
    
    subroutine remove_cell_col(iCol, iCellBefore)
    
        integer, intent(in) :: iCol, iCellBefore
        
        type(Particle), pointer :: courant => null()
        type(Particle), pointer :: suivant => null(), precedent => null()
    
        precedent => cellsBegin(iCellBefore)%particle
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
            
            if (.not. associated(suivant%next)) exit
            
            courant => suivant
        
        end do
            
    end subroutine remove_cell_col   
    
    subroutine add_cell_col(iCol, iCellAfter)
    
        integer, intent(in) :: iCol, iCellAfter
    
        type(Particle), pointer :: nouveau => null()
        type(Particle), pointer :: suivant => null(), precedent => null()           
          
        
        precedent => cellsBegin(iCellAfter)%particle
        
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
            
    end  subroutine add_cell_col
    
! -----------------------------------------------------------------------------
! VOISINS :
! -----------------------------------------------------------------------------
    
    function cell_coord_to_ind(coord)
    
        integer, dimension(dim), intent(in) :: coord
        
        integer :: cell_coord_to_ind
        
        cell_coord_to_ind = coord(1) + cell_coordMax(1)*(coord(2)-1) + &
            cell_coordMax(1)*cell_coordMax(2)*(coord(3)-1)
    
    end function cell_coord_to_ind
    
    function cell_neigh_coord_to_ind(coord) ! changer de nom ?
    
        integer, dimension(dim), intent(in) :: coord
        
        integer :: cell_neigh_coord_to_ind
        
        cell_neigh_coord_to_ind = coord(1) + cell_neigh_coordMax(1)*(coord(2)-1)&
            + cell_neigh_coordMax(1)*cell_neigh_coordMax(2)*(coord(3)-1)
    
    end function cell_neigh_coord_to_ind
    
    function cell_period(coord)
    
        integer, dimension(dim), intent(in) :: coord
        
        integer, dimension(dim) :: cell_period
        
        cell_period(:) = coord(:)
        
        where (cell_period(:) < 1)
            cell_period(:) = cell_period(:) + cell_coordMax(:)
        end where
        
        where (cell_period(:) > cell_coordMax(:))
            cell_period(:) = cell_period(:) - cell_coordMax(:)
        end where
    
    end function cell_period
    
    subroutine ini_cell_neighs()
    
        integer :: i, j, k, ind
        integer :: neigh_i, neigh_j, neigh_k, neigh_ind
        integer, dimension(dim) :: coord, neigh_coord
        
        do i = 1, cell_iMax
        do j = 1, cell_jMax
        do k = 1, cell_kMax
            
            ind = cell_coord_to_ind([i, j, k])

            do neigh_i = 1, cell_neigh_iMax
            do neigh_j = 1, cell_neigh_jMax
            do neigh_k = 1, cell_neigh_kMax
            
                neigh_coord(:) = [neigh_i, neigh_j, neigh_k]                
                neigh_ind = cell_neigh_coord_to_ind(neigh_coord(:))          
                neigh_coord(:) = neigh_coord(:) - cell_neigh_coordMax(:) + 1
                    ! Par rapport au centre (i, j, k)
                
                coord(:) = [i, j, k] + neigh_coord(:)
                
                cell_neighs(ind, neigh_ind) = cell_coord_to_ind( cell_period(&
                    coord(:)) )
                    
            end do
            end do
            end do
        
        end do
        end do
        end do
            
    end subroutine ini_cell_neighs

end module mod_neighbours
