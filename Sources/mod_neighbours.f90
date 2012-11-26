module mod_neighbours

use data_cell
use data_neighbours
use data_particles

implicit none

contains
    
    ! Assignation : particule -> cellule
    
    function col_to_cell(iCol)
    
        integer, intent(in) :: iCol
        integer, dimension(Dim) :: cell_coord
        integer :: col_to_cell
    
        cell_coord = int( x(:, iCol)/cell_Lsize(:) ) + 1
        col_to_cell = cell_coord(1) + cell_coordMax(1) * (cell_coord(2)-1) + &
            cell_coordMax(1)*cell_coordMax(2)*(cell_coord(3)-1)
    
    end function col_to_cell
    
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
    
    ! Repérage d'une particule dans une cellule :
    
    subroutine locate_cell_iCol(iCol, iCell, cell_iCol)
    
        integer, intent(in) :: iCol
        integer, intent(out) :: iCell, cell_iCol

        type(Particle), pointer :: courant => null(), suivant => null()
        
        iCell = col_to_cell(iCol)
        
        courant => cellsBegin(iCell)%particle
        cell_iCol = 0
        
        do
        
            suivant => courant%next
        
            if ( courant%iCol == iCol ) then
            
                return                
                
            end if
            
            courant => suivant            
            
            cell_iCol = cell_iCol + 1
            
            if (.not. associated(suivant%next)) exit
        
        end do
    
    end subroutine locate_cell_iCol
    
    ! Mise à jour des TV
    
    subroutine update_cell(iCol, iCellBefore)
    
        integer, intent(in) :: iCol, iCellBefore
        
        integer :: iCellAfter
    
        type(Particle), pointer :: courant => null(), nouveau => null()
        type(Particle), pointer :: suivant => null(), precedent => null()
        
         iCellAfter = col_to_cell(iCol)
        
        if ( iCellBefore /= iCellAfter ) then
            
            ! remove
        
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
        
                courant => suivant
                
                if (.not. associated(suivant%next)) exit
            
            end do
            
            ! add        
            
            precedent => cellsBegin(iCellAfter)%particle
            
            do
            
                suivant => precedent%next
                
                if (suivant%iCol == 0) then
                
                    allocate(nouveau)
                    nouveau%next => precedent%next
                    precedent%next => nouveau
                    nouveau%iCol = iCol
                    
                end if
                
                precedent => suivant
                
                if (.not. associated(suivant%next)) exit
                
            end do
        
        end if
                   
   
    end subroutine update_cell
    
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
    
    subroutine cell_neighs_ini()
    
        integer :: i, j, k, ind
        integer :: neigh_i, neigh_j, neigh_k, neigh_ind, neigh_center_ind
        integer, dimension(dim) :: coord, neigh_coord
        integer, dimension(product(cell_coordMax)) :: local_neighs
        
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
                
                if ( neigh_coord(1)==0 .and. neigh_coord(2)==0 .and. &
                    neigh_coord(3)==0 ) then            
                     neigh_center_ind = neigh_ind
                end if
                
                coord(:) = [i, j, k] + neigh_coord(:)
                
                local_neighs(neigh_ind) = cell_coord_to_ind( cell_period(&
                    coord(:)) )
                    
            end do
            end do
            end do
 
            do neigh_ind = 1, neigh_center_ind-1
                cell_neighs(ind, neigh_ind) = local_neighs(neigh_ind)
            end do
            
            do neigh_ind = neigh_center_ind+1, cell_neigh_iMax*cell_neigh_jMax*&
                cell_neigh_kMax
                cell_neighs(ind, neigh_ind-1) = local_neighs(neigh_ind)
            end do
        
        end do
        end do
        end do
            
    end subroutine cell_neighs_ini

end module mod_neighbours
