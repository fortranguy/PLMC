module mod_neighbours

use data_constants
use data_cell
use data_particles
use data_potentiel
use data_neighbours

implicit none

contains
    
    
    
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
        
        do i = 1, cell_coordMax(1)
        do j = 1, cell_coordMax(2)
        do k = 1, cell_coordMax(3)
            
            ind = cell_coord_to_ind([i, j, k])

            do neigh_i = 1, cell_neigh_coordMax(1)
            do neigh_j = 1, cell_neigh_coordMax(2)
            do neigh_k = 1, cell_neigh_coordMax(3)
            
                neigh_coord(:) = [neigh_i, neigh_j, neigh_k]                
                neigh_ind = cell_neigh_coord_to_ind(neigh_coord(:))          
                neigh_coord(:) = neigh_coord(:) - cell_neigh_coordMax(:) + 1
                    ! Par rapport au centre (i, j, k)
                
                coord(:) = [i, j, k] + neigh_coord(:)
                
                cell_neighs(neigh_ind, ind) = cell_coord_to_ind( cell_period(&
                    coord(:)) )
                    
            end do
            end do
            end do
        
        end do
        end do
        end do
            
    end subroutine ini_cell_neighs

end module mod_neighbours
