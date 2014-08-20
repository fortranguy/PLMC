integer, parameter :: NnearCell_layer = NnearCell_dim(1) * NnearCell_dim(2)
integer, parameter :: NnearCell = NnearCell_layer * NnearCell_dim(3)

class HS

this%Volume = Lsize(1)*Lsize(2)*(Height - this%sigma)
density = Ncol_avg / this%Volume

        iNearCell_min = 1
        iNearCell_max = NnearCell
        NtotalCell_layer = this%sameCells%get_NtotalCell_layer()
        
        if (iTotalCell <= NtotalCell_layer) then
            iNearCell_min = NnearCell_layer + 1
        else if ((this%sameCells%get_NtotalCell()-iTotalCell) < NtotalCell_layer) then
            iNearCell_max = (Ndim-1) * NnearCell_layer
        end if
    
        do iNearCell = iNearCell_min, iNearCell_max

class mix

        iNearCell_min = 1
        iNearCell_max = NnearCell
        NtotalCell_layer = neighCells%get_NtotalCell_layer()
        
        if (iTotalCell <= NtotalCell_layer) then
            iNearCell_min = NnearCell_layer + 1
        else if ((neighCells%get_NtotalCell()-iTotalCell) < NtotalCell_layer) then
            iNearCell_max = (Ndim-1) * NnearCell_layer
        end if
    
        do iNearCell = iNearCell_min, iNearCell_max
        
        
class neigh

integer :: NtotalCell, NtotalCell_layer
        procedure :: get_NtotalCell => NeighbourCells_get_NtotalCell ?
        procedure :: get_NtotalCell_layer => NeighbourCells_get_NtotalCell_layer ?
this%NtotalCell_layer = this%NtotalCell_dim(1)  * this%NtotalCell_dim(2)

pure subroutine NeighbourCells_init_near_among_total(this)
    
        class(NeighbourCells), intent(inout) :: this
    
        integer :: iTotalCell, jTotalCell, kTotalCell, totalCell_index
        integer :: iNearCell, jNearCell, kNearCell, nearCell_index
        integer, dimension(Ndim) :: totalCell_coord, nearCell_coord
        integer :: kNearCell_min, kNearCell_max
        
        do kTotalCell = 1, this%NtotalCell_dim(3)

            kNearCell_min = 1
            kNearCell_max = NnearCell_dim(3)
        
            if (kTotalCell == 1) then
                kNearCell_min = NnearCell_dim(3) - 1
                kNearCell_max = NnearCell_dim(3)
            else if (kTotalCell == this%NtotalCell_dim(3)) then
                kNearCell_min = 1
                kNearCell_max = 2
            end if
        
            do jTotalCell = 1, this%NtotalCell_dim(2)
            do iTotalCell = 1, this%NtotalCell_dim(1)

                totalCell_index = index_from_coord([iTotalCell, jTotalCell, kTotalCell], &
                                                   this%NtotalCell_dim)

                do kNearCell = kNearCell_min, kNearCell_max
                do jNearCell = 1, NnearCell_dim(2)
                do iNearCell = 1, NnearCell_dim(1)
                
                    nearCell_coord(:) = [iNearCell, jNearCell, kNearCell]
                    nearCell_index = index_from_coord(nearCell_coord, NnearCell_dim)
                    nearCell_coord(:) = nearCell_coord(:) - NnearCell_dim(:) + 1
                        ! with respect to the center (?) [iTotalCell, jTotalCell, kTotalCell]
                    
                    totalCell_coord(:) = [iTotalCell, jTotalCell, kTotalCell] + nearCell_coord(:)
                    totalCell_coord(:) = coord_PBC(totalCell_coord, this%NtotalCell_dim)
  
                    this%near_among_total(nearCell_index, totalCell_index) = &
                        index_from_coord(totalCell_coord, this%NtotalCell_dim)
                        
                end do
                end do
                end do
            
            end do
            end do
            
        end do
            
    end subroutine NeighbourCells_init_near_among_total
    
module algo

C

if (xNew(3)+this%get_sigma()/2._DP > Height .or. xNew(3)-this%get_sigma()/2._DP < 0._DP) then
            this_obs%move_Nreject = this_obs%move_Nreject + 1
            return
        end if

call random_number(xRand)
            xTest(1:2) = Lsize(1:2) * xRand(1:2)
            xTest(3) = (Height - this%get_sigma()) * xRand(3) + this%get_sigma()/2._DP

GC
xNew(1:2) = Lsize(1:2) * xRand(1:2)
        xNew(3) = (Height - this%get_sigma()) * xRand(3) + this%get_sigma()/2._DP
        this%Volume (GC)
