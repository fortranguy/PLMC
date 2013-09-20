    subroutine HardSpheres_move(this, this_obs, other, mix, mix_Epot)
    
        class(HardSpheres), intent(inout) :: this
        class(Observables) :: this_obs
        class(Spheres), intent(inout) :: other
        class(MixingPotential), intent(in) :: mix        
        real(DP), intent(inout) :: mix_Epot
        
        real(DP) :: random
        integer :: iOld
        real(DP), dimension(Ndim) :: xOld, xRand, xNew
        logical :: overlap
        integer :: this_iCellOld, this_iCellNew
        integer :: mix_iCellOld, mix_iCellNew
        real(DP) :: mix_deltaEpot
        real(DP) :: mix_EpotNew, mix_EpotOld
        
        call random_number(random)
        iOld = int(random*this%Ncol) + 1
        xOld(:) = this%positions(:, iOld)

        call random_number(xRand)
        xNew(:) = xOld(:) + (xRand(:)-0.5_DP)*this%move_delta(:)
        xNew(:) = modulo(xNew(:), Lsize(:))
        
        if (this%Ncol >= other%Ncol) then        
            this_iCellNew = this%sameCells%index_from_position(xNew)
            call this%Epot_neighCells(iOld, xNew, this_iCellNew, overlap)            
        else        
            mix_iCellNew = this%mixCells%index_from_position(xNew)
            call mix%Epot_neighCells(xNew, mix_iCellNew, this%mixCells, other%positions, overlap, &
                                     mix_EpotNew)        
        end if
        
        if (.not. overlap) then
        
            if (this%Ncol >= other%Ncol) then        
                mix_iCellNew = this%mixCells%index_from_position(xNew)
                call mix%Epot_neighCells(xNew, mix_iCellNew, this%mixCells, other%positions, overlap, &
                                         mix_EpotNew)
            else                
                this_iCellNew = this%sameCells%index_from_position(xNew)
                call this%Epot_neighCells(iOld, xNew, this_iCellNew, overlap)
            end if
                        
            if (.not. overlap) then
    
                this_iCellOld = this%sameCells%index_from_position(xOld)
                call this%Epot_neighCells(iOld, xOld, this_iCellOld, overlap)
                    
                mix_iCellOld = this%mixCells%index_from_position(xOld)
                call mix%Epot_neighCells(xOld, mix_iCellOld, this%mixCells, other%positions, overlap, &
                                         mix_EpotOld)
                
                mix_deltaEpot = mix_EpotNew - mix_EpotOld
                
                call random_number(random)
                if (random < exp(-mix_deltaEpot/Temperature)) then
                
                    this%positions(:, iOld) = xNew(:)
                    mix_Epot = mix_Epot + mix_deltaEpot
                    
                    if (this_iCellOld /= this_iCellNew) then                
                        call this%sameCells%remove_col_from_cell(iOld, this_iCellOld)
                        call this%sameCells%add_col_to_cell(iOld, this_iCellNew)
                    end if
                    
                    if (mix_iCellOld /= mix_iCellNew) then                
                        call other%mixCells%remove_col_from_cell(iOld, mix_iCellOld)
                        call other%mixCells%add_col_to_cell(iOld, mix_iCellNew)
                    end if
                    
                else
                    this_obs%move_Nreject = this_obs%move_Nreject + 1
                end if
         
            else
                this_obs%move_Nreject = this_obs%move_Nreject + 1
            end if            
            
        else        
            this_obs%move_Nreject = this_obs%move_Nreject + 1
        end if
    
    end subroutine HardSpheres_move
