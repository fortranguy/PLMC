    !> Particle switch
    
    subroutine before_switch_energy(this, iCol, other, mix, indicesOld, EpotsOld)
        
        class(HardSpheres), intent(in) :: this, other
        class(MixingPotential), intent(in) :: mix
        integer, intent(in) :: iCol
        integer, dimension(:), intent(out) :: indicesOld
        real(DP), dimension(:), intent(out) :: EpotsOld
        
        real(DP), dimension(Ndim) :: xOld, mCol
        logical :: overlap
        
        xOld(:) = this%positions(:, iCol)
        
        indicesOld(1) = this%sameCells%index_from_position(xOld)
        select type (this)
            type is (DipolarSpheres)
                mCol(:) = this%orientations(:, iCol)
                EpotsOld(1) = this%Epot_real_solo(iCol, xOld, mCol) - &
                              this%deltaEpot_reci_exchange(xOld, -mCol)
            type is (InteractingSpheres)
                call this%Epot_neighCells(iCol, xOld, indicesOld(1), overlap, EpotsOld(1))
            type is (HardSpheres)
                EpotsOld(1) = 0._DP
        end select
        
        indicesOld(2) = this%mixCells%index_from_position(xOld)
        call mix%Epot_neighCells(0, xOld, indicesOld(2), this%mixCells, other%positions, overlap, &
                                 EpotsOld(2))
        
    end subroutine before_switch_energy
    
    subroutine after_switch_energy(this, this_iCol, other, other_iCol, mix, iCellNew, mix_iCellNew, &
                                   overlap, EpotNew, mix_EpotNew)
        
        class(HardSpheres), intent(in) :: this, other
        class(MixingPotential), intent(in) :: mix
        integer, intent(in) :: this_iCol, other_iCol
        integer, intent(ou) :: iCellNew, mix_iCellNew
        logical, intent(out) :: overlap
        real(DP), intent(out) :: EpotNew, mix_EpotNew
        
        real(DP), dimension(Ndim) :: xNew, mCol
        
        xNew(:) = other%positions(:, other_iCol)
        
        if (this%get_Ncol() >= other%get_Ncol()) then ! optimisation : more chance to overlap
            iCellNew = this%sameCells%index_from_position(xNew)
            call this%Epot_neighCells(this_iCol, xNew, iCellNew, overlap, EpotNew)
        else
            mix_iCellNew = this%mixCells%index_from_position(xNew)
            call mix%Epot_neighCells(other_iCol, xNew, mix_iCellNew, this%mixCells, other%positions, &
                                     overlap, mix_EpotNew)
        end if
        
        if (.not. overlap) then
        
            if (this%get_Ncol() >= other%get_Ncol()) then
                mix_iCellNew = this%mixCells%index_from_position(xNew)
                call mix%Epot_neighCells(other_iCol, xNew, mix_iCellNew, this%mixCells, &
                                         other%positions, overlap, mix_EpotNew)
            else
                iCellNew = this%sameCells%index_from_position(xNew)
                call this%Epot_neighCells(this_iCol, xNew, iCellNew, overlap, EpotNew)
            end if
            
            if (.not. overlap) then
            
                select type (this)
                    type is (DipolarSpheres)
                        mCol(:) = this%orientations(:, iCol)
                        EpotNew = this%Epot_real_solo(this_iCol, xNew, mCol) + &
                                       this%deltaEpot_reci_exchange(xNew, +mCol)
                end select
            
            end if
                
        end if
    
    end subroutine after_switch_energy
    
    subroutine after_switch_update(this, this_iCol, this_iCellOld, mix_iCellOld, this_iCellNew, &
                                   mix_iCellNew, other, other_iCol, mix)
        
        class(HardSpheres), intent(inout) :: this, other
        class(MixingPotential), intent(in) :: mix
        integer, intent(in) :: this_iOld, other_iCol
        integer, intent(in) :: this_iCellOld, mix_iCellOld
        integer, intent(in) :: this_iCellNew, mix_iCellNew
        
        real(DP), dimension(Ndim) :: xOld, xNew, mCol
        
        select type (this)
            type is (DipolarSpheres)
                xOld(:) = this%positions(:, this_iCol)
                mCol(:) = this%orientations(:, this_iCol)
                call this%reci_update_structure_exchange(xOld, -mCol)
        end select
    
        call this%sameCells%remove_col_from_cell(this_iCol, this_iCellOld)
        call other%mixCells%remove_col_from_cell(this_iCol, mix_iCellOld)
        
        this%positions(:, this_iCol) = other%positions(:, other_iCol)
        select type (this)
            type is (DipolarSpheres)
                xNew(:) = this%positions(:, this_iCol)
                mCol(:) = this%orientations(:, this_iCol)
                call this%reci_update_structure_exchange(xNew, +mCol)
        end select
        
        call this%sameCells%add_col_to_cell(this_iCol, this_iCellNew)
        call other%mixCells%add_col_to_cell(this_iCol, mix_iCellNew)
        
    end subroutine after_switch_update
    
    subroutine switch(type1, type1_obs, type2, type2_obs, mix, mix_Epot)
    
        class(HardSpheres), intent(inout) :: type1, type2
        class(Observables), intent(inout) :: type1_obs, type2_obs
        class(MixingPotential), intent(in) :: mix
        real(DP), intent(inout) :: mix_Epot
        
        real(DP) :: random
        integer :: type1_iCol, type2_iCol
        integer, dimension(2) :: type1_indicesOld, type2_indicesOld ! (1): same, (2): other
        integer, dimension(2) :: type1_indicesNew, type2_indicesNew
        logical :: overlap
        real(DP) :: deltaEpot, type1_deltaEpot, type2_deltaEpot
        real(DP) :: type1_mix_deltaEpot, type2_mix_deltaEpot
        real(DP), dimension(2) :: type1_EpotsOld, type2_EpotsOld ! (1): same, (2): mix
        real(DP), dimension(2) :: type1_EpotsNew, type2_EpotsNew
        
        type1_obs%Nswitch = type1_obs%Nswitch + 1
        
        if (type1%get_Ncol()==0 .or. type2%get_Ncol()==0) then
            type1_obs%switch_Nreject = type1_obs%switch_Nreject + 1
            return
        end if
        
        ! Old : before switch
        call random_number(random)
        type1_iCol = int(random*type1%get_Ncol()) + 1
        call before_switch_energy(type1, type1_iCol, type2, mix, type1_indicesOld, type1_EpotsOld)
        
        call random_number(random)
        type2_iCol = int(random*type2%get_Ncol()) + 1
        call before_switch_energy(type2, type2_iCol, type1, mix, type2_indicesOld, type2_EpotsOld)
        
        ! New : after switch
        call after_switch_energy(type1, type1_iCol, type2, type2_iCol, mix, type1_indicesNew, &
                                 overlap, type1_EpotsNew)
        
        if (.not. overlap) then
        
            call after_switch_energy(type2, type2_iCol, type1, type1_iCol, mix, type2_indicesNew, &
                                     overlap, type2_EpotNew, type2_EpotsNew)

            type1_deltaEpot = type1_EpotsNew(1) - type1_EpotsOld(1)
            type1_mix_deltaEpot = type1_EpotsNew(2) - type1_EpotsOld(2)
            type2_deltaEpot = type2_EpotsNew(1) - type2_EpotsOld(1)
            type2_mix_deltaEpot = type2_EpotsNew(2) - type2_EpotsOld(2)
            deltaEpot = type1_deltaEpot + type1_mix_deltaEpot + type2_deltaEpot + type2_mix_deltaEpot
            
            call random_number(random)
            if (random < exp(-deltaEpot/Temperature)) then
            
                call after_switch_update(type1, type1_iCol, type1_indicesOld, &
                                         type1_indicesNew, type2, mix)
                call after_switch_update(type2, type2_iCol, type2_indicesOld, &
                                         type2_indicesNew, type1, mix)
                                         
                type1_obs%Epot = type1_obs%Epot + type1_deltaEpot
                type2_obs%Epot = type2_obs%Epot + type2_deltaEpot
                mix_Epot = mix_Epot + type1_mix_deltaEpot + type2_mix_deltaEpot
            
            end if
        
        end if
        
    end subroutine switch
