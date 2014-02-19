    !> Particle switch
    
    subroutine before_switch_energy(this, iOld, other, mix, EpotOld, mix_EpotOld)
        
        class(HardSpheres), intent(in) :: this, other
        class(MixingPotential), intent(in) :: mix
        integer, intent(in) :: iOld
        real(DP), intent(out) :: EpotOld, mix_EpotOld
        
        xOld(:) = this%positions(:, iOld)
        
        iCellOld = this%sameCells%index_from_position(xOld)
        select type (this)
            type is (DipolarSpheres)
                mOld(:) = this%orientations(:, iOld)
                EpotOld = this%Epot_real_solo(iOld, xOld, mOld) - &
                          this%deltaEpot_reci_exchange(xOld, -mOld) - &
                          this%Epot_self_solo(mOld) - &
                          this%deltaEpot_bound_exchange(-mOld)
            type is (InteractingSpheres)
                call this%Epot_neighCells(iOld, xOld, iCellOld, overlap, EpotOld)
            type is (HardSpheres)
                EpotOld = 0._DP
        end select
        
        mix_iCellOld = this%mixCells%index_from_position(xOld)
        call mix%Epot_neighCells(xOld, mix_iCellOld, this%mixCells, other%positions, overlap, &
                                 mix_EpotOld)
        
    end subroutine before_switch_energy
    
    subroutine after_switch_energy(this, iNew, other, mix, overlap, EpotNew, mix_EpotNew)
        
        class(HardSpheres), intent(in) :: this, other
        class(MixingPotential), intent(in) :: mix
        integer, intent(in) :: iNew
        logical, intent(out) :: overlap
        real(DP), intent(out) :: EpotNew, mix_EpotNew
        
        xNew(:) = this%positions(:, iNew)
        
        if (this%get_Ncol() >= other%get_Ncol()) then ! optimisation : more chance to overlap
            iCellNew = this%sameCells%index_from_position(xNew)
            call this%Epot_neighCells(iNew, xNew, iCellNew, overlap, EpotNew)
        else
            mix_iCellNew = this%mixCells%index_from_position(xNew)
            call mix%Epot_neighCells(xNew, mix_iCellNew, this%mixCells, other%positions, overlap, &
                                     mix_EpotNew)
        end if
        
        if (.not. overlap) then
        
            if (this%get_Ncol() >= other%get_Ncol()) then
                mix_iCellNew = this%mixCells%index_from_position(xNew)
                call mix%Epot_neighCells(xNew, mix_iCellNew, this%mixCells, other%positions, overlap, &
                                         mix_EpotNew)
            else
                this_iCellNew = this%sameCells%index_from_position(xNew)
                call this%Epot_neighCells(iNew, xNew, this_iCellNew, overlap, this_EpotNew)
            end if
            
            if (.not. overlap) then
            
                select type (this)
                    type is (DipolarSpheres)
                        mNew(:) = random_surface()
                        EpotNew = this%Epot_real_solo(iNew, xNew, mNew) + &
                                       this%deltaEpot_reci_exchange(xNew, +mNew) - &
                                       this%Epot_self_solo(mNew) + &
                                       this%deltaEpot_bound_exchange(+mNew)
                end select
            
                EpotNew = EpotNew + mix_EpotNew
                
        end if
    
    end subroutine after_switch_energy
    
    subroutine after_switch_update(this, iOld, other, mix)
        
        class(HardSpheres), intent(inout) :: this, other
        class(MixingPotential), intent(in) :: mix
        integer, intent(in) :: iOld
        
        select type (this)
            type is (DipolarSpheres)
                xOld(:) = this%positions(:, iOld)
                mOld(:) = this%orientations(:, iOld)
                call this%reci_update_structure_exchange(xOld, -mOld)
                call this%update_totalMoment_exchange(-mOld)
        end select
    
        call this%sameCells%remove_col_from_cell(iOld, this_iCellOld)
        call other%mixCells%remove_col_from_cell(iOld, mix_iCellOld)
        
        this%positions(:, iNew) = other%positions(:, iOld)
        select type (this)
            type is (DipolarSpheres)
                xNew(:) = this%positions(:, iNew)
                mOld(:) = this%orientations(:, iNew)
                call this%reci_update_structure_exchange(xNew, +mOld)
                call this%update_totalMoment_exchange(+mOld)
        end select
        
        call this%sameCells%add_col_to_cell(iOld, this_iCellNew)
        call other%mixCells%add_col_to_cell(iOld, mix_iCellNew)
        
    end subroutine after_switch_update
    
    subroutine switch(type1, type1_obs, type2, type2_obs, mix, mix_Epot)
    
        class(HardSpheres), intent(inout) :: type1, type2
        class(Observables), intent(inout) :: type1_obs, type2_obs
        class(MixingPotential), intent(in) :: mix
        real(DP), intent(inout) :: mix_Epot
        
        real(DP) :: random
        integer :: type1_iOld, type2_iOld
        real(DP), dimension(Ndim) :: type1_xOld, type2_xOld, xRand, type1_xNew, type2_xNew
        logical :: overlap
        integer :: type1_iCellOld, type1_iCellNew
        integer :: mix_iCellOld, mix_iCellNew
        real(DP) :: deltaEpot
        real(DP) :: type1_deltaEpot, mix_deltaEpot
        real(DP) :: type1_EpotNew, type1_EpotOld
        real(DP) :: mix_EpotNew, mix_EpotOld
        
        real(DP), dimension(Ndim) :: mCol
        real(DP) :: type1_EpotNew_real, type1_EpotOld_real
        
        type1_obs%Nswitch = type1_obs%Nswitch + 1
        
        if (type1%get_Ncol()==0 .or. type2%get_Ncol()==0) then
            type1_obs%switch_Nreject = type1_obs%switch_Nreject + 1
            return
        end if
        
        ! Old : before switch
        call random_number(random)
        type1_iOld = int(random*type1%get_Ncol()) + 1        
        call before_switch_energy(type1, type1_iOld, type2, type1_EpotOld)
        
        call random_number(random)
        type2_iOld = int(random*type2%get_Ncol()) + 1        
        call before_switch_energy(type2, type2_iOld, type1, type2_EpotOld)
        
        ! New : after switch
        call after_switch_energy(type1, type2_iOld, type2, mix, overlap, type1_EpotNew)
        
        if (.not. overlap) then
        
            call after_switch_energy(type2, type1_iOld, type1, mix, overlap, type2_EpotNew)
            
            deltaEpot = type1_EpotNew - type1_EpotOld + type2_EpotNew - type2_EpotOld
            
            call random_number(random)
            if (random < exp(-deltaEpot/Temperature)) then
            
                call after_switch_update(type1, type1_iOld, type2, mix)
                call after_switch_update(type2, type2_iOld, type1, mix)
            
            end if
        
        end if
        
    end subroutine switch
