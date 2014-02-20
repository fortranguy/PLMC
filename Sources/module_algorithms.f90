module module_algorithms

use data_precisions, only : DP
use data_box, only : Ndim, Lsize
use data_monteCarlo, only : Temperature
use module_physics, only : random_surface, markov_surface
use class_hardSpheres
use class_interactingSpheres
use class_dipolarSpheres
use class_mixingPotential
use class_observables

implicit none
private
public move, widom, switch, rotate

contains

    !> Particle move
    
    subroutine move(this, this_obs, other, mix, mix_Epot)
    
        class(HardSpheres), intent(inout) :: this
        class(Observables), intent(inout) :: this_obs
        class(HardSpheres), intent(inout) :: other
        class(MixingPotential), intent(in) :: mix
        real(DP), intent(inout) :: mix_Epot
        
        real(DP) :: random
        integer :: iOld
        real(DP), dimension(Ndim) :: xOld, xRand, xNew
        logical :: overlap
        integer :: this_iCellOld, this_iCellNew
        integer :: mix_iCellOld, mix_iCellNew
        real(DP) :: deltaEpot
        real(DP) :: this_deltaEpot, mix_deltaEpot
        real(DP) :: this_EpotNew, this_EpotOld
        real(DP) :: mix_EpotNew, mix_EpotOld
        
        real(DP), dimension(Ndim) :: mCol
        real(DP) :: this_EpotNew_real, this_EpotOld_real
        
        this_obs%move_Nhit = this_obs%move_Nhit + 1
        
        call random_number(random)
        iOld = int(random*this%get_Ncol()) + 1
        xOld(:) = this%positions(:, iOld)
        
        ! Random new position
        call random_number(xRand)
        xNew(:) = xOld(:) + (xRand(:)-0.5_DP)*this%move_delta(:)
        xNew(:) = modulo(xNew(:), Lsize(:))
        
        if (this%get_Ncol() >= other%get_Ncol()) then
            this_iCellNew = this%sameCells%index_from_position(xNew)
            call this%Epot_neighCells(iOld, xNew, this_iCellNew, overlap, this_EpotNew)
        else
            mix_iCellNew = this%mixCells%index_from_position(xNew)
            call mix%Epot_neighCells(0, xNew, mix_iCellNew, this%mixCells, other%positions, overlap, &
                                     mix_EpotNew)
        end if
        
        if (.not. overlap) then
        
            if (this%get_Ncol() >= other%get_Ncol()) then
                mix_iCellNew = this%mixCells%index_from_position(xNew)
                call mix%Epot_neighCells(0, xNew, mix_iCellNew, this%mixCells, other%positions, &
                                         overlap, mix_EpotNew)
            else
                this_iCellNew = this%sameCells%index_from_position(xNew)
                call this%Epot_neighCells(iOld, xNew, this_iCellNew, overlap, this_EpotNew)
            end if
                        
            if (.not. overlap) then
    
                this_iCellOld = this%sameCells%index_from_position(xOld)
                select type (this)
                    type is (DipolarSpheres)
                        mCol(:) = this%orientations(:, iOld)
                        this_EpotNew_real = this%Epot_real_solo(iOld, xNew, mCol)
                        this_EpotOld_real = this%Epot_real_solo(iOld, xOld, mCol)
                        this_deltaEpot = (this_EpotNew_real-this_EpotOld_real) + &
                                         this%deltaEpot_reci_move(xOld, xNew, mCol)
                    class default
                        call this%Epot_neighCells(iOld, xOld, this_iCellOld, overlap, this_EpotOld)
                        this_deltaEpot = this_EpotNew - this_EpotOld
                end select
                    
                mix_iCellOld = this%mixCells%index_from_position(xOld)
                call mix%Epot_neighCells(0, xOld, mix_iCellOld, this%mixCells, other%positions, &
                                         overlap, mix_EpotOld)
                
                mix_deltaEpot = mix_EpotNew - mix_EpotOld

                deltaEpot = this_deltaEpot + mix_deltaEpot
                
                call random_number(random)
                if (random < exp(-deltaEpot/Temperature)) then
                
                    select type (this)
                        type is (DipolarSpheres)
                            call this%reci_update_structure_move(xOld, xNew, mCol)
                    end select
                
                    this%positions(:, iOld) = xNew(:)
                    this_obs%Epot = this_obs%Epot + this_deltaEpot
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
    
    end subroutine move
    
    !> Widom's method

    subroutine widom(this, this_obs, other, mix)
        
        class(HardSpheres), intent(in) :: this
        class(Observables), intent(inout) :: this_obs
        class(HardSpheres), intent(in) :: other
        class(MixingPotential), intent(in) :: mix
        
        integer :: iWidom
        real(DP) :: widTestSum
        real(DP), dimension(Ndim) :: xRand, xTest
        integer :: this_iCellTest, mix_iCellTest
        logical :: overlap
        real(DP) :: EpotTest
        real(DP) :: this_EpotTest, mix_EpotTest
        
        real(DP), dimension(Ndim) :: mTest
        
        widTestSum = 0._DP
        
        do iWidom = 1, this%get_Nwidom()
            
            call random_number(xRand)
            xTest(:) = Lsize(:) * xRand(:)

            if (this%get_Ncol() >= other%get_Ncol()) then
                this_iCellTest = this%sameCells%index_from_position(xTest)
                call this%Epot_neighCells(0, xTest, this_iCellTest, overlap, this_EpotTest)
            else
                mix_iCellTest = this%mixCells%index_from_position(xTest)
                call mix%Epot_neighCells(0, xTest, mix_iCellTest, this%mixCells, other%positions, &
                                         overlap, mix_EpotTest)
            end if
            
            if (.not. overlap) then
            
                if (this%get_Ncol() >= other%get_Ncol()) then
                    mix_iCellTest = this%mixCells%index_from_position(xTest)
                    call mix%Epot_neighCells(0, xTest, mix_iCellTest, this%mixCells, other%positions, &
                                             overlap, mix_EpotTest)
                else
                    this_iCellTest = this%sameCells%index_from_position(xTest)
                    call this%Epot_neighCells(0, xTest, this_iCellTest, overlap, this_EpotTest)
                end if
                
                if (.not. overlap) then
                
                    select type (this)
                        type is (DipolarSpheres)
                            mTest(:) = random_surface()
                            this_EpotTest = this%Epot_real_solo(0, xTest, mTest) + &
                                            this%deltaEpot_reci_exchange(xTest, +mTest) - &
                                            this%Epot_self_solo(mTest) + &
                                            this%deltaEpot_bound_exchange(+mTest)
                    end select
                
                    EpotTest = this_EpotTest + mix_EpotTest
                    widTestSum = widTestSum + exp(-EpotTest/Temperature)
                    
                end if
                
            end if
            
        end do
        
        this_obs%activ = widTestSum/real(this%get_Nwidom(), DP)
        
    end subroutine widom
    
    !> Particle switch
    
    subroutine before_switch_energy(this, this_iCol, other, other_iCol, mix, indicesOld, xOld, &
                                    EpotsOld)
        
        class(HardSpheres), intent(in) :: this, other
        class(MixingPotential), intent(in) :: mix
        integer, intent(in) :: this_iCol, other_iCol
        integer, dimension(:), intent(out) :: indicesOld
        real(DP), dimension(:), intent(out) :: xOld, EpotsOld
        
        real(DP), dimension(Ndim) :: mCol
        logical :: overlap
        
        xOld(:) = this%positions(:, this_iCol)
        
        indicesOld(1) = this%sameCells%index_from_position(xOld)
        select type (this)
            type is (DipolarSpheres)
                mCol(:) = this%orientations(:, this_iCol)
                EpotsOld(1) = this%Epot_real_solo(this_iCol, xOld, mCol) ! Epot_reci: 
                                                                         ! cf. after_switch_energy
            type is (InteractingSpheres)
                call this%Epot_neighCells(this_iCol, xOld, indicesOld(1), overlap, EpotsOld(1))
            type is (HardSpheres)
                EpotsOld(1) = 0._DP
        end select
        
        indicesOld(2) = other%mixCells%index_from_position(xOld)
        call mix%Epot_neighCells(other_iCol, xOld, indicesOld(2), this%mixCells, other%positions, &
                                 overlap, EpotsOld(2))
        
    end subroutine before_switch_energy
    
    subroutine after_switch_energy(this, this_iCol, xOld, other, other_iCol, mix, overlap, &
                                   indicesNew, xNew, EpotsNew)
        
        class(HardSpheres), intent(in) :: this, other
        class(MixingPotential), intent(in) :: mix
        integer, intent(in) :: this_iCol, other_iCol
        real(DP), dimension(:), intent(in) :: xOld
        integer, dimension(:), intent(out) :: indicesNew
        logical, intent(out) :: overlap
        real(DP), dimension(:), intent(out) :: xNew, EpotsNew
        
        real(DP), dimension(Ndim) :: mCol
        
        xNew(:) = other%positions(:, other_iCol)
        
        if (this%get_Ncol() >= other%get_Ncol()) then ! optimisation : more chance to overlap
            indicesNew(1) = this%sameCells%index_from_position(xNew)
            call this%Epot_neighCells(this_iCol, xNew, indicesNew(1), overlap, EpotsNew(1))
        else
            indicesNew(2) = other%mixCells%index_from_position(xNew)
            call mix%Epot_neighCells(other_iCol, xNew, indicesNew(2), this%mixCells, other%positions, &
                                     overlap, EpotsNew(2))
        end if
        
        if (.not. overlap) then
        
            if (this%get_Ncol() >= other%get_Ncol()) then
                indicesNew(2) = other%mixCells%index_from_position(xNew)
                call mix%Epot_neighCells(other_iCol, xNew, indicesNew(2), this%mixCells, &
                                         other%positions, overlap, EpotsNew(2))
            else
                indicesNew(1) = this%sameCells%index_from_position(xNew)
                call this%Epot_neighCells(this_iCol, xNew, indicesNew(1), overlap, EpotsNew(1))
            end if
            
            if (.not. overlap) then
            
                select type (this)
                    type is (DipolarSpheres)
                        mCol(:) = this%orientations(:, this_iCol)
                        EpotsNew(1) = this%Epot_real_solo(this_iCol, xNew, mCol) + &
                                      this%deltaEpot_reci_move(xOld, xNew, mCol)
                end select
            
            end if
                
        end if
    
    end subroutine after_switch_energy
    
    subroutine after_switch_update(this, iCol, indicesOld, indicesNew, xOld, xNew, other)
        
        class(HardSpheres), intent(inout) :: this, other
        integer, intent(in) :: iCol
        integer, dimension(:), intent(in) :: indicesOld, indicesNew
        real(DP), dimension(:), intent(in) :: xOld, xNew
        
        real(DP), dimension(Ndim) :: mCol
        
        this%positions(:, iCol) = xNew(:)
        
        select type (this)
            type is (DipolarSpheres)
                mCol(:) = this%orientations(:, iCol)
                call this%reci_update_structure_exchange(xOld, -mCol)
        end select
        
        select type (this)
            type is (DipolarSpheres)
                call this%reci_update_structure_exchange(xNew, +mCol)
        end select
        
        if (indicesOld(1) /= indicesNew(1)) then
            call this%sameCells%remove_col_from_cell(iCol, indicesOld(1))
            call this%sameCells%add_col_to_cell(iCol, indicesNew(1))
        end if
        if (indicesOld(2) /= indicesNew(2)) then
            call other%mixCells%remove_col_from_cell(iCol, indicesOld(2))
            call other%mixCells%add_col_to_cell(iCol, indicesNew(2))
        end if
        
    end subroutine after_switch_update
    
    subroutine switch(type1, type1_obs, type2, type2_obs, mix, mix_Epot, switch_Nreject)
    
        class(HardSpheres), intent(inout) :: type1, type2
        class(Observables), intent(inout) :: type1_obs, type2_obs
        class(MixingPotential), intent(in) :: mix
        real(DP), intent(inout) :: mix_Epot
        integer, intent(inout) :: switch_Nreject
        
        real(DP) :: random
        integer :: type1_iCol, type2_iCol
        integer, dimension(2) :: type1_indicesOld, type1_indicesNew ! (1): same, (2): other
        integer, dimension(2) :: type2_indicesOld , type2_indicesNew
        logical :: overlap
        real(DP), dimension(Ndim) :: type1_xOld, type1_xNew
        real(DP), dimension(Ndim) :: type2_xOld, type2_xNew
        real(DP) :: deltaEpot, type1_deltaEpot, type2_deltaEpot
        real(DP) :: type1_mix_deltaEpot, type2_mix_deltaEpot
        real(DP), dimension(2) :: type1_EpotsOld, type1_EpotsNew ! (1): same, (2): mix
        real(DP), dimension(2) :: type2_EpotsOld, type2_EpotsNew
        
        if (type1%get_Ncol()==0 .or. type2%get_Ncol()==0) then
            switch_Nreject = switch_Nreject + 1
            return
        end if
        
        ! Old : before switch
        call random_number(random)
        type1_iCol = int(random*type1%get_Ncol()) + 1
        call before_switch_energy(type1, type1_iCol, type2, type2_iCol, mix, type1_indicesOld, &
                                  type1_xOld, type1_EpotsOld)
        
        call random_number(random)
        type2_iCol = int(random*type2%get_Ncol()) + 1
        call before_switch_energy(type2, type2_iCol, type1, type1_iCol, mix, type2_indicesOld, &
                                  type2_xOld, type2_EpotsOld)
        
        ! New : after switch
        call after_switch_energy(type1, type1_iCol, type1_xOld, type2, type2_iCol, mix, overlap, &
                                 type1_indicesNew, type1_xNew, type1_EpotsNew)
        
        if (.not. overlap) then
        
            call after_switch_energy(type2, type2_iCol, type2_xOld, type1, type1_iCol, mix, overlap, &
                                     type2_indicesNew, type2_xNew, type2_EpotsNew)
            
            if (.not. overlap) then

                type1_deltaEpot = type1_EpotsNew(1) - type1_EpotsOld(1)
                type1_mix_deltaEpot = type1_EpotsNew(2) - type1_EpotsOld(2)
                type2_deltaEpot = type2_EpotsNew(1) - type2_EpotsOld(1)
                type2_mix_deltaEpot = type2_EpotsNew(2) - type2_EpotsOld(2)
                deltaEpot = type1_deltaEpot + type1_mix_deltaEpot + type2_deltaEpot + &
                            type2_mix_deltaEpot
                
                call random_number(random)
                if (random < exp(-deltaEpot/Temperature)) then
                
                    call after_switch_update(type1, type1_iCol, type1_indicesOld, type1_indicesNew, &
                                             type1_xOld, type1_xNew, type2)
                    call after_switch_update(type2, type2_iCol, type2_indicesOld, type2_indicesNew, &
                                             type2_xOld, type2_xNew, type1)
                                             
                    type1_obs%Epot = type1_obs%Epot + type1_deltaEpot
                    type2_obs%Epot = type2_obs%Epot + type2_deltaEpot
                    mix_Epot = mix_Epot + type1_mix_deltaEpot + type2_mix_deltaEpot
                    
                else
                
                    switch_Nreject = switch_Nreject + 1
                    
                end if
                
            else
            
                switch_Nreject = switch_Nreject + 1
            
            end if
            
        else
        
            switch_Nreject = switch_Nreject + 1
            
        end if
        
    end subroutine switch
    
    !> Dipole rotation
    
    subroutine rotate(this, obs)
    
        class(DipolarSpheres), intent(inout) :: this
        class(MoreObservables), intent(inout) :: obs
        
        real(DP) :: random
        integer :: iOld
        real(DP), dimension(Ndim) :: xCol
        real(DP), dimension(Ndim) :: mOld, mNew
        real(DP) :: deltaEpot
        real(DP) :: deltaEpot_real, deltaEpot_self
        real(DP) :: real_EpotNew, real_EpotOld
        
        obs%rotate_Nhit = obs%rotate_Nhit + 1

        call random_number(random)
        iOld = int(random*this%get_Ncol()) + 1
        
        xCol(:) = this%positions(:, iOld)
        mOld(:) = this%orientations(:, iOld)
        mNew(:) = mOld(:)
        call markov_surface(mNew, this%rotate_delta)
        
        real_EpotNew = this%Epot_real_solo(iOld, xCol, mNew)
        real_EpotOld = this%Epot_real_solo(iOld, xCol, mOld)
        deltaEpot_real = real_EpotNew - real_EpotOld
        
        deltaEpot_self = this%Epot_self_solo(mNew) - this%Epot_self_solo(mOld)
        
        deltaEpot = deltaEpot_real + this%deltaEpot_reci_rotate(xCol, mOld, mNew) - deltaEpot_self + &
                    this%deltaEpot_bound_rotate(mOld, mNew)
        
        call random_number(random)
        if (random < exp(-deltaEpot/Temperature)) then
        
            call this%reci_update_structure_rotate(xCol, mOld, mNew)
            call this%update_totalMoment_rotate(mOld, mNew)
            this%orientations(:, iOld) = mNew(:)
            
            obs%Epot = obs%Epot + deltaEpot
            
        else
            obs%rotate_Nreject = obs%rotate_Nreject + 1
        end if
    
    end subroutine rotate

end module module_algorithms
