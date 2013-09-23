module module_algorithms

use data_precisions, only : DP
use data_box, only : Ndim, Lsize
use data_monteCarlo, only : Temperature
use class_observables
use class_hardSpheres
use class_interactingSpheres
use class_dipolarSpheres
use class_mixingPotential


implicit none

contains

    !> Particle move
    
    subroutine move(this, this_obs, other, mix, mix_Epot)
    
        class(HardSpheres), intent(inout) :: this
        class(Observables) :: this_obs
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
        
        call random_number(random)
        iOld = int(random*this%Ncol) + 1
        xOld(:) = this%positions(:, iOld)
        
        ! Random new position
        call random_number(xRand)
        xNew(:) = xOld(:) + (xRand(:)-0.5_DP)*this%move_delta(:)
        xNew(:) = modulo(xNew(:), Lsize(:))
        
        select type (this)
            type is (DipolarSpheres)
                mCol(:) = this%orientations(:, iOld)
        end select
        
        if (this%Ncol >= other%Ncol) then        
            this_iCellNew = this%sameCells%index_from_position(xNew)
            call this%Epot_neighCells(iOld, xNew, this_iCellNew, overlap, this_EpotNew)
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
                call this%Epot_neighCells(iOld, xNew, this_iCellNew, overlap, this_EpotNew)
            end if
                        
            if (.not. overlap) then
    
                this_iCellOld = this%sameCells%index_from_position(xOld)                
                select type (this)
                    type is (DipolarSpheres)
                        this_EpotNew_real = this%Epot_real_solo(iOld, xNew, mCol)
                        this_EpotOld_real = this%Epot_real_solo(iOld, xOld, mCol)
                        this_deltaEpot = (this_EpotNew_real-this_EpotOld_real) + &
                                         this%deltaEpot_reci_move(xOld, xNew, mCol)
                    class default
                        call this%Epot_neighCells(iOld, xOld, this_iCellOld, overlap, this_EpotOld)
                        this_deltaEpot = this_EpotNew - this_EpotOld
                end select
                    
                mix_iCellOld = this%mixCells%index_from_position(xOld)
                call mix%Epot_neighCells(xOld, mix_iCellOld, this%mixCells, other%positions, overlap, &
                                         mix_EpotOld)
                
                mix_deltaEpot = mix_EpotNew - mix_EpotOld

                deltaEpot = this_deltaEpot + mix_deltaEpot
                
                call random_number(random)
                if (random < exp(-deltaEpot/Temperature)) then
                
                    select type (this)
                        type is (DipolarSpheres)
                            call this%deltaEpot_reci_move_update_structure(xOld, xNew, mCol)
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

    subroutine widom(this, other, mix, activ)
        
        class(HardSpheres), intent(in) :: this
        class(HardSpheres), intent(in) :: other
        class(MixingPotential), intent(in) :: mix
        real(DP), intent(inOut) :: activ 
        
        integer :: iWidom
        real(DP) :: widTestSum
        real(DP), dimension(Ndim) :: xRand, xTest
        integer :: this_iCellTest, mix_iCellTest
        logical :: overlap
        real(DP) :: EpotTest, this_EpotTest, mix_EpotTest
        
        widTestSum = 0._DP
        
        do iWidom = 1, this%Nwidom           
            
            call random_number(xRand)
            xTest(:) = Lsize(:) * xRand(:)

            if (this%Ncol >= other%Ncol) then
                this_iCellTest = this%sameCells%index_from_position(xTest)
                call this%Epot_neighCells(0, xTest, this_iCellTest, overlap, this_EpotTest)
            else
                mix_iCellTest = this%mixCells%index_from_position(xTest)
                call mix%Epot_neighCells(xTest, mix_iCellTest, this%mixCells, other%positions, &
                                         overlap, mix_EpotTest)
            end if
            
            if (.not. overlap) then
            
                if (this%Ncol >= other%Ncol) then
                    mix_iCellTest = this%mixCells%index_from_position(xTest)
                    call mix%Epot_neighCells(xTest, mix_iCellTest, this%mixCells, other%positions, &
                                            overlap, mix_EpotTest)
                else
                    this_iCellTest = this%sameCells%index_from_position(xTest)
                    call this%Epot_neighCells(0, xTest, this_iCellTest, overlap, this_EpotTest)
                end if
                
                if (.not. overlap) then
                
                    EpotTest = this_EpotTest + mix_EpotTest
                    widTestSum = widTestSum + exp(-EpotTest/Temperature)
                    
                end if
                
            end if
            
        end do
        
        activ = widTestSum/real(this%Nwidom, DP)
        
    end subroutine widom

end module module_algorithms
