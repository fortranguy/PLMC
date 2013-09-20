module module_algorithms

use data_precisions, only : DP
use data_box, only : Ndim, Lsize
use data_monteCarlo, only : Temperature
use class_observables
use class_spheres
use class_mixingPotential
use class_hardSpheres
use class_interactingSpheres
use class_dipolarSpheres
use class_observables


implicit none

contains

    subroutine polymorph(sph)
    
        class(Spheres), intent(in) :: sph
        
        select type (sph)
        
            type is (hardSpheres)
            
                write(*, *) "hard : ", sph%name
        
            type is (interactingSpheres)
            
                write(*, *) "inter : ", sph%name
            
            type is (dipolarSpheres)
            
                write(*, *) "dipol : ", sph%name
        
        end select 

    end subroutine polymorph
    
    !> Particle move
    
    subroutine move(this, this_obs, other, mix, mix_Epot)
    
        class(Spheres), intent(inout) :: this
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
        real(DP) :: deltaEpot
        real(DP) :: this_deltaEpot, mix_deltaEpot
        real(DP) :: this_EpotNew, this_EpotOld
        real(DP) :: mix_EpotNew, mix_EpotOld
        
        call random_number(random)
        iOld = int(random*this%Ncol) + 1
        xOld(:) = this%positions(:, iOld)
        
        ! Random new position
        call random_number(xRand)
        xNew(:) = xOld(:) + (xRand(:)-0.5_DP)*this%move_delta(:)
        xNew(:) = modulo(xNew(:), Lsize(:))
        
        if (this%Ncol >= other%Ncol) then        
            this_iCellNew = this%sameCells%index_from_position(xNew)
            select type (this)
                type is (hardSpheres)
                    call this%Epot_neighCells(iOld, xNew, this_iCellNew, overlap)
                type is (interactingSpheres)
                    call this%Epot_neighCells(iOld, xNew, this_iCellNew, overlap, this_EpotNew)                
            end select
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
                select type (this)
                    type is (hardSpheres)
                        call this%Epot_neighCells(iOld, xNew, this_iCellNew, overlap)
                    type is (interactingSpheres)
                        call this%Epot_neighCells(iOld, xNew, this_iCellNew, overlap, this_EpotNew)
                end select
            end if
                        
            if (.not. overlap) then
    
                this_iCellOld = this%sameCells%index_from_position(xOld)
                select type (this)
                    type is (hardSpheres)
                        call this%Epot_neighCells(iOld, xOld, this_iCellOld, overlap)
                        this_deltaEpot = 0._DP
                    type is (interactingSpheres)
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

end module module_algorithms
