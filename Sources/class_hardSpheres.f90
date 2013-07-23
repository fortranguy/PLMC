!> \brief Description of the Hard Spheres class

module class_hardSpheres

use, intrinsic :: iso_fortran_env, only : output_unit
use data_precisions, only : DP
use data_box, only : Ndim, Lsize
use data_particles, only : hard_radius, hard_rMin, hard_Ncol
use data_potential, only : hard_rCut
use data_monteCarlo, only : Temperature, hard_deltaX, hard_rejectFix, hard_Nadapt, hard_Nwidom
use data_neighbourCells, only : NnearCell, hard_cell_size
use data_distribution, only : hard_snap_factor
use module_physics, only : dist_PBC
use class_observables
use class_neighbourCells
use class_mixingPotential
use class_spheres

implicit none

private

    type, extends(Spheres), public :: HardSpheres
    
        private
    
        ! Potential
        real(DP) :: Epot
        
    contains

        !> Construction and destruction of the class
        procedure :: construct => HardSpheres_construct
        procedure :: destroy => HardSpheres_destroy
        
        !> Print a report of the component in a file
        procedure :: printReport => HardSpheres_printReport
              
        !> Potential energy
        procedure :: Epot_print => HardSpheres_Epot_print
        procedure, private :: Epot_neighCells => HardSpheres_Epot_neighCells
        procedure :: Epot_conf => HardSpheres_Epot_conf
        procedure :: consistTest => HardSpheres_consistTest
        
        !> Monte-Carlo
        procedure :: move => HardSpheres_move
        procedure :: widom => HardSpheres_widom
        
    end type HardSpheres
    
contains

    subroutine HardSpheres_construct(this, mix_cell_size, mix_rCut)
    
        class(HardSpheres), intent(out) :: this
        real(DP), dimension(:), intent(in) :: mix_cell_size
        real(DP), intent(in) :: mix_rCut
        
        this%name = "hardS"
        write(output_unit, *) this%name, " class construction"
    
        ! Particles
        this%radius = hard_radius
        this%rMin = hard_rMin
        this%Ncol = hard_Ncol
        allocate(this%positions(Ndim, this%Ncol))
        
        ! Snapshot
        this%snap_factor = hard_snap_factor
        
        ! Monte-Carlo
        this%deltaX = hard_deltaX
        this%deltaXSave = this%deltaX
        this%rejectFix = hard_rejectFix
        this%Nadapt = hard_Nadapt
        this%Nwidom = hard_Nwidom
                
        ! Potential
        this%rCut = hard_rCut
        this%Epot = 0._DP
        
        ! Neighbour Cells
        call this%sameCells%construct(hard_cell_size, this%rCut) !< same kind
        call this%mixCells%construct(mix_cell_size, mix_rCut) !< other kind
    
    end subroutine HardSpheres_construct
    
    subroutine HardSpheres_destroy(this)
    
        class(HardSpheres), intent(inout) :: this
        
        write(output_unit, *) this%name, " class destruction"
        
        if (allocated(this%positions)) then
            deallocate(this%positions)
        end if
        
        call this%sameCells%destroy()
        call this%mixCells%destroy()
    
    end subroutine HardSpheres_destroy
    
    !> Report
    
    subroutine HardSpheres_printReport(this, report_unit)
    
        class(HardSpheres), intent(in) :: this
        integer, intent(in) :: report_unit    
        
        write(report_unit, *) "Data :"
        
        write(report_unit ,*) "    Ncol = ", this%Ncol
        write(report_unit ,*) "    Nwidom = ", this%Nwidom
        write(report_unit ,*) "    Nadapt = ", this%Nadapt
        
        write(report_unit, *) "    rCut = ", this%rCut
        
        write(report_unit, *) "    same_NtotalCell_dim(:) = ", this%sameCells%getNtotalCell_dim()
        write(report_unit, *) "    same_cell_size(:) = ", this%sameCells%getCell_size()        
        write(report_unit, *) "    mix_NtotalCell_dim(:) = ", this%mixCells%getNtotalCell_dim()
        write(report_unit, *) "    mix_cell_size(:) = ", this%mixCells%getCell_size()
        
    end subroutine HardSpheres_printReport
    
    !> Print the potential : dummy
    
    subroutine HardSpheres_Epot_print(this, Epot_unit)
    
        class(HardSpheres), intent(in) :: this
        integer, intent(in) :: Epot_unit

        write(Epot_unit, *) this%rCut, this%Epot
    
    end subroutine HardSpheres_Epot_print
    
    subroutine HardSpheres_Epot_neighCells(this, iCol, xCol, iTotalCell, overlap)
        
        class(HardSpheres), intent(in) :: this        
        integer, intent(in) :: iCol, iTotalCell
        real(DP), dimension(:), intent(in) :: xCol
        logical, intent(out) :: overlap
    
        integer :: iNearCell,  nearCell_index
        real(DP) :: r_ij
    
        type(Link), pointer :: current => null(), next => null()
        
        overlap = .false.
    
        do iNearCell = 1, NnearCell
        
            nearCell_index = this%sameCells%nearCells_among_totalCells(iNearCell, iTotalCell)
            current => this%sameCells%beginCells(nearCell_index)%particle%next            
            if (.not. associated(current%next)) cycle
            
            do
            
                next => current%next
            
                if (current%iCol /= iCol) then
                
                    r_ij = dist_PBC(xCol(:), this%positions(:, current%iCol))
                    if (r_ij < this%rMin) then
                        overlap = .true.
                        return
                    end if
       
                end if
                
                if (.not. associated(next%next)) exit
                
                current => next
            
            end do            
            
        end do
    
    end subroutine HardSpheres_Epot_neighCells
    
    !> Particle move
    
    subroutine HardSpheres_move(this, iOld, other, mix, same_obs, mix_Epot)
    
        class(HardSpheres), intent(inout) :: this
        integer, intent(in) :: iOld
        class(Spheres), intent(inout) :: other
        class(MixingPotential), intent(in) :: mix
        class(Observables) :: same_obs
        real(DP), intent(inout) :: mix_Epot
        
        real(DP), dimension(Ndim) :: xRand
        logical :: overlap
        real(DP), dimension(Ndim) :: xNew
        integer :: same_iCellOld, same_iCellNew
        integer :: mix_iCellOld, mix_iCellNew
        real(DP) :: mix_deltaEpot
        real(DP) :: mix_eNew, mix_eOld
        real(DP) :: random
        
        ! Random new position
        call random_number(xRand)
        xNew(:) = this%positions(:, iOld) + (xRand(:)-0.5_DP)*this%deltaX(:)
        xNew(:) = modulo(xNew(:), Lsize(:))
        
        same_iCellNew = this%sameCells%index_from_position(xNew)
        call this%Epot_neighCells(iOld, xNew, same_iCellNew, overlap)
        
        if (.not. overlap) then
        
            mix_iCellNew = this%mixCells%index_from_position(xNew)
            call mix%Epot_neighCells(xNew, mix_iCellNew, this%mixCells, other%positions, overlap, &
                                     mix_eNew)
                        
            if (.not. overlap) then
    
                same_iCellOld = this%sameCells%index_from_position(this%positions(:, iOld))
                call this%Epot_neighCells(iOld, this%positions(:, iOld), same_iCellOld, overlap)
                    
                mix_iCellOld = this%mixCells%index_from_position(this%positions(:, iOld))
                call mix%Epot_neighCells(this%positions(:, iOld), mix_iCellOld, this%mixCells, &
                                         other%positions, overlap, mix_eOld)
                
                mix_deltaEpot = mix_eNew - mix_eOld
                
                call random_number(random)
                if (random < exp(-mix_deltaEpot/Temperature)) then
                
                    this%positions(:, iOld) = xNew(:)
                    same_obs%Epot = same_obs%Epot + 0._DP
                    mix_Epot = mix_Epot + mix_deltaEpot
                    
                    if (same_iCellOld /= same_iCellNew) then                
                        call this%sameCells%remove_col_from_cell(iOld, same_iCellOld)
                        call this%sameCells%add_col_to_cell(iOld, same_iCellNew)
                    end if
                    
                    if (mix_iCellOld /= mix_iCellNew) then                
                        call other%mixCells%remove_col_from_cell(iOld, mix_iCellOld)
                        call other%mixCells%add_col_to_cell(iOld, mix_iCellNew)
                    end if
                    
                else
                    same_obs%Nreject = same_obs%Nreject + 1
                end if
         
            else
                same_obs%Nreject = same_obs%Nreject + 1
            end if            
            
        else        
            same_obs%Nreject = same_obs%Nreject + 1
        end if
    
    end subroutine HardSpheres_move
    
    !> Widom's method : with other type ?

    subroutine HardSpheres_widom(this, other_positions, mix, activ)
        
        class(HardSpheres), intent(in) :: this
        real(DP), dimension(:, :), intent(in) :: other_positions
        class(MixingPotential), intent(in) :: mix
        real(DP), intent(inOut) :: activ 
        
        integer :: iWidom
        real(DP) :: widTestSum
        real(DP), dimension(Ndim) :: xRand, xTest
        integer :: same_iCellTest, mix_iCellTest
        logical :: overlap
        real(DP) :: enTest, mix_enTest
        
        widTestSum = 0._DP
        
        do iWidom = 1, this%Nwidom           
            
            call random_number(xRand)
            xTest(:) = Lsize(:) * xRand(:)
            same_iCellTest = this%sameCells%index_from_position(xTest)
            call this%Epot_neighCells(0, xTest, same_iCellTest, overlap)
            
            if (.not. overlap) then
            
                mix_iCellTest = this%mixCells%index_from_position(xTest)
                call mix%Epot_neighCells(xTest, mix_iCellTest, this%mixCells, other_positions, &
                                         overlap, mix_enTest)
                
                if (.not. overlap) then
                
                    enTest = 0._DP + mix_enTest
                    widTestSum = widTestSum + exp(-enTest/Temperature)
                    
                end if
                
            end if
            
        end do
        
        activ = widTestSum/real(this%Nwidom, DP)
        
    end subroutine HardSpheres_widom
    
    !> Total potential energy : dummy
    
    pure function HardSpheres_Epot_conf(this) result(Epot_conf)
    
        class(HardSpheres), intent(in) :: this        
        real(DP) :: Epot_conf
    
        Epot_conf = this%Epot
        
    end function HardSpheres_Epot_conf
    
    !> Consistency test : dummy
    
    subroutine HardSpheres_consistTest(this, Epot, report_unit)
    
        class(HardSpheres), intent(in) :: this
        real(DP), intent(in) :: Epot
        integer, intent(in) :: report_unit
        
        real(DP) :: Epot_conf
        real(DP) :: difference
    
        Epot_conf = this%Epot_conf()
        difference = abs(Epot_conf-Epot)
        
        write(report_unit, *) "Consistency test:"
        write(report_unit, *) "    Epot = ", Epot
        write(report_unit, *) "    Epot_conf = ", Epot_conf
        write(report_unit, *) "    absolute difference = ", difference
        
        if (difference /= 0._DP) then
            write(report_unit, *) "    WARNING !"
        else
            write(report_unit, *) "    OK !"
        end if
    
    end subroutine HardSpheres_consistTest

end module class_hardSpheres
