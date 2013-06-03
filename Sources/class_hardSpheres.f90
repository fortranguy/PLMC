!> \brief Description of the Hard Spheres class

module class_hardSpheres

use data_constants
use data_cell
use data_particles
use data_potentiel
use data_mc
use data_neighbours
use data_distrib
use mod_physics
use class_neighbours
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
        procedure :: report => HardSpheres_report
              
        !> Potential energy
        procedure :: Epot_print => HardSpheres_Epot_print
        procedure :: Epot_neigh => HardSpheres_Epot_neigh
        procedure :: Epot_conf => HardSpheres_Epot_conf
        procedure :: consistTest => HardSpheres_consistTest
        
        !> Monte-Carlo
        procedure :: move => HardSpheres_move
        procedure :: widom => HardSpheres_widom
        
    end type HardSpheres
    
contains

    subroutine HardSpheres_construct(this, shared_cell_Lsize, shared_rCut)
    
        class(HardSpheres), intent(out) :: this
        real(DP), dimension(:), intent(in) :: shared_cell_Lsize
        real(DP), intent(in) :: shared_rCut
        
        this%name = "hardS"
    
        ! Particles
        this%radius = hard_radius
        this%rMin = hard_rMin
        this%Ncol = hard_Ncol
        allocate(this%positions(Dim, this%Ncol))
        
        ! Snapshot
        this%snap_factor = hard_snap_factor
        
        ! Monte-Carlo
        this%deltaX = hard_deltaX
        this%deltaXSave = this%deltaX
        this%rejFix = hard_rejFix
        this%Nadapt = hard_Nadapt
        this%Nwidom = hard_Nwidom
                
        ! Potential
        this%rCut = hard_rCut
        this%Epot = 0._DP
        
        ! Neighbours : same kind
        call this%same%construct(hard_cell_Lsize, this%rCut)
        call this%same%alloc_cells()
        call this%same%ini_cell_neighs()
        ! Neighbours : other kind
        call this%mix%construct(shared_cell_Lsize, shared_rCut)
        call this%mix%alloc_cells()
        call this%mix%ini_cell_neighs()
    
    end subroutine HardSpheres_construct
    
    subroutine HardSpheres_destroy(this)
    
        class(HardSpheres), intent(inout) :: this
        
        if (allocated(this%positions)) then
            deallocate(this%positions)
        end if
        
        call this%same%destroy()
        call this%mix%destroy()
    
    end subroutine HardSpheres_destroy
    
    !> Report
    
    subroutine HardSpheres_report(this, report_unit)
    
        class(HardSpheres), intent(in) :: this
        integer, intent(in) :: report_unit    
        
        write(report_unit, *) "Data :"
        
        write(report_unit ,*) "    Ncol = ", this%Ncol
        write(report_unit ,*) "    Nwidom = ", this%Nwidom
        write(report_unit ,*) "    Nadapt = ", this%Nadapt
        
        write(report_unit, *) "    rCut = ", this%rCut
        
        write(report_unit, *) "    same_cell_coordMax(:) = ", this%same%cell_coordMax(:)
        write(report_unit, *) "    same_cell_Lsize(:) = ", this%same%cell_Lsize(:)        
        write(report_unit, *) "    mix_cell_coordMax(:) = ", this%mix%cell_coordMax(:)
        write(report_unit, *) "    mix_cell_Lsize(:) = ", this%mix%cell_Lsize(:)
        
    end subroutine HardSpheres_report
    
    !> Print the potential : dummy
    
    subroutine HardSpheres_Epot_print(this, Epot_unit)
    
        class(HardSpheres), intent(in) :: this
        integer, intent(in) :: Epot_unit

        write(Epot_unit, *) this%rCut, this%Epot
    
    end subroutine HardSpheres_Epot_print
    
    subroutine HardSpheres_Epot_neigh(this, iCol, xCol, iCell, overlap)
        
        class(HardSpheres), intent(in) :: this        
        integer, intent(in) :: iCol, iCell
        real(DP), dimension(:), intent(in) :: xCol
        logical, intent(out) :: overlap
    
        integer :: iNeigh,  iCell_neigh
        real(DP) :: r
    
        type(Link), pointer :: current => null(), next => null()
        
        overlap = .false.
    
        do iNeigh = 1, cell_neighs_nb
        
            iCell_neigh = this%same%cell_neighs(iNeigh, iCell)
            current => this%same%cellsBegin(iCell_neigh)%particle%next            
            if (.not. associated(current%next)) cycle
            
            do
            
                next => current%next
            
                if (current%iCol /= iCol) then
                
                    r = dist(xCol(:), this%positions(:, current%iCol))
                    if (r < this%rMin) then
                        overlap = .true.
                        return
                    end if
       
                end if
                
                if (.not. associated(next%next)) exit
                
                current => next
            
            end do            
            
        end do
    
    end subroutine HardSpheres_Epot_neigh
    
    !> Particle move
    
    subroutine HardSpheres_move(this, iOld, other, mix, same_Epot, mix_Epot, Nrej)
    
        class(HardSpheres), intent(inout) :: this
        integer, intent(in) :: iOld
        class(Spheres), intent(inout) :: other
        class(MixingPotential), intent(in) :: mix
        real(DP), intent(inout) :: same_Epot, mix_Epot
        integer, intent(inout) :: Nrej
        
        real(DP), dimension(Dim) :: xRand
        logical :: overlap
        real(DP), dimension(Dim) :: xNew
        integer :: same_iCellOld, same_iCellNew
        integer :: mix_iCellOld, mix_iCellNew
        real(DP) :: mix_dEpot
        real(DP) :: mix_eNew, mix_eOld
        real(DP) :: rand
        
        ! Random new position
        call random_number(xRand)
        xNew(:) = this%positions(:, iOld) + (xRand(:)-0.5_DP)*this%deltaX(:)
        xNew(:) = modulo(xNew(:), Lsize(:))
        
        same_iCellNew = this%same%position_to_cell(xNew)
        call this%Epot_neigh(iOld, xNew, same_iCellNew, overlap)
        
        if (.not. overlap) then
        
            mix_iCellNew = this%mix%position_to_cell(xNew)
            call mix%Epot_neigh(xNew, mix_iCellNew, this%mix, other%positions, overlap, mix_eNew)
                        
            if (.not. overlap) then
    
                same_iCellOld = this%same%position_to_cell(this%positions(:, iOld))
                call this%Epot_neigh(iOld, this%positions(:, iOld), same_iCellOld, overlap)
                    
                mix_iCellOld = this%mix%position_to_cell(this%positions(:, iOld))
                call mix%Epot_neigh(this%positions(:, iOld), mix_iCellOld, this%mix, other%positions, &
                                    overlap, mix_eOld)
                
                mix_dEpot = mix_eNew - mix_eOld
                
                call random_number(rand)
                if (rand < exp(-mix_dEpot/Tstar)) then
                
                    this%positions(:, iOld) = xNew(:)
                    same_Epot = same_Epot + 0._DP
                    mix_Epot = mix_Epot + mix_dEpot
                    
                    if (same_iCellOld /= same_iCellNew) then                
                        call this%same%remove_cell_col(iOld, same_iCellOld)
                        call this%same%add_cell_col(iOld, same_iCellNew)
                    end if
                    
                    if (mix_iCellOld /= mix_iCellNew) then                
                        call other%mix%remove_cell_col(iOld, mix_iCellOld)
                        call other%mix%add_cell_col(iOld, mix_iCellNew)
                    end if
                    
                else
                    Nrej = Nrej + 1
                end if
         
            else
                Nrej = Nrej + 1                
            end if            
            
        else        
            Nrej = Nrej + 1            
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
        real(DP), dimension(Dim) :: xRand, xTest
        integer :: same_iCellTest, mix_iCellTest
        logical :: overlap
        real(DP) :: enTest, mix_enTest
        
        widTestSum = 0._DP
        
        do iWidom = 1, this%Nwidom           
            
            call random_number(xRand)
            xTest(:) = Lsize(:) * xRand(:)    
            same_iCellTest = this%same%position_to_cell(xTest)
            call this%Epot_neigh(0, xTest, same_iCellTest, overlap) 
            
            if (.not. overlap) then
            
                mix_iCellTest = this%mix%position_to_cell(xTest)
                call mix%Epot_neigh(xTest, mix_iCellTest, this%mix, other_positions, overlap, &
                                    mix_enTest)
                
                if (.not. overlap) then
                
                    enTest = 0._DP + mix_enTest
                    widTestSum = widTestSum + exp(-enTest/Tstar)
                    
                end if
                
            end if
            
        end do
        
        activ = widTestSum/real(this%Nwidom, DP)
        
    end subroutine HardSpheres_widom
    
    !> Total potential energy : dummy
    
    function HardSpheres_Epot_conf(this) result(Epot_conf)
    
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
