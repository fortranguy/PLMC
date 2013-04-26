!> \brief Description of the Hard Spheres class

module class_hardSpheres

use data_constants
use data_cell
use data_particles
use data_potentiel
use data_mc
use data_neighbours
use mod_physics
use class_neighbours
use class_mixingPotential
use class_spheres

implicit none

private

    type, extends(Spheres), public :: HardSpheres
    
        ! Potential
        real(DP) :: ePot
        
    contains

        !> Construction and destruction of the class
        procedure :: construct => HardSpheres_construct
        procedure :: destroy => HardSpheres_destroy
        
        !> Print a report of the component in a file
        procedure :: report => HardSpheres_report
              
        !> Potential energy
        procedure :: ePot_neigh => HardSpheres_ePot_neigh
        procedure :: ePot_conf => HardSpheres_ePot_conf
        procedure :: consistTest => HardSpheres_consistTest
        
        !> Monte-Carlo
        procedure :: move => HardSpheres_move
        procedure :: widom => HardSpheres_widom
        
    end type HardSpheres
    
contains

    subroutine HardSpheres_construct(this)
    
        class(HardSpheres), intent(out) :: this
        
        this%name = "hardS"
    
        ! Particles
        this%radius = hard_radius
        this%rMin = hard_rMin
        this%Ncol = hard_Ncol
        allocate(this%X(Dim, this%Ncol))
        
        ! Monte-Carlo
        this%dx = hard_dx
        this%dx_save = hard_dx
        this%Nadapt = hard_Nadapt
        this%Nwidom = hard_Nwidom
                
        ! Potential
        this%ePot = 0._DP
        this%rCut = hard_rCut
        
        ! Neighbours : same kind
        call this%same%construct(this%rCut)
        call this%same%alloc_cells()
        call this%same%ini_cell_neighs()
        ! Neighbours : other kind
        call this%mix%construct(mix_rCut)
        call this%mix%alloc_cells()
        call this%mix%ini_cell_neighs()
    
    end subroutine HardSpheres_construct
    
    subroutine HardSpheres_destroy(this)
    
        class(HardSpheres), intent(inout) :: this
        
        if (allocated(this%X)) then
            deallocate(this%X)
        end if
        
        call this%same%destroy()
        call this%mix%destroy()
    
    end subroutine HardSpheres_destroy
    
    !> Report
    
    subroutine HardSpheres_report(this, report_unit)
    
        class(HardSpheres), intent(in) :: this
        integer, intent(in) :: report_unit    
        
        write(report_unit, *) "Simulation MC_C :"
        write(report_unit ,*) "    Ncol = ", this%Ncol
        write(report_unit ,*) "    Nwidom = ", this%Nwidom
        write(report_unit, *) "    rCut = ", this%rCut
        write(report_unit, *) "    cell_coordMax(:) = ", &
        	this%same%cell_coordMax(:)
        write(report_unit, *) "    cell_Lsize(:) = ", this%same%cell_Lsize(:)
        
    end subroutine HardSpheres_report
    
    subroutine HardSpheres_ePot_neigh(this, iCol, xCol, iCell, overlap)
        
        class(HardSpheres), intent(in) :: this        
        integer, intent(in) :: iCol, iCell
        real(DP), dimension(Dim), intent(in) :: xCol
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
                
                    r = dist(xCol(:), this%X(:, current%iCol))
                    if (r < this%rMin) then
                        overlap = .true.
                        return
                    end if
       
                end if
                
                if (.not. associated(next%next)) exit
                
                current => next
            
            end do            
            
        end do
    
    end subroutine HardSpheres_ePot_neigh
    
    !> Particle move
    
    subroutine HardSpheres_move(this, other, mix, same_ePot, mix_ePot, Nrej)
    
        class(HardSpheres), intent(inout) :: this
        class(Spheres), intent(inout) :: other
        class(MixingPotential), intent(in) :: mix
        real(DP), intent(inout) :: same_ePot, mix_ePot
        integer, intent(inout) :: Nrej
        
        logical :: overlap
        integer :: iOld
        real(DP) :: rand
        real(DP), dimension(Dim) :: xRand, xNew
        integer :: same_iCellOld, same_iCellNew
        integer :: mix_iCellOld, mix_iCellNew
        real(DP) :: mix_dEpot
        real(DP) :: mix_eNew, mix_eOld
        
        call random_number(rand)
        iOld = int(rand*this%Ncol) + 1
        
        call random_number(xRand)
        xNew(:) = this%X(:, iOld) + (xRand(:)-0.5_DP)*this%dx(:)
        xNew(:) = modulo(xNew(:), Lsize(:))
        same_iCellNew = this%same%position_to_cell(xNew)
        call this%ePot_neigh(iOld, xNew, same_iCellNew, overlap)
        
        if (.not. overlap) then
        
            mix_iCellNew = this%mix%position_to_cell(xNew)
            call mix%ePot_neigh(xNew, mix_iCellNew, this%mix, other%X, &
                overlap, mix_eNew)
                        
            if (.not. overlap) then
    
                same_iCellOld = this%same%position_to_cell(this%X(:, iOld))
                call this%ePot_neigh(iOld, this%X(:, iOld), same_iCellOld, &
                    overlap)
                    
                mix_iCellOld = this%mix%position_to_cell(this%X(:, iOld))
                call mix%ePot_neigh(this%X(:, iOld), mix_iCellOld, &
                    this%mix, other%X, overlap, mix_eOld)
                
                mix_dEpot = mix_eNew - mix_eOld
                    
                call random_number(rand)
                if (rand < exp(-mix_dEpot/Tstar)) then
                
                    this%X(:, iOld) = xNew(:)
                    same_ePot = same_ePot + 0._DP
                    mix_ePot = mix_ePot + mix_dEpot
                    
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
    
    !> Widom's method

    subroutine HardSpheres_widom(this, activ)
        
        class(HardSpheres), intent(in) :: this
        real(DP), intent(inOut) :: activ 
        
        integer :: iWid
        real(DP) :: widTestSum
        real(DP), dimension(Dim) :: xRand, xTest
        integer :: iCellTest
        logical :: overlap        
        
        widTestSum = 0._DP
        
        do iWid = 1, this%Nwidom           
            
            call random_number(xRand)
            xTest(:) = Lsize(:) * xRand(:)    
            iCellTest = this%same%position_to_cell(xTest)
            call this%ePot_neigh(0, xTest, iCellTest, overlap) 
            
            if (.not. overlap) then
                widTestSum = widTestSum + 1._DP
            end if
            
        end do
        
        activ = widTestSum/real(this%Nwidom, DP)
        
    end subroutine HardSpheres_widom
    
    !> Total potential energy : dummy
    
    function HardSpheres_ePot_conf(this) result(ePot_conf)
    
        class(HardSpheres), intent(in) :: this
        
        real(DP) :: ePot_conf
    
        ePot_conf = this%ePot
        
    end function HardSpheres_ePot_conf
    
    !> Consistency test : dummy
    
    subroutine HardSpheres_consistTest(this, ePot, report_unit)
    
        class(HardSpheres), intent(in) :: this
        real(DP), intent(in) :: ePot
        integer, intent(in) :: report_unit
        
        real(DP) :: ePot_conf
    
        ePot_conf = this%ePot_conf()
        write(report_unit, *) "Consistency test:"
        write(report_unit, *) "    ePot = ", ePot
        write(report_unit, *) "    ePot_conf = ", ePot_conf
        write(report_unit, *) "    absolute difference = ", &
            abs(ePot_conf-ePot)
    
    end subroutine HardSpheres_consistTest

end module class_hardSpheres
