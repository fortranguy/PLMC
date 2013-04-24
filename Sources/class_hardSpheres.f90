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
        
    contains

        !> Construction and destruction of the class
        procedure :: construct => HardSpheres_construct
        procedure :: destroy => HardSpheres_destroy
        
        !> Print a report of the component in a file
        procedure :: report => HardSpheres_report
              
        !> Potential energy
        procedure :: ePot_neigh => HardSpheres_ePot_neigh
        
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
        this%rCut = hard_rCut
        
        ! Neighbours : same kind
        call this%same%construct(this%rCut)
        call this%same%alloc_cells()
        call this%same%ini_cell_neighs()
        ! Neighbours : other kind
        call this%other%construct(mix_rCut)
        call this%other%alloc_cells()
        call this%other%ini_cell_neighs()
    
    end subroutine HardSpheres_construct
    
    subroutine HardSpheres_destroy(this)
    
        class(HardSpheres), intent(inout) :: this
        
        if (allocated(this%X)) then
            deallocate(this%X)
        end if
        
        call this%same%destroy()
        call this%other%destroy()
    
    end subroutine HardSpheres_destroy
    
    !> Report
    
    subroutine HardSpheres_report(this, unitReport)
    
        class(HardSpheres), intent(in) :: this
        integer, intent(in) :: unitReport    
        
        write(unitReport, *) "Simulation MC_C :"
        write(unitReport ,*) "    Ncol = ", this%Ncol
        write(unitReport ,*) "    Nwidom = ", this%Nwidom
        write(unitReport, *) "    rCut = ", this%rCut
        write(unitReport, *) "    cell_coordMax(:) = ", &
        	this%same%cell_coordMax(:)
        write(unitReport, *) "    cell_Lsize(:) = ", this%same%cell_Lsize(:)
        
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
    
    subroutine HardSpheres_move(this, mix, other_X, Nrej)
    
        class(HardSpheres), intent(inout) :: this
        class(MixingPotential), intent(in) :: mix
        real(DP), dimension(:, :) :: other_X
        integer, intent(inout) :: Nrej
        
        logical :: overlap
        integer :: iOld
        real(DP) :: rand
        real(DP), dimension(Dim) :: xRand, xNew
        integer :: same_iCellOld, same_iCellNew
        integer :: other_iCellOld, other_iCellNew
        real(DP) :: other_eNew, other_eOld
        
        call random_number(rand)
        iOld = int(rand*this%Ncol) + 1
        
        call random_number(xRand)
        xNew(:) = this%X(:, iOld) + (xRand(:)-0.5_DP)*this%dx(:)
        xNew(:) = modulo(xNew(:), Lsize(:))
        same_iCellNew = this%same%position_to_cell(xNew)
        call this%ePot_neigh(iOld, xNew, same_iCellNew, overlap)
        
        if (.not. overlap) then
        
            other_iCellNew = this%other%position_to_cell(xNew)
            call mix%ePot_neigh(xNew, other_iCellNew, this%other, other_X, &
                overlap, other_eNew)
        
            if (.not. overlap) then
        
                same_iCellOld = this%same%position_to_cell(this%X(:, iOld))
                call this%ePot_neigh(iOld, this%X(:, iOld), same_iCellOld, &
                    overlap)
                    
                other_iCellOld = this%other%position_to_cell(this%X(:, iOld))
                call mix%ePot_neigh(this%X(:, iOld), other_iCellOld, &
                    this%other, other_X, overlap, other_eOld)
            
                this%X(:, iOld) = xNew(:)
                    
                if ( same_iCellOld /= same_iCellNew ) then                
                    call this%same%remove_cell_col(iOld, same_iCellOld)
                    call this%same%add_cell_col(iOld, same_iCellNew)
                end if
                
            else
                Nrej = Nrej + 1
            end if
            
        else
        
            Nrej = Nrej + 1
            
        end if
    
    end subroutine HardSpheres_move
    
    !> Widom's method

    subroutine HardSpheres_widom(this, activExInv)
        
        class(HardSpheres), intent(in) :: this
        real(DP), intent(inOut) :: activExInv 
        
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
        
        activExInv = widTestSum/real(this%Nwidom, DP)
        
    end subroutine HardSpheres_widom

end module class_hardSpheres
