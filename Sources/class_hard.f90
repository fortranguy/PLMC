!> \brief Description of the Hard Sphere class

module class_hard

use data_constants
use data_cell
use data_particles
use data_potentiel
use data_mc
use data_neighbours
use class_neighbours
use mod_physics
use class_spheres

implicit none

private

    type, extends(Spheres), public :: Hard
        
    contains

        !> Construction and destruction of the class
        procedure :: construct => Hard_construct
        procedure :: destroy => Hard_destroy
        
        !> Print a report of the component in a file
        procedure :: report => Hard_report
              
        !> Potential energy
        procedure :: ePot_neigh => Hard_ePot_neigh
        
        !> Monte-Carlo
        procedure :: move => Hard_move
        procedure :: widom => Hard_widom
        
    end type Hard
    
contains

    subroutine Hard_construct(this)
    
        class(Hard), intent(out) :: this
    
        ! Particles
        this%radius = hard_radius
        this%rMin = hard_rMin
        this%Ncol = hard_Ncol
        allocate(this%X(Dim, hard_Ncol))
        
        ! Monte-Carlo
        this%dx = hard_dx
        this%Nwidom = hard_Nwidom
                
        ! Potential
        this%rCut = hard_rCut
        
        ! Neighbours        
        call this%same%construct(hard_rCut)
        call this%same%alloc_cells()
        call this%same%ini_cell_neighs()
    
    end subroutine Hard_construct
    
    subroutine Hard_destroy(this)
    
        class(Hard), intent(inout) :: this
        
        deallocate(this%X)
        call this%same%destroy()
    
    end subroutine Hard_destroy
    
    !> Report
    
    subroutine Hard_report(this, unitReport)
    
        class(Hard), intent(in) :: this
        integer, intent(in) :: unitReport    
        
        write(unitReport, *) "Simulation MC_C :"
        write(unitReport ,*) "    Lsize(:) = ", Lsize(:)
        write(unitReport ,*) "    Vol = ", product(Lsize)
        write(unitReport ,*) "    Ncol = ", this%Ncol
        write(unitReport ,*) "    Nwidom = ", this%Nwidom
        write(unitReport, *) "    Nstep = ", Nstep
        write(unitReport, *) "    Ntherm = ", Ntherm
        write(unitReport, *) "    Nmove = ", Nmove
        write(unitReport, *) "    rCut = ", this%rCut
        write(unitReport, *) "    cell_coordMax(:) = ", &
        	this%same%cell_coordMax(:)
        write(unitReport, *) "    cell_Lsize(:) = ", this%same%cell_Lsize(:)
        
    end subroutine Hard_report
    
    subroutine Hard_ePot_neigh(this, iCol, xCol, iCell, overlap)
        
        class(Hard), intent(in) :: this        
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
    
    end subroutine Hard_ePot_neigh
    
    !> Particle move
    
    subroutine Hard_move(this, Nrejects)
    
        class(Hard), intent(inout) :: this
        integer, intent(inout) :: Nrejects
        
        logical :: overlap
        integer :: iOld
        real(DP) :: rand
        real(DP), dimension(Dim) :: xRand, xNew
        integer :: iCellBefore, iCellAfter
        
        call random_number(rand)
        iOld = int(rand*this%Ncol) + 1
        
        call random_number(xRand)
        xNew(:) = this%X(:, iOld) + (xRand(:)-0.5_DP)*this%dx(:)
        xNew(:) = modulo(xNew(:), Lsize(:))
        iCellAfter = this%same%position_to_cell(xNew)
        call this%ePot_neigh(iOld, xNew, iCellAfter, overlap)
        
        if (.not. overlap) then
        
            iCellBefore = this%same%position_to_cell(this%X(:, iOld))
            call this%ePot_neigh(iOld, this%X(:, iOld), iCellBefore, overlap)
        
            this%X(:, iOld) = xNew(:)
                
            if ( iCellBefore /= iCellAfter ) then                
                call this%same%remove_cell_col(iOld, iCellBefore)
                call this%same%add_cell_col(iOld, iCellAfter)
            end if
            
        else
        
            Nrejects = Nrejects + 1
            
        end if
    
    end subroutine Hard_move
    
    !> Widom's method

    subroutine Hard_widom(this, activExInv)
        
        class(Hard), intent(in) :: this
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
        
    end subroutine Hard_widom

end module class_hard
