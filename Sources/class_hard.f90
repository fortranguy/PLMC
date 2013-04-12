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
public :: inter_constructor

    type, extends(Spheres), public :: Hard
        
    contains

        !> Destructor of the class
        procedure :: destructor => Hard_destructor
        
        !> Print a report of the component in a file
        procedure :: report => Hard_report
              
        !> Potential energy
        procedure :: ePot_neigh => Hard_ePot_neigh
        
        !> Monte-Carlo
        procedure :: mcMove => Hard_mcMove
        procedure :: widom => Hard_widom
        
    end type Hard
    
contains

    function inter_constructor()
    
        type(Hard) :: inter_constructor
    
        ! Construction                

        inter_constructor%radius = inter_radius
        inter_constructor%rMin = inter_rMin
        inter_constructor%Ncol = inter_Ncol
        allocate(inter_constructor%X(Dim, inter_Ncol))
        
        inter_constructor%dx = inter_dx
        
        inter_constructor%rCut = inter_rCut
        
        !	Neighbours        
        inter_constructor%same = neigh_constructor(inter_rCut)
        call inter_constructor%same%alloc_cells()
        call inter_constructor%same%ini_cell_neighs()
    
    end function inter_constructor
    
    subroutine Hard_destructor(this)
    
        class(Hard), intent(inout) :: this
        
        deallocate(this%X)
        call this%same%destructor()
    
    end subroutine Hard_destructor
    
    !> Report
    
    subroutine Hard_report(this, nWidom, unitReport)
    
        class(Hard), intent(in) :: this
        integer, intent(in) :: nWidom
        integer, intent(in) :: unitReport    
        
        write(unitReport, *) "Simulation MC_C :"
        write(unitReport ,*) "    Lsize(:) = ", Lsize(:)
        write(unitReport ,*) "    Vol = ", product(Lsize)
        write(unitReport ,*) "    Ncol = ", this%Ncol
        write(unitReport ,*) "    nWidom = ", nWidom
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
    
    subroutine Hard_mcMove(this, Nrejects)
    
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
    
    end subroutine Hard_mcMove
    
    !> Widom's method

    subroutine Hard_widom(this, nWidom, activExInv)
        
        class(Hard), intent(in) :: this
        integer, intent(in) :: nWidom
        real(DP), intent(inOut) :: activExInv 
        
        integer :: iWid
        real(DP) :: widTestSum
        real(DP), dimension(Dim) :: xRand, xTest
        integer :: iCellTest
        logical :: overlap        
        
        widTestSum = 0._DP
        
        do iWid = 1, nWidom           
            
            call random_number(xRand)
            xTest(:) = Lsize(:) * xRand(:)    
            iCellTest = this%same%position_to_cell(xTest)
            call this%ePot_neigh(0, xTest, iCellTest, overlap) 
            
            if (.not. overlap) then
                widTestSum = widTestSum + 1._DP
            end if
            
        end do
        
        activExInv = widTestSum/real(nWidom, DP)
        
    end subroutine Hard_widom

end module class_hard
