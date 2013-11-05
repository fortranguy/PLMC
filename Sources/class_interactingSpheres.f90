!> \brief Description of the InteractingSpheres class

module class_interactingSpheres

use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
use data_precisions, only : DP, consist_tiny
use data_box, only : Ndim, Lsize
use data_particles, only : inter_rMin, inter_Ncol
use data_potential, only : inter_rCut, inter_dr, inter_epsilon, inter_alpha
use data_monteCarlo, only : inter_move_delta, inter_move_rejectFix, inter_Nwidom
use data_neighbourCells, only : NnearCell
use data_distribution, only : inter_snap_factor
use module_physics, only : set_discrete_length, dist_PBC
use class_neighbourCells
use class_hardSpheres

implicit none

private

    type, extends(HardSpheres), public :: InteractingSpheres

        private

        ! Potential :
        real(DP)  :: dr !< discretisation step
        integer :: iMin !< minimum index of tabulation : minimum distance
        integer :: iCut !< maximum index of tabulation : until potential cut
        real(DP) :: epsilon !< factor in Yukawa
        real(DP) :: alpha !< coefficient in Yukawa
        real(DP), dimension(:), allocatable :: Epot_tab !< tabulation
        
    contains

        !> Construction and destruction of the class
        procedure :: construct => InteractingSpheres_construct
        procedure :: destroy => InteractingSpheres_destroy
        
        !> Print a report of the component in a file
        procedure :: print_report => InteractingSpheres_print_report
        
        !> Potential energy
        procedure, private :: Epot_init => InteractingSpheres_Epot_init
        procedure :: Epot_print => InteractingSpheres_Epot_print
        procedure :: Epot_pair => InteractingSpheres_Epot_pair
        procedure :: Epot_neighCells => InteractingSpheres_Epot_neighCells
        procedure :: Epot_conf => InteractingSpheres_Epot_conf
        procedure :: test_consist => InteractingSpheres_test_consist

    end type InteractingSpheres
    
contains

    subroutine InteractingSpheres_construct(this)
    
        class(InteractingSpheres), intent(out) :: this
        
        real(DP), dimension(Ndim) :: cell_size
        
        this%name = "inter"
        write(output_unit, *) this%name, " class construction"
    
        ! Particles
        this%rMin = inter_rMin
        this%radius = this%rMin/2._DP
        this%Ncol = inter_Ncol
        allocate(this%positions(Ndim, this%Ncol))
        
        ! Snapshot
        this%snap_factor = inter_snap_factor
        
        ! Monte-Carlo
        this%move_delta = inter_move_delta
        this%move_deltaSave = this%move_delta
        this%move_rejectFix = inter_move_rejectFix
        this%Nwidom = inter_Nwidom
        
        ! Potential
        this%rCut = inter_rCut
        this%dr = inter_dr
        call set_discrete_length(this%rMin, this%dr)
        this%iMin = int(this%rMin/this%dr)
        this%iCut = int(this%rCut/this%dr) + 1
        this%epsilon = inter_epsilon
        this%alpha = inter_alpha        
        allocate(this%Epot_tab(this%iMin:this%iCut))
        call this%Epot_init()
        
        ! Neighbour Cells
        cell_size(:) = this%rCut
        call this%sameCells%construct(cell_size, this%rCut)
    
    end subroutine InteractingSpheres_construct
    
    subroutine InteractingSpheres_destroy(this)
    
        class(InteractingSpheres), intent(inout) :: this
        
        write(output_unit, *) this%name, " class destruction"
        
        if (allocated(this%positions)) then
            deallocate(this%positions)
        end if
        
        if (allocated(this%Epot_tab)) then
            deallocate(this%Epot_tab)
        endif
        
        call this%sameCells%destroy()
        call this%mixCells%destroy()
    
    end subroutine InteractingSpheres_destroy
    
    !> Report
    
    subroutine InteractingSpheres_print_report(this, report_unit)
    
        class(InteractingSpheres), intent(in) :: this
        integer, intent(in) :: report_unit
        
        call this%HardSpheres%print_report(report_unit)
        
        write(report_unit, *) "    epsilon = ", this%epsilon
        write(report_unit, *) "    alpha = ", this%alpha
        write(report_unit, *) "    dr = ", this%dr
        
    end subroutine InteractingSpheres_print_report
    
    !> Potential energy
    !> Tabulation of Yukawa potential    
    !> \f[ \epsilon \frac{e^{-\alpha (r-r_{min})}}{r} \f]
    
    pure subroutine InteractingSpheres_Epot_init(this)
    
        class(InteractingSpheres), intent(inout) :: this

        integer :: i
        real(DP) :: r_i
       
        ! cut
        do i = this%iMin, this%iCut       
            r_i = real(i, DP)*this%dr
            this%Epot_tab(i) = this%epsilon * exp(-this%alpha*(r_i-this%rMin)) / r_i
        end do
        
        ! shift        
        this%Epot_tab(:) = this%Epot_tab(:) - this%Epot_tab(this%iCut)

    end subroutine InteractingSpheres_Epot_init
    
    !> Print the tabulated potential
    
    subroutine InteractingSpheres_Epot_print(this, Epot_unit)

        class(InteractingSpheres), intent(in) :: this
        integer, intent(in) :: Epot_unit

        integer :: i
        real(DP) :: r_i

        do i = this%iMin, this%iCut
            r_i = real(i, DP)*this%dr
            write(Epot_unit, *) r_i, this%Epot_tab(i)
        end do

    end subroutine InteractingSpheres_Epot_print

    pure function InteractingSpheres_Epot_pair(this, r) result(Epot_pair)
        
        class(InteractingSpheres), intent(in) :: this
        real(DP), intent(in) :: r
        real(DP) :: Epot_pair
        
        integer :: i
        real(DP) :: r_i
       
        if (r < this%rCut) then
       
            i = int(r/this%dr)
            r_i = real(i, DP)*this%dr
            Epot_pair = this%Epot_tab(i) + (r-r_i)/this%dr * (this%Epot_tab(i+1)-this%Epot_tab(i))
           
        else
       
            Epot_pair = 0._DP
           
        end if
        
    end function InteractingSpheres_Epot_pair
    
    subroutine InteractingSpheres_Epot_neighCells(this, iCol, xCol, iTotalCell, overlap, energ)
        
        class(InteractingSpheres), intent(in) :: this        
        integer, intent(in) :: iCol, iTotalCell
        real(DP), dimension(:), intent(in) :: xCol
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energ
    
        integer :: iNearCell,  nearCell_index
        real(DP) :: r_ij
    
        type(Node), pointer :: current => null(), next => null()
        
        overlap = .false.
        energ = 0._DP
    
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
                    energ = energ + this%Epot_pair(r_ij)
       
                end if
                
                if (.not. associated(next%next)) exit
                
                current => next
            
            end do            
            
        end do
    
    end subroutine InteractingSpheres_Epot_neighCells

    !> Total potential energy
    
    pure function InteractingSpheres_Epot_conf(this) result(Epot_conf)
    
        class(InteractingSpheres), intent(in) :: this
        real(DP) :: Epot_conf
        
        integer :: iCol, jCol
        real(DP) :: r_ij
    
        Epot_conf = 0._DP
        
        do jCol = 1, this%Ncol
            do iCol = 1, this%Ncol
                if (iCol /= jCol) then
                
                    r_ij = dist_PBC(this%positions(:, iCol), this%positions(:, jCol))
                    Epot_conf = Epot_conf + this%Epot_pair(r_ij)
                    
                end if
            end do
        end do
        
        Epot_conf = 0.5_DP*Epot_conf
    
    end function InteractingSpheres_Epot_conf
    
    !> Consistency test 
    
    subroutine InteractingSpheres_test_consist(this, Epot, report_unit)
    
        class(InteractingSpheres), intent(in) :: this
        real(DP), intent(in) :: Epot
        integer, intent(in) :: report_unit
        
        real(DP) :: Epot_conf
        real(DP) :: difference
        
        Epot_conf = this%Epot_conf()
        difference = abs((Epot_conf-Epot)/Epot_conf)
        
        write(report_unit, *) "Consistency test:"
        write(report_unit, *) "    Epot = ", Epot
        write(report_unit, *) "    Epot_conf = ", Epot_conf
        write(report_unit, *) "    relative difference = ", difference
        
        if (difference > consist_tiny) then
            write(report_unit, *) "    WARNING !"
        else
            write(report_unit, *) "    OK !"
        end if
    
    end subroutine InteractingSpheres_test_consist

end module class_interactingSpheres
