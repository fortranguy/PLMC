!> \brief Description of the InteractingSpheres class

module class_interactingSpheres

use data_precisions, only : DP, consist_tiny
use data_cell, only : Ndim, Lsize
use data_particles, only : inter_radius, inter_rMin, inter_Ncol
use data_potential, only : inter_rCut, inter_dr, inter_epsilon, inter_alpha
use data_mc, only : Temperature, inter_deltaX, inter_rejectFix, inter_Nadapt, inter_Nwidom
use data_neighbours, only : cell_neighs_nb, inter_cell_Lsize
use data_distrib, only : inter_snap_factor
use mod_physics, only : dist
use class_observables
use class_neighbourCells
use class_mixingPotential
use class_spheres

implicit none

private

    type, extends(Spheres), public :: InteractingSpheres

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
        procedure :: PrintReport => InteractingSpheres_printReport
        
        !> Potential energy
        procedure, private :: Epot_init => InteractingSpheres_Epot_init
        procedure :: Epot_print => InteractingSpheres_Epot_print
        procedure :: Epot_pair => InteractingSpheres_Epot_pair
        procedure, private :: Epot_neigh => InteractingSpheres_Epot_neigh
        procedure :: Epot_conf => InteractingSpheres_Epot_conf
        procedure :: consistTest => InteractingSpheres_consistTest
        
        !> Monte-Carlo
        procedure :: move => InteractingSpheres_move
        procedure :: widom => InteractingSpheres_widom
        
    end type InteractingSpheres
    
contains

    subroutine InteractingSpheres_construct(this, mix_cell_Lsize, mix_rCut)
    
        class(InteractingSpheres), intent(out) :: this
        real(DP), dimension(:), intent(in) :: mix_cell_Lsize
        real(DP), intent(in) :: mix_rCut
        
        this%name = "inter"
    
        ! Particles
        this%radius = inter_radius
        this%rMin = inter_rMin
        this%Ncol = inter_Ncol
        allocate(this%positions(Ndim, this%Ncol))
        
        ! Snapshot
        this%snap_factor = inter_snap_factor
        
        ! Monte-Carlo
        this%deltaX = inter_deltaX
        this%deltaXSave = this%deltaX
        this%rejectFix = inter_rejectFix
        this%Nadapt = inter_Nadapt
        this%Nwidom = inter_Nwidom
        
        ! Potential
        this%rCut = inter_rCut
        this%dr = inter_dr
        this%iMin = int(this%rMin/this%dr)
        this%iCut = int(this%rCut/this%dr)
        this%epsilon = inter_epsilon
        this%alpha = inter_alpha        
        allocate(this%Epot_tab(this%iMin:this%iCut))
        call this%Epot_init()
        
        ! Neighbours : same kind
        call this%same%construct(inter_cell_Lsize, this%rCut)
        call this%same%alloc_cells()
        call this%same%cell_neighs_init()
        ! Neighbours : other kind
        call this%mix%construct(mix_cell_Lsize, mix_rCut)
        call this%mix%alloc_cells()
        call this%mix%cell_neighs_init()
    
    end subroutine InteractingSpheres_construct
    
    subroutine InteractingSpheres_destroy(this)
    
        class(InteractingSpheres), intent(inout) :: this
        
        if (allocated(this%positions)) then
            deallocate(this%positions)
        end if
        
        if (allocated(this%Epot_tab)) then
            deallocate(this%Epot_tab)
        endif
        
        call this%same%destroy()
        call this%mix%destroy()
    
    end subroutine InteractingSpheres_destroy
    
    !> Report
    
    subroutine InteractingSpheres_printReport(this, report_unit)
    
        class(InteractingSpheres), intent(in) :: this
        integer, intent(in) :: report_unit    
        
        write(report_unit, *) "Data :"
        
        write(report_unit ,*) "    Ncol = ", this%Ncol
        write(report_unit ,*) "    Nwidom = ", this%Nwidom
        write(report_unit ,*) "    Nadapt = ", this%Nadapt
        
        write(report_unit, *) "    epsilon = ", this%epsilon
        write(report_unit, *) "    alpha = ", this%alpha
        write(report_unit, *) "    rCut = ", this%rCut
        write(report_unit, *) "    dr = ", this%dr
        
        write(report_unit, *) "    same_cell_coordMax(:) = ", this%same%cell_coordMax(:)
        write(report_unit, *) "    same_cell_Lsize(:) = ", this%same%cell_Lsize(:)
        write(report_unit, *) "    mix_cell_coordMax(:) = ", this%mix%cell_coordMax(:)
        write(report_unit, *) "    mix_cell_Lsize(:) = ", this%mix%cell_Lsize(:)
        
    end subroutine InteractingSpheres_printReport
    
    !> Potential energy
    !> Tabulation of Yukawa potential    
    !> \f[ \epsilon \frac{e^{-\alpha (r-r_{min})}}{r} \f]
    
    subroutine InteractingSpheres_Epot_init(this)
    
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
    
    subroutine InteractingSpheres_Epot_neigh(this, iCol, xCol, iCell, overlap, energ)
        
        class(InteractingSpheres), intent(in) :: this        
        integer, intent(in) :: iCol, iCell
        real(DP), dimension(:), intent(in) :: xCol
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energ
    
        integer :: iNeigh,  iCell_neigh
        real(DP) :: r_ij
    
        type(Link), pointer :: current => null(), next => null()
        
        overlap = .false.
        energ = 0._DP
    
        do iNeigh = 1, cell_neighs_nb
        
            iCell_neigh = this%same%cell_neighs(iNeigh, iCell)
            current => this%same%cellsBegin(iCell_neigh)%particle%next            
            if (.not. associated(current%next)) cycle
            
            do
            
                next => current%next
            
                if (current%iCol /= iCol) then
                
                    r_ij = dist(xCol(:), this%positions(:, current%iCol))
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
    
    end subroutine InteractingSpheres_Epot_neigh
    
    !> Particle move
    
    subroutine InteractingSpheres_move(this, iOld, other, mix, same_obs, mix_Epot)
    
        class(InteractingSpheres), intent(inout) :: this
        integer, intent(in) :: iOld
        class(Spheres), intent(inout) :: other
        class(MixingPotential), intent(in) :: mix
        class(Observables) :: same_obs
        real(DP), intent(inout) :: mix_Epot
        
        logical :: overlap
        real(DP) :: rand
        real(DP), dimension(Ndim) :: xRand, xNew
        integer :: same_iCellOld, same_iCellNew
        integer :: mix_iCellOld, mix_iCellNew
        real(DP) :: deltaEpot
        real(DP) :: same_deltaEpot, mix_deltaEpot
        real(DP) :: same_eNew, same_eOld
        real(DP) :: mix_eNew, mix_eOld
        
        call random_number(xRand)
        xNew(:) = this%positions(:, iOld) + (xRand(:)-0.5_DP)*this%deltaX(:)
        xNew(:) = modulo(xNew(:), Lsize(:))
        same_iCellNew = this%same%position_to_cell(xNew)
        call this%Epot_neigh(iOld, xNew, same_iCellNew, overlap, same_eNew)
        
        if (.not. overlap) then
        
            mix_iCellNew = this%mix%position_to_cell(xNew)
            call mix%Epot_neigh(xNew, mix_iCellNew, this%mix, other%positions, overlap, mix_eNew)
                        
            if (.not. overlap) then
    
                same_iCellOld = this%same%position_to_cell(this%positions(:, iOld))
                call this%Epot_neigh(iOld, this%positions(:, iOld), same_iCellOld, overlap, same_eOld)
                same_deltaEpot = same_eNew - same_eOld
                    
                mix_iCellOld = this%mix%position_to_cell(this%positions(:, iOld))
                call mix%Epot_neigh(this%positions(:, iOld), mix_iCellOld, this%mix, other%positions, &
                                    overlap, mix_eOld)
                mix_deltaEpot = mix_eNew - mix_eOld
                
                deltaEpot = same_deltaEpot + mix_deltaEpot
                
                call random_number(rand)
                if (rand < exp(-deltaEpot/Temperature)) then
                
                    this%positions(:, iOld) = xNew(:)
                    same_obs%Epot = same_obs%Epot + same_deltaEpot
                    mix_Epot = mix_Epot + mix_deltaEpot
                    
                    if (same_iCellOld /= same_iCellNew) then
                        call this%same%remove_col_from_cell(iOld, same_iCellOld)
                        call this%same%add_col_to_cell(iOld, same_iCellNew)
                    end if
                    
                    if (mix_iCellOld /= mix_iCellNew) then
                        call other%mix%remove_col_from_cell(iOld, mix_iCellOld)
                        call other%mix%add_col_to_cell(iOld, mix_iCellNew)
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
    
    end subroutine InteractingSpheres_move
    
    !> Widom's method : with other type ?

    subroutine InteractingSpheres_widom(this, other_positions, mix, activ)
        
        class(InteractingSpheres), intent(in) :: this
        real(DP), dimension(:, :), intent(in) :: other_positions
        class(MixingPotential), intent(in) :: mix
        real(DP), intent(inOut) :: activ
        
        integer :: iWidom
        real(DP) :: widTestSum
        real(DP), dimension(Ndim) :: xRand, xTest
        integer :: same_iCellTest, mix_iCellTest
        logical :: overlap        
        real(DP) :: enTest, same_enTest, mix_enTest
        
        widTestSum = 0._DP
        
        do iWidom = 1, this%Nwidom
            
            call random_number(xRand)
            xTest(:) = Lsize(:) * xRand(:)    
            same_iCellTest = this%same%position_to_cell(xTest)
            call this%Epot_neigh(0, xTest, same_iCellTest, overlap, same_enTest) 
            
            if (.not. overlap) then
                
                mix_iCellTest = this%mix%position_to_cell(xTest)
                call mix%Epot_neigh(xTest, mix_iCellTest, this%mix, other_positions, overlap, &
                                    mix_enTest)
                
                if (.not. overlap) then
                
                    enTest = same_enTest + mix_enTest
                    widTestSum = widTestSum + exp(-enTest/Temperature)
                    
                end if
                
            end if
            
        end do
        
        activ = widTestSum/real(this%Nwidom, DP)
        
    end subroutine InteractingSpheres_widom

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
                
                    r_ij = dist(this%positions(:, iCol), this%positions(:, jCol))
                    Epot_conf = Epot_conf + this%Epot_pair(r_ij)
                    
                end if
            end do
        end do
        
        Epot_conf = 0.5_DP*Epot_conf
    
    end function InteractingSpheres_Epot_conf
    
    !> Consistency test 
    
    subroutine InteractingSpheres_consistTest(this, Epot, report_unit)
    
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
    
    end subroutine InteractingSpheres_consistTest

end module class_interactingSpheres
