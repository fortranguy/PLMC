!> \brief Description of the InteractingSpheres class

module class_interactingSpheres

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
        procedure :: report => InteractingSpheres_report
        
        !> Potential energy
        procedure :: Epot_init => InteractingSpheres_Epot_init
        procedure :: Epot_print => InteractingSpheres_Epot_print
        procedure :: Epot_pair => InteractingSpheres_Epot_pair
        procedure :: Epot_neigh => InteractingSpheres_Epot_neigh
        procedure :: Epot_conf => InteractingSpheres_Epot_conf
        procedure :: consistTest => InteractingSpheres_consistTest
        
        !> Monte-Carlo
        procedure :: move => InteractingSpheres_move
        procedure :: widom => InteractingSpheres_widom
        
    end type InteractingSpheres
    
contains

    subroutine InteractingSpheres_construct(this, shared_rCut)
    
        class(InteractingSpheres), intent(out) :: this
        real(DP), intent(in) :: shared_rCut
        
        this%name = "inter"
    
        ! Particles
        this%radius = inter_radius
        this%rMin = inter_rMin
        this%Ncol = inter_Ncol
        allocate(this%X(Dim, this%Ncol))
        
        ! Monte-Carlo
        this%dx = inter_dx
        this%dx_save = inter_dx
        this%rejFix = inter_rejFix
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
        call this%same%construct(this%rCut)
        call this%same%alloc_cells()
        call this%same%ini_cell_neighs()
        ! Neighbours : other kind
        call this%mix%construct(shared_rCut)
        call this%mix%alloc_cells()
        call this%mix%ini_cell_neighs()
    
    end subroutine InteractingSpheres_construct
    
    subroutine InteractingSpheres_destroy(this)
    
        class(InteractingSpheres), intent(inout) :: this
        
        if (allocated(this%X)) then
            deallocate(this%X)
        end if
        
        if (allocated(this%Epot_tab)) then
            deallocate(this%Epot_tab)
        endif
        
        call this%same%destroy()
        call this%mix%destroy()
    
    end subroutine InteractingSpheres_destroy
    
    !> Report
    
    subroutine InteractingSpheres_report(this, report_unit)
    
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
        
    end subroutine InteractingSpheres_report
    
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

    function InteractingSpheres_Epot_pair(this, r) result(Epot_pair)
        
        class(InteractingSpheres), intent(in) :: this
        real(DP), intent(in) :: r
        
        integer :: i
        real(DP) :: r_i, Epot_pair
       
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
        real(DP) :: r
    
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
                
                    r = dist(xCol(:), this%X(:, current%iCol))
                    if (r < this%rMin) then
                        overlap = .true.
                        return
                    end if
                    energ = energ + this%Epot_pair(r)
       
                end if
                
                if (.not. associated(next%next)) exit
                
                current => next
            
            end do            
            
        end do
    
    end subroutine InteractingSpheres_Epot_neigh
    
    !> Particle move
    
    subroutine InteractingSpheres_move(this, other, mix, same_Epot, mix_Epot, Nrej)
    
        class(InteractingSpheres), intent(inout) :: this
        class(Spheres), intent(inout) :: other
        class(MixingPotential), intent(in) :: mix        
        real(DP), intent(inout) :: same_Epot, mix_Epot
        integer, intent(inout) :: Nrej
        
        logical :: overlap
        integer :: iOld
        real(DP) :: rand
        real(DP), dimension(Dim) :: xRand, xNew
        integer :: same_iCellOld, same_iCellNew
        integer :: mix_iCellOld, mix_iCellNew
        real(DP) :: dEpot
        real(DP) :: same_dEpot, mix_dEpot
        real(DP) :: same_eNew, same_eOld
        real(DP) :: mix_eNew, mix_eOld
        
        call random_number(rand)
        iOld = int(rand*this%Ncol) + 1
        
        call random_number(xRand)
        xNew(:) = this%X(:, iOld) + (xRand(:)-0.5_DP)*this%dx(:)
        xNew(:) = modulo(xNew(:), Lsize(:))
        same_iCellNew = this%same%position_to_cell(xNew)
        call this%Epot_neigh(iOld, xNew, same_iCellNew, overlap, same_eNew)
        
        if (.not. overlap) then
        
            mix_iCellNew = this%mix%position_to_cell(xNew)
            call mix%Epot_neigh(xNew, mix_iCellNew, this%mix, other%X, overlap, mix_eNew)
                        
            if (.not. overlap) then
    
                same_iCellOld = this%same%position_to_cell(this%X(:, iOld))
                call this%Epot_neigh(iOld, this%X(:, iOld), same_iCellOld, overlap, same_eOld)                    
                same_dEpot = same_eNew - same_eOld
                    
                mix_iCellOld = this%mix%position_to_cell(this%X(:, iOld))
                call mix%Epot_neigh(this%X(:, iOld), mix_iCellOld, this%mix, other%X, overlap, mix_eOld)
                mix_dEpot = mix_eNew - mix_eOld
                
                dEpot = same_dEpot + mix_dEpot
                
                call random_number(rand)
                if (rand < exp(-dEpot/Tstar)) then
                
                    this%X(:, iOld) = xNew(:)
                    same_Epot = same_Epot + same_dEpot
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
    
    end subroutine InteractingSpheres_move
    
    !> Widom's method

    subroutine InteractingSpheres_widom(this, activ)
        
        class(InteractingSpheres), intent(in) :: this
        real(DP), intent(inOut) :: activ 
        
        integer :: iWid
        real(DP) :: widTestSum
        real(DP), dimension(Dim) :: xRand, xTest
        integer :: iCellTest
        logical :: overlap        
        real(DP) :: enTest
        
        widTestSum = 0._DP
        
        do iWid = 1, this%Nwidom
            
            call random_number(xRand)
            xTest(:) = Lsize(:) * xRand(:)    
            iCellTest = this%same%position_to_cell(xTest)
            call this%Epot_neigh(0, xTest, iCellTest, overlap, enTest) 
            
            if (.not. overlap) then
                widTestSum = widTestSum + exp(-enTest/Tstar)
            end if
            
        end do
        
        activ = widTestSum/real(this%Nwidom, DP)
        
    end subroutine InteractingSpheres_widom

    !> Total potential energy
    
    function InteractingSpheres_Epot_conf(this) result(Epot_conf)
    
        class(InteractingSpheres), intent(in) :: this
        
        integer :: iCol, jCol
        real(DP) :: r_ij
        real(DP) :: Epot_conf
    
        Epot_conf = 0._DP
        
        do jCol = 1, this%Ncol
            do iCol = 1, this%Ncol
                if (iCol /= jCol) then
                
                    r_ij = dist(this%X(:, iCol), this%X(:, jCol))
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
    
        Epot_conf = this%Epot_conf()
        write(report_unit, *) "Consistency test:"
        write(report_unit, *) "    Epot = ", Epot
        write(report_unit, *) "    Epot_conf = ", Epot_conf
        write(report_unit, *) "    relative difference = ", abs((Epot_conf-Epot)/Epot_conf)
    
    end subroutine InteractingSpheres_consistTest

end module class_interactingSpheres
