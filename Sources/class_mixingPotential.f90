!> \brief Description of the  Mixing Potential class

module class_mixingPotential

use iso_fortran_env
use data_constants
use data_particles
use data_mc
use data_potentiel
use data_neighbours
use mod_physics
use class_neighbours

implicit none

private

    type, public :: MixingPotential

        private
        
        character(len=5) :: name
        
        real(DP) :: rMin !< minimum distance between two particles
        real(DP) :: rCut !< short-range cut
        real(DP) :: dr !< discretisation step
        integer :: iMin !< minimum index of tabulation : minimum distance
        integer :: iCut !< maximum index of tabulation : until potential cut
        real(DP) :: epsilon !< factor in Yukawa
        real(DP) :: alpha !< coefficient in Yukawa
        real(DP), dimension(:), allocatable :: ePot_tab !< tabulation

    contains

        procedure :: construct => MixingPotential_construct
        procedure :: destroy => MixingPotential_destroy

        procedure :: report => MixingPotential_report

        procedure ::getRmin => MixingPotential_getRmin
        
        procedure :: overlapTest => MixingPotential_overlapTest

        procedure :: ePot_init => MixingPotential_ePot_init
        procedure :: ePot_pair => MixingPotential_ePot_pair
        procedure :: ePot_neigh => MixingPotential_ePot_neigh
        procedure :: ePot_total => MixingPotential_ePot_total

    end type

contains

    subroutine MixingPotential_construct(this)

        class(MixingPotential), intent(out) :: this
        
        this%name = "[mix]"
        
        ! Particles
        this%rMin = mix_rMin
        
        ! MixingPotential
        this%rCut = mix_rCut
        this%dr = mix_dr
        this%iMin = int(this%rMin/this%dr)
        this%iCut = int(this%rCut/this%dr)
        this%epsilon = mix_epsilon
        this%alpha = mix_alpha
        allocate(this%ePot_tab(this%iMin:this%iCut))
        call this%ePot_init()

    end subroutine MixingPotential_construct

    subroutine MixingPotential_destroy(this)
    
        class(MixingPotential), intent(inout) :: this
        
        if (allocated(this%ePot_tab)) then
            deallocate(this%ePot_tab)
        end if
    
    end subroutine MixingPotential_destroy
    
    !> Report
    
    subroutine MixingPotential_report(this, unitReport)
    
        class(MixingPotential), intent(in) :: this
        integer, intent(in) :: unitReport    
        
        write(unitReport, *) "Simulation MC_C :"
        write(unitReport, *) "    epsilon = ", this%epsilon
        write(unitReport, *) "    alpha = ", this%alpha
        write(unitReport, *) "    rCut = ", this%rCut
        write(unitReport, *) "    dr = ", this%dr
        
    end subroutine MixingPotential_report
    
    !> Accessor : rMin
    
    function MixingPotential_getRmin(this) result(getRmin)
    
        class(MixingPotential), intent(in) :: this
        
        real(DP) :: getRmin
        
        getRmin = this%rMin
    
    end function MixingPotential_getRmin
    
    !> Overlapt test
    
    subroutine MixingPotential_overlapTest(this, type1_X, type2_X)
    
        class(MixingPotential), intent(in) :: this
        real(DP), dimension(:, :), intent(in) :: type1_X, type2_X
        
        integer :: Ncol1, Ncol2
        integer :: iCol1, iCol2
        real(DP) :: r_mix
        
        Ncol1 = size(type1_X, 2)
        Ncol2 = size(type2_X, 2)
        
        do iCol1 = 1, Ncol1
            do iCol2 = 1, Ncol2
                    
                r_mix = dist(type1_X(:, iCol1), type2_X(:, iCol2))
                if (r_mix < this%rMin) then
                    write(output_unit, *) this%name, " :    Overlap !", &
                        iCol1, iCol2
                    write(output_unit, *) "    r_mix = ", r_mix
                    stop
                end if

            end do
        end do

        write(output_unit, *) this%name, " :    Overlap test : OK !"
    
    end subroutine MixingPotential_overlapTest

    !> MixingPotential energy
    !> Tabulation of Yukawa potential
    !> \f[ \epsilon \frac{e^{-\alpha (r-r_{min})}}{r} \f]
    
    subroutine MixingPotential_ePot_init(this)
    
        class(MixingPotential), intent(inout) :: this

        integer :: i
        real(DP) :: r_i
       
        ! cut
        do i = this%iMin, this%iCut       
            r_i = real(i, DP)*this%dr
            this%ePot_tab(i) = this%epsilon * exp(-this%alpha*(r_i-this%rMin))&
            /r_i
        end do
        
        ! shift        
        this%ePot_tab(:) = this%ePot_tab(:) - this%epsilon * &
            exp(-this%alpha*(this%rCut-this%rMin)) / this%rCut

    end subroutine MixingPotential_ePot_init

    function MixingPotential_ePot_pair(this, r) result(ePot_pair)
        
        class(MixingPotential), intent(in) :: this
        real(DP), intent(in) :: r
        
        integer :: i
        real(DP) :: r_i, ePot_pair
       
        if (r < this%rCut) then
       
            i = int(r/this%dr)
            r_i = real(i, DP)*this%dr
            ePot_pair = this%ePot_tab(i) + (r-r_i)/this%dr * &
                (this%ePot_tab(i+1)-this%ePot_tab(i))
           
        else
       
            ePot_pair = 0._DP
           
        end if
        
    end function MixingPotential_ePot_pair
    
    subroutine MixingPotential_ePot_neigh(this, xCol, iCell, neigh, other_X, &
        overlap, energ)
        
        class(MixingPotential), intent(in) :: this
        real(DP), dimension(:), intent(in) :: xCol !< type A
        integer, intent(in) :: iCell !< type A in mix grid
        type(Neighbours), intent(in) :: neigh
        real(DP), dimension(:, :), intent(in) :: other_X
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energ
    
        integer :: iNeigh,  iCell_neigh
        real(DP) :: r
    
        type(Link), pointer :: current => null(), next => null()
        
        overlap = .false.
        energ = 0._DP
    
        do iNeigh = 1, cell_neighs_nb
        
            iCell_neigh = neigh%cell_neighs(iNeigh, iCell)
            current => neigh%cellsBegin(iCell_neigh)%particle%next            
            if (.not. associated(current%next)) cycle
            
            do
            
                next => current%next
                
                r = dist(xCol(:), other_X(:, current%iCol))
                if (r < this%rMin) then
                    overlap = .true.
                    return
                end if
                energ = energ + this%ePot_pair(r)
                
                if (.not. associated(next%next)) exit
                
                current => next
            
            end do            
            
        end do
    
    end subroutine MixingPotential_ePot_neigh
    
    !> Total potential energy
    
    function MixingPotential_ePot_total(this, type1_X, type2_X) &
        result(ePot_total)
    
        class(MixingPotential), intent(in) :: this
        real(DP), dimension(:, :), intent(in) :: type1_X, type2_X
        
        integer :: Ncol1, Ncol2
        integer :: iCol1, iCol2
        real(DP) :: r_mix
        real(DP) :: ePot_total
        
        Ncol1 = size(type1_X, 2)
        Ncol2 = size(type2_X, 2)
        
        ePot_total = 0._DP
        
        do iCol1 = 1, Ncol1
            do iCol2 = 1, Ncol2
                
                r_mix = dist(type1_X(:, iCol1), type2_X(:, iCol2))
                ePot_total = ePot_total + this%ePot_pair(r_mix)

            end do
        end do
    
    end function MixingPotential_ePot_total

end module class_mixingPotential