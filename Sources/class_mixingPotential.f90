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
        real(DP), dimension(:), allocatable :: Epot_tab !< tabulation

    contains

        procedure :: construct => MixingPotential_construct
        procedure :: destroy => MixingPotential_destroy

        procedure :: report => MixingPotential_report

        procedure :: getRmin => MixingPotential_getRmin
        procedure :: getRcut => MixingPotential_getRcut
        
        procedure :: overlapTest => MixingPotential_overlapTest

        procedure :: Epot_init => MixingPotential_Epot_init
        procedure :: Epot_print => MixingPotential_Epot_print
        procedure :: Epot_pair => MixingPotential_Epot_pair
        procedure :: Epot_neigh => MixingPotential_Epot_neigh
        procedure :: Epot_conf => MixingPotential_Epot_conf

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
        allocate(this%Epot_tab(this%iMin:this%iCut))
        call this%Epot_init()

    end subroutine MixingPotential_construct

    subroutine MixingPotential_destroy(this)
    
        class(MixingPotential), intent(inout) :: this
        
        if (allocated(this%Epot_tab)) then
            deallocate(this%Epot_tab)
        end if
    
    end subroutine MixingPotential_destroy
    
    !> Report
    
    subroutine MixingPotential_report(this, report_unit)
    
        class(MixingPotential), intent(in) :: this
        integer, intent(in) :: report_unit    
        
        write(report_unit, *) "Data :"
        write(report_unit, *) "    epsilon = ", this%epsilon
        write(report_unit, *) "    alpha = ", this%alpha
        write(report_unit, *) "    rCut = ", this%rCut
        write(report_unit, *) "    dr = ", this%dr
        
    end subroutine MixingPotential_report
    
    !> Accessor : rMin
    
    function MixingPotential_getRmin(this) result(getRmin)
    
        class(MixingPotential), intent(in) :: this
        
        real(DP) :: getRmin
        
        getRmin = this%rMin
    
    end function MixingPotential_getRmin
    
    !> Accessor : rCut
    
    function MixingPotential_getRcut(this) result(getRcut)
    
        class(MixingPotential), intent(in) :: this
        
        real(DP) :: getRcut
        
        getRcut = this%rCut
    
    end function MixingPotential_getRcut
    
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
                    write(error_unit, *) this%name, " :    Overlap !", iCol1, iCol2
                    write(error_unit, *) "    r_mix = ", r_mix
                    stop
                end if

            end do
        end do

        write(output_unit, *) this%name, " :    Overlap test : OK !"
    
    end subroutine MixingPotential_overlapTest

    !> MixingPotential energy
    !> Tabulation of Yukawa potential
    !> \f[ \epsilon \frac{e^{-\alpha (r-r_{min})}}{r} \f]
    
    subroutine MixingPotential_Epot_init(this)
    
        class(MixingPotential), intent(inout) :: this

        integer :: i
        real(DP) :: r_i
       
        ! cut
        do i = this%iMin, this%iCut       
            r_i = real(i, DP)*this%dr
            this%Epot_tab(i) = this%epsilon * exp(-this%alpha*(r_i-this%rMin)) / r_i
        end do
        
        write(*, *) "mix"
        write(*, *) "cutSimple", this%Epot_tab(this%iCut)
        write(*, *) "cutComplx", this%epsilon * exp(-this%alpha*(this%rCut-this%rMin)) / this%rCut
        
        ! shift        
        this%Epot_tab(:) = this%Epot_tab(:) - this%Epot_tab(this%iCut)

    end subroutine MixingPotential_Epot_init
    
    !> Print the tabulated potential
    
    subroutine MixingPotential_Epot_print(this, Epot_unit)

	    class(MixingPotential), intent(in) :: this
	    integer, intent(in) :: Epot_unit
	
        integer :: i
        real(DP) :: r_i
	
	    do i = this%iMin, this%iCut
		    r_i = real(i, DP)*this%dr
		    write(Epot_unit, *) r_i, this%Epot_tab(i)
	    end do

    end subroutine MixingPotential_Epot_print

    function MixingPotential_Epot_pair(this, r) result(Epot_pair)
        
        class(MixingPotential), intent(in) :: this
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
        
    end function MixingPotential_Epot_pair
    
    subroutine MixingPotential_Epot_neigh(this, xCol, iCell, neigh, other_X, overlap, energ)
        
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
                energ = energ + this%Epot_pair(r)
                
                if (.not. associated(next%next)) exit
                
                current => next
            
            end do            
            
        end do
    
    end subroutine MixingPotential_Epot_neigh
    
    !> Total potential energy
    
    function MixingPotential_Epot_conf(this, type1_X, type2_X) result(Epot_conf)
    
        class(MixingPotential), intent(in) :: this
        real(DP), dimension(:, :), intent(in) :: type1_X, type2_X
        
        integer :: Ncol1, Ncol2
        integer :: iCol1, iCol2
        real(DP) :: r_mix
        real(DP) :: Epot_conf
        
        Ncol1 = size(type1_X, 2)
        Ncol2 = size(type2_X, 2)
        
        Epot_conf = 0._DP
        
        do iCol1 = 1, Ncol1
            do iCol2 = 1, Ncol2
                
                r_mix = dist(type1_X(:, iCol1), type2_X(:, iCol2))
                Epot_conf = Epot_conf + this%Epot_pair(r_mix)

            end do
        end do
    
    end function MixingPotential_Epot_conf

end module class_mixingPotential
