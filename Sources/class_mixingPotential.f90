!> \brief Description of the  Mixing Potential class

module class_mixingPotential

use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
use data_precisions, only : DP
use data_cell, only : Ndim
use data_particles, only : mix_rMin
use data_potential, only : mix_rCut, mix_dr, mix_epsilon, mix_alpha
use data_neighbourCells, only : NnearCell, mix_cell_size
use mod_physics, only : dist
use class_neighbourCells

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
        
        real(DP), dimension(Ndim) :: cell_size

    contains

        procedure :: construct => MixingPotential_construct
        procedure :: destroy => MixingPotential_destroy

        procedure :: PrintReport => MixingPotential_printReport

        procedure :: getRmin => MixingPotential_getRmin
        procedure :: getRcut => MixingPotential_getRcut
        procedure :: getCell_size => MixingPotential_getCell_size
        
        procedure :: overlapTest => MixingPotential_overlapTest

        procedure, private :: Epot_init => MixingPotential_Epot_init
        procedure :: Epot_print => MixingPotential_Epot_print
        procedure :: Epot_pair => MixingPotential_Epot_pair
        procedure :: Epot_neighCells => MixingPotential_Epot_neighCells
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
        
        ! Neighbours        
        this%cell_size(:) = mix_cell_size(:)


    end subroutine MixingPotential_construct

    subroutine MixingPotential_destroy(this)
    
        class(MixingPotential), intent(inout) :: this
        
        if (allocated(this%Epot_tab)) then
            deallocate(this%Epot_tab)
        end if
    
    end subroutine MixingPotential_destroy
    
    !> Report
    
    subroutine MixingPotential_printReport(this, report_unit)
    
        class(MixingPotential), intent(in) :: this
        integer, intent(in) :: report_unit    
        
        write(report_unit, *) "Data :"
        write(report_unit, *) "    epsilon = ", this%epsilon
        write(report_unit, *) "    alpha = ", this%alpha
        write(report_unit, *) "    rCut = ", this%rCut
        write(report_unit, *) "    dr = ", this%dr
        
    end subroutine MixingPotential_printReport
    
    !> Accessor : rMin
    
    pure function MixingPotential_getRmin(this) result(getRmin)
    
        class(MixingPotential), intent(in) :: this        
        real(DP) :: getRmin
        
        getRmin = this%rMin
    
    end function MixingPotential_getRmin
    
    !> Accessor : rCut
    
    pure function MixingPotential_getRcut(this) result(getRcut)
    
        class(MixingPotential), intent(in) :: this        
        real(DP) :: getRcut
        
        getRcut = this%rCut
    
    end function MixingPotential_getRcut
    
    !> Accessor : cell_size
    
    pure function MixingPotential_getCell_size(this) result(getCell_size)
    
        class(MixingPotential), intent(in) :: this        
        real(DP), dimension(Ndim) :: getCell_size
        
        getCell_size(:) = this%cell_size(:)
    
   end function MixingPotential_getCell_size
    
    !> Overlapt test
    
    subroutine MixingPotential_overlapTest(this, type1_positions, type2_positions)
    
        class(MixingPotential), intent(in) :: this
        real(DP), dimension(:, :), intent(in) :: type1_positions, type2_positions
        
        integer :: Ncol1, Ncol2
        integer :: iCol1, iCol2
        real(DP) :: r_mix
        
        Ncol1 = size(type1_positions, 2)
        Ncol2 = size(type2_positions, 2)
        
        do iCol1 = 1, Ncol1
            do iCol2 = 1, Ncol2
                    
                r_mix = dist(type1_positions(:, iCol1), type2_positions(:, iCol2))
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

    pure function MixingPotential_Epot_pair(this, r) result(Epot_pair)
        
        class(MixingPotential), intent(in) :: this
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
        
    end function MixingPotential_Epot_pair
    
    subroutine MixingPotential_Epot_neighCells(this, xCol, iTotalCell, neigh, other_positions, &
                                               overlap, energ)
        
        class(MixingPotential), intent(in) :: this
        real(DP), dimension(:), intent(in) :: xCol !< type A
        integer, intent(in) :: iTotalCell !< type A in mix grid
        type(NeighbourCells), intent(in) :: neigh
        real(DP), dimension(:, :), contiguous, intent(in) :: other_positions
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energ
    
        integer :: iNearCell,  nearCell_index
        real(DP) :: r
    
        type(Link), pointer :: current => null(), next => null()
        
        overlap = .false.
        energ = 0._DP
    
        do iNearCell = 1, NnearCell
        
            nearCell_index = neigh%nearCells_among_totalCells(iNearCell, iTotalCell)
            current => neigh%beginCells(nearCell_index)%particle%next            
            if (.not. associated(current%next)) cycle
            
            do
            
                next => current%next
                
                r = dist(xCol(:), other_positions(:, current%iCol))
                if (r < this%rMin) then
                    overlap = .true.
                    return
                end if
                energ = energ + this%Epot_pair(r)
                
                if (.not. associated(next%next)) exit
                
                current => next
            
            end do            
            
        end do
    
    end subroutine MixingPotential_Epot_neighCells
    
    !> Total potential energy
    
    pure function MixingPotential_Epot_conf(this, type1_positions, type2_positions) result(Epot_conf)
    
        class(MixingPotential), intent(in) :: this
        real(DP), dimension(:, :), intent(in) :: type1_positions, type2_positions
        real(DP) :: Epot_conf
        
        integer :: Ncol1, Ncol2
        integer :: iCol1, iCol2
        real(DP) :: r_mix
        
        Ncol1 = size(type1_positions, 2)
        Ncol2 = size(type2_positions, 2)
        
        Epot_conf = 0._DP
        
        do iCol1 = 1, Ncol1
            do iCol2 = 1, Ncol2
                
                r_mix = dist(type1_positions(:, iCol1), type2_positions(:, iCol2))
                Epot_conf = Epot_conf + this%Epot_pair(r_mix)

            end do
        end do
    
    end function MixingPotential_Epot_conf

end module class_mixingPotential
