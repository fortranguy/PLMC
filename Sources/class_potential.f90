!> \brief Description of the Potential class

module class_potential

use data_constants
use data_particles
use data_potentiel

implicit none

private

	type, public :: Potential
	
		private
		
		real(DP) :: rMin !< minimum distance between two particles
		real(DP) :: rCut !< short-range cut
		real(DP) :: dr !< discretisation step
        integer :: iMin !< minimum index of tabulation : minimum distance
        integer :: iCut !< maximum index of tabulation : until potential cut
        real(DP) :: epsilon !< factor in Yukawa
        real(DP) :: alpha !< coefficient in Yukawa
        real(DP), dimension(:), allocatable :: ePot_tab !< tabulation
	
	contains
	
		procedure :: construct => Potential_construct
		procedure :: destroy => Potential_destroy
	
		procedure :: ePot_init => Potential_ePot_init
        procedure :: ePot => Potential_ePot
        procedure :: ePot_neigh => Potential_ePot_neigh
        !procedure :: ePot_total => Potential_ePot_total
	
	end type
	
contains

	subroutine Potential_construct(this)
	
		class(Potential), intent(out) :: this
		
		! Particles
        this%rMin = mix_rMin
		
		! Potential
        this%rCut = mix_rCut
        this%dr = mix_dr
        this%iMin = int(this%rMin/this%dr)
        this%iCut = int(this%rCut/this%dr)
        this%epsilon = mix_epsilon
        this%alpha = mix_alpha
        allocate(this%ePot_tab(this%iMin:this%iCut))
        call this%ePot_init()
	
	end subroutine Potential_construct
	
	subroutine Potential_destroy(this)
    
        class(interactingSpheres), intent(inout) :: this
        
        if (allocated(this%ePot_tab)) then
	        deallocate(this%ePot_tab)
        end if
    
    end subroutine Potential_destroy

	!> Potential energy
    !> Tabulation of Yukawa potential
    !> \f[ \epsilon \frac{e^{-\alpha (r-r_{min})}}{r} \f]
    
    subroutine Potential_ePot_init(this)
    
        class(Potential), intent(inout) :: this

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

    end subroutine Potential_ePot_init

    function Potential_ePot(this, r) result(ePot)
        
        class(Potential), intent(in) :: this
        real(DP), intent(in) :: r
        
        integer :: i
        real(DP) :: r_i, ePot
       
        if (r < this%rCut) then
       
            i = int(r/this%dr)
            r_i = real(i, DP)*this%dr
            ePot = this%ePot_tab(i) + (r-r_i)/this%dr * &
                (this%ePot_tab(i+1)-this%ePot_tab(i))
           
        else
       
            ePot = 0._DP
           
        end if
        
    end function Potential_ePot
    
    subroutine Potential_ePot_neigh(this, xCol, iCell, other_neigh, other_X, 
    	overlap, energ)
        
        class(Potential), intent(in) :: this
        real(DP), dimension(Dim), intent(in) :: xCol !< type A
        integer, intent(in) :: iCell < ! type A in mix grid
        type(Neighbours), intent(in) :: other_neigh
        real(DP), dimension(Dim, :) :: other_X
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energ
    
        integer :: iNeigh,  iCell_neigh
        real(DP) :: r
    
        type(Link), pointer :: current => null(), next => null()
        
        overlap = .false.
        energ = 0._DP
    
        do iNeigh = 1, cell_neighs_nb
        
            iCell_neigh = other%cell_neighs(iNeigh, iCell)
            current => other%cellsBegin(iCell_neigh)%particle%next            
            if (.not. associated(current%next)) cycle
            
            do
            
                next => current%next
                
                r = dist(xCol(:), other_X(:, current%iCol))
                if (r < this%rMin) then
                    overlap = .true.
                    return
                end if
                energ = energ + this%ePot(r)
                
                if (.not. associated(next%next)) exit
                
                current => next
            
            end do            
            
        end do
    
    end subroutine Potential_ePot_neigh
    
    !> Total potential energy

end module class_potential
