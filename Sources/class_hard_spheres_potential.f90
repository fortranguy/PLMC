module class_hard_spheres_potential

use data_precisions, only: DP
use data_potential, only: hard_rMin_factor
use data_neighbour_cells, only: NnearCell
use module_types_micro, only: Box_Dimensions, Node, Particle_Index
use module_physics_micro, only: dist_PBC
use class_neighbour_cells
use class_hard_spheres

implicit none

private

    type, public :: Hard_Spheres_Potential
    
        private
        real(DP) :: min_distance
        real(DP) :: range_cut
        
    contains
    
        procedure :: construct => Hard_Spheres_Potential_construct
        procedure :: write => Hard_Spheres_Potential_write
        procedure :: get_range_cut => Hard_Spheres_Potential_get_range_cut
        procedure :: neighCells => Hard_Spheres_Potential_neighCells
        procedure :: conf => Hard_Spheres_Potential_conf
        
    end type Hard_Spheres_Potential
    
contains

    subroutine Hard_Spheres_Potential_construct(this, diameter)    
        class(Hard_Spheres_Potential), intent(inout) :: this
        real(DP), intent(in) :: diameter
        
        this%min_distance = hard_rMin_factor * diameter ! careful for Dipoles !
        this%range_cut = this%min_distance
        
    end subroutine Hard_Spheres_Potential_construct
    
    subroutine Hard_Spheres_Potential_write(this, unit)    
        class(Hard_Spheres_Potential), intent(in) :: this
        integer, intent(in) :: unit

        write(unit, *) this%range_cut, 0._DP
    
    end subroutine Hard_Spheres_Potential_write
    
    pure function Hard_Spheres_Potential_get_range_cut(this) result(get_range_cut)
        class(Hard_Spheres_Potential), intent(in) :: this
        real(DP) :: get_range_cut
        
        get_range_cut = this%range_cut
    end function Hard_Spheres_Potential_get_range_cut
    
    
    subroutine Hard_Spheres_write_report(this, unit)
        class(Hard_Spheres_Potential), intent(in) :: this
        integer, intent(in) :: unit
    
        write(unit, *) "    range_cut = ", this%range_cut
    end subroutine Hard_Spheres_write_report
    
    subroutine Hard_Spheres_Potential_neighCells(this, Box_size, this_spheres, this_cells, particle, &
                                                 overlap, energ)
        
        class(Hard_Spheres_Potential), intent(in) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        class(Hard_Spheres), intent(in) :: this_spheres
        class(Neighbour_Cells), intent(in) :: this_cells
        type(Particle_Index), intent(in) :: particle
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energ
    
        integer :: iNearCell,  nearCell_index
        real(DP) :: r_ij
    
        type(Node), pointer :: current => null(), next => null()
        
        overlap = .false.
        energ = 0._DP
    
        do iNearCell = 1, NnearCell
        
            nearCell_index = this_cells%near_among_total(iNearCell, particle%same_iCell)
            call this_cells%point_to_begin(current, nearCell_index)
            if (.not. associated(current%next)) cycle
            
            do
            
                next => current%next
            
                if (current%number /= particle%number) then
                    r_ij = dist_PBC(Box_size, particle%position, &
                                    this_spheres%get_position(current%number))
                    if (r_ij < this%min_distance) then
                        overlap = .true.
                        return
                    end if
                end if
                
                if (.not. associated(next%next)) exit
                
                current => next
            
            end do
            
        end do
    
    end subroutine Hard_Spheres_Potential_neighCells
    
    !> Total potential energy: dummy
    
    pure function Hard_Spheres_Potential_conf(this) result(conf)
    
        class(Hard_Spheres_Potential), intent(in) :: this
        real(DP) :: conf
    
        conf = 0._DP
        
    end function Hard_Spheres_Potential_conf

end module class_hard_spheres_potential
