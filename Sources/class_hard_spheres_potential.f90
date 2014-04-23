module class_hard_spheres_potential

use data_precisions, only: DP
use module_types, only: Box_Dimensions, Node, Particle_Index
use data_neighbour_cells, only: NnearCell
use class_neighbour_cells
use class_hard_spheres

implicit none

private

    type, public :: Hard_Spheres_Potential
        private
        real(DP) :: min_distance
        real(DP) :: range_cut
    contains
        procedure :: set_domain => Hard_Spheres_Potential_set_domain
        procedure :: write => Hard_Spheres_Potential_write
        procedure :: neighCells => Hard_Spheres_Potential_neighCells
        procedure :: conf => Hard_Spheres_Potential_conf
    end type Hard_Spheres_Potential
    
contains

    subroutine Hard_Spheres_Potential_set_domain(this, diameter)
    
        class(Hard_Spheres_Potential), intent(inout) :: this
        real(DP), intent(in) :: diameter
        
        this%min_distance = hard_rMin_factor * diameter
        this%range_cut = this%min_distance
        
    end subroutine Hard_Spheres_Potential_set_domain
    
    subroutine Hard_Spheres_Potential_write(this, unit)
    
        class(Hard_Spheres_Potential), intent(in) :: this
        integer, intent(in) :: unit

        write(unit, *) this%range_cut, 0._DP
    
    end subroutine Hard_Spheres_Potential_write
    
    subroutine Hard_Spheres_write_report(this, unit)
        class(Hard_Spheres_Potential), intent(in) :: this
        integer, intent(in) :: unit
    
        write(report_unit, *) "    range_cut = ", this%range_cut
    subroutine Hard_Spheres_write_report
    
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
                    r_ij = dist_PBC(Box_size, particle%position, 
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
    
    pure function Hard_Spheres_Potential_conf(this, Box) result(conf)
    
        class(Hard_Spheres_Potential), intent(in) :: this
        type(Box_Dimensions), intent(in) :: Box
        real(DP) :: conf
        
        real(DP) :: volume_dummy
        volume_dummy = product(Box%size)
    
        conf = this%num_particles * 0._DP
        
    end function Hard_Spheres_Potential_conf

end module class_hard_spheres_potential
