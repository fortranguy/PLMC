module class_hard_spheres_potential

use data_precisions, only: DP
use data_box, only: Ndim
use data_neighbour_cells, only: NnearCell
use json_module, only: json_file
use module_types_micro, only: Box_Dimensions, Node, Particle_Index
use module_physics_micro, only: PBC_distance
use class_neighbour_cells, only: Neighbour_Cells
use class_hard_spheres, only: Hard_Spheres

implicit none

private

    type, public :: Hard_Spheres_Potential_Energy
    
        private
        real(DP) :: min_distance
        real(DP) :: range_cut
        
    contains
    
        procedure :: construct => Hard_Spheres_Potential_Energy_construct
        procedure :: write => Hard_Spheres_Potential_Energy_write
        procedure :: get_min_distance => Hard_Spheres_Potential_Energy_get_min_distance
        procedure :: get_range_cut => Hard_Spheres_Potential_Energy_get_range_cut
        procedure :: neighCells => Hard_Spheres_Potential_Energy_neighCells
        procedure :: total => Hard_Spheres_Potential_Energy_total
        procedure, private :: solo => Hard_Spheres_Potential_Energy_solo
        procedure, private :: pair => Hard_Spheres_Potential_Energy_pair
        
    end type Hard_Spheres_Potential_Energy
    
    type, extends(Hard_Spheres_Potential_Energy), public :: Between_Hard_Spheres_Potential_Energy
    
    contains
    
        procedure :: neighCells => Between_Hard_Spheres_Potential_Energy_neighCells
        procedure :: total => Between_Hard_Spheres_Potential_Energy_total
        procedure, private :: solo => Between_Hard_Spheres_Potential_Energy_solo
    
    end type Between_Hard_Spheres_Potential_Energy
    
contains

    subroutine Hard_Spheres_Potential_Energy_construct(this, json, type_name, diameter)  
    
        class(Hard_Spheres_Potential_Energy), intent(inout) :: this
        type(json_file), intent(inout) :: json
        character(len=*), intent(in) :: type_name
        real(DP), intent(in) ::  diameter
        
        character(len=4096) :: data_name
        logical :: found
        real(DP) :: min_distance_factor
        
        data_name = "Potential Energy."//type_name//".minimum distance factor"
        call json%get(data_name, min_distance_factor, found)
        call test_data_found(data_name, found)
        
        this%min_distance = min_distance_factor * diameter
        this%range_cut = this%min_distance
        
    end subroutine Hard_Spheres_Potential_Energy_construct
    
    subroutine Hard_Spheres_Potential_Energy_write(this, unit)    
    
        class(Hard_Spheres_Potential_Energy), intent(in) :: this
        integer, intent(in) :: unit

        write(unit, *) this%range_cut, this%pair(this%range_cut)
    
    end subroutine Hard_Spheres_Potential_Energy_write
    
    pure function Hard_Spheres_Potential_Energy_get_min_distance(this) result(get_min_distance)
    
        class(Hard_Spheres_Potential_Energy), intent(in) :: this
        real(DP) :: get_min_distance
        
        get_min_distance = this%min_distance
        
    end function Hard_Spheres_Potential_Energy_get_min_distance
    
    pure function Hard_Spheres_Potential_Energy_get_range_cut(this) result(get_range_cut)
    
        class(Hard_Spheres_Potential_Energy), intent(in) :: this
        real(DP) :: get_range_cut
        
        get_range_cut = this%range_cut
        
    end function Hard_Spheres_Potential_Energy_get_range_cut    
    
    subroutine Hard_Spheres_write_report(this, unit)
    
        class(Hard_Spheres_Potential_Energy), intent(in) :: this
        integer, intent(in) :: unit
    
        write(unit, *) "    range_cut = ", this%range_cut
        
    end subroutine Hard_Spheres_write_report
    
    subroutine Hard_Spheres_Potential_Energy_neighCells(this, Box_size, spheres, this_cells, &
                                                        particle, overlap, energ)
        
        class(Hard_Spheres_Potential_Energy), intent(in) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        class(Hard_Spheres), intent(in) :: spheres
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
            current => this_cells%beginCells(nearCell_index)%particle%next
            if (.not. associated(current%next)) cycle
            
            do
            
                next => current%next
            
                if (current%number /= particle%number) then
                    r_ij = PBC_distance(Box_size, particle%position, &
                                    spheres%get_position(current%number))
                    if (r_ij < this%min_distance) then
                        overlap = .true.
                        return
                    end if
                end if
                
                if (.not. associated(next%next)) exit
                
                current => next
            
            end do
            
        end do
    
    end subroutine Hard_Spheres_Potential_Energy_neighCells
    
    subroutine Between_Hard_Spheres_Potential_Energy_neighCells(this, Box_size, spheres, &
                                                                this_cells, particle, overlap, energ)
        
        class(Between_Hard_Spheres_Potential_Energy), intent(in) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        type(Hard_Spheres), intent(in) :: spheres
        type(Neighbour_Cells), intent(in) :: this_cells
        type(Particle_Index), intent(in) :: particle
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energ
    
        integer :: iNearCell,  nearCell_index
        real(DP) :: r
    
        type(Node), pointer :: current => null(), next => null()
        
        overlap = .false.
        energ = 0._DP
        
        do iNearCell = 1, NnearCell
        
            nearCell_index = this_cells%near_among_total(iNearCell, particle%mix_iCell)
            current => this_cells%beginCells(nearCell_index)%particle%next
            if (.not. associated(current%next)) cycle
            
            do
            
                next => current%next

                if (current%number /= particle%other_number) then
                    r = PBC_distance(Box_size, particle%position, &
                                               spheres%get_position(current%number))
                    if (r < this%min_distance) then
                        overlap = .true.
                        return
                    end if
                    energ = energ + this%pair(r)
                end if
                
                if (.not. associated(next%next)) exit
                
                current => next
            
            end do
            
        end do
    
    end subroutine Between_Hard_Spheres_Potential_Energy_neighCells
    
    !> Total potential_energy energy: dummy
    
    pure function Hard_Spheres_Potential_Energy_total(this, Box_size, spheres, other_spheres) &
        result(total)
    
        class(Hard_Spheres_Potential_Energy), intent(in) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        class(Hard_Spheres), intent(in) :: spheres
        class(Hard_Spheres), intent(in), optional :: other_spheres
        real(DP) :: total
        
        integer :: i_particle
        type(Particle_Index) :: particle
        
        total = 0._DP
        
        if (.not.present(other_spheres)) then            
            do i_particle = 1, spheres%get_num_particles()
                particle%number = i_particle
                particle%position(:) = spheres%get_position(particle%number)
                total = total + this%solo(Box_size, spheres, particle)
            end do            
        end if
        
    end function Hard_Spheres_Potential_Energy_total
    
    pure function Hard_Spheres_Potential_Energy_solo(this, Box_size, spheres, particle) result(solo)
    
        class(Hard_Spheres_Potential_Energy), intent(in) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        class(Hard_Spheres), intent(in) :: spheres
        type(Particle_Index), intent(in) :: particle
        real(DP) :: solo
        
        integer :: j_particle
        real(DP), dimension(Ndim) :: position_j
        real(DP) :: distance_ij
        
        solo = 0._DP
        
        do j_particle = 1, spheres%get_num_particles()
            if (j_particle /= particle%number) then
            
                position_j(:) = spheres%get_position(j_particle)
                distance_ij = PBC_distance(Box_size, particle%position, position_j)
                
                solo = solo + this%pair(distance_ij)
            
            end if            
        end do
    
    end function Hard_Spheres_Potential_Energy_solo
    
    pure function Between_Hard_Spheres_Potential_Energy_total(this, Box_size, spheres, other_spheres) &
        result(total)
    
        class(Between_Hard_Spheres_Potential_Energy), intent(in) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        class(Hard_Spheres), intent(in) :: spheres
        class(Hard_Spheres), intent(in), optional :: other_spheres
        real(DP) :: total
        
        integer :: i_particle
        type(Particle_Index) :: type1_particle
    
        total = 0._DP
        
        if (present(other_spheres)) then        
            do i_particle = 1, spheres%get_num_particles()            
                type1_particle%number = i_particle
                type1_particle%position(:) = spheres%get_position(type1_particle%number)
                total = total + this%solo(Box_size, other_spheres, type1_particle)                
            end do            
        end if
        
    end function Between_Hard_Spheres_Potential_Energy_total
    
    pure function Between_Hard_Spheres_Potential_Energy_solo(this, Box_size, spheres, particle) &
        result(solo)
    
        class(Between_Hard_Spheres_Potential_Energy), intent(in) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        class(Hard_Spheres), intent(in) :: spheres
        type(Particle_Index), intent(in) :: particle
        real(DP) :: solo
        
        integer :: j_particle
        real(DP), dimension(Ndim) :: position_j
        real(DP) :: distance_ij
        
        solo = 0._DP
        
        do j_particle = 1, spheres%get_num_particles()
            
            position_j(:) = spheres%get_position(j_particle)
            distance_ij = PBC_distance(Box_size, particle%position, position_j)
            
            solo = solo + this%pair(distance_ij)
                      
        end do
    
    end function Between_Hard_Spheres_Potential_Energy_solo
    
    !> Dummy
    
    pure function Hard_Spheres_Potential_Energy_pair(this, distance) result(pair)
    
        class(Hard_Spheres_Potential_Energy), intent(in) :: this
        real(DP), intent(in) :: distance
        real(DP) :: pair
        
        if (distance >= this%min_distance) then
            pair = 0._DP
        end if
    
    end function Hard_Spheres_Potential_Energy_pair    

end module class_hard_spheres_potential
