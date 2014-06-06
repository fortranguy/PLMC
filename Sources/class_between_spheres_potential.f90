!> \brief Description of the Between Spheres Potential_Energy class

module class_between_spheres_potential

use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
use data_precisions, only: DP, real_zero
use data_box, only: Ndim
use data_neighbour_cells, only: NnearCell
use json_module, only: json_file
use module_types_micro, only: Node, Particle_Index
use module_physics_micro, only: set_discrete_length, PBC_distance, potential_energy_yukawa
use module_data, only: test_data_found
use class_hard_spheres, only: Hard_Spheres
use class_neighbour_cells, only: Neighbour_Cells


implicit none

private

    type, public :: Between_Hard_Spheres_Potential_Energy

        private
        
        character(len=5) :: name
        
        real(DP) :: non_additivity
        real(DP) :: diameter
        real(DP) :: min_distance
        real(DP) :: range_cut
        real(DP) :: delta
        integer :: i_min_distance
        integer :: i_range_cut
        real(DP) :: epsilon !< factor in Yukawa
        real(DP) :: alpha !< coefficient in Yukawa
        real(DP), dimension(:), allocatable :: tabulation
        
        real(DP), dimension(Ndim) :: cell_size

    contains
    
        procedure :: construct => Between_Hard_Spheres_Potential_Energy_construct
        procedure :: destroy => Between_Hard_Spheres_Potential_Energy_destroy

        procedure :: get_diameter => Between_Hard_Spheres_Potential_Energy_get_diameter
        procedure :: get_range_cut => Between_Hard_Spheres_Potential_Energy_get_range_cut
        procedure :: set_cell_size => Between_Hard_Spheres_Potential_Energy_set_cell_size
        procedure :: get_cell_size => Between_Hard_Spheres_Potential_Energy_get_cell_size
        
        procedure :: test_overlap => Between_Hard_Spheres_Potential_Energy_test_overlap
        procedure, private :: set_tabulation => Between_Hard_Spheres_Potential_Energy_set_tabulation
        procedure :: write => Between_Hard_Spheres_Potential_Energy_write
        procedure, private :: pair => Between_Hard_Spheres_Potential_Energy_pair
        procedure :: neighCells => Between_Hard_Spheres_Potential_Energy_neighCells
        procedure :: conf => Between_Hard_Spheres_Potential_Energy_conf

    end type

contains

    subroutine Between_Hard_Spheres_Potential_Energy_construct(this, json, type1_diameter, type2_diameter)

        class(Between_Hard_Spheres_Potential_Energy), intent(out) :: this
        type(json_file), intent(inout) :: json
        real(DP), intent(in) :: type1_diameter, type2_diameter
        
        character(len=4096) :: data_name
        logical :: found
        
        real(DP) :: min_distance_factor
        
        this%name = "[mix]"
        write(output_unit, *) this%name, " class construction"
        
        data_name = "Particles.Between Spheres.non addivity"
        call json%get(data_name, this%non_additivity, found)
        call test_data_found(data_name, found)        
        this%diameter = (type1_diameter + type2_diameter)/2._DP + this%non_additivity
        
        data_name = "Potential Energy.Between Spheres.minimum distance factor"
        call json%get(data_name, min_distance_factor, found)
        call test_data_found(data_name, found)        
        this%min_distance = min_distance_factor * this%diameter
        
        data_name = "Potential Energy.Between Spheres.range cut"
        call json%get(data_name, this%range_cut, found)
        call test_data_found(data_name, found)        
        if (this%range_cut < this%min_distance) then
            this%range_cut = this%min_distance
        end if
        
        data_name = "Potential Energy.Between Spheres.delta"
        call json%get(data_name, this%delta, found)
        call test_data_found(data_name, found)        
        call set_discrete_length(this%min_distance, this%delta)
        this%i_min_distance = int(this%min_distance/this%delta)
        this%i_range_cut = int(this%range_cut/this%delta) + 1
        
        data_name = "Potential Energy.Between Spheres.Yukawa.epsilon"
        call json%get(data_name, this%epsilon, found)
        call test_data_found(data_name, found)
        
        data_name = "Potential Energy.Between Spheres.Yukawa.alpha"
        call json%get(data_name, this%alpha, found)
        call test_data_found(data_name, found)
        
        if (allocated(this%tabulation)) deallocate(this%tabulation)
        allocate(this%tabulation(this%i_min_distance:this%i_range_cut))
        call this%set_tabulation()

    end subroutine Between_Hard_Spheres_Potential_Energy_construct

    subroutine Between_Hard_Spheres_Potential_Energy_destroy(this)
    
        class(Between_Hard_Spheres_Potential_Energy), intent(inout) :: this
        
        write(output_unit, *) this%name, " class destruction"
        
        if (allocated(this%tabulation)) deallocate(this%tabulation)
    
    end subroutine Between_Hard_Spheres_Potential_Energy_destroy

    !> Accessors & Mutators

    pure function Between_Hard_Spheres_Potential_Energy_get_diameter(this) result(get_diameter)
        class(Between_Hard_Spheres_Potential_Energy), intent(in) :: this
        real(DP) :: get_diameter
        get_diameter = this%diameter
    end function Between_Hard_Spheres_Potential_Energy_get_diameter
    
    pure function Between_Hard_Spheres_Potential_Energy_get_range_cut(this) result(get_range_cut)
        class(Between_Hard_Spheres_Potential_Energy), intent(in) :: this
        real(DP) :: get_range_cut
        get_range_cut = this%range_cut
    end function Between_Hard_Spheres_Potential_Energy_get_range_cut
    
    pure subroutine Between_Hard_Spheres_Potential_Energy_set_cell_size(this)
        class(Between_Hard_Spheres_Potential_Energy), intent(inout) :: this
        this%cell_size(:) = this%range_cut
    end subroutine Between_Hard_Spheres_Potential_Energy_set_cell_size
    
    pure function Between_Hard_Spheres_Potential_Energy_get_cell_size(this) result(get_cell_size)
        class(Between_Hard_Spheres_Potential_Energy), intent(in) :: this
        real(DP), dimension(Ndim) :: get_cell_size
        get_cell_size(:) = this%cell_size(:)
   end function Between_Hard_Spheres_Potential_Energy_get_cell_size
    
    !> Overlapt test
    
    subroutine Between_Hard_Spheres_Potential_Energy_test_overlap(this, Box_size, type1, type2)
    
        class(Between_Hard_Spheres_Potential_Energy), intent(in) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        class(Hard_Spheres), intent(in) :: type1, type2
        
        integer :: type1_i_particle, type2_i_particle
        real(DP) :: r_mix
        real(DP), dimension(Ndim) :: type1_xCol, type2_xCol
        
        do type1_i_particle = 1, type1%get_num_particles()
            do type2_i_particle = 1, type2%get_num_particles()
                    
                type1_xCol(:) = type1%get_position(type1_i_particle)
                type2_xCol(:) = type2%get_position(type2_i_particle)
                r_mix = PBC_distance(Box_size, type1_xCol, type2_xCol)
                if (r_mix < this%min_distance) then
                    write(error_unit, *) this%name, ":    Overlap !", type1_i_particle, type2_i_particle
                    write(error_unit, *) "    r_mix = ", r_mix
                    error stop
                end if

            end do
        end do

        write(output_unit, *) this%name, ":    Overlap test: OK !"
    
    end subroutine Between_Hard_Spheres_Potential_Energy_test_overlap
    
    !> Tabulation of Yukawa potential_energy
    !> \f[ \epsilon \frac{e^{-\alpha (r-r_{min})}}{r} \f]
    
    pure subroutine Between_Hard_Spheres_Potential_Energy_set_tabulation(this)
    
        class(Between_Hard_Spheres_Potential_Energy), intent(inout) :: this

        integer :: i
        real(DP) :: r_i
       
        ! cut
        do i = this%i_min_distance, this%i_range_cut
            r_i = real(i, DP)*this%delta
            this%tabulation(i) = potential_energy_yukawa(this%epsilon, this%alpha, this%min_distance, r_i)
        end do
        
        ! shift
        this%tabulation(:) = this%tabulation(:) - this%tabulation(this%i_range_cut)

    end subroutine Between_Hard_Spheres_Potential_Energy_set_tabulation
    
    !> Write the tabulated potential_energy
    
    subroutine Between_Hard_Spheres_Potential_Energy_write(this, potential_energy_unit)

        class(Between_Hard_Spheres_Potential_Energy), intent(in) :: this
        integer, intent(in) :: potential_energy_unit

        integer :: i
        real(DP) :: r_i

        do i = this%i_min_distance, this%i_range_cut
            r_i = real(i, DP)*this%delta
            write(potential_energy_unit, *) r_i, this%tabulation(i)
        end do

    end subroutine Between_Hard_Spheres_Potential_Energy_write
    
    !> Total potential_energy energy
    
    pure function Between_Hard_Spheres_Potential_Energy_conf(this, Box_size, type1, type2) result(conf)
    
        class(Between_Hard_Spheres_Potential_Energy), intent(in) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        class(Hard_Spheres), intent(in) :: type1, type2
        real(DP) :: conf

        integer :: type1_i_particle, type2_i_particle
        real(DP) :: r_mix
        real(DP), dimension(Ndim) :: type1_xCol, type2_xCol

        conf = 0._DP
        
        if (this%epsilon < real_zero) then
            return
        end if
        
        do type1_i_particle = 1, type1%get_num_particles()
            do type2_i_particle = 1, type2%get_num_particles()
                
                type1_xCol(:) = type1%get_position(type1_i_particle)
                type2_xCol(:) = type2%get_position(type2_i_particle)
                r_mix = PBC_distance(Box_size, type1_xCol, type2_xCol)
                conf = conf + this%pair(r_mix)

            end do
        end do
    
    end function Between_Hard_Spheres_Potential_Energy_conf
    
    subroutine Between_Hard_Spheres_Potential_Energy_neighCells(this, Box_size, particle, this_cells, other, overlap, &
                                                energ)
        
        class(Between_Hard_Spheres_Potential_Energy), intent(in) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        type(Particle_Index), intent(in) :: particle
        type(Neighbour_Cells), intent(in) :: this_cells
        type(Hard_Spheres), intent(in) :: other
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
                    r = PBC_distance(Box_size, particle%position, other%get_position(current%number))
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
    
    pure function Between_Hard_Spheres_Potential_Energy_pair(this, r) result(pair)
        
        class(Between_Hard_Spheres_Potential_Energy), intent(in) :: this
        real(DP), intent(in) :: r
        real(DP) :: pair
        
        integer :: i
        real(DP) :: r_i
       
        if (r < this%range_cut) then
            i = int(r/this%delta)
            r_i = real(i, DP)*this%delta
            pair = this%tabulation(i) + (r-r_i)/this%delta * (this%tabulation(i+1)-this%tabulation(i))
        else
            pair = 0._DP
        end if
        
    end function Between_Hard_Spheres_Potential_Energy_pair

end module class_between_spheres_potential
