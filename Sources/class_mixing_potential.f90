!> \brief Description of the  Mixing Potential class

module class_mixing_potential

use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
use data_precisions, only: DP, real_zero
use data_box, only: Ndim
use data_neighbour_cells, only: NnearCell
use json_module, only: json_file
use module_types_micro, only: Node, Particle_Index
use module_physics_micro, only: set_discrete_length, PBC_distance, Epot_yukawa
use module_data, only: test_data_found
use class_neighbour_cells
use class_hard_spheres

implicit none

private

    type, public :: Mixing_Potential

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
        real(DP), dimension(:), allocatable :: tabulation !< tabulation
        
        real(DP), dimension(Ndim) :: cell_size

    contains
    
        procedure :: construct => Mixing_Potential_construct
        procedure :: destroy => Mixing_Potential_destroy
        procedure :: write_report => Mixing_Potential_write_report

        procedure :: get_diameter => Mixing_Potential_get_diameter
        procedure :: get_range_cut => Mixing_Potential_get_range_cut
        procedure :: set_cell_size => Mixing_Potential_set_cell_size
        procedure :: get_cell_size => Mixing_Potential_get_cell_size
        
        procedure :: test_overlap => Mixing_Potential_test_overlap
        procedure, private :: set_tabulation => Mixing_Potential_set_tabulation
        procedure :: write => Mixing_Potential_write
        procedure, private :: Epot_pair => Mixing_Potential_Epot_pair
        procedure :: Epot_neighCells => Mixing_Potential_Epot_neighCells
        procedure :: Epot_conf => Mixing_Potential_Epot_conf

    end type

contains

    subroutine Mixing_Potential_construct(this, json, type1_diameter, type2_diameter)

        class(Mixing_Potential), intent(out) :: this
        type(json_file), intent(inout) :: json
        real(DP), intent(in) :: type1_diameter, type2_diameter
        
        character(len=4096) :: data_name
        logical :: found
        
        real(DP) :: min_distance_factor
        
        this%name = "[mix]"
        write(output_unit, *) this%name, " class construction"
        
        data_name = "Particles.Mixing.non addivity"
        call json%get(data_name, this%non_additivity, found)
        call test_data_found(data_name, found)        
        this%diameter = (type1_diameter + type2_diameter)/2._DP + this%non_additivity
        
        data_name = "Potential.Mixing.minimum distance factor"
        call json%get(data_name, min_distance_factor, found)
        call test_data_found(data_name, found)        
        this%min_distance = min_distance_factor * this%diameter
        
        data_name = "Potential.Mixing.range cut"
        call json%get(data_name, this%range_cut, found)
        call test_data_found(data_name, found)        
        if (this%range_cut < this%min_distance) then
            this%range_cut = this%min_distance
        end if
        
        data_name = "Potential.Mixing.delta"
        call json%get(data_name, this%delta, found)
        call test_data_found(data_name, found)        
        call set_discrete_length(this%min_distance, this%delta)
        this%i_min_distance = int(this%min_distance/this%delta)
        this%i_range_cut = int(this%range_cut/this%delta) + 1
        
        data_name = "Potential.Mixing.Yukawa.epsilon"
        call json%get(data_name, this%epsilon, found)
        call test_data_found(data_name, found)
        
        data_name = "Potential.Mixing.Yukawa.alpha"
        call json%get(data_name, this%alpha, found)
        call test_data_found(data_name, found)
        
        if (allocated(this%tabulation)) deallocate(this%tabulation)
        allocate(this%tabulation(this%i_min_distance:this%i_range_cut))
        call this%set_tabulation()

    end subroutine Mixing_Potential_construct

    subroutine Mixing_Potential_destroy(this)
    
        class(Mixing_Potential), intent(inout) :: this
        
        write(output_unit, *) this%name, " class destruction"
        
        if (allocated(this%tabulation)) deallocate(this%tabulation)
    
    end subroutine Mixing_Potential_destroy
    
    !> Report
    
    subroutine Mixing_Potential_write_report(this, report_unit)
    
        class(Mixing_Potential), intent(in) :: this
        integer, intent(in) :: report_unit
        
        write(report_unit, *) "Data: "
        
        write(report_unit, *) "    non_additivity = ", this%non_additivity
        
        write(report_unit, *) "    epsilon = ", this%epsilon
        write(report_unit, *) "    alpha = ", this%alpha
        write(report_unit, *) "    range_cut = ", this%range_cut
        write(report_unit, *) "    delta = ", this%delta
        
    end subroutine Mixing_Potential_write_report

    !> Accessors & Mutators

    pure function Mixing_Potential_get_diameter(this) result(get_diameter)
        class(Mixing_Potential), intent(in) :: this
        real(DP) :: get_diameter
        get_diameter = this%diameter
    end function Mixing_Potential_get_diameter
    
    pure function Mixing_Potential_get_range_cut(this) result(get_range_cut)
        class(Mixing_Potential), intent(in) :: this
        real(DP) :: get_range_cut
        get_range_cut = this%range_cut
    end function Mixing_Potential_get_range_cut
    
    pure subroutine Mixing_Potential_set_cell_size(this)
        class(Mixing_Potential), intent(inout) :: this
        this%cell_size(:) = this%range_cut
    end subroutine Mixing_Potential_set_cell_size
    
    pure function Mixing_Potential_get_cell_size(this) result(get_cell_size)
        class(Mixing_Potential), intent(in) :: this
        real(DP), dimension(Ndim) :: get_cell_size
        get_cell_size(:) = this%cell_size(:)
   end function Mixing_Potential_get_cell_size
    
    !> Overlapt test
    
    subroutine Mixing_Potential_test_overlap(this, Box_size, type1, type2)
    
        class(Mixing_Potential), intent(in) :: this
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
    
    end subroutine Mixing_Potential_test_overlap
    
    !> Tabulation of Yukawa potential
    !> \f[ \epsilon \frac{e^{-\alpha (r-r_{min})}}{r} \f]
    
    pure subroutine Mixing_Potential_set_tabulation(this)
    
        class(Mixing_Potential), intent(inout) :: this

        integer :: i
        real(DP) :: r_i
       
        ! cut
        do i = this%i_min_distance, this%i_range_cut
            r_i = real(i, DP)*this%delta
            this%tabulation(i) = Epot_yukawa(this%epsilon, this%alpha, this%min_distance, r_i)
        end do
        
        ! shift
        this%tabulation(:) = this%tabulation(:) - this%tabulation(this%i_range_cut)

    end subroutine Mixing_Potential_set_tabulation
    
    !> Write the tabulated potential
    
    subroutine Mixing_Potential_write(this, Epot_unit)

        class(Mixing_Potential), intent(in) :: this
        integer, intent(in) :: Epot_unit

        integer :: i
        real(DP) :: r_i

        do i = this%i_min_distance, this%i_range_cut
            r_i = real(i, DP)*this%delta
            write(Epot_unit, *) r_i, this%tabulation(i)
        end do

    end subroutine Mixing_Potential_write
    
    !> Total potential energy
    
    pure function Mixing_Potential_Epot_conf(this, Box_size, type1, type2) result(Epot_conf)
    
        class(Mixing_Potential), intent(in) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        class(Hard_Spheres), intent(in) :: type1, type2
        real(DP) :: Epot_conf

        integer :: type1_i_particle, type2_i_particle
        real(DP) :: r_mix
        real(DP), dimension(Ndim) :: type1_xCol, type2_xCol

        Epot_conf = 0._DP
        
        if (this%epsilon < real_zero) then
            return
        end if
        
        do type1_i_particle = 1, type1%get_num_particles()
            do type2_i_particle = 1, type2%get_num_particles()
                
                type1_xCol(:) = type1%get_position(type1_i_particle)
                type2_xCol(:) = type2%get_position(type2_i_particle)
                r_mix = PBC_distance(Box_size, type1_xCol, type2_xCol)
                Epot_conf = Epot_conf + this%Epot_pair(r_mix)

            end do
        end do
    
    end function Mixing_Potential_Epot_conf
    
    subroutine Mixing_Potential_Epot_neighCells(this, Box_size, particle, this_cells, other, overlap, &
                                                energ)
        
        class(Mixing_Potential), intent(in) :: this
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
                    energ = energ + this%Epot_pair(r)
                end if
                
                if (.not. associated(next%next)) exit
                
                current => next
            
            end do
            
        end do
    
    end subroutine Mixing_Potential_Epot_neighCells
    
    pure function Mixing_Potential_Epot_pair(this, r) result(Epot_pair)
        
        class(Mixing_Potential), intent(in) :: this
        real(DP), intent(in) :: r
        real(DP) :: Epot_pair
        
        integer :: i
        real(DP) :: r_i
       
        if (r < this%range_cut) then
            i = int(r/this%delta)
            r_i = real(i, DP)*this%delta
            Epot_pair = this%tabulation(i) + (r-r_i)/this%delta * (this%tabulation(i+1)-this%tabulation(i))
        else
            Epot_pair = 0._DP
        end if
        
    end function Mixing_Potential_Epot_pair

end module class_mixing_potential
