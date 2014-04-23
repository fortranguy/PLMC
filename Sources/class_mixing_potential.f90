!> \brief Description of the  Mixing Potential class

module class_mixing_potential

use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
use data_precisions, only: DP, real_zero
use data_box, only: Ndim
use data_particles, only: mix_delta
use data_potential, only: mix_rMin_factor, mix_rCut, mix_dr, mix_epsilon, mix_alpha
use data_neighbour_cells, only: NnearCell
use module_types, only: Node, Particle_Index
use module_physics_micro, only: set_discrete_length, dist_PBC, Epot_yukawa
use class_neighbour_cells
use class_hard_spheres

implicit none

private

    type, public :: Mixing_Potential

        private
        
        character(len=5) :: name
        
        real(DP) :: delta !< demixing length
        real(DP) :: min_distance
        real(DP) :: rMin !< minimum distance between two particles
        real(DP) :: rCut !< short-range cut
        real(DP) :: dr !< discretisation step
        integer :: iMin !< minimum index of tabulation: minimum distance
        integer :: iCut !< maximum index of tabulation: until potential cut
        real(DP) :: epsilon !< factor in Yukawa
        real(DP) :: alpha !< coefficient in Yukawa
        real(DP), dimension(:), allocatable :: Epot_tab !< tabulation
        
        real(DP), dimension(Ndim) :: cell_size

    contains
    
        procedure :: construct => Mixing_Potential_construct
        procedure :: destroy => Mixing_Potential_destroy
        procedure :: write_report => Mixing_Potential_write_report

        procedure :: get_min_distance => Mixing_Potential_get_min_distance
        procedure :: get_rCut => Mixing_Potential_get_rCut
        procedure :: set_cell_size => Mixing_Potential_set_cell_size
        procedure :: get_cell_size => Mixing_Potential_get_cell_size
        
        procedure :: test_overlap => Mixing_Potential_test_overlap
        procedure, private :: set_Epot_tab => Mixing_Potential_set_Epot_tab
        procedure :: set_Epot => Mixing_Potential_set_Epot
        procedure :: write_Epot => Mixing_Potential_write_Epot
        procedure, private :: Epot_pair => Mixing_Potential_Epot_pair
        procedure :: Epot_neighCells => Mixing_Potential_Epot_neighCells
        procedure :: Epot_conf => Mixing_Potential_Epot_conf

    end type

contains

    subroutine Mixing_Potential_construct(this, type1_diameter, type2_diameter)

        class(Mixing_Potential), intent(out) :: this
        real(DP), intent(in) :: type1_diameter, type2_diameter
        
        this%name = "[mix]"
        write(output_unit, *) this%name, " class construction"
        
        this%delta = mix_delta
        this%min_distance = (type1_diameter + type2_diameter)/2._DP + this%delta

    end subroutine Mixing_Potential_construct

    subroutine Mixing_Potential_destroy(this)
    
        class(Mixing_Potential), intent(inout) :: this
        
        write(output_unit, *) this%name, " class destruction"
        
        if (allocated(this%Epot_tab)) deallocate(this%Epot_tab)
    
    end subroutine Mixing_Potential_destroy
    
    !> Report
    
    subroutine Mixing_Potential_write_report(this, report_unit)
    
        class(Mixing_Potential), intent(in) :: this
        integer, intent(in) :: report_unit
        
        write(report_unit, *) "Data: "
        
        write(report_unit, *) "    delta = ", this%delta
        
        write(report_unit, *) "    epsilon = ", this%epsilon
        write(report_unit, *) "    alpha = ", this%alpha
        write(report_unit, *) "    rCut = ", this%rCut
        write(report_unit, *) "    dr = ", this%dr
        
    end subroutine Mixing_Potential_write_report

    !> Accessors & Mutators

    pure function Mixing_Potential_get_min_distance(this) result(get_min_distance)
        class(Mixing_Potential), intent(in) :: this
        real(DP) :: get_min_distance
        get_min_distance = this%min_distance
    end function Mixing_Potential_get_min_distance
    
    pure function Mixing_Potential_get_rCut(this) result(get_rCut)
        class(Mixing_Potential), intent(in) :: this
        real(DP) :: get_rCut
        get_rCut = this%rCut
    end function Mixing_Potential_get_rCut
    
    pure subroutine Mixing_Potential_set_cell_size(this)
        class(Mixing_Potential), intent(inout) :: this
        this%cell_size(:) = this%rCut
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
                    
                type1_xCol(:) = type1%all_positions(:, type1_i_particle)
                type2_xCol(:) = type2%all_positions(:, type2_i_particle)
                r_mix = dist_PBC(Box_size, type1_xCol, type2_xCol)
                if (r_mix < this%rMin) then
                    write(error_unit, *) this%name, ":    Overlap !", type1_i_particle, type2_i_particle
                    write(error_unit, *) "    r_mix = ", r_mix
                    error stop
                end if

            end do
        end do

        write(output_unit, *) this%name, ":    Overlap test: OK !"
    
    end subroutine Mixing_Potential_test_overlap

    !> Mixing_Potential energy
    !> Tabulation of Yukawa potential
    !> \f[ \epsilon \frac{e^{-\alpha (r-r_{min})}}{r} \f]
    
    pure subroutine Mixing_Potential_set_Epot_tab(this)
    
        class(Mixing_Potential), intent(inout) :: this

        integer :: i
        real(DP) :: r_i
       
        ! cut
        do i = this%iMin, this%iCut
            r_i = real(i, DP)*this%dr
            this%Epot_tab(i) = Epot_yukawa(this%epsilon, this%alpha, this%rMin, r_i)
        end do
        
        ! shift
        this%Epot_tab(:) = this%Epot_tab(:) - this%Epot_tab(this%iCut)

    end subroutine Mixing_Potential_set_Epot_tab
    
    subroutine Mixing_Potential_set_Epot(this)
    
        class(Mixing_Potential), intent(inout) :: this

        this%rMin = mix_rMin_factor * this%min_distance
        this%rCut = mix_rCut
        
        if (this%rCut < this%rMin) then
            this%rCut = this%rMin
        end if
        
        this%dr = mix_dr
        call set_discrete_length(this%rMin, this%dr)
        this%iMin = int(this%rMin/this%dr)
        this%iCut = int(this%rCut/this%dr) + 1
        this%epsilon = mix_epsilon
        this%alpha = mix_alpha
        
        if (allocated(this%Epot_tab)) deallocate(this%Epot_tab)
        allocate(this%Epot_tab(this%iMin:this%iCut))
        call this%set_Epot_tab()
        
    end subroutine Mixing_Potential_set_Epot
    
    !> Write the tabulated potential
    
    subroutine Mixing_Potential_write_Epot(this, Epot_unit)

        class(Mixing_Potential), intent(in) :: this
        integer, intent(in) :: Epot_unit

        integer :: i
        real(DP) :: r_i

        do i = this%iMin, this%iCut
            r_i = real(i, DP)*this%dr
            write(Epot_unit, *) r_i, this%Epot_tab(i)
        end do

    end subroutine Mixing_Potential_write_Epot

    pure function Mixing_Potential_Epot_pair(this, r) result(Epot_pair)
        
        class(Mixing_Potential), intent(in) :: this
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
        
    end function Mixing_Potential_Epot_pair
    
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
            call this_cells%point_to_begin(current, nearCell_index)
            if (.not. associated(current%next)) cycle
            
            do
            
                next => current%next

                if (current%number /= particle%other_number) then
                    r = dist_PBC(Box_size, particle%position, other%get_position(current%number))
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
    
    end subroutine Mixing_Potential_Epot_neighCells
    
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
                
                type1_xCol(:) = type1%all_positions(:, type1_i_particle)
                type2_xCol(:) = type2%all_positions(:, type2_i_particle)
                r_mix = dist_PBC(Box_size, type1_xCol, type2_xCol)
                Epot_conf = Epot_conf + this%Epot_pair(r_mix)

            end do
        end do
    
    end function Mixing_Potential_Epot_conf

end module class_mixing_potential
