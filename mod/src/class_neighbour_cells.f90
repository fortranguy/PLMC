!> \brief Description of the neighbours class

module class_neighbour_cells

use, intrinsic :: iso_fortran_env, only: DP => REAL64, error_unit
use data_precisions, only: real_zero
use data_box, only: num_dimensions
use data_neighbour_cells, only: num_near_cells_dim, num_near_cells_layer, num_near_cells
use module_types_micro, only: Node, Linked_List
use module_geometry, only: geometry
use module_physics_micro, only: index_from_coord, coord_PBC
use class_hard_spheres, only: Hard_Spheres

implicit none
private

    type, public :: Neighbour_Cells
    
        private
        
        integer :: num_total_cell, num_total_cell_layer
        integer, dimension(num_dimensions) :: num_total_cell_dim
        real(DP), dimension(num_dimensions) :: cell_size
        integer, dimension(:, :), allocatable :: near_among_total
        type(Linked_List), dimension(:), allocatable :: begin_cells
        type(Linked_List), dimension(:), allocatable :: current_cells, next_cells
        
    contains
    
        procedure :: construct => Neighbour_Cells_construct
        procedure :: destroy => Neighbour_Cells_destroy
        
        procedure :: alloc_nodes => Neighbour_Cells_alloc_nodes
        procedure :: alloc_cells => Neighbour_Cells_alloc_cells
        procedure :: dealloc_nodes => Neighbour_Cells_dealloc_nodes
        procedure :: dealloc_cells => Neighbour_Cells_dealloc_cells        
        procedure :: check_cellsSize => Neighbour_Cells_check_cellsSize
        procedure :: index_from_position => Neighbour_Cells_index_from_position
        procedure :: all_cols_to_cells => Neighbour_Cells_all_cols_to_cells

        procedure :: near_cell_bounds => Neighbour_Cells_near_cell_bounds
        procedure :: point_to_begin => Neighbour_Cells_point_to_begin
        procedure :: remove_col_from_cell => Neighbour_Cells_remove_col_from_cell
        procedure :: add_col_to_cell => Neighbour_Cells_add_col_to_cell
        procedure, private :: init_near_among_total => &
                              Neighbour_Cells_init_near_among_total
        
    end type Neighbour_Cells
    
contains

    subroutine Neighbour_Cells_construct(this, Box_size, proposed_cell_size, cutoff)
    
        class(Neighbour_Cells), intent(out) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        real(DP), dimension(:), intent(in) :: proposed_cell_size
        real(DP), intent(in) :: cutoff
        
        this%num_total_cell_dim(:) = floor(Box_size(:)/proposed_cell_size(:))
        this%num_total_cell = product(this%num_total_cell_dim)
        this%num_total_cell_layer = product(this%num_total_cell_dim(1:2))
        this%cell_size(:) = Box_size(:)/real(this%num_total_cell_dim(:), DP)
        
        allocate(this%near_among_total(num_near_cells, this%num_total_cell))
            
        call this%check_CellsSize(Box_size, cutoff)
        call this%alloc_cells()
        call this%init_near_among_total()
    
    end subroutine Neighbour_Cells_construct
    
    subroutine Neighbour_Cells_destroy(this)
    
        class(Neighbour_Cells), intent(inout) :: this
        
        if (allocated(this%near_among_total)) deallocate(this%near_among_total)
        call this%dealloc_cells()
        
    end subroutine Neighbour_Cells_destroy
    
    !> Linked-list allocation
    
    subroutine Neighbour_Cells_alloc_nodes(this)
    
        class(Neighbour_Cells), intent(inout) :: this
    
        integer :: i_cell
    
        do i_cell = 1, this%num_total_cell

            allocate(this%begin_cells(i_cell)%particle)
            this%current_cells(i_cell)%particle => this%begin_cells(i_cell)%particle
            this%current_cells(i_cell)%particle%number = 0
            
            allocate(this%next_cells(i_cell)%particle)
            this%next_cells(i_cell)%particle%number = 0
            this%current_cells(i_cell)%particle%next => this%next_cells(i_cell)%particle
            this%current_cells(i_cell)%particle => this%next_cells(i_cell)%particle
    
        end do
    
    end subroutine Neighbour_Cells_alloc_nodes
    
    subroutine Neighbour_Cells_alloc_cells(this)
    
        class(Neighbour_Cells), intent(inout) :: this

        allocate(this%begin_cells(this%num_total_cell))
        allocate(this%current_cells(this%num_total_cell))
        allocate(this%next_cells(this%num_total_cell))
        
        call this%alloc_nodes()
        
    end subroutine Neighbour_Cells_alloc_cells
    
    ! Linked-list deallocation
    
    recursive pure subroutine free_link(current)

        type(Node), pointer :: current
        
        if (associated(current%next)) then
            call free_link(current%next)
        end if
        deallocate(current)
        
    end subroutine free_link
    
    pure subroutine Neighbour_Cells_dealloc_nodes(this)
    
        class(Neighbour_Cells), intent(inout) :: this
    
        integer :: i_cell

        do i_cell = 1, this%num_total_cell
            if (associated(this%begin_cells(i_cell)%particle)) then
                call free_link(this%begin_cells(i_cell)%particle)
            end if
        end do
    
    end subroutine Neighbour_Cells_dealloc_nodes
    
    pure subroutine Neighbour_Cells_dealloc_cells(this)
    
        class(Neighbour_Cells), intent(inout) :: this
        
        call this%dealloc_nodes()
        
        if (allocated(this%begin_cells)) deallocate(this%begin_cells)
        if (allocated(this%current_cells)) deallocate(this%current_cells)
        if (allocated(this%next_cells)) deallocate(this%next_cells)
    
    end subroutine Neighbour_Cells_dealloc_cells
    
    ! Neighbour_Cells cells size check
    
    subroutine Neighbour_Cells_check_cellsSize(this, Box_size, cutoff)
    
        class(Neighbour_Cells), intent(in) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        real(DP), intent(in) :: cutoff
        
        integer :: i_dim
        real(DP) :: Box_size_mod_cell_size
        
        do i_dim = 1, num_dimensions
        
            if (this%cell_size(i_dim) < cutoff) then
                write(error_unit, *) "    cutoff too big in the dimension", i_dim, ": "
                write(error_unit, *) "    ", this%cell_size(i_dim), "<", cutoff
                error stop
            end if
            
            if (this%num_total_cell_dim(i_dim) < num_near_cells_dim(i_dim)) then
                write(error_unit, *) "    Too few cells in the dimension", i_dim, ": "
                write(error_unit, *) "    ", this%num_total_cell_dim(i_dim), "<", &
                                             num_near_cells_dim(i_dim)
                stop
            end if
            
            Box_size_mod_cell_size = modulo(Box_size(i_dim), this%cell_size(i_dim))
            if (Box_size_mod_cell_size > real_zero .and. &
            abs(Box_size_mod_cell_size - this%cell_size(i_dim)) > real_zero) then
                write(error_unit, *) "    Cell size is not a divisor of the system size"
                write(error_unit, *) "    in the dimension", i_dim, ": "
                write(error_unit, *) "    Box_size", Box_size(i_dim)
                write(error_unit, *) "    cell_size", this%cell_size(i_dim)
                write(error_unit, *) "    modulo(Box_size, cell_size) = ", Box_size_mod_cell_size
                stop
            end if
            
        end do
        
    end subroutine Neighbour_Cells_check_cellsSize
    
    ! Assignment: particle -> cell
    
    pure function Neighbour_Cells_index_from_position(this, position) result(index_from_position)
    
        class(Neighbour_Cells), intent(in) :: this
        real(DP), dimension(:), intent(in) :: position
        integer :: index_from_position
        
        integer, dimension(num_dimensions) :: cell_coord
    
        cell_coord(:) = int(position(:)/this%cell_size(:)) + 1
        index_from_position = cell_coord(1) + this%num_total_cell_dim(1)*(cell_coord(2)-1) + &
                             this%num_total_cell_dim(1)*this%num_total_cell_dim(2)*(cell_coord(3)-1)
    
    end function Neighbour_Cells_index_from_position
    
    pure subroutine Neighbour_Cells_all_cols_to_cells(this, num_particles, spheres)
    
        class(Neighbour_Cells), intent(inout) :: this
        integer, intent(in) :: num_particles
        class(Hard_Spheres), intent(in) :: spheres
    
        integer :: i_particle
        integer :: i_cell
    
        do i_particle = 1, num_particles
    
            i_cell = this%index_from_position(spheres%get_position(i_particle))
            this%current_cells(i_cell)%particle%number = i_particle
            
            allocate(this%next_cells(i_cell)%particle)
            this%next_cells(i_cell)%particle%number = 0
            this%current_cells(i_cell)%particle%next => this%next_cells(i_cell)%particle
            this%current_cells(i_cell)%particle => this%next_cells(i_cell)%particle
            
        end do
        
        do i_cell = 1, this%num_total_cell
            this%current_cells(i_cell)%particle%next => null()
        end do
        
    end subroutine Neighbour_Cells_all_cols_to_cells

    pure function Neighbour_Cells_near_cell_bounds(this, i_total_cell) result(near_cell_bounds)

        class(Neighbour_Cells), intent(in) :: this
        integer, intent(in) :: i_total_cell
        integer, dimension(2) :: near_cell_bounds

        near_cell_bounds(1) = 1
        near_cell_bounds(2) = num_near_cells

         if (geometry%slab) then
            if (i_total_cell <= this%num_total_cell_layer) then
                near_cell_bounds(1) = num_near_cells_layer + 1
            else if ((this%num_total_cell-i_total_cell) < this%num_total_cell_layer) then
                near_cell_bounds(2) = (2) * num_near_cells_layer
            end if
        end if     

    end function Neighbour_Cells_near_cell_bounds
    
    subroutine Neighbour_Cells_point_to_begin(this, current, i_near_cell, i_total_cell)
    
        class(Neighbour_Cells), intent(in) :: this
        type(Node), pointer, intent(out) :: current
        integer, intent(in) :: i_near_cell, i_total_cell
        
        integer :: i_cell
                
        i_cell = this%near_among_total(i_near_cell, i_total_cell)        
        current => this%begin_cells(i_cell)%particle%next
        
    end subroutine Neighbour_Cells_point_to_begin
    
    ! Neighbour cells update
    
    subroutine Neighbour_Cells_remove_col_from_cell(this, i_particle, i_cellOld)
    
        class(Neighbour_Cells), intent(inout) :: this
        integer, intent(in) :: i_particle, i_cellOld
        
        type(Node), pointer :: current => null()
        type(Node), pointer :: next => null(), previous => null()
    
        previous => this%begin_cells(i_cellOld)%particle
        current => previous%next
        
        do
        
            next => current%next
        
            if (current%number == i_particle) then
                previous%next => current%next
                deallocate(current)
                current => next
                exit
            else
                previous => current
            end if
            
            current => next
        
        end do
            
    end subroutine Neighbour_Cells_remove_col_from_cell
    
    subroutine Neighbour_Cells_add_col_to_cell(this, i_particle, i_cellNew)
    
        class(Neighbour_Cells), intent(inout) :: this
        integer, intent(in) :: i_particle, i_cellNew
    
        type(Node), pointer :: new => null()
        type(Node), pointer :: next => null(), previous => null()
        
        previous => this%begin_cells(i_cellNew)%particle
        next => previous%next

        allocate(new)
        new%next => previous%next
        previous%next => new
        new%number = i_particle

    end subroutine Neighbour_Cells_add_col_to_cell
    
    ! Neighbour cells initialisation
    
    pure subroutine Neighbour_Cells_init_near_among_total(this)
    
        class(Neighbour_Cells), intent(inout) :: this
    
        integer :: i1_total_cell, i2_total_cell, i3_total_cell, i_total_cell
        integer :: i1_near_cell, i2_near_cell, i3_near_cell, i_near_cell
        integer :: i3_near_cell_min, i3_near_cell_max
        integer, dimension(num_dimensions) :: total_cell_coord, near_cell_coord
        
        do i3_total_cell = 1, this%num_total_cell_dim(3)

            i3_near_cell_min = 1
            i3_near_cell_max = num_near_cells_dim(3)

            if (geometry%slab) then
                if (i3_total_cell == 1) then
                    i3_near_cell_min = num_near_cells_dim(3) - 1
                    i3_near_cell_max = num_near_cells_dim(3)
                else if (i3_total_cell == this%num_total_cell_dim(3)) then
                    i3_near_cell_min = 1
                    i3_near_cell_max = 2
                end if
            end if
            
        do i2_total_cell = 1, this%num_total_cell_dim(2)
        do i1_total_cell = 1, this%num_total_cell_dim(1)
            
            i_total_cell = index_from_coord([i1_total_cell, i2_total_cell, i3_total_cell], &
                                            this%num_total_cell_dim)

            do i3_near_cell = i3_near_cell_min,  i3_near_cell_max
            do i2_near_cell = 1, num_near_cells_dim(2)
            do i1_near_cell = 1, num_near_cells_dim(1)
            
                near_cell_coord(:) = [i1_near_cell, i2_near_cell, i3_near_cell]
                i_near_cell = index_from_coord(near_cell_coord, num_near_cells_dim)
                near_cell_coord(:) = near_cell_coord(:) - num_near_cells_dim(:) + 1
                    ! symmetry: to the center
                
                total_cell_coord(:) = [i1_total_cell, i2_total_cell, i3_total_cell] + &
                                      near_cell_coord(:)
                total_cell_coord(:) = coord_PBC(total_cell_coord, this%num_total_cell_dim)
                
                this%near_among_total(i_near_cell, i_total_cell) = &
                    index_from_coord(total_cell_coord, this%num_total_cell_dim)
                    
            end do
            end do
            end do
        
        end do
        end do
        end do
            
    end subroutine Neighbour_Cells_init_near_among_total

end module class_neighbour_cells
