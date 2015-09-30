module class_visitable_cells

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_precisions, only: real_zero
use data_geometry, only: num_dimensions
use procedures_errors, only: error_exit
use procedures_checks, only: check_positive
use class_periodic_box, only: Abstract_Periodic_Box
use class_particles_positions, only: Abstract_Particles_Positions
use types_particle, only: Concrete_Particle
use class_visitable_list, only: Abstract_Visitable_List
use procedures_visitable_cells, only: pbc_3d_index, local_reindex
use class_pair_potential, only: Abstract_Pair_Potential

implicit none

private

    integer, parameter :: nums_local_cells(num_dimensions) = 3

    type, abstract, public :: Abstract_Visitable_Cells
    private
        integer :: nums(num_dimensions)
        real(DP) :: size(num_dimensions)
        integer, dimension(num_dimensions) :: global_lbounds, global_ubounds
        logical :: skip_bottom_layer(nums_local_cells(1), nums_local_cells(2), nums_local_cells(3))
        logical :: skip_top_layer(nums_local_cells(1), nums_local_cells(2), nums_local_cells(3))
        class(Abstract_Visitable_List), allocatable :: visitable_lists(:, :, :)
        integer, allocatable :: neighbours(:, :, :, :, :, :, :)
        class(Abstract_Particles_Positions), pointer :: positions
        class(Abstract_Periodic_Box), pointer :: periodic_box
        class(Abstract_Pair_Potential), pointer :: pair_potential
    contains
        procedure :: construct => Abstract_Visitable_Cells_construct
        procedure :: destroy => Abstract_Visitable_Cells_destroy
        procedure :: visit => Abstract_Visitable_Cells_visit
        procedure :: move => Abstract_Visitable_Cells_move
        procedure :: add => Abstract_Visitable_Cells_add
        procedure :: remove => Abstract_Visitable_Cells_remove
        procedure, private :: set_nums => Abstract_Visitable_Cells_set_nums
        procedure(Abstract_Visitable_Cells_check_nums), private, deferred :: check_nums
        procedure, private :: set_division => Abstract_Visitable_Cells_set_division
        procedure, private :: check_division => Abstract_Visitable_Cells_check_division
        procedure, private :: construct_visitable_lists => &
            Abstract_Visitable_Cells_construct_visitable_lists
        procedure, private :: set_neighbours => Abstract_Visitable_Cells_set_neighbours
        procedure, private :: index => Abstract_Visitable_Cells_index
        procedure(Abstract_Visitable_Cells_set_skip_layers), private, deferred :: set_skip_layers
        procedure, private :: skip_local => Abstract_Visitable_Cells_skip_local
        procedure, private :: fill => Abstract_Visitable_Cells_fill
    end type Abstract_Visitable_Cells

    abstract interface

        subroutine Abstract_Visitable_Cells_check_nums(this)
        import :: Abstract_Visitable_Cells
            class(Abstract_Visitable_Cells), intent(in) :: this
        end subroutine Abstract_Visitable_Cells_check_nums

        pure subroutine Abstract_Visitable_Cells_set_skip_layers(this)
        import :: Abstract_Visitable_Cells
            class(Abstract_Visitable_Cells), intent(inout) :: this
        end subroutine Abstract_Visitable_Cells_set_skip_layers

    end interface

    type, extends(Abstract_Visitable_Cells), public :: Null_Visitable_Cells
    contains
        procedure :: construct => Null_Visitable_Cells_construct
        procedure :: destroy => Null_Visitable_Cells_destroy
        procedure :: visit => Null_Visitable_Cells_visit
        procedure :: move => Null_Visitable_Cells_move
        procedure :: add => Null_Visitable_Cells_add
        procedure :: remove => Null_Visitable_Cells_remove
        procedure, private :: check_nums => Null_Visitable_Cells_check_nums
        procedure, private :: set_skip_layers => Null_Visitable_Cells_set_skip_layers
        procedure, private :: fill => Null_Visitable_Cells_fill
    end type Null_Visitable_Cells

    type, extends(Abstract_Visitable_Cells), public :: XYZ_PBC_Visitable_Cells
    contains
        procedure, private :: check_nums => XYZ_PBC_Visitable_Cells_check_nums
        procedure, private :: set_skip_layers => XYZ_PBC_Visitable_Cells_set_skip_layers
    end type XYZ_PBC_Visitable_Cells

    type, extends(Abstract_Visitable_Cells), public :: XY_PBC_Visitable_Cells
    private

    contains
        procedure, private :: check_nums => XY_PBC_Visitable_Cells_check_nums
        procedure, private :: set_skip_layers => XY_PBC_Visitable_Cells_set_skip_layers
    end type XY_PBC_Visitable_Cells

contains

!implementation Abstract_Visitable_Cells

    subroutine Abstract_Visitable_Cells_construct(this, mold, periodic_box, positions, &
            pair_potential)
        class(Abstract_Visitable_Cells), intent(out) :: this
        class(Abstract_Visitable_List), intent(in) :: mold
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Particles_Positions), target, intent(in) :: positions
        class(Abstract_Pair_Potential), target, intent(in) :: pair_potential

        this%periodic_box => periodic_box
        this%positions => positions
        this%pair_potential => pair_potential
        call this%set_nums()
        call this%set_division()
        this%global_lbounds = -this%nums/2
        this%global_ubounds = this%global_lbounds + this%nums - 1
        allocate(this%visitable_lists(this%global_lbounds(1):this%global_ubounds(1), &
                                      this%global_lbounds(2):this%global_ubounds(2), &
                                      this%global_lbounds(3):this%global_ubounds(3)), &
                                      mold=mold)
        call this%construct_visitable_lists(periodic_box)

        call this%set_skip_layers()
        allocate(this%neighbours(3, nums_local_cells(1), nums_local_cells(2), &
                                    nums_local_cells(3), &
                                    this%global_lbounds(1):this%global_ubounds(1), &
                                    this%global_lbounds(2):this%global_ubounds(2), &
                                    this%global_lbounds(3):this%global_ubounds(3)))
        call this%set_neighbours()
        call this%fill()
    end subroutine Abstract_Visitable_Cells_construct

    subroutine Abstract_Visitable_Cells_set_nums(this)
        class(Abstract_Visitable_Cells), intent(inout) :: this

        call check_positive("Abstract_Visitable_Cells", "this%pair_potential%get_max_distance()", &
            this%pair_potential%get_max_distance())
        this%nums = floor(this%periodic_box%get_size()/this%pair_potential%get_max_distance())
        call check_positive("Abstract_Visitable_Cells", "this%nums", this%nums)
        call this%check_nums()
    end subroutine Abstract_Visitable_Cells_set_nums

    subroutine Abstract_Visitable_Cells_set_division(this)
        class(Abstract_Visitable_Cells), intent(inout) :: this

        this%size = this%periodic_box%get_size() / real(this%nums, DP)
        call check_positive("Abstract_Visitable_Cells", "this%size", this%size)
        call this%check_division()
    end subroutine Abstract_Visitable_Cells_set_division

    subroutine Abstract_Visitable_Cells_check_division(this)
        class(Abstract_Visitable_Cells), intent(in) :: this

        real(DP) :: box_mod_cell(num_dimensions)

        box_mod_cell = modulo(this%periodic_box%get_size(), this%size)
        if (any(box_mod_cell > real_zero .and. abs(box_mod_cell - this%size) > real_zero)) then
            call error_exit("Abstract_Visitable_Cells:"//&
                            "this%size is not a divisor of periodic_box%get_size()")
        end if
    end subroutine Abstract_Visitable_Cells_check_division

    subroutine Abstract_Visitable_Cells_construct_visitable_lists(this, periodic_box)
        class(Abstract_Visitable_Cells), intent(inout) :: this
        class(Abstract_Periodic_Box), intent(in) :: periodic_box

        integer :: global_i1, global_i2, global_i3

        do global_i3 = this%global_lbounds(3), this%global_ubounds(3)
        do global_i2 = this%global_lbounds(2), this%global_ubounds(2)
        do global_i1 = this%global_lbounds(1), this%global_ubounds(1)
            call this%visitable_lists(global_i1, global_i2, global_i3)%construct(periodic_box)
        end do
        end do
        end do
    end subroutine Abstract_Visitable_Cells_construct_visitable_lists

    subroutine Abstract_Visitable_Cells_set_neighbours(this)
        class(Abstract_Visitable_Cells), intent(inout) :: this

        integer :: global_i1, global_i2, global_i3
        logical :: bottom_layer, top_layer
        integer :: local_i1, local_i2, local_i3
        integer :: i_cell(num_dimensions)

        this%neighbours = 0
        do global_i3 = this%global_lbounds(3), this%global_ubounds(3)
            bottom_layer= (global_i3 == this%global_lbounds(3))
            top_layer = (global_i3 == this%global_ubounds(3))
        do global_i2 = this%global_lbounds(2), this%global_ubounds(2)
        do global_i1 = this%global_lbounds(1), this%global_ubounds(1)
            do local_i3 = 1, nums_local_cells(3)
            do local_i2 = 1, nums_local_cells(2)
            do local_i1 = 1, nums_local_cells(1)
                if (this%skip_local(bottom_layer, top_layer, local_i1, local_i2, local_i3)) cycle
                i_cell = [global_i1, global_i2, global_i3] + &
                    local_reindex([local_i1, local_i2, local_i3], nums_local_cells)
                i_cell = pbc_3d_index(i_cell, this%nums)
                this%neighbours(:, local_i1, local_i2, local_i3, &
                    global_i1, global_i2, global_i3) = i_cell
            end do
            end do
            end do
        end do
        end do
        end do
    end subroutine Abstract_Visitable_Cells_set_neighbours

    pure function Abstract_Visitable_Cells_skip_local(this, bottom_layer, top_layer, &
        local_i1, local_i2, local_i3) result(skip_local)
        class(Abstract_Visitable_Cells), intent(in) :: this
        logical, intent(in) :: bottom_layer, top_layer
        integer, intent(in) :: local_i1, local_i2, local_i3
        logical :: skip_local

        skip_local = (bottom_layer .and. this%skip_bottom_layer(local_i1, local_i2, local_i3)) &
            .or. (top_layer .and. this%skip_top_layer(local_i1, local_i2, local_i3))
    end function Abstract_Visitable_Cells_skip_local

    subroutine Abstract_Visitable_Cells_fill(this)
        class(Abstract_Visitable_Cells), intent(inout) :: this

        type(Concrete_Particle) :: particle
        integer :: i_particle

        do i_particle = 1, this%positions%get_num()
            particle%i = i_particle
            particle%position = this%positions%get(particle%i)
            call this%add(particle)
        end do
    end subroutine Abstract_Visitable_Cells_fill

    subroutine Abstract_Visitable_Cells_destroy(this)
        class(Abstract_Visitable_Cells), intent(inout) :: this

        integer :: global_i1, global_i2, global_i3

        this%pair_potential => null()
        this%positions => null()

        do global_i3 = this%global_ubounds(3), this%global_lbounds(3), -1
        do global_i2 = this%global_ubounds(2), this%global_lbounds(2), -1
        do global_i1 = this%global_ubounds(1), this%global_lbounds(1), -1
            call this%visitable_lists(global_i1, global_i2, global_i3)%destroy()
        end do
        end do
        end do

        if (allocated(this%neighbours)) deallocate(this%neighbours)
        if (allocated(this%visitable_lists)) deallocate(this%visitable_lists)
    end subroutine Abstract_Visitable_Cells_destroy

    subroutine Abstract_Visitable_Cells_visit(this, overlap, energy, particle)
        class(Abstract_Visitable_Cells), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy
        type(Concrete_Particle), intent(in) :: particle

        real(DP) :: energy_i
        integer, dimension(num_dimensions) :: i_cell, i_local_cell
        logical :: bottom_layer, top_layer
        integer :: local_i1, local_i2, local_i3

        i_cell = this%index(particle%position)
        bottom_layer= (i_cell(3) == this%global_lbounds(3))
        top_layer = (i_cell(3) == this%global_ubounds(3))
        energy = 0._DP
        do local_i3 = 1, nums_local_cells(3)
        do local_i2 = 1, nums_local_cells(2)
        do local_i1 = 1, nums_local_cells(1)
            if (this%skip_local(bottom_layer, top_layer, local_i1, local_i2, local_i3)) cycle
            i_local_cell = this%neighbours(:, local_i1, local_i2, local_i3, &
                i_cell(1), i_cell(2), i_cell(3))
            call this%visitable_lists(i_local_cell(1), i_local_cell(2), &
                i_local_cell(3))%visit(overlap, energy_i, particle, this%pair_potential)
            if (overlap) return
            energy = energy + energy_i
        end do
        end do
        end do
    end subroutine Abstract_Visitable_Cells_visit

    subroutine Abstract_Visitable_Cells_move(this, from, to)
        class(Abstract_Visitable_Cells), intent(inout) :: this
        type(Concrete_Particle), intent(in) :: from, to

        integer, dimension(num_dimensions) :: from_i_cell, to_i_cell

        from_i_cell = this%index(from%position)
        to_i_cell = this%index(to%position)
        if (any(from_i_cell /= to_i_cell)) then
            call this%visitable_lists(from_i_cell(1), from_i_cell(2), &
                from_i_cell(3))%remove(from%i)
            call this%visitable_lists(to_i_cell(1), to_i_cell(2), &
                to_i_cell(3))%add(to)
        else
            call this%visitable_lists(from_i_cell(1), from_i_cell(2), &
                from_i_cell(3))%set(from%i, to)
        end if
    end subroutine Abstract_Visitable_Cells_move

    subroutine Abstract_Visitable_Cells_add(this, particle)
        class(Abstract_Visitable_Cells), intent(inout) :: this
        type(Concrete_Particle), intent(in) :: particle

        integer :: i_cell(num_dimensions)

        i_cell = this%index(particle%position)
        call this%visitable_lists(i_cell(1), i_cell(2), i_cell(3))%add(particle)
    end subroutine Abstract_Visitable_Cells_add

    subroutine Abstract_Visitable_Cells_remove(this, particle)
        class(Abstract_Visitable_Cells), intent(inout) :: this
        type(Concrete_Particle), intent(in) :: particle

        integer :: i_cell(num_dimensions)

        i_cell = this%index(particle%position)
        call this%visitable_lists(i_cell(1), i_cell(2), i_cell(3))%remove(particle%i)
        if (particle%i < this%positions%get_num()) then
            call this%visitable_lists(i_cell(1), i_cell(2), &
                i_cell(3))%set(this%positions%get_num(), particle)
        end if
    end subroutine Abstract_Visitable_Cells_remove

    pure function Abstract_Visitable_Cells_index(this, position) result(index)
        class(Abstract_Visitable_Cells), intent(in) :: this
        real(DP), intent(in) :: position(:)
        integer :: index(num_dimensions)

        where (mod(this%nums, 2) == 0)
            index = floor(position/this%size)
        elsewhere
            index = nint(position/this%size)
        end where
    end function Abstract_Visitable_Cells_index

!end implementation Abstract_Visitable_Cells

!implementation Null_Visitable_Cells

    subroutine Null_Visitable_Cells_construct(this, mold, periodic_box, positions, &
            pair_potential)
        class(Null_Visitable_Cells), intent(out) :: this
        class(Abstract_Visitable_List), intent(in) :: mold
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Particles_Positions), target, intent(in) :: positions
        class(Abstract_Pair_Potential), target, intent(in) :: pair_potential
    end subroutine Null_Visitable_Cells_construct

    subroutine Null_Visitable_Cells_fill(this)
        class(Null_Visitable_Cells), intent(inout) :: this
    end subroutine Null_Visitable_Cells_fill

    subroutine Null_Visitable_Cells_destroy(this)
        class(Null_Visitable_Cells), intent(inout) :: this
    end subroutine Null_Visitable_Cells_destroy

    subroutine Null_Visitable_Cells_visit(this, overlap, energy, particle)
        class(Null_Visitable_Cells), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy
        type(Concrete_Particle), intent(in) :: particle
        overlap = .false.
        energy = 0._DP
    end subroutine Null_Visitable_Cells_visit

    subroutine Null_Visitable_Cells_move(this, from, to)
        class(Null_Visitable_Cells), intent(inout) :: this
        type(Concrete_Particle), intent(in) :: from, to
    end subroutine Null_Visitable_Cells_move

    subroutine Null_Visitable_Cells_add(this, particle)
        class(Null_Visitable_Cells), intent(inout) :: this
        type(Concrete_Particle), intent(in) :: particle
    end subroutine Null_Visitable_Cells_add

    subroutine Null_Visitable_Cells_remove(this, particle)
        class(Null_Visitable_Cells), intent(inout) :: this
        type(Concrete_Particle), intent(in) :: particle
    end subroutine Null_Visitable_Cells_remove

    subroutine Null_Visitable_Cells_check_nums(this)
        class(Null_Visitable_Cells), intent(in) :: this
    end subroutine Null_Visitable_Cells_check_nums

    pure subroutine Null_Visitable_Cells_set_skip_layers(this)
        class(Null_Visitable_Cells), intent(inout) :: this
    end subroutine Null_Visitable_Cells_set_skip_layers

!end implementation Null_Visitable_Cells

!implementation XYZ_PBC_Visitable_Cells

    subroutine XYZ_PBC_Visitable_Cells_check_nums(this)
        class(XYZ_PBC_Visitable_Cells), intent(in) :: this

        if (any(this%nums < nums_local_cells)) then
            call error_exit("XYZ_PBC_Visitable_Cells: this%nums is too small.")
        end if
    end subroutine XYZ_PBC_Visitable_Cells_check_nums

    pure subroutine XYZ_PBC_Visitable_Cells_set_skip_layers(this)
        class(XYZ_PBC_Visitable_Cells), intent(inout) :: this

        this%skip_top_layer = .false.
        this%skip_bottom_layer = .false.
    end subroutine XYZ_PBC_Visitable_Cells_set_skip_layers

!end implementation XYZ_PBC_Visitable_Cells

!implementation XY_PBC_Visitable_Cells

    subroutine XY_PBC_Visitable_Cells_check_nums(this)
        class(XY_PBC_Visitable_Cells), intent(in) :: this

        if (any(this%nums(1:2) < nums_local_cells(1:2))) then
            call error_exit("XY_PBC_Visitable_Cells: this%nums is too small.")
        end if
    end subroutine XY_PBC_Visitable_Cells_check_nums

    pure subroutine XY_PBC_Visitable_Cells_set_skip_layers(this)
        class(XY_PBC_Visitable_Cells), intent(inout) :: this

        integer :: i_local_cell(num_dimensions)
        integer :: local_i1, local_i2, local_i3

        do local_i3 = 1, nums_local_cells(3)
            do local_i2 = 1, nums_local_cells(2)
                do local_i1 = 1, nums_local_cells(1)
                    i_local_cell = local_reindex([local_i1, local_i2, local_i3], nums_local_cells)
                    if (i_local_cell(3) == 1) then
                        this%skip_top_layer(local_i1, local_i2, local_i3) = .true.
                    else
                        this%skip_top_layer(local_i1, local_i2, local_i3) = .false.
                    end if
                    if (i_local_cell(3) == -1) then
                        this%skip_bottom_layer(local_i1, local_i2, local_i3) = .true.
                    else
                        this%skip_bottom_layer(local_i1, local_i2, local_i3) = .false.
                    end if
                end do
            end do
        end do
    end subroutine XY_PBC_Visitable_Cells_set_skip_layers

!end implementation XY_PBC_Visitable_Cells

end module class_visitable_cells
