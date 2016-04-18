module class_visitable_cells

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions, real_zero
use data_cells, only: nums_local_cells
use procedures_errors, only: error_exit
use class_periodic_box, only: Abstract_Periodic_Box
use class_component_coordinates, only: Abstract_Component_Coordinates
use types_temporary_particle, only: Concrete_Temporary_Particle
use class_visitable_list, only: Abstract_Visitable_List
use class_pair_potential, only: Abstract_Pair_Potential

implicit none

private

    type, abstract, public :: Abstract_Visitable_Cells
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
        class(Abstract_Component_Coordinates), pointer :: positions => null()
        class(Abstract_Pair_Potential), pointer :: pair_potential => null()
        integer :: nums(num_dimensions) = 0
        real(DP) :: size(num_dimensions) = 0._DP
        integer, dimension(num_dimensions) :: global_lbounds = 0, global_ubounds = 0
        logical :: skip_bottom_layer(-nums_local_cells(3)/2:nums_local_cells(3)/2) = .false.
        logical :: skip_top_layer(-nums_local_cells(3)/2:nums_local_cells(3)/2) = .false.
        class(Abstract_Visitable_List), allocatable :: visitable_lists(:, :, :), list_mold
        integer, allocatable :: neighbours(:, :, :, :, :, :, :)
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: visit => Abstract_visit
        procedure :: move => Abstract_move
        procedure :: add => Abstract_add
        procedure :: remove => Abstract_remove
        procedure(Abstract_set_skip_layers), private, deferred :: set_skip_layers
        procedure, private :: create_space => Abstract_create_space
        procedure(Abstract_check_nums), private, deferred :: check_nums
        procedure, private :: check_size => Abstract_check_size
        procedure, private :: construct_visitable_lists => &
            Abstract_construct_visitable_lists
        procedure, private :: set_neighbours => Abstract_set_neighbours
        procedure, private :: index => Abstract_index
        procedure, private :: skip_local => Abstract_skip_local
        procedure, private :: fill_with_particles => Abstract_fill_with_particles
    end type Abstract_Visitable_Cells

    abstract interface

        subroutine Abstract_check_nums(this)
        import :: Abstract_Visitable_Cells
            class(Abstract_Visitable_Cells), intent(in) :: this
        end subroutine Abstract_check_nums

        pure subroutine Abstract_set_skip_layers(this)
        import :: Abstract_Visitable_Cells
            class(Abstract_Visitable_Cells), intent(inout) :: this
        end subroutine Abstract_set_skip_layers

    end interface

    type, extends(Abstract_Visitable_Cells), public :: Null_Visitable_Cells
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: visit => Null_visit
        procedure :: move => Null_move
        procedure :: add => Null_add
        procedure :: remove => Null_remove
        procedure, private :: check_nums => Null_check_nums
        procedure, private :: set_skip_layers => Null_set_skip_layers
        procedure, private :: fill_with_particles => Null_fill_with_particles
    end type Null_Visitable_Cells

    type, extends(Abstract_Visitable_Cells), public :: XYZ_PBC_Visitable_Cells
    contains
        procedure, private :: check_nums => XYZ_check_nums
        procedure, private :: set_skip_layers => XYZ_set_skip_layers
    end type XYZ_PBC_Visitable_Cells

    type, extends(Abstract_Visitable_Cells), public :: XY_PBC_Visitable_Cells
    private

    contains
        procedure, private :: check_nums => XY_check_nums
        procedure, private :: set_skip_layers => XY_set_skip_layers
    end type XY_PBC_Visitable_Cells

contains

!implementation Abstract_Visitable_Cells

    subroutine Abstract_construct(this, periodic_box, positions, pair_potential, list_mold)
        class(Abstract_Visitable_Cells), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Component_Coordinates), target, intent(in) :: positions
        class(Abstract_Pair_Potential), target, intent(in) :: pair_potential
        class(Abstract_Visitable_List), intent(in) :: list_mold

        this%periodic_box => periodic_box
        this%positions => positions
        this%pair_potential => pair_potential
        allocate(this%list_mold, mold=list_mold)
        call this%set_skip_layers()

        call this%create_space()
        call this%fill_with_particles()
    end subroutine Abstract_construct

    subroutine Abstract_create_space(this)
        class(Abstract_Visitable_Cells), intent(inout) :: this

        this%nums = floor(this%periodic_box%get_size()/this%pair_potential%get_max_distance())
        call this%check_nums()
        this%size = this%periodic_box%get_size() / real(this%nums, DP)
        call this%check_size()

        this%global_lbounds = -this%nums/2
        this%global_ubounds = this%global_lbounds + this%nums - 1
        allocate(this%visitable_lists(this%global_lbounds(1):this%global_ubounds(1), &
                                      this%global_lbounds(2):this%global_ubounds(2), &
                                      this%global_lbounds(3):this%global_ubounds(3)), &
                                      mold=this%list_mold)
        call this%construct_visitable_lists()

        allocate(this%neighbours(3, -nums_local_cells(1)/2:nums_local_cells(1)/2, &
                                    -nums_local_cells(2)/2:nums_local_cells(2)/2, &
                                    -nums_local_cells(3)/2:nums_local_cells(3)/2, &
                                    this%global_lbounds(1):this%global_ubounds(1), &
                                    this%global_lbounds(2):this%global_ubounds(2), &
                                    this%global_lbounds(3):this%global_ubounds(3)))
        call this%set_neighbours()
    end subroutine Abstract_create_space

    subroutine Abstract_check_size(this)
        class(Abstract_Visitable_Cells), intent(in) :: this

        real(DP) :: box_mod_cell(num_dimensions)

        box_mod_cell = modulo(this%periodic_box%get_size(), this%size)
        if (any(box_mod_cell > real_zero .and. abs(box_mod_cell - this%size) > real_zero)) then
            call error_exit("Abstract_Visitable_Cells: check_size: "//&
                            "this%size is not a divisor of periodic_box%get_size()")
        end if
    end subroutine Abstract_check_size

    subroutine Abstract_construct_visitable_lists(this)
        class(Abstract_Visitable_Cells), intent(inout) :: this

        integer :: global_i1, global_i2, global_i3

        do global_i3 = this%global_lbounds(3), this%global_ubounds(3)
        do global_i2 = this%global_lbounds(2), this%global_ubounds(2)
        do global_i1 = this%global_lbounds(1), this%global_ubounds(1)
            call this%visitable_lists(global_i1, global_i2, global_i3)%construct(this%periodic_box,&
                this%positions)
        end do
        end do
        end do
    end subroutine Abstract_construct_visitable_lists

    subroutine Abstract_set_neighbours(this)
        class(Abstract_Visitable_Cells), intent(inout) :: this

        integer :: global_i1, global_i2, global_i3
        logical :: bottom_layer, top_layer
        integer :: local_i1, local_i2, local_i3
        integer :: ijk_cell(num_dimensions)

        this%neighbours = 0
        do global_i3 = this%global_lbounds(3), this%global_ubounds(3)
            bottom_layer= (global_i3 == this%global_lbounds(3))
            top_layer = (global_i3 == this%global_ubounds(3))
        do global_i2 = this%global_lbounds(2), this%global_ubounds(2)
        do global_i1 = this%global_lbounds(1), this%global_ubounds(1)
            do local_i3 = -nums_local_cells(3)/2, nums_local_cells(3)/2
                if (this%skip_local(bottom_layer, top_layer, local_i3)) cycle
            do local_i2 = -nums_local_cells(2)/2, nums_local_cells(2)/2
            do local_i1 = -nums_local_cells(1)/2, nums_local_cells(1)/2
                ijk_cell = [global_i1, global_i2, global_i3] + [local_i1, local_i2, local_i3]
                ijk_cell = pbc_3d_index(ijk_cell, this%nums)
                this%neighbours(:, local_i1, local_i2, local_i3, &
                    global_i1, global_i2, global_i3) = ijk_cell
            end do
            end do
            end do
        end do
        end do
        end do
    end subroutine Abstract_set_neighbours

    pure function pbc_3d_index(ijk_cell, nums_cells)
        integer, intent(in) :: ijk_cell(:), nums_cells(:)
        integer :: pbc_3d_index(3)

        pbc_3d_index = modulo(ijk_cell + nums_cells/2, nums_cells) - nums_cells/2
    end function pbc_3d_index

    pure function Abstract_skip_local(this, bottom_layer, top_layer, local_i3) result(skip_local)
        class(Abstract_Visitable_Cells), intent(in) :: this
        logical, intent(in) :: bottom_layer, top_layer
        integer, intent(in) :: local_i3
        logical :: skip_local

        skip_local = (bottom_layer .and. this%skip_bottom_layer(local_i3)) .or. &
            (top_layer .and. this%skip_top_layer(local_i3))
    end function Abstract_skip_local

    subroutine Abstract_fill_with_particles(this)
        class(Abstract_Visitable_Cells), intent(inout) :: this

        type(Concrete_Temporary_Particle) :: particle
        integer :: i_particle

        do i_particle = 1, this%positions%get_num()
            particle%i = i_particle
            particle%position = this%positions%get(particle%i)
            call this%add(particle)
        end do
    end subroutine Abstract_fill_with_particles

    subroutine Abstract_destroy(this)
        class(Abstract_Visitable_Cells), intent(inout) :: this

        integer :: global_i1, global_i2, global_i3

        do global_i3 = this%global_ubounds(3), this%global_lbounds(3), -1
        do global_i2 = this%global_ubounds(2), this%global_lbounds(2), -1
        do global_i1 = this%global_ubounds(1), this%global_lbounds(1), -1
            call this%visitable_lists(global_i1, global_i2, global_i3)%destroy()
        end do
        end do
        end do
        if (allocated(this%neighbours)) deallocate(this%neighbours)
        if (allocated(this%visitable_lists)) deallocate(this%visitable_lists)
        if (allocated(this%list_mold)) deallocate(this%list_mold)
        this%pair_potential => null()
        this%positions => null()
        this%periodic_box => null()
    end subroutine Abstract_destroy

    subroutine Abstract_visit(this, overlap, energy, particle, i_exclude)
        class(Abstract_Visitable_Cells), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy
        type(Concrete_Temporary_Particle), intent(in) :: particle
        integer, intent(in) :: i_exclude

        real(DP) :: energy_i
        integer, dimension(num_dimensions) :: ijk_cell, ijk_local_cell
        logical :: bottom_layer, top_layer
        integer :: local_i1, local_i2, local_i3

        ijk_cell = this%index(particle%position)
        bottom_layer = (ijk_cell(3) == this%global_lbounds(3))
        top_layer = (ijk_cell(3) == this%global_ubounds(3))
        energy = 0._DP
        do local_i3 = -nums_local_cells(3)/2, nums_local_cells(3)/2
            if (this%skip_local(bottom_layer, top_layer, local_i3)) cycle
        do local_i2 = -nums_local_cells(2)/2, nums_local_cells(2)/2
        do local_i1 = -nums_local_cells(1)/2, nums_local_cells(1)/2
            ijk_local_cell = this%neighbours(:, local_i1, local_i2, local_i3, &
                ijk_cell(1), ijk_cell(2), ijk_cell(3))
            call this%visitable_lists(ijk_local_cell(1), ijk_local_cell(2), ijk_local_cell(3))%&
                visit(overlap, energy_i, particle, this%pair_potential, i_exclude)
            if (overlap) return
            energy = energy + energy_i
        end do
        end do
        end do
    end subroutine Abstract_visit

    subroutine Abstract_move(this, to_position, from)
        class(Abstract_Visitable_Cells), intent(inout) :: this
        real(DP), intent(in) :: to_position(:)
        type(Concrete_Temporary_Particle), intent(in) :: from

        integer, dimension(num_dimensions) :: from_ijk_cell, to_ijk_cell

        from_ijk_cell = this%index(from%position)
        to_ijk_cell = this%index(to_position)
        if (any(from_ijk_cell /= to_ijk_cell)) then
            call this%visitable_lists(from_ijk_cell(1), from_ijk_cell(2), from_ijk_cell(3))%&
                remove(from%i)
            call this%visitable_lists(to_ijk_cell(1), to_ijk_cell(2), to_ijk_cell(3))%add(from%i)
        end if
    end subroutine Abstract_move

    subroutine Abstract_add(this, particle)
        class(Abstract_Visitable_Cells), intent(inout) :: this
        type(Concrete_Temporary_Particle), intent(in) :: particle

        integer :: ijk_cell(num_dimensions)

        ijk_cell = this%index(particle%position)
        call this%visitable_lists(ijk_cell(1), ijk_cell(2), ijk_cell(3))%add(particle%i)
    end subroutine Abstract_add

    subroutine Abstract_remove(this, particle)
        class(Abstract_Visitable_Cells), intent(inout) :: this
        type(Concrete_Temporary_Particle), intent(in) :: particle

        integer :: ijk_cell(num_dimensions)

        ijk_cell = this%index(particle%position)
        call this%visitable_lists(ijk_cell(1), ijk_cell(2), ijk_cell(3))%remove(particle%i)
        if (particle%i < this%positions%get_num()) then
            call this%visitable_lists(ijk_cell(1), ijk_cell(2), &
                ijk_cell(3))%set(this%positions%get_num(), particle%i)
        end if
    end subroutine Abstract_remove

    pure function Abstract_index(this, position) result(index)
        class(Abstract_Visitable_Cells), intent(in) :: this
        real(DP), intent(in) :: position(:)
        integer :: index(num_dimensions)

        where (mod(this%nums, 2) == 0)
            index = floor(position/this%size)
        elsewhere
            index = nint(position/this%size)
        end where
    end function Abstract_index

!end implementation Abstract_Visitable_Cells

!implementation Null_Visitable_Cells

    subroutine Null_construct(this, periodic_box, positions, pair_potential, list_mold)
        class(Null_Visitable_Cells), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Component_Coordinates), target, intent(in) :: positions
        class(Abstract_Pair_Potential), target, intent(in) :: pair_potential
        class(Abstract_Visitable_List), intent(in) :: list_mold
    end subroutine Null_construct

    subroutine Null_fill_with_particles(this)
        class(Null_Visitable_Cells), intent(inout) :: this
    end subroutine Null_fill_with_particles

    subroutine Null_destroy(this)
        class(Null_Visitable_Cells), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_visit(this, overlap, energy, particle, i_exclude)
        class(Null_Visitable_Cells), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy
        type(Concrete_Temporary_Particle), intent(in) :: particle
        integer, intent(in) :: i_exclude
        overlap = .false.
        energy = 0._DP
    end subroutine Null_visit

    subroutine Null_move(this, to_position, from)
        class(Null_Visitable_Cells), intent(inout) :: this
        real(DP), intent(in) :: to_position(:)
        type(Concrete_Temporary_Particle), intent(in) :: from
    end subroutine Null_move

    subroutine Null_add(this, particle)
        class(Null_Visitable_Cells), intent(inout) :: this
        type(Concrete_Temporary_Particle), intent(in) :: particle
    end subroutine Null_add

    subroutine Null_remove(this, particle)
        class(Null_Visitable_Cells), intent(inout) :: this
        type(Concrete_Temporary_Particle), intent(in) :: particle
    end subroutine Null_remove

    subroutine Null_check_nums(this)
        class(Null_Visitable_Cells), intent(in) :: this
    end subroutine Null_check_nums

    pure subroutine Null_set_skip_layers(this)
        class(Null_Visitable_Cells), intent(inout) :: this
    end subroutine Null_set_skip_layers

!end implementation Null_Visitable_Cells

!implementation XYZ_PBC_Visitable_Cells

    pure subroutine XYZ_set_skip_layers(this)
        class(XYZ_PBC_Visitable_Cells), intent(inout) :: this

        this%skip_top_layer = .false.
        this%skip_bottom_layer = .false.
    end subroutine XYZ_set_skip_layers

    subroutine XYZ_check_nums(this)
        class(XYZ_PBC_Visitable_Cells), intent(in) :: this

        if (any(this%nums < nums_local_cells)) then
            call error_exit("XYZ_PBC_Visitable_Cells: this%nums is too small.")
        end if
    end subroutine XYZ_check_nums

!end implementation XYZ_PBC_Visitable_Cells

!implementation XY_PBC_Visitable_Cells

    pure subroutine XY_set_skip_layers(this)
        class(XY_PBC_Visitable_Cells), intent(inout) :: this

        integer :: local_i3

        do local_i3 = -nums_local_cells(3)/2, nums_local_cells(3)/2
            if (local_i3 == 1) then
                this%skip_top_layer(local_i3) = .true.
            else
                this%skip_top_layer(local_i3) = .false.
            end if
            if (local_i3 == -1) then
                this%skip_bottom_layer(local_i3) = .true.
            else
                this%skip_bottom_layer(local_i3) = .false.
            end if
        end do
    end subroutine XY_set_skip_layers

    subroutine XY_check_nums(this)
        class(XY_PBC_Visitable_Cells), intent(in) :: this

        if (any(this%nums(1:2) < nums_local_cells(1:2))) then
            call error_exit("XY_PBC_Visitable_Cells: this%nums is too small.")
        end if
    end subroutine XY_check_nums

!end implementation XY_PBC_Visitable_Cells

end module class_visitable_cells
