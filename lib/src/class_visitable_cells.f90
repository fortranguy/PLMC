module class_visitable_cells

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_precisions, only: real_zero
use data_geometry, only: num_dimensions
use procedures_errors, only: error_exit
use procedures_checks, only: check_positive
use class_periodic_box, only: Abstract_Periodic_Box
use class_positions, only: Abstract_Positions
use module_particles, only: Concrete_Particle
use class_pair_potential, only: Abstract_Pair_Potential
use class_visitable_list, only: Abstract_Visitable_List
use procedures_cells, only: pbc_3d_index

implicit none

private

    integer, parameter :: num_local_cells(num_dimensions) = 3

    type, abstract, public :: Abstract_Visitable_Cells
    private
        integer :: nums(num_dimensions)
        real(DP) :: size(num_dimensions)
        integer, dimension(num_dimensions) :: global_lbounds, global_ubounds
        integer, dimension(num_dimensions) :: local_lbounds, local_ubounds
        class(Abstract_Visitable_List), allocatable :: visitable_lists(:, :, :)
        integer, allocatable :: neighbours(:, :, :, :, :, :, :)
        class(Abstract_Positions), pointer :: positions => null()
    contains
        procedure :: construct => Abstract_Visitable_Cells_construct
        procedure :: destroy => Abstract_Visitable_Cells_destroy
        procedure(Abstract_Visitable_Cells_visit), deferred :: visit
        procedure :: move => Abstract_Visitable_Cells_move
        procedure :: add => Abstract_Visitable_Cells_add
        procedure :: remove => Abstract_Visitable_Cells_remove
        procedure, private :: check_size => Abstract_Visitable_Cells_check_size
        procedure, private :: construct_visitable_lists => &
            Abstract_Visitable_Cells_construct_visitable_lists
        procedure(Abstract_Visitable_Cells_set_neighbours), private, deferred :: set_neighbours
        procedure, private :: fill => Abstract_Visitable_Cells_fill
        procedure, private :: index => Abstract_Visitable_Cells_index
    end type Abstract_Visitable_Cells

    abstract interface

        subroutine Abstract_Visitable_Cells_set_neighbours(this)
        import :: Abstract_Visitable_Cells
            class(Abstract_Visitable_Cells), intent(inout) :: this
        end subroutine Abstract_Visitable_Cells_set_neighbours

        subroutine Abstract_Visitable_Cells_visit(this, particle, pair_potential, overlap, energy)
        import :: DP, Abstract_Pair_Potential, Concrete_Particle, Abstract_Visitable_Cells
            class(Abstract_Visitable_Cells), intent(in) :: this
            type(Concrete_Particle), intent(in) :: particle
            class(Abstract_Pair_Potential), intent(in) :: pair_potential
            logical, intent(out) :: overlap
            real(DP), intent(out) :: energy
        end subroutine Abstract_Visitable_Cells_visit

    end interface

    type, extends(Abstract_Visitable_Cells), public :: PBC_3D_Visitable_Cells
    contains
        procedure :: set_neighbours => PBC_3D_Visitable_Cells_set_neighbours
        procedure :: visit => PBC_3D_Visitable_Cells_visit
    end type

contains

!implementation Abstract_Visitable_Cells

    subroutine Abstract_Visitable_Cells_construct(this, mold, periodic_box, positions, &
            min_cell_edge)
        class(Abstract_Visitable_Cells), intent(out) :: this
        class(Abstract_Visitable_List), intent(in) :: mold
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Positions), target, intent(in) :: positions
        real(DP), intent(in) :: min_cell_edge

        this%nums = floor(periodic_box%get_size()/min_cell_edge)
        call check_positive("Abstract_Visitable_Cells", "this%nums", this%nums)
        this%size = periodic_box%get_size() / real(this%nums, DP)
        call check_positive("Abstract_Visitable_Cells", "this%size", this%size)
        call this%check_size(periodic_box)

        this%global_lbounds = -this%nums/2
        this%global_ubounds = this%global_lbounds + this%nums - 1
        allocate(this%visitable_lists(this%global_lbounds(1):this%global_ubounds(1), &
                                      this%global_lbounds(2):this%global_ubounds(2), &
                                      this%global_lbounds(3):this%global_ubounds(3)), &
                                      mold=mold)
        call this%construct_visitable_lists(periodic_box, positions)

        this%local_lbounds = -num_local_cells/2
        this%local_ubounds = this%local_lbounds + num_local_cells - 1
        allocate(this%neighbours(3, this%local_lbounds(1):this%local_ubounds(1), &
                                    this%local_lbounds(2):this%local_ubounds(2), &
                                    this%local_lbounds(3):this%local_ubounds(3), &
                                    this%global_lbounds(1):this%global_ubounds(1), &
                                    this%global_lbounds(2):this%global_ubounds(2), &
                                    this%global_lbounds(3):this%global_ubounds(3)))
        call this%set_neighbours()

        this%positions => positions
        call this%fill()
    end subroutine Abstract_Visitable_Cells_construct

    subroutine Abstract_Visitable_Cells_construct_visitable_lists(this, periodic_box, positions)
        class(Abstract_Visitable_Cells), intent(inout) :: this
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Positions), intent(in) :: positions

        integer :: global_i1, global_i2, global_i3

        do global_i3 = this%global_lbounds(3), this%global_ubounds(3)
        do global_i2 = this%global_lbounds(2), this%global_ubounds(2)
        do global_i1 = this%global_lbounds(1), this%global_ubounds(1)
            call this%visitable_lists(global_i1, global_i2, global_i3)%construct(periodic_box, &
                positions)
        end do
        end do
        end do
    end subroutine Abstract_Visitable_Cells_construct_visitable_lists

    subroutine Abstract_Visitable_Cells_check_size(this, periodic_box)
        class(Abstract_Visitable_Cells), intent(in) :: this
        class(Abstract_Periodic_Box), intent(in) :: periodic_box

        real(DP) :: box_mod_cell(num_dimensions)

        if (any(this%nums < num_local_cells)) then
            call error_exit("Abstract_Visitable_Cells: this%nums is too small.")
            !go back to 3 then
        end if
        box_mod_cell = modulo(periodic_box%get_size(), this%size)
        if (any(box_mod_cell > real_zero .and. abs(box_mod_cell - this%size) > real_zero)) then
            call error_exit("Abstract_Visitable_Cells:"//&
                            "this%size size is not a divisor of periodic_box%get_size()")
        end if
    end subroutine Abstract_Visitable_Cells_check_size

    subroutine Abstract_Visitable_Cells_destroy(this)
        class(Abstract_Visitable_Cells), intent(inout) :: this

        integer :: global_i1, global_i2, global_i3

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

    subroutine Abstract_Visitable_Cells_move(this, from, to)
        class(Abstract_Visitable_Cells), intent(inout) :: this
        type(Concrete_Particle), intent(in) :: from, to

        integer, dimension(num_dimensions) :: from_i_cell, to_i_cell

        from_i_cell = this%index(from%position)
        to_i_cell = this%index(to%position)
        if (any(from_i_cell /= to_i_cell)) then
            call this%visitable_lists(from_i_cell(1), from_i_cell(2), &
                from_i_cell(3))%deallocate(from%i)
            call this%visitable_lists(to_i_cell(1), to_i_cell(2), &
                to_i_cell(3))%allocate(to%i)
        end if
    end subroutine Abstract_Visitable_Cells_move

    subroutine Abstract_Visitable_Cells_add(this, particle)
        class(Abstract_Visitable_Cells), intent(inout) :: this
        type(Concrete_Particle), intent(in) :: particle

        integer :: i_cell(num_dimensions)

        i_cell = this%index(particle%position)
        call this%visitable_lists(i_cell(1), i_cell(2), i_cell(3))%allocate(particle%i)
    end subroutine Abstract_Visitable_Cells_add

    subroutine Abstract_Visitable_Cells_remove(this, particle)
        class(Abstract_Visitable_Cells), intent(inout) :: this
        type(Concrete_Particle), intent(in) :: particle

        integer :: i_cell(num_dimensions)

        i_cell = this%index(particle%position)
        call this%visitable_lists(i_cell(1), i_cell(2), i_cell(3))%deallocate(particle%i)
        if (particle%i < this%positions%get_num()) then
            call this%visitable_lists(i_cell(1), i_cell(2), &
                i_cell(3))%overwrite(this%positions%get_num(), particle%i)
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

!implementation PBC_3D_Visitable_Cells

    subroutine PBC_3D_Visitable_Cells_set_neighbours(this)
        class(PBC_3D_Visitable_Cells), intent(inout) :: this

        integer :: global_i1, global_i2, global_i3
        integer :: local_i1, local_i2, local_i3
        integer :: i_cell(num_dimensions)

        this%neighbours = 0
        do global_i3 = this%global_lbounds(3), this%global_ubounds(3)
        do global_i2 = this%global_lbounds(2), this%global_ubounds(2)
        do global_i1 = this%global_lbounds(1), this%global_ubounds(1)
            do local_i3 = this%local_lbounds(3), this%local_ubounds(3)
            do local_i2 = this%local_lbounds(2), this%local_ubounds(2)
            do local_i1 = this%local_lbounds(1), this%local_ubounds(1)
                i_cell = [global_i1, global_i2, global_i3] + [local_i1, local_i2, local_i3]
                i_cell = pbc_3d_index(i_cell, this%nums)
                this%neighbours(:, local_i1, local_i2, local_i3, &
                    global_i1, global_i2, global_i3) = i_cell
            end do
            end do
            end do
        end do
        end do
        end do
    end subroutine PBC_3D_Visitable_Cells_set_neighbours

    subroutine PBC_3D_Visitable_Cells_visit(this, particle, pair_potential, overlap, energy)
        class(PBC_3D_Visitable_Cells), intent(in) :: this
        type(Concrete_Particle), intent(in) :: particle
        class(Abstract_Pair_Potential), intent(in) :: pair_potential
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy

        real(DP) :: energy_i
        integer, dimension(num_dimensions) :: i_cell, i_local_cell
        integer :: local_i1, local_i2, local_i3

        i_cell = this%index(particle%position)
        energy = 0._DP
        do local_i3 = this%local_lbounds(3), this%local_ubounds(3)
        do local_i2 = this%local_lbounds(2), this%local_ubounds(2)
        do local_i1 = this%local_lbounds(1), this%local_ubounds(1)
            i_local_cell = this%neighbours(:, local_i1, local_i2, local_i3, &
                i_cell(1), i_cell(2), i_cell(3))
            call this%visitable_lists(i_local_cell(1), i_local_cell(2), &
                i_local_cell(3))%visit(particle, pair_potential, overlap, energy_i)
            if (overlap) return
            energy = energy + energy_i
        end do
        end do
        end do
    end subroutine PBC_3D_Visitable_Cells_visit

!implementation PBC_3D_Visitable_Cells

end module class_visitable_cells
