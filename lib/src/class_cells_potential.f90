module class_cells_potential

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions
use module_particles, only: Concrete_Particle
use class_periodic_box, only: Abstract_Periodic_Box
use class_positions, only: Abstract_Positions
use module_particles, only: Concrete_Particle
use module_cells, only: Abstract_Cells

implicit none

private

    type, abstract, public :: Abstract_Cells_Potential
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box
        class(Abstract_Positions), pointer :: positions
        class(Abstract_Cells), pointer :: cells
        integer, allocatable :: neighbours(:, :, :, :, :, :)
    contains
        procedure :: construct => Abstract_Cells_Potential_construct
        procedure, private :: fill => Abstract_Cells_Potential_fill
        procedure :: destroy => Abstract_Cells_Potential_destroy
        procedure :: visit => Abstract_Cells_Potential_visit
        procedure :: move => Abstract_Cells_Potential_move
        procedure :: add => Abstract_Cells_Potential_add
        procedure :: remove => Abstract_Cells_Potential_remove
    end type Abstract_Cells_Potential

    type, extends(Abstract_Cells_Potential), public :: Concrete_Neighbour_Cells

    end type Concrete_Neighbour_Cells

contains

    subroutine Abstract_Cells_Potential_construct(this, periodic_box, positions, cells)
        class(Abstract_Cells_Potential), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Positions), target, intent(in) :: positions
        class(Abstract_Cells), target, intent(in) :: cells

        this%periodic_box => periodic_box
        this%positions => positions
        this%cells => cells
    end subroutine Abstract_Cells_Potential_construct

    subroutine Abstract_Cells_Potential_fill(this)
        class(Abstract_Cells_Potential), intent(out) :: this

        type(Concrete_Particle) :: particle
        integer :: i_particle

        do i_particle = 1, this%positions%get_num()
            particle%i = i_particle
            particle%position = this%positions%get(particle%i)
            call this%add(particle)
        end do
    end subroutine Abstract_Cells_Potential_fill

    subroutine Abstract_Cells_Potential_destroy(this)
        class(Abstract_Cells_Potential), intent(inout) :: this

        this%cells => null()
        this%positions => null()
        this%periodic_box => null()
    end subroutine Abstract_Cells_Potential_destroy

    subroutine Abstract_Cells_Potential_visit(this, same_type, particle, overlap, energy)
        class(Abstract_Cells_Potential), intent(in) :: this
        logical, intent(in) :: same_type
        type(Concrete_Particle), intent(in) :: particle
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy

        real(DP) :: energy_i
        integer, dimension(num_dimensions) :: i_cell, i_local_cell
        integer :: local_i1, local_i2, local_i3

        i_cell = nint(particle%position/this%cells%size)
        energy = 0._DP
        do local_i3 = this%cells%local_lbounds(3), this%cells%global_ubounds(3)
        do local_i2 = this%cells%local_lbounds(2), this%cells%global_ubounds(2)
        do local_i1 = this%cells%local_lbounds(1), this%cells%global_ubounds(1)
            i_local_cell = this%cells%neighbours(:, local_i1, local_i2, local_i3, &
                                                 i_cell(1), i_cell(2), i_cell(3))
            call this%cells%visitables_lists(i_local_cell(1), i_local_cell(2), &
                i_local_cell(3))%visit(same_type, particle, overlap, energy_i)
            if (overlap) return
            energy = energy + energy_i
        end do
        end do
        end do
        !call this%linked_list%visit(particle, energy)
    end subroutine Abstract_Cells_Potential_visit

    subroutine Abstract_Cells_Potential_move(this, from, to)
        class(Abstract_Cells_Potential), intent(inout) :: this
        type(Concrete_Particle), intent(in) :: from, to

        integer, dimension(num_dimensions) :: from_i_cell, to_i_cell

        from_i_cell = nint(from%position/this%cells%size)
        to_i_cell = nint(to%position/this%cells%size)
        if (any(from_i_cell /= to_i_cell)) then
            call this%cells%visitables_lists(from_i_cell(1), from_i_cell(2), &
                                             from_i_cell(3))%deallocate(from%i)
            call this%cells%visitables_lists(to_i_cell(1), to_i_cell(2), &
                                             to_i_cell(3))%allocate(to%i)
        end if
    end subroutine Abstract_Cells_Potential_move

    subroutine Abstract_Cells_Potential_add(this, particle)
        class(Abstract_Cells_Potential), intent(inout) :: this
        type(Concrete_Particle), intent(in) :: particle

        integer :: i_cell(num_dimensions)

        i_cell = nint(particle%position/this%cells%size)
        call this%cells%visitables_lists(i_cell(1), i_cell(2), i_cell(3))%allocate(particle%i)
    end subroutine Abstract_Cells_Potential_add

    subroutine Abstract_Cells_Potential_remove(this, particle)
        class(Abstract_Cells_Potential), intent(inout) :: this
        type(Concrete_Particle), intent(in) :: particle

        integer :: i_cell(num_dimensions)

        i_cell = nint(particle%position/this%cells%size)
        call this%cells%visitables_lists(i_cell(1), i_cell(2), i_cell(3))%deallocate(particle%i)
        if (particle%i < this%positions%get_num()) then
            call this%cells%visitables_lists(i_cell(1), i_cell(2), i_cell(3))%overwrite(&
                this%positions%get_num(), particle%i)
        end if
    end subroutine Abstract_Cells_Potential_remove

end module class_cells_potential
