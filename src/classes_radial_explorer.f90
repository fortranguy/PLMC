module classes_radial_explorer

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use data_strings, only: max_line_length
use classes_number_to_string, only: Concrete_Number_to_String
use procedures_errors, only: error_exit
use procedures_checks, only: check_positive, check_in_range, check_string_not_empty
use procedures_elementary_geometry, only: sphere_surface
use classes_periodic_box, only: Abstract_Periodic_Box
use types_raw_coordinates, only: Concrete_Raw_Coordinates
use types_component_coordinates_reader_selector, only: Component_Coordinates_Reader_Selector
use procedures_complete_coordinates_reader, only: complete_coordinates_read => read

implicit none

private

    !Calculate the radial distribution function
    !> \[
    !>      g(r) = \frac{\mathrm{d}N}{\mathrm{d}r} \frac{1}{\rho S(r)}
    !> \]
    !> with \( S(r) \) the area of a sphere.
    type, abstract, public :: Abstract_Radial_Explorer
    private
        integer :: num_components = 0
        real(DP) :: density_sum = 0._DP
        type(Component_Coordinates_Reader_Selector) :: selector
        real(DP) :: max_distance, delta_distance
        real(DP), allocatable :: bins_function(:)
        character(len=:), allocatable :: filename
    contains
        procedure :: destroy => Abstract_destroy
        procedure(Abstract_read_and_fill), deferred :: read_and_fill
        procedure :: write => Abstract_write
        procedure, private :: set_core => Abstract_set_core
    end type Abstract_Radial_Explorer

    abstract interface

        subroutine Abstract_read_and_fill(this, periodic_box, snap_filename)
        import :: Abstract_Periodic_Box, Abstract_Radial_Explorer
            class(Abstract_Radial_Explorer), intent(inout) :: this
            class(Abstract_Periodic_Box), intent(inout) :: periodic_box
        character(len=*), intent(in) :: snap_filename
        end subroutine Abstract_read_and_fill

    end interface

    type, extends(Abstract_Radial_Explorer), public :: Intra_Radial_Explorer
    private
        integer :: i_component = 0
    contains
        procedure :: construct => Intra_construct
        procedure :: read_and_fill => Intra_read_and_fill
    end type Intra_Radial_Explorer

    type, extends(Abstract_Radial_Explorer), public :: Inter_Radial_Explorer
    private
        integer :: ij_components(2) = 0
    contains
        procedure :: construct => Inter_construct
        procedure :: read_and_fill => Inter_read_and_fill
    end type Inter_Radial_Explorer

contains

!implementation Abstract_Radial_Explorer

    subroutine Abstract_set_core(this, num_components, selector, max_distance, delta_distance, &
        filename)
        class(Abstract_Radial_Explorer), intent(out) :: this
        integer, intent(in) :: num_components
        type(Component_Coordinates_Reader_Selector), intent(in) :: selector
        real(DP), intent(in) :: max_distance, delta_distance
        character(len=*), intent(in) :: filename

        call check_positive("Abstract_Radial_Explorer: construct", "num_components", num_components)
        this%num_components = num_components
        this%selector%read_positions = selector%read_positions
        this%selector%read_orientations = selector%read_orientations
        call check_positive("Abstract_Radial_Explorer: construct", "max_distance", max_distance)
        this%max_distance = max_distance
        call check_positive("Abstract_Radial_Explorer: construct", "delta_distance", delta_distance)
        this%delta_distance = delta_distance
        call check_string_not_empty("Abstract_Radial_Explorer: construct", filename)
        this%filename = filename
        allocate(this%bins_function(nint(this%max_distance / this%delta_distance)))
        this%bins_function = 0._DP
    end subroutine Abstract_set_core

    subroutine Abstract_destroy(this)
        class(Abstract_Radial_Explorer), intent(inout) :: this

        if (allocated(this%bins_function)) deallocate(this%bins_function)
        if (allocated(this%filename)) deallocate(this%filename)
    end subroutine Abstract_destroy

    subroutine Abstract_write(this, num_snaps)
        class(Abstract_Radial_Explorer), intent(in) :: this
        integer, intent(in) :: num_snaps

        integer :: bins_unit
        real(DP) :: bins_function(size(this%bins_function))
        integer :: i_bin
        real(DP) :: distance_i
        real(DP) :: density

        bins_function = this%bins_function / real(num_snaps, DP)
        density = this%density_sum / real(num_snaps)

        open(newunit=bins_unit, recl=max_line_length, file=this%filename, action="write")
        do i_bin = 1, size(bins_function)
            distance_i = real(i_bin, DP) * this%delta_distance
            bins_function(i_bin) = bins_function(i_bin) / this%delta_distance / density / &
                sphere_surface(distance_i)
            write(bins_unit, *) distance_i, bins_function(i_bin)
        end do
        close(bins_unit)
    end subroutine Abstract_write

!end implementation Abstract_Radial_Explorer

!implementation Intra_Radial_Explorer

    subroutine Intra_construct(this, num_components, i_component, selector, max_distance, &
        delta_distance, filename)
        class(Intra_Radial_Explorer), intent(out) :: this
        integer, intent(in) :: num_components, i_component
        type(Component_Coordinates_Reader_Selector), intent(in) :: selector
        real(DP), intent(in) :: max_distance, delta_distance
        character(len=*), intent(in) :: filename

        call this%set_core(num_components, selector, max_distance, delta_distance, filename)
        call check_in_range("Intra_Radial_Explorer: construct", this%num_components, "i_component",&
            i_component)
        this%i_component = i_component
    end subroutine Intra_construct

    subroutine Intra_read_and_fill(this, periodic_box, snap_filename)
        class(Intra_Radial_Explorer), intent(inout) :: this
        class(Abstract_Periodic_Box), intent(inout) :: periodic_box
        character(len=*), intent(in) :: snap_filename

        real(DP) :: bins_snap(size(this%bins_function))
        integer :: i_bin
        real(DP) :: box_size(num_dimensions)
        type(Concrete_Raw_Coordinates) :: raw_coordinates
        integer :: num_particles, i_particle, j_particle
        real(DP) :: distance_ij

        call complete_coordinates_read(box_size, raw_coordinates, this%num_components, &
            this%i_component, snap_filename, this%selector)
        call periodic_box%set(box_size)
        num_particles = size(raw_coordinates%positions, 2)

        bins_snap = 0._DP
        do j_particle = 1, num_particles
            do i_particle = 1, j_particle-1
                distance_ij = periodic_box%distance(raw_coordinates%positions(:, i_particle), &
                    raw_coordinates%positions(:, j_particle))
                if (distance_ij < this%max_distance) then
                    i_bin = nint(distance_ij / this%delta_distance)
                    bins_snap(i_bin) = bins_snap(i_bin) + 1._DP
                end if
            end do
        end do

        if (num_particles > 0) then
            this%bins_function = this%bins_function + 2._DP * bins_snap / real(num_particles, DP)
            this%density_sum = this%density_sum + real(num_particles, DP) / &
                product(periodic_box%get_size())
        end if
    end subroutine Intra_read_and_fill

!end implementation Intra_Radial_Explorer

!implementation Inter_Radial_Explorer

    subroutine Inter_construct(this, num_components, ij_components, selector, max_distance, &
        delta_distance, filename)
        class(Inter_Radial_Explorer), intent(out) :: this
        integer, intent(in) :: num_components, ij_components(:)
        type(Component_Coordinates_Reader_Selector), intent(in) :: selector
        real(DP), intent(in) :: max_distance, delta_distance
        character(len=*), intent(in) :: filename

        integer :: i_partner
        type(Concrete_Number_to_String) :: string

        call this%set_core(num_components, selector, max_distance, delta_distance, filename)
        do i_partner = 1, size(this%ij_components)
            call check_in_range("Inter_Radial_Explorer: construct", this%num_components, &
                "ij_components("//string%get(i_partner)//")", ij_components(i_partner))
        end do
        if (ij_components(1) == ij_components(2)) then
            call error_exit("Inter_Radial_Explorer: construct: i_component = j_component.")
        end if
        this%ij_components = ij_components
    end subroutine Inter_construct

    subroutine Inter_read_and_fill(this, periodic_box, snap_filename)
        class(Inter_Radial_Explorer), intent(inout) :: this
        class(Abstract_Periodic_Box), intent(inout) :: periodic_box
        character(len=*), intent(in) :: snap_filename

        real(DP) :: bins_snap(size(this%bins_function))
        integer :: i_bin
        real(DP) :: box_size(num_dimensions)
        type(Concrete_Raw_Coordinates) :: raw_coordinates(2)
        integer :: nums_particles(2), i_particle, j_particle, i_parter
        real(DP) :: distance_ij

        do i_parter = 1, size(raw_coordinates)
            call complete_coordinates_read(box_size, raw_coordinates(i_parter), this%&
                num_components, this%ij_components(i_parter), snap_filename, this%selector)
            nums_particles(i_parter) = size(raw_coordinates(i_parter)%positions, 2)
        end do
        call periodic_box%set(box_size)

        bins_snap = 0._DP
        do j_particle = 1, nums_particles(2)
            do i_particle = 1, nums_particles(1)
                distance_ij = periodic_box%distance(raw_coordinates(1)%positions(:, i_particle), &
                    raw_coordinates(2)%positions(:, j_particle))
                if (distance_ij < this%max_distance) then
                    i_bin = nint(distance_ij / this%delta_distance)
                    bins_snap(i_bin) = bins_snap(i_bin) + 1._DP
                end if
            end do
        end do

        if (all(nums_particles > 0)) then
            this%bins_function = this%bins_function + bins_snap / real(nums_particles(2), DP)
            this%density_sum = this%density_sum + real(nums_particles(1), DP) / &
                product(periodic_box%get_size())
        end if
    end subroutine Inter_read_and_fill

!end implementation Inter_Radial_Explorer

end module classes_radial_explorer
