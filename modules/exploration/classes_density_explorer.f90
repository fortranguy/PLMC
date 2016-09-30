module classes_density_explorer

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use data_constants, only: num_dimensions
use data_strings, only: max_line_length
use procedures_checks, only: check_positive
use procedures_elementary_statistics, only: average, standard_deviation
use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain
use procedures_boxes_factory, only: boxes_destroy => destroy

implicit none

private

    type, abstract, public :: Abstract_Density_Explorer
    private
        class(Abstract_Parallelepiped_Domain), allocatable :: parallelepiped_domains(:)
    contains
        procedure(Abstract_destroy), deferred :: destroy
        procedure(Abstract_fill), deferred :: fill
        procedure(Abstract_write), deferred :: write
    end type Abstract_Density_Explorer

    abstract interface

        subroutine Abstract_destroy(this)
        import :: Abstract_Density_Explorer
            class(Abstract_Density_Explorer), intent(inout) :: this
        end subroutine Abstract_destroy

        subroutine Abstract_fill(this, i_snap, positions)
        import :: DP, Abstract_Density_Explorer
            class(Abstract_Density_Explorer), intent(inout) :: this
            integer, intent(in) :: i_snap
            real(DP), intent(in) :: positions(:, :)
        end subroutine Abstract_fill

        subroutine Abstract_write(this)
            import :: Abstract_Density_Explorer
            class(Abstract_Density_Explorer), intent(in) :: this
        end subroutine Abstract_write

    end interface

    type, extends(Abstract_Density_Explorer), public :: Plain_Density_Explorer
    private
        real(DP), allocatable :: density(:)
    contains
        procedure :: construct => Plain_construct
        procedure :: destroy => Plain_destroy
        procedure :: fill => Plain_fill
        procedure :: write => Plain_write
    end type Plain_Density_Explorer

    type, extends(Abstract_Density_Explorer), public :: Z_Density_Explorer
    private
        real(DP), allocatable :: density(:, :)
        real(DP) :: delta_z
        integer :: i_z_min, i_z_max
        character(len=:), allocatable :: filename
    contains
        procedure :: construct => Z_construct
        procedure :: destroy => Z_destroy
        procedure :: fill => Z_fill
        procedure :: write => Z_write
    end type Z_Density_Explorer

    type, extends(Abstract_Density_Explorer), public :: XZ_Density_Explorer
    private
        real(DP), allocatable :: density(:, :, :)
        real(DP) :: delta_x, delta_z
        integer :: i_x_min, i_x_max, i_z_min, i_z_max
        character(len=:), allocatable :: filename
    contains
        procedure :: construct => XZ_construct
        procedure :: destroy => XZ_destroy
        procedure :: fill => XZ_fill
        procedure :: write => XZ_write
    end type XZ_Density_Explorer

contains

!implementation Plain_Density_Explorer

    subroutine Plain_construct(this, parallelepiped_domains, num_snaps)
        class(Plain_Density_Explorer), intent(out) :: this
        class(Abstract_Parallelepiped_Domain), intent(in) :: parallelepiped_domains(:)
        integer, intent(in) :: num_snaps

        allocate(this%parallelepiped_domains, source=parallelepiped_domains)
        allocate(this%density(num_snaps))
        this%density = 0._DP
    end subroutine

    subroutine Plain_destroy(this)
        class(Plain_Density_Explorer), intent(inout) :: this

        if (allocated(this%density)) deallocate(this%density)
        call boxes_destroy(this%parallelepiped_domains)
    end subroutine Plain_destroy

    subroutine Plain_fill(this, i_snap, positions)
        class(Plain_Density_Explorer), intent(inout) :: this
        integer, intent(in) :: i_snap
        real(DP), intent(in) :: positions(:, :)

        integer :: i_particle

        do i_particle = 1, size(positions, 2)
            if (this%parallelepiped_domains(1)%is_inside(positions(:, i_particle))) then
                this%density(i_snap) = this%density(i_snap) + 1._DP
            end if
        end do
        this%density(i_snap) = this%density(i_snap) / &
            product(this%parallelepiped_domains(1)%get_size())
    end subroutine Plain_fill

    !> \[ \rho = \left\langle \frac{N}{V} \right\rangle \]
    subroutine Plain_write(this)
        class(Plain_Density_Explorer), intent(in) :: this

        real(DP) :: density_average, density_std_dev

        density_average = average(this%density)
        write(output_unit, *) "average", density_average
        density_std_dev = standard_deviation(this%density)
        write(output_unit, *) "standard_deviation", density_std_dev
    end subroutine Plain_write

!end implementation Plain_Density_Explorer

!implementation Z_Density_Explorer

    subroutine Z_construct(this, parallelepiped_domains, delta_z, num_snaps, filename)
        class(Z_Density_Explorer), intent(out) :: this
        class(Abstract_Parallelepiped_Domain), intent(in) :: parallelepiped_domains(:)
        real(DP), intent(in) :: delta_z
        integer, intent(in) :: num_snaps
        character(len=*), intent(in) :: filename

        real(DP), dimension(num_dimensions) :: domain_origin, domain_size

        allocate(this%parallelepiped_domains, source=parallelepiped_domains)
        call check_positive("Z_Density_Explorer: construct", "delta_z", delta_z)
        this%delta_z = delta_z
        domain_origin = this%parallelepiped_domains(1)%get_origin()
        domain_size = this%parallelepiped_domains(1)%get_size()
        this%i_z_min = nint((domain_origin(3) - domain_size(3)/2) / this%delta_z)
        this%i_z_max = nint((domain_origin(3) + domain_size(3)/2) / this%delta_z)
        allocate(this%density(this%i_z_min:this%i_z_max, num_snaps))
            !! @note swap dimensions? cf. [[Z_write]]
        this%density = 0._DP
        this%filename = filename
    end subroutine Z_construct

    subroutine Z_destroy(this)
        class(Z_Density_Explorer), intent(inout) :: this

        if (allocated(this%filename)) deallocate(this%filename)
        if (allocated(this%density)) deallocate(this%density)
        call boxes_destroy(this%parallelepiped_domains)
    end subroutine Z_destroy

    subroutine Z_fill(this, i_snap, positions)
        class(Z_Density_Explorer), intent(inout) :: this
        integer, intent(in) :: i_snap
        real(DP), intent(in) :: positions(:, :)

        real(DP) :: domain_size(num_dimensions)
        integer :: i_particle, i_z

        do i_particle = 1, size(positions, 2)
            if (this%parallelepiped_domains(1)%is_inside(positions(:, i_particle))) then
                i_z = nint(positions(3, i_particle) / this%delta_z)
                this%density(i_z, i_snap) = this%density(i_z, i_snap) + 1._DP
            end if
        end do
        domain_size = this%parallelepiped_domains(1)%get_size()
        this%density(:, i_snap) = this%density(:, i_snap) / product(domain_size(1:2)) / this%delta_z
    end subroutine Z_fill

    !> \[
    !>      \rho(z) = \left\langle
    !>          \frac{1}{L_x L_y}\frac{\mathrm{d} N}{\mathrm{d} z}(z)
    !>      \right\rangle
    !> \]
    subroutine Z_write(this)
        class(Z_Density_Explorer), intent(in) :: this

        real(DP), dimension(this%i_z_min:this%i_z_max) :: density_average, density_std_dev
        integer :: i_z
        integer :: density_unit

        open(newunit=density_unit, recl=max_line_length, file=this%filename, action="write")
        write(density_unit, *) "# z    average    standard_deviation"
        do i_z = lbound(density_average, 1), ubound(density_average, 1)
            density_average(i_z) = average(this%density(i_z, :))
            density_std_dev(i_z) = standard_deviation(this%density(i_z, :))
            write(density_unit, *) real(i_z, DP) * this%delta_z, density_average(i_z), &
                density_std_dev(i_z)
        end do
        close(density_unit)
    end subroutine Z_write

!end implementation Z_Density_Explorer

!implementation XZ_Density_Explorer

    subroutine XZ_construct(this, parallelepiped_domains, delta_xz, num_snaps, filename)
        class(XZ_Density_Explorer), intent(out) :: this
        class(Abstract_Parallelepiped_Domain), intent(in) :: parallelepiped_domains(:)
        real(DP), intent(in) :: delta_xz
        integer, intent(in) :: num_snaps
        character(len=*), intent(in) :: filename

        real(DP), dimension(num_dimensions) :: domain_origin, domain_size

        allocate(this%parallelepiped_domains, source=parallelepiped_domains)
        call check_positive("XZ_Density_Explorer: construct", "delta_xz", delta_xz)
        this%delta_x = delta_xz
        this%delta_z = delta_xz

        domain_origin = this%parallelepiped_domains(1)%get_origin()
        domain_size = this%parallelepiped_domains(1)%get_size()
        this%i_x_min = nint((domain_origin(1) - domain_size(1)/2) / this%delta_x)
        this%i_x_max = nint((domain_origin(1) + domain_size(1)/2) / this%delta_x)
        this%i_z_min = nint((domain_origin(3) - domain_size(3)/2) / this%delta_z)
        this%i_z_max = nint((domain_origin(3) + domain_size(3)/2) / this%delta_z)
        allocate(this%density(this%i_x_min:this%i_x_max, this%i_z_min:this%i_z_max, num_snaps))
        this%density = 0
        this%filename = filename
    end subroutine XZ_construct

    subroutine XZ_destroy(this)
        class(XZ_Density_Explorer), intent(inout) :: this

        if (allocated(this%filename)) deallocate(this%filename)
        if (allocated(this%density)) deallocate(this%density)
        call boxes_destroy(this%parallelepiped_domains)
    end subroutine XZ_destroy

    subroutine XZ_fill(this, i_snap, positions)
        class(XZ_Density_Explorer), intent(inout) :: this
        integer, intent(in) :: i_snap
        real(DP), intent(in) :: positions(:, :)

        real(DP) :: domain_size(num_dimensions)
        integer :: i_particle, i_x, i_z

        do i_particle = 1, size(positions, 2)
            if (this%parallelepiped_domains(1)%is_inside(positions(:, i_particle))) then
                i_x = nint(positions(1, i_particle) / this%delta_x)
                i_z = nint(positions(3, i_particle) / this%delta_z)
                this%density(i_x, i_z, i_snap) = this%density(i_x, i_z, i_snap) + 1._DP
            end if
        end do
        domain_size = this%parallelepiped_domains(1)%get_size()
        this%density(:, :, i_snap) = this%density(:, :, i_snap) / domain_size(2) / &
            (this%delta_x * this%delta_z)
    end subroutine XZ_fill

    subroutine XZ_write(this)
        class(XZ_Density_Explorer), intent(in) :: this

        real(DP), dimension(this%i_x_min:this%i_x_max, this%i_z_min:this%i_z_max) :: &
            density_average, density_std_dev
        integer :: i_x, i_z
        integer :: density_unit

        open(newunit=density_unit, recl=max_line_length, file=this%filename, action="write")
        write(density_unit, *) "# x    z    average    standard_deviation"
        do i_z = lbound(density_average, 2), ubound(density_average, 2)
            do i_x = lbound(density_average, 1), ubound(density_average, 1)
                density_average(i_x, i_z) = average(this%density(i_x, i_z, :))
                density_std_dev(i_x, i_z) = standard_deviation(this%density(i_x, i_z, :))
                write(density_unit, *) real(i_x, DP) * this%delta_x, real(i_z, DP) * this%delta_z, &
                    density_average(i_x, i_z), density_std_dev(i_x, i_z)
            end do
            write(density_unit, *)
        end do
        close(density_unit)
    end subroutine XZ_write

!end implementation XZ_Density_Explorer

end module classes_density_explorer
