module classes_density_explorator

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain
use procedures_box_factory, only: box_destroy => destroy

implicit none

private

    type, abstract, public :: Abstract_Density_Explorator
    private
        class(Abstract_Parallelepiped_Domain), allocatable :: parallelepiped_domain
        !character(len=:), allocatable :: snap_filename
    contains
        procedure(Abstract_destroy), deferred :: destroy
        procedure(Abstract_fill), deferred :: fill
        procedure(Abstract_write), deferred :: write
    end type Abstract_Density_Explorator

    abstract interface

        subroutine Abstract_destroy(this)
        import :: Abstract_Density_Explorator
            class(Abstract_Density_Explorator), intent(inout) :: this
        end subroutine Abstract_destroy

        subroutine Abstract_fill(this, i_snap, positions)
        import :: DP, Abstract_Density_Explorator
            class(Abstract_Density_Explorator), intent(inout) :: this
            integer, intent(in) :: i_snap
            real(DP), intent(in) :: positions(:, :)
        end subroutine Abstract_fill

        subroutine Abstract_write(this)
            import :: Abstract_Density_Explorator
            class(Abstract_Density_Explorator), intent(in) :: this
        end subroutine Abstract_write

    end interface

    type, extends(Abstract_Density_Explorator), public :: Plain_Density_Explorator
    private
        real(DP), allocatable :: nums_inside(:)
    contains
        procedure :: construct => Plain_construct
        procedure :: destroy => Plain_destroy
        procedure :: fill => Plain_fill
        procedure :: write => Plain_write
    end type Plain_Density_Explorator

contains

!implementation Plain_Density_Explorator

    subroutine Plain_construct(this, parallelepiped_domain, num_snaps)
        class(Plain_Density_Explorator), intent(out) :: this
        class(Abstract_Parallelepiped_Domain), intent(in) :: parallelepiped_domain
        integer, intent(in) :: num_snaps

        allocate(this%parallelepiped_domain, source=parallelepiped_domain)
        allocate(this%nums_inside(num_snaps))
        this%nums_inside = 0._DP
    end subroutine

    subroutine Plain_destroy(this)
        class(Plain_Density_Explorator), intent(inout) :: this

        if (allocated(this%nums_inside)) deallocate(this%nums_inside)
        call box_destroy(this%parallelepiped_domain)
    end subroutine Plain_destroy

    subroutine Plain_fill(this, i_snap, positions)
        class(Plain_Density_Explorator), intent(inout) :: this
        integer, intent(in) :: i_snap
        real(DP), intent(in) :: positions(:, :)

        integer :: i_particle

        do i_particle = 1, size(positions, 2)
            if (this%parallelepiped_domain%is_inside(positions(:, i_particle))) then
                this%nums_inside(i_snap) = this%nums_inside(i_snap) + 1._DP
            end if
        end do
    end subroutine Plain_fill

    subroutine Plain_write(this)
        class(Plain_Density_Explorator), intent(in) :: this

        real(DP) :: avg_num, rms_num

        avg_num = sum(this%nums_inside) / size(this%nums_inside)
        write(output_unit, *) "domain.density.avg", avg_num / &
            product(this%parallelepiped_domain%get_size())
        rms_num = sqrt(sum((this%nums_inside - avg_num)**2) / (size(this%nums_inside) - 1))
        write(output_unit, *) "domain.density.rms", rms_num / &
            product(this%parallelepiped_domain%get_size())
    end subroutine Plain_write

!end implementation Plain_Density_Explorator

end module classes_density_explorator
