module procedures_des_surf_mixture_factory

use procedures_errors, only: error_exit
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_permittivity, only: Abstract_Permittivity
use procedures_environment_inquirers, only: periodicity_is_xyz, periodicity_is_xy
use classes_mixture_total_moment, only: Abstract_Mixture_Total_Moment, &
    Concrete_Mixture_Total_Moment, Null_Mixture_Total_Moment
use classes_des_surf_mixture, only: Abstract_DES_Surf_Mixture, Spherical_DES_Surf_Mixture, &
    Rectangular_DES_Surf_Mixture, Null_DES_Surf_Mixture

implicit none

private
public :: create, destroy

contains

    subroutine create(mixture, periodic_boxes, permittivity, total_moments)
        class(Abstract_DES_Surf_Mixture), allocatable, intent(out) :: mixture(:)
        class(Abstract_Periodic_Box), intent(in) :: periodic_boxes(:)
        class(Abstract_Permittivity), intent(in) :: permittivity
        class(Abstract_Mixture_Total_Moment), intent(in) :: total_moments(:)

        integer :: i_box

        select type(total_moments)
            type is (Concrete_Mixture_Total_Moment)
                if (all(periodicity_is_xyz(periodic_boxes))) then
                    allocate(Spherical_DES_Surf_Mixture :: mixture(size(total_moments)))
                else if (all(periodicity_is_xy(periodic_boxes))) then
                    allocate(Rectangular_DES_Surf_Mixture :: mixture(size(total_moments)))
                else
                    call error_exit("procedures_des_surf_mixture_factory: create: "//&
                            "box periodicity is unknown.")
                end if
            type is (Null_Mixture_Total_Moment)
                allocate(Null_DES_Surf_Mixture :: mixture(size(total_moments)))
            class default
                call error_exit("procedures_des_surf_mixture_factory: create: //"&
                    "total_moments type unknown.")
        end select

        do i_box = 1, size(mixture)
            call mixture(i_box)%construct(periodic_boxes(i_box), permittivity, total_moments(i_box))
        end do
    end subroutine create

    subroutine destroy(mixture)
        class(Abstract_DES_Surf_Mixture), allocatable, intent(inout) :: mixture(:)

        integer :: i_box

        if (allocated(mixture)) then
            do i_box = size(mixture), 1, -1
                call mixture(i_box)%destroy()
            end do
            deallocate(mixture)
        end if
    end subroutine destroy

end module procedures_des_surf_mixture_factory
