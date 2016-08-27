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

    subroutine create(mixture, periodic_box, permittivity, total_moment)
        class(Abstract_DES_Surf_Mixture), allocatable, intent(out) :: mixture
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Permittivity), intent(in) :: permittivity
        class(Abstract_Mixture_Total_Moment), intent(in) :: total_moment

        select type(total_moment)
            type is (Concrete_Mixture_Total_Moment)
                if (periodicity_is_xyz(periodic_box)) then
                    allocate(Spherical_DES_Surf_Mixture :: mixture)
                else if (periodicity_is_xy(periodic_box)) then
                    allocate(Rectangular_DES_Surf_Mixture :: mixture)
                else
                    call error_exit("procedures_des_surf_mixture_factory: create: "//&
                            "box periodicity is unknown.")
                end if
            type is (Null_Mixture_Total_Moment)
                allocate(Null_DES_Surf_Mixture :: mixture)
            class default
                call error_exit("procedures_des_surf_mixture_factory: create: //"&
                    "total_moment type unknown.")
        end select
        call mixture%construct(periodic_box, permittivity, total_moment)
    end subroutine create

    subroutine destroy(mixture)
        class(Abstract_DES_Surf_Mixture), allocatable, intent(inout) :: mixture

        if (allocated(mixture)) then
            call mixture%destroy()
            deallocate(mixture)
        end if
    end subroutine destroy

end module procedures_des_surf_mixture_factory
