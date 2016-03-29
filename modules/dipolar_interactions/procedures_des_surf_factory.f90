module procedures_des_surf_factory

use procedures_errors, only: error_exit
use class_periodic_box, only: Abstract_Periodic_Box, XYZ_Periodic_Box, XY_Periodic_Box
use class_permittivity, only: Abstract_Permittivity
use class_mixture_total_moment, only: Abstract_Mixture_Total_Moment, &
    Concrete_Mixture_Total_Moment, Null_Mixture_Total_Moment
use types_component_wrapper, only: Component_Wrapper
use class_des_surf_mixture, only: Abstract_DES_Surf_Mixture, Spheric_DES_Surf_Mixture, &
    Rectangular_DES_Surf_Mixture, Null_DES_Surf_Mixture

implicit none

private
public :: des_surf_create, des_surf_destroy

interface des_surf_create
    module procedure :: create_mixture
end interface des_surf_create

interface des_surf_destroy
    module procedure :: destroy_mixture
end interface des_surf_destroy

contains

    subroutine create_mixture(mixture, periodic_box, permittivity, total_moment)
        class(Abstract_DES_Surf_Mixture), allocatable, intent(out) :: mixture
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Permittivity), intent(in) :: permittivity
        class(Abstract_Mixture_Total_Moment), intent(in) :: total_moment

        select type(total_moment)
            type is (Concrete_Mixture_Total_Moment)
                select type(periodic_box)
                    type is (XYZ_Periodic_Box)
                        allocate(Spheric_DES_Surf_Mixture :: mixture)
                    type is (XY_Periodic_Box)
                        allocate(Rectangular_DES_Surf_Mixture :: mixture)
                    class default
                        call error_exit("create_mixture: periodic_box type unknown.")
                end select
            type is (Null_Mixture_Total_Moment)
                allocate(Null_DES_Surf_Mixture :: mixture)
            class default
                call error_exit("create_mixture: total_moment type unknown.")
        end select
        call mixture%construct(periodic_box, permittivity, total_moment)
    end subroutine create_mixture

    subroutine destroy_mixture(mixture)
        class(Abstract_DES_Surf_Mixture), allocatable, intent(inout) :: mixture

        if (allocated(mixture)) then
            call mixture%destroy()
            deallocate(mixture)
        end if
    end subroutine destroy_mixture

end module procedures_des_surf_factory
