module procedures_ewald_surf_factory

use procedures_errors, only: error_exit
use class_periodic_box, only: Abstract_Periodic_Box, XYZ_Periodic_Box, XY_Periodic_Box
use class_permittivity, only: Abstract_Permittivity
use class_mixture_total_moment, only: Abstract_Mixture_Total_Moment, &
    Concrete_Mixture_Total_Moment, Null_Mixture_Total_Moment
use types_component_wrapper, only: Component_Wrapper
use class_ewald_surf_mixture, only: Abstract_Ewald_Surf_Mixture, Spheric_Ewald_Surf_Mixture, &
    Rectangular_Ewald_Surf_Mixture, Null_Ewald_Surf_Mixture

implicit none

private
public :: ewald_surf_create, ewald_surf_destroy

interface ewald_surf_create
    module procedure :: create_surf_mixture
end interface ewald_surf_create

interface ewald_surf_destroy
    module procedure :: destroy_surf_mixture
end interface ewald_surf_destroy

contains

    subroutine create_surf_mixture(surf_mixture, periodic_box, permittivity, total_moment)
        class(Abstract_Ewald_Surf_Mixture), allocatable, intent(out) :: surf_mixture
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Permittivity), intent(in) :: permittivity
        class(Abstract_Mixture_Total_Moment), intent(in) :: total_moment

        select type(total_moment)
            type is (Concrete_Mixture_Total_Moment)
                select type(periodic_box)
                    type is (XYZ_Periodic_Box)
                        allocate(Spheric_Ewald_Surf_Mixture :: surf_mixture)
                    type is (XY_Periodic_Box)
                        allocate(Rectangular_Ewald_Surf_Mixture :: surf_mixture)
                    class default
                        call error_exit("create_surf_mixture: periodic_box type unknown.")
                end select
            type is (Null_Mixture_Total_Moment)
                allocate(Null_Ewald_Surf_Mixture :: surf_mixture)
            class default
                call error_exit("create_surf_mixture: total_moment type unknown.")
        end select
        call surf_mixture%construct(periodic_box, permittivity, total_moment)
    end subroutine create_surf_mixture

    subroutine destroy_surf_mixture(surf_mixture)
        class(Abstract_Ewald_Surf_Mixture), allocatable, intent(inout) :: surf_mixture

        if (allocated(surf_mixture)) then
            call surf_mixture%destroy()
            deallocate(surf_mixture)
        end if
    end subroutine destroy_surf_mixture

end module procedures_ewald_surf_factory
