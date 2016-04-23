module procedures_real_writer_factory

use classes_real_writer, only: Abstract_Real_Writer, Concrete_Real_Writer, Null_Real_Writer

implicit none

private
public :: create, destroy

interface create
    module procedure :: create_dipolar_mixture_energy
end interface create

interface destroy
    module procedure :: destroy_dipolar_mixture_energy
end interface destroy

contains

    subroutine create_dipolar_mixture_energy(dipolar_mixture_energy, dipoles_exist, filename)
        class(Abstract_Real_Writer), allocatable, intent(out) :: dipolar_mixture_energy
        logical, intent(in) :: dipoles_exist
        character(len=*), intent(in) :: filename

        if (dipoles_exist) then
            allocate(Concrete_Real_Writer :: dipolar_mixture_energy)
        else
            allocate(Null_Real_Writer :: dipolar_mixture_energy)
        end if
        call dipolar_mixture_energy%construct(filename)
    end subroutine create_dipolar_mixture_energy

    subroutine destroy_dipolar_mixture_energy(dipolar_mixture_energy)
        class(Abstract_Real_Writer), allocatable, intent(inout) :: dipolar_mixture_energy

        if (allocated(dipolar_mixture_energy)) then
            call dipolar_mixture_energy%destroy()
            deallocate(dipolar_mixture_energy)
        end if
    end subroutine destroy_dipolar_mixture_energy

end module procedures_real_writer_factory
