module procedures_exploring_writers_factory

use classes_widom_method, only: Abstract_Widom_Method
use procedures_real_writer_factory, only: real_writer_create => create, &
    real_writer_destroy => destroy
use procedures_line_writer_factory, only: line_writer_create => create, &
    line_writer_destroy => destroy
use types_exploring_writers_wrapper, only: Exploring_Writers_Wrapper
use procedures_property_inquirers, only: measure_chemical_potentials

implicit none

private
public :: create, destroy

contains

    subroutine create(writers, widom_method, num_components)
        type(Exploring_Writers_Wrapper), intent(out) :: writers
        class(Abstract_Widom_Method), intent(in) :: widom_method
        integer, intent(in) :: num_components

        logical :: selector(num_components)

        selector = measure_chemical_potentials(widom_method)
        call line_writer_create(writers%widom_successes, selector, "widom_successes.out")
        call line_writer_create(writers%inv_pow_activities, selector, "inv_pow_activities.out")
    end subroutine create

    subroutine destroy(writers)
        type(Exploring_Writers_Wrapper), intent(inout) :: writers

        call line_writer_destroy(writers%inv_pow_activities)
        call line_writer_destroy(writers%widom_successes)
    end subroutine destroy

end module procedures_exploring_writers_factory
