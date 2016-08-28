module classes_maximum_box_compression_explorer

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_errors, only: error_exit
use types_environment_wrapper, only: Environment_Wrapper
use types_component_wrapper, only: Component_Wrapper
use types_short_interactions_wrapper, only: Short_Interactions_Wrapper
use classes_maximum_box_compression, only: Abstract_Maximum_Box_Compression
use types_exploring_observables_wrapper, only: Exploring_Observables_Wrapper
use classes_exploring_algorithm, only: Abstract_Exploring_Algorithm
use procedures_plmc_visit, only: visit_short

implicit none

private

    type, extends(Abstract_Exploring_Algorithm), abstract, public :: &
        Abstract_Maximum_Box_Compression_Explorer
    private
        type(Environment_Wrapper), pointer :: environment => null()
        type(Component_Wrapper), pointer :: components(:) => null()
        type(Short_Interactions_Wrapper), pointer :: short_interactions => null()
        class(Abstract_Maximum_Box_Compression), allocatable :: maximum_box_compression
        real(DP) :: min_distance
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: try => Abstract_try
        procedure, private :: set_min_distance => Abstract_set_min_distance
    end type Abstract_Maximum_Box_Compression_Explorer

    type, extends(Abstract_Maximum_Box_Compression_Explorer), public :: &
        Concrete_Maximum_Box_Compression_Explorer
    end type Concrete_Maximum_Box_Compression_Explorer

    type, extends(Abstract_Maximum_Box_Compression_Explorer), public :: &
        Null_Maximum_Box_Compression_Explorer
    contains
        procedure :: Null_construct
        procedure :: Null_destroy
        procedure :: Null_try
    end type Null_Maximum_Box_Compression_Explorer

contains

!implementation Abstract_Maximum_Box_Compression_Explorer

    subroutine Abstract_construct(this, environment, components, short_interactions, &
        maximum_box_compression)
        class(Abstract_Maximum_Box_Compression_Explorer), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        type(Component_Wrapper), target, intent(in) :: components(:)
        type(Short_Interactions_Wrapper), target, intent(in) :: short_interactions
        class(Abstract_Maximum_Box_Compression), intent(in) :: maximum_box_compression

        this%environment => environment
        this%components => components
        this%short_interactions => short_interactions
        call this%set_min_distance()
        allocate(this%maximum_box_compression, source=maximum_box_compression)
    end subroutine Abstract_construct

    subroutine Abstract_set_min_distance(this)
        class(Abstract_Maximum_Box_Compression_Explorer), intent(inout) :: this

        integer :: i_component, j_component
        real(DP) :: min_distance_ij

        this%min_distance = this%environment%periodic_box%get_max_distance()
        do j_component = 1, size(this%components)
            do i_component = 1, j_component
                min_distance_ij = this%short_interactions%components_pairs(j_component)%&
                    line(i_component)%potential%get_min_distance()
                if (min_distance_ij < this%min_distance) this%min_distance = min_distance_ij
            end do
        end do
    end subroutine Abstract_set_min_distance

    subroutine Abstract_destroy(this)
        class(Abstract_Maximum_Box_Compression_Explorer), intent(inout) :: this

        if (allocated(this%maximum_box_compression)) deallocate(this%maximum_box_compression)
        this%short_interactions => null()
        this%components => null()
        this%environment => null()
    end subroutine Abstract_destroy

    subroutine Abstract_try(this, observables)
        class(Abstract_Maximum_Box_Compression_Explorer), intent(in) :: this
        type(Exploring_Observables_Wrapper), intent(inout) :: observables

        logical :: overlap
        real(DP) :: min_distance_ratio, max_distance_ratio

        max_distance_ratio = this%environment%periodic_box%get_max_distance() / this%min_distance
        call visit_short(overlap, min_distance_ratio, max_distance_ratio, this%components, this%&
            short_interactions)
        if (overlap) call error_exit("Abstract_Maximum_Box_Compression_Explorer: try: "//&
            "visit_short: overlap")
        observables%maximum_box_compression_delta = this%maximum_box_compression%&
            get_delta(min_distance_ratio)
    end subroutine Abstract_try

!end implementation Abstract_Maximum_Box_Compression_Explorer

!implementation Null_Maximum_Box_Compression_Explorer

    subroutine Null_construct(this, environment, components, short_interactions)
        class(Null_Maximum_Box_Compression_Explorer), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        type(Component_Wrapper), target, intent(in) :: components(:)
        type(Short_Interactions_Wrapper), target, intent(in) :: short_interactions
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Maximum_Box_Compression_Explorer), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_try(this, min_distance_ratio)
        class(Null_Maximum_Box_Compression_Explorer), intent(in) :: this
        real(DP), intent(out) :: min_distance_ratio
        min_distance_ratio = 0._DP
    end subroutine Null_try

!implementation Null_Maximum_Box_Compression_Explorer

end module classes_maximum_box_compression_explorer