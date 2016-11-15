module classes_dipoles_neighbourhood

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_checks, only: check_positive
use procedures_dipolar_interactions_micro, only: dipolar_energy_is_negative

implicit none

private

    type, abstract, public :: Abstract_Dipolar_Neighbourhood
    private
        real(DP) :: max_distance = 0._DP
    contains
        procedure :: set => Abstract_set
        procedure :: meet => Abstract_meet
    end type Abstract_Dipolar_Neighbourhood

    type, extends(Abstract_Dipolar_Neighbourhood), public :: Concrete_Dipolar_Neighbourhood

    end type Concrete_Dipolar_Neighbourhood

    type, extends(Abstract_Dipolar_Neighbourhood), public :: Null_Dipolar_Neighbourhood
    contains
        procedure :: set => Null_set
        procedure :: are_neighbour => Null_are_neighbour
    end type Null_Dipolar_Neighbourhood

    type, public :: Dipolar_Neighbourhood_Wrapper
        class(Abstract_Dipolar_Neighbourhood), allocatable :: neighbourhood
    end type Dipolar_Neighbourhood_Wrapper

    type, public :: Dipolar_Neighbourhood_Line
        type(Dipolar_Neighbourhood_Wrapper), allocatable :: line(:)
    end type Dipolar_Neighbourhood_Line

contains

!implementation Abstract_Dipolar_Neighbourhood

    subroutine Abstract_set(this, max_distance)
        class(Abstract_Dipolar_Neighbourhood), intent(inout) :: this
        real(DP), intent(in) :: max_distance

        call check_positive("Abstract_Dipolar_Neighbourhood: set", "max_distance", max_distance)
        this%max_distance = max_distance
    end subroutine Abstract_set

    !> @todo Make \( \vec{r}_ij \cdot \mu_i > 0 \) stricter and symmetric?
    pure subroutine Abstract_meet(this, overlap, ij_are_neighbour, min_distance, vector_ij, &
        orientation_i, orientation_j)
        class(Abstract_Dipolar_Neighbourhood), intent(in) :: this
        logical, intent(out) :: overlap, ij_are_neighbour
        real(DP), intent(in) :: min_distance
        real(DP), dimension(:), intent(in) :: vector_ij,  orientation_i, orientation_j

        real(DP) :: distance_ij

        distance_ij = norm2(vector_ij)
        if (distance_ij < min_distance) then
            overlap = .true.
            return
        else
            overlap = .false.
        end if

        ij_are_neighbour = distance_ij < this%max_distance .and. &
            dipolar_energy_is_negative(vector_ij, orientation_i, orientation_j) .and. &
            dot_product(vector_ij, orientation_i) > 0._DP
    end subroutine Abstract_meet

!end implementation Abstract_Dipolar_Neighbourhood

!implementation Null_Dipolar_Neighbourhood

    subroutine Null_set(this, max_distance)
        class(Null_Dipolar_Neighbourhood), intent(inout) :: this
        real(DP), intent(in) :: max_distance
    end subroutine Null_set

    pure subroutine Null_are_neighbour(this, overlap, ij_are_neighbour, min_distance, vector_ij,&
        orientation_i, orientation_j)
        class(Null_Dipolar_Neighbourhood), intent(in) :: this
        logical, intent(out) :: overlap, ij_are_neighbour
        real(DP), intent(in) :: min_distance
        real(DP), dimension(:), intent(in) :: vector_ij,  orientation_i, orientation_j
        overlap = .false.; ij_are_neighbour = .false.
    end subroutine Null_are_neighbour

!end implementation Null_Dipolar_Neighbourhood

end module classes_dipoles_neighbourhood
