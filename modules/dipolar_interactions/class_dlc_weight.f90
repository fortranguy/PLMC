module class_dlc_weight

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions, PI
use class_periodic_box, only: Abstract_Periodic_Box
use class_reciprocal_lattice, only: Abstract_Reciprocal_Lattice
use class_permittivity, only: Abstract_Permittivity

implicit none

private

    type, abstract, public :: Abstract_DLC_Weight
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
        integer :: reci_numbers(num_dimensions) = 0
        real(DP) :: permittivity = 0._DP
        real(DP), dimension(:, :), allocatable :: weight
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: reset => Abstract_set
        procedure :: get => Abstract_get
        procedure, private :: set => Abstract_set
    end type Abstract_DLC_Weight

    type, extends(Abstract_DLC_Weight), public :: Concrete_DLC_Weight

    end type Concrete_DLC_Weight

    type, extends(Abstract_DLC_Weight), public :: Null_DLC_Weight
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: reset => Null_set
        procedure :: get => Null_get
    end type Null_DLC_Weight

contains

!implementation Abstract_DLC_Weight

    subroutine Abstract_construct(this, periodic_box, reciprocal_lattice, permittivity)
        class(Abstract_DLC_Weight), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Reciprocal_Lattice), intent(in) :: reciprocal_lattice
        class(Abstract_Permittivity), intent(in) :: permittivity

        this%periodic_box => periodic_box
        this%reci_numbers = reciprocal_lattice%get_numbers()
        this%permittivity = permittivity%get()
        allocate(this%weight(0:this%reci_numbers(1), 0:this%reci_numbers(2)))
        call this%set()
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_DLC_Weight), intent(inout) :: this

        if (allocated(this%weight)) deallocate(this%weight)
        this%periodic_box => null()
    end subroutine Abstract_destroy

    subroutine Abstract_set(this)
        class(Abstract_DLC_Weight), intent(inout) :: this

        real(DP) :: box_size(num_dimensions), wave_vector(2), wave_norm
        integer :: n_1, n_2

        box_size = this%periodic_box%get_size()
        do n_2 = 0, this%reci_numbers(2)
            wave_vector(2) = 2._DP*PI * real(n_2, DP) / box_size(2)
        do n_1 = 0, this%reci_numbers(1)
            wave_vector(1) = 2._DP*PI * real(n_1, DP) / box_size(1)

            if (n_1**2 + n_2**2 > this%reci_numbers(1)**2) cycle

            if (n_1 == 0 .and. n_2 == 0) then
                this%weight(n_1, n_2) = 0._DP
            else
                wave_norm = norm2(wave_vector)
                this%weight(n_1, n_2) = 1._DP / (2._DP*this%permittivity*product(box_size(1:2))) / &
                    (wave_norm * (exp(wave_norm*box_size(3)) - 1._DP))
            end if

        end do
        end do
    end subroutine Abstract_set

    !> \[
    !>      w(\vec{k}_{1:2}) =
    !>          \begin{cases}
    !>              0   &   \text{if } \vec{k}_{1:2} = \vec{0} \\
    !>              \frac{1}{2\epsilon S} \frac{1}{k_{1:2} (e^{k_{1:2} L_3} - 1)}   & \text{else}
    !>          \end{cases}
    !> \]
    pure real(DP) function Abstract_get(this, n_1, n_2) result(weight)
        class(Abstract_DLC_Weight), intent(in) :: this
        integer, intent(in) :: n_1, n_2

        weight = this%weight(abs(n_1), abs(n_2))
    end function Abstract_get

!end implementation Abstract_DLC_Weight

!implementation Null_DLC_Weight

    subroutine Null_construct(this, periodic_box, reciprocal_lattice, permittivity)
        class(Null_DLC_Weight), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Reciprocal_Lattice), intent(in) :: reciprocal_lattice
        class(Abstract_Permittivity), intent(in) :: permittivity
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_DLC_Weight), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_set(this)
        class(Null_DLC_Weight), intent(inout) :: this
    end subroutine Null_set

    pure real(DP) function Null_get(this, n_1, n_2) result(weight)
        class(Null_DLC_Weight), intent(in) :: this
        integer, intent(in) :: n_1, n_2
        weight = 0._DP
    end function Null_get

!end implementation Null_DLC_Weight

end module class_dlc_weight
