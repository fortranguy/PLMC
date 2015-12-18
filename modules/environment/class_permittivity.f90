module class_permittivity

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_checks, only: check_positive

implicit none

private

    type, abstract, public :: Abstract_Permittivity
    private
        real(DP) :: permittivity
    contains
        procedure :: set => Abstract_set
        procedure :: get => Abstract_get
    end type Abstract_Permittivity

    type, extends(Abstract_Permittivity), public :: Concrete_Permittivity

    end type Concrete_Permittivity

    type, extends(Abstract_Permittivity), public :: Null_Permittivity
    contains
        procedure :: set => Null_set
        procedure :: get => Null_get
    end type Null_Permittivity

contains

!implementation Abstract_Permittivity

    subroutine Abstract_set(this, permittivity)
        class(Abstract_Permittivity), intent(inout) :: this
        real(DP), intent(in) :: permittivity

        call check_positive("Abstract_set", "permittivity", permittivity)
        this%permittivity = permittivity
    end subroutine Abstract_set

    pure real(DP) function Abstract_get(this) result(permittivity)
        class(Abstract_Permittivity), intent(in) :: this

        permittivity = this%permittivity
    end function Abstract_get

!end implementation Abstract_Permittivity

!implementation Null_Permittivity

    subroutine Null_set(this, permittivity)
        class(Null_Permittivity), intent(inout) :: this
        real(DP), intent(in) :: permittivity
    end subroutine Null_set

    pure real(DP) function Null_get(this) result(permittivity)
        class(Null_Permittivity), intent(in) :: this
        permittivity = 0._DP
    end function Null_get

!end implementation Null_Permittivity

end module class_permittivity
