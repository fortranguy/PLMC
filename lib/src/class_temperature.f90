module class_temperature

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_checks, only: check_positive

implicit none

    type, abstract, public :: Abstract_Temperature
    private
        real(DP) :: temperature
    contains
        procedure :: set => Abstract_Temperature_set
        procedure :: get => Abstract_Temperature_get
    end type Abstract_Temperature
    
    type, extends(Abstract_Temperature), public :: Concrete_Temperature
    
    end type Concrete_Temperature
    
contains

    subroutine Abstract_Temperature_set(this, temperature)
        class(Abstract_Temperature), intent(inout) :: this
        real(DP), intent(in) :: temperature
        
        call check_positive("Abstract_Temperature", "temperature", temperature)
        this%temperature = temperature
    end subroutine Abstract_Temperature_set
    
    pure function Abstract_Temperature_get(this) result(temperature)
        class(Abstract_Temperature), intent(in) :: this
        real(DP) :: temperature
        
        temperature = this%temperature
    end function Abstract_Temperature_get

end module class_temperature
