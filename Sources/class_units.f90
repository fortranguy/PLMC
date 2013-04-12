!> \brief Description of the Units class

module unitsCounter

    integer :: iUnit = 10
    
end firstUnit

module unitsCounter

use unitsCounter

implicit none

private

    type, public :: Units
    
        integer :: obs
        integer :: obsTherm
        integer :: dx
        
        integer :: snapIni
        integer :: snapFin  
        integer :: snapShots              
        
        integer :: report        
    
    contains
    
        procedure :: init => Units_init
    
    end type Units
    
contains

    subroutine Units_init(this)
    
        class(Units), intent(out) :: this
        
        this%obs = iUnit
        this%obsTherm = iUnit
        this%dx = iUnit
        
        this%snapIni = iUnit
        this%snapFin = iUnit
        this%snapShots = iUnit
        
        this%report = iUnit
    
    end subroutine Units_init

end module class_Units
