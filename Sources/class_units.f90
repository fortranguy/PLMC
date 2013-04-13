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
        
        this%obs = iUnit ;      iUnit = iUnit + 1
        this%obsTherm = iUnit ; iUnit = iUnit + 1
        this%dx = iUnit ;       iUnit = iUnit + 1
        
        this%snapIni = iUnit ;  iUnit = iUnit + 1
        this%snapFin = iUnit ;  iUnit = iUnit + 1
        this%snapShots = iUnit ;iUnit = iUnit + 1
        
        this%report = iUnit ;   iUnit = iUnit + 1
    
    end subroutine Units_init

end module class_Units
