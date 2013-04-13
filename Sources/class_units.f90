!> \brief Description of the Units class

module unitsCounter

    integer :: iUnit = 10

end module unitsCounter

module class_Units

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
    
        procedure :: open => Units_open
        procedure :: close => Units_close
    
    end type Units
    
contains

    subroutine Units_open(this)
    
        class(Units), intent(out) :: this
        
        this%obs = iUnit ;      iUnit = iUnit + 1
        this%obsTherm = iUnit ; iUnit = iUnit + 1
        this%dx = iUnit ;       iUnit = iUnit + 1
        
        this%snapIni = iUnit ;  iUnit = iUnit + 1
        this%snapFin = iUnit ;  iUnit = iUnit + 1
        this%snapShots = iUnit ;iUnit = iUnit + 1
        
        this%report = iUnit ;   iUnit = iUnit + 1
        
        open(unit=this%obs, recl=4096, file="obs.out", status='new', &
            action='write')
        open(unit=this%obsTherm, recl=4096, file="obsTherm.out", status='new', &
            action='write')
        open(unit=this%dx, recl=4096, file="dx.out", status='new', &
            action='write')
        
        open(unit=this%snapIni, recl=4096, file="snapShotIni.out", &
            status='new', action='write')
        open(unit=this%snapFin, recl=4096, file="snapShotFin.out", &
            status='new', action='write')
        open(unit=this%snapShots, recl=4096, file="snap.shot", status='new', &
            action='write')
            
        open(unit=this%report, recl=4096, file="report.out", status='new', &
            action='write')
        
    end subroutine Units_open
    
    subroutine Units_close(this)
    
        class(Units), intent(inout) :: this
        
        close(this%obs)
        close(this%obsTherm)
        close(this%dx)
        
        close(this%snapIni)
        close(this%snapFin)
        close(this%snapShots)
        
        close(this%report)
    
    end subroutine Units_close

end module class_Units