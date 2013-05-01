!> \brief Description of the Units class

module class_Units

implicit none

private

    type, public :: Units
    
        integer :: obs
        integer :: obsTherm
        integer :: dx

        integer :: configOld
        integer :: snapIni
        integer :: snapFin
        integer :: snapShots
        
        integer :: report
        integer :: Epot
    
    contains
    
        procedure :: open => Units_open
        procedure :: close => Units_close
    
    end type Units
    
contains

    subroutine Units_open(this, name)
    
        class(Units), intent(out) :: this
        character(len=*), intent(in) :: name
        
        open(newunit=this%obs, recl=4096, file=name//"_obs.out", status='new', action='write')
        open(newunit=this%obsTherm, recl=4096, file=name//"_obsTherm.out", status='new', action='write')
        open(newunit=this%dx, recl=4096, file=name//"_dx.out", status='new', action='write')
        
        open(newunit=this%snapIni, recl=4096, file=name//"_snapIni.out", status='new', action='write')
        open(newunit=this%snapFin, recl=4096, file=name//"_snapFin.out", status='new', action='write')
        open(newunit=this%snapShots, recl=4096, file=name//"_snap.shots", status='new', action='write')
        
        open(newunit=this%report, recl=4096, file=name//"_report.out", status='new', action='write')
        open(newunit=this%Epot, recl=4096, file=name//"_Epot.out", status='new', action='write')
        
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
        close(this%Epot)
    
    end subroutine Units_close

end module class_Units
