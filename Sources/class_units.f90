!> \brief Description of the Units class

module class_Units

implicit none

private

    type, public :: Units
    
        integer :: obs
        integer :: obsTherm
        integer :: dx

        integer :: snapIni_X
        integer :: snapFin_X
        integer :: snapShots_X
        
        integer :: report
        integer :: Epot
    
    contains
    
        procedure :: open => Units_open
        procedure :: close => Units_close
    
    end type Units
    
    type, extends(Units), public :: MoreUnits
    
        integer :: dm

        integer :: snapIni_M
        integer :: snapFin_M
        integer :: snapShots_M
    
    contains
    
        procedure :: openMore => MoreUnits_openMore
        procedure :: closeMore => MoreUnits_closeMore
    
    end type MoreUnits
    
contains

    subroutine Units_open(this, name)
    
        class(Units), intent(out) :: this
        character(len=*), intent(in) :: name
        
        open(newunit=this%obs, recl=4096, file=name//"_obs.out", status='new', action='write')
        open(newunit=this%obsTherm, recl=4096, file=name//"_obsTherm.out", status='new', action='write')
        open(newunit=this%dx, recl=4096, file=name//"_dx.out", status='new', action='write')
        
        open(newunit=this%snapIni_X, recl=4096, file=name//"_snapIni_X.out", status='new', action='write')
        open(newunit=this%snapFin_X, recl=4096, file=name//"_snapFin_X.out", status='new', action='write')
        open(newunit=this%snapShots_X, recl=4096, file=name//"_snap_X.shots", status='new', action='write')
        
        open(newunit=this%report, recl=4096, file=name//"_report.out", status='new', action='write')
        open(newunit=this%Epot, recl=4096, file=name//"_Epot.out", status='new', action='write')
        
    end subroutine Units_open
    
    subroutine MoreUnits_openMore(this, name)
    
        class(MoreUnits), intent(inout) :: this
        character(len=*), intent(in) :: name
        
        open(newunit=this%dm, recl=4096, file=name//"_dm.out", status='new', action='write')
        
        open(newunit=this%snapIni_M, recl=4096, file=name//"_snapIni_M.out", status='new', action='write')
        open(newunit=this%snapFin_M, recl=4096, file=name//"_snapFin_M.out", status='new', action='write')
        open(newunit=this%snapShots_M, recl=4096, file=name//"_snap_M.shots", status='new', action='write')
    
    end subroutine MoreUnits_openMore
    
    subroutine Units_close(this)
    
        class(Units), intent(inout) :: this
        
        close(this%obs)
        close(this%obsTherm)
        close(this%dx)

        close(this%snapIni_X)
        close(this%snapFin_X)
        close(this%snapShots_X)
        
        close(this%report)
        close(this%Epot)
    
    end subroutine Units_close
    
    subroutine MoreUnits_closeMore(this)
    
        class(MoreUnits), intent(inout) :: this
        
        close(this%dm)
        
        close(this%snapIni_M)
        close(this%snapFin_M)
        close(this%snapShots_M)
    
    end subroutine MoreUnits_closeMore

end module class_Units
