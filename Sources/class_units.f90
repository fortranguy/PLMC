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
    
    end type MoreUnits
    
contains

    subroutine Units_open(this, name)
    
        class(Units), intent(out) :: this
        character(len=*), intent(in) :: name
        
        open(newunit=this%obs, recl=4096, file=name//"_obs.out", status='new', action='write')
        open(newunit=this%obsTherm, recl=4096, file=name//"_obsTherm.out", status='new', &
             action='write')
        open(newunit=this%dx, recl=4096, file=name//"_dx.out", status='new', action='write')
        
        open(newunit=this%snapIni_X, recl=4096, file=name//"_snapIni_X.out", status='new', &
             action='write')
        open(newunit=this%snapFin_X, recl=4096, file=name//"_snapFin_X.out", status='new', &
             action='write')
        open(newunit=this%snapShots_X, recl=4096, file=name//"_snap_X.shots", status='new', &
             action='write')
        
        open(newunit=this%report, recl=4096, file=name//"_report.txt", status='new', action='write')
        open(newunit=this%Epot, recl=4096, file=name//"_Epot.out", status='new', action='write')
        
        select type (this)
        
            type is (Units)
            
            class is (MoreUnits)
            
                open(newunit=this%dm, recl=4096, file=name//"_dm.out", status='new', action='write')
                
                open(newunit=this%snapIni_M, recl=4096, file=name//"_snapIni_M.out", status='new', &
                     action='write')
                open(newunit=this%snapFin_M, recl=4096, file=name//"_snapFin_M.out", status='new', &
                     action='write')
                open(newunit=this%snapShots_M, recl=4096, file=name//"_snap_M.shots", status='new', &
                     action='write')
                
            class default
                
                stop "Units_open : expected type for Units object !"
                
        end select
        
    end subroutine Units_open
    
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
        
        select type (this)
        
            type is (Units)
            
            class is (MoreUnits)
            
                close(this%dm)
        
                close(this%snapIni_M)
                close(this%snapFin_M)
                close(this%snapShots_M)
                    
            class default
                
                stop "Units_close : expected type for Units object !"
                
        end select
    
    end subroutine Units_close

end module class_Units
