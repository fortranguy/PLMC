!> \brief Description of the Units class

module class_Units

implicit none

private

    type, public :: Units
    
        integer :: obsEqb
        integer :: obsTherm
        integer :: dx

        integer :: snapIni_positions
        integer :: snapFin_positions
        integer :: snapShots_positions
        
        integer :: report
        integer :: Epot
    
    contains
    
        procedure :: open => Units_open
        procedure :: close => Units_close
    
    end type Units
    
    type, extends(Units), public :: MoreUnits
    
        integer :: dm

        integer :: snapIni_moments
        integer :: snapFin_moments
        integer :: snapShots_moments
        
        integer :: structure_moduli
    
    end type MoreUnits
    
contains

    subroutine Units_open(this, name)
    
        class(Units), intent(out) :: this
        character(len=*), intent(in) :: name
        
        open(newunit=this%obsEqb, recl=4096, file=name//"_obsEqb.out", status='new', action='write')
        open(newunit=this%obsTherm, recl=4096, file=name//"_obsTherm.out", status='new', &
             action='write')
        open(newunit=this%dx, recl=4096, file=name//"_dx.out", status='new', action='write')
        
        open(newunit=this%snapIni_positions, recl=4096, file=name//"_snapIni_positions.out", &
             status='new', action='write')
        open(newunit=this%snapFin_positions, recl=4096, file=name//"_snapFin_positions.out", &
             status='new', action='write')
        open(newunit=this%snapShots_positions, recl=4096, file=name//"_snap_positions.shots", &
             status='new', action='write')
        
        open(newunit=this%report, recl=4096, file=name//"_report.txt", status='new', action='write')
        open(newunit=this%Epot, recl=4096, file=name//"_Epot.out", status='new', action='write')
        
        select type (this)
        
            type is (Units)
            
            class is (MoreUnits)
            
                open(newunit=this%dm, recl=4096, file=name//"_dm.out", status='new', action='write')
                
                open(newunit=this%snapIni_moments, recl=4096, file=name//"_snapIni_moments.out", &
                     status='new', action='write')
                open(newunit=this%snapFin_moments, recl=4096, file=name//"_snapFin_moments.out", &
                     status='new', action='write')
                open(newunit=this%snapShots_moments, recl=4096, file=name//"_snap_moments.shots", &
                     status='new', action='write')
                     
                 open(newunit=this%structure_moduli, recl=4096, file=name//"_structure_moduli.out", &
                    status='new', action='write')
                
            class default
                
                stop "Units_open : expected type for Units object !"
                
        end select
        
    end subroutine Units_open
    
    subroutine Units_close(this)
    
        class(Units), intent(inout) :: this
        
        close(this%obsEqb)
        close(this%obsTherm)
        close(this%dx)

        close(this%snapIni_positions)
        close(this%snapFin_positions)
        close(this%snapShots_positions)
        
        close(this%report)
        close(this%Epot)
        
        select type (this)
        
            type is (Units)
            
            class is (MoreUnits)
            
                close(this%dm)
        
                close(this%snapIni_moments)
                close(this%snapFin_moments)
                close(this%snapShots_moments)
                
                close(this%structure_moduli)
                    
            class default
                
                stop "Units_close : expected type for Units object !"
                
        end select
    
    end subroutine Units_close

end module class_Units
