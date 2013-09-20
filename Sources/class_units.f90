!> \brief Description of the Units class

module class_Units

implicit none
private

    type, public :: Units
    
        integer :: obsThermal
        integer :: obsEquilib
        integer :: move_delta

        integer :: snapIni_positions
        integer :: snapFin_positions
        integer :: snap_positions
        
        integer :: report
        integer :: Epot
    
    contains
    
        procedure :: open => Units_open
        procedure :: close => Units_close
    
    end type Units
    
    type, extends(Units), public :: MoreUnits
    
        integer :: rotate_delta

        integer :: snapIni_orientations
        integer :: snapFin_orientations
        integer :: snap_orientations

        integer :: waveVectors
        integer :: structure_modulus
        
        integer :: totalMoment_modulus
    
    end type MoreUnits
    
contains

    subroutine Units_open(this, name)
    
        class(Units), intent(out) :: this
        character(len=*), intent(in) :: name
        
        open(newunit=this%obsThermal, recl=4096, file=name//"_obsThermal.out", status='new', &
             action='write')
        open(newunit=this%obsEquilib, recl=4096, file=name//"_obsEquilib.out", status='new', &
             action='write')
        open(newunit=this%move_delta, recl=4096, file=name//"_move_delta.out", status='new', &
            action='write')
        
        open(newunit=this%snapIni_positions, recl=4096, file=name//"_snapIni_positions.out", &
             status='new', action='write')
        open(newunit=this%snapFin_positions, recl=4096, file=name//"_snapFin_positions.out", &
             status='new', action='write')
        open(newunit=this%snap_positions, recl=4096, file=name//"_snap_positions.shots", &
             status='new', action='write')
        
        open(newunit=this%report, recl=4096, file=name//"_report.txt", status='new', action='write')
        open(newunit=this%Epot, recl=4096, file=name//"_Epot.tmp", status='new', action='write')
        
        select type (this)

            type is (MoreUnits)
            
                open(newunit=this%rotate_delta, recl=4096, file=name//"_rotate_delta.out", &
                    status='new', action='write')
                
                open(newunit=this%snapIni_orientations, recl=4096, &
                     file=name//"_snapIni_orientations.out", status='new', action='write')
                open(newunit=this%snapFin_orientations, recl=4096, &
                     file=name//"_snapFin_orientations.out", status='new', action='write')
                open(newunit=this%snap_orientations, recl=4096, &
                     file=name//"_snap_orientations.shots", status='new', action='write')

                open(newunit=this%structure_modulus, recl=4096, file=name//"_structure_modulus.out", &
                     status='new', action='write')
                open(newunit=this%waveVectors, recl=4096, file=name//"_waveVectors.tmp", &
                     status='new', action='write')
                     
                open(newunit=this%totalMoment_modulus, recl=4096, &
                     file=name//"_totalMoment_modulus.out", status='new', action='write')
                
        end select
        
    end subroutine Units_open
    
    subroutine Units_close(this)
    
        class(Units), intent(inout) :: this
        
        close(this%obsThermal)
        close(this%obsEquilib)
        close(this%move_delta)

        close(this%snapIni_positions)
        close(this%snapFin_positions)
        close(this%snap_positions)
        
        close(this%report)
        close(this%Epot)
        
        select type (this)
        
            type is (Units)
            
            class is (MoreUnits)
            
                close(this%rotate_delta)
        
                close(this%snapIni_orientations)
                close(this%snapFin_orientations)
                close(this%snap_orientations)

                close(this%waveVectors)
                close(this%structure_modulus)
                
                close(this%totalMoment_modulus)
                    
            class default
                
                stop "Units_close : expected type for Units object !"
                
        end select
    
    end subroutine Units_close

end module class_Units
