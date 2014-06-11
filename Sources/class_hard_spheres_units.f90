!> \brief Description of the Hard_Spheres_Units class

module class_hard_spheres_units

implicit none
private

    type, public :: Hard_Spheres_Units
    
        integer :: observables_thermalisation
        integer :: observables_equilibrium
        integer :: move_delta

        integer :: snapIni_positions
        integer :: snapFin_positions
        integer :: snap_positions
        
        integer :: report
        integer :: potential_energy
    
    contains
    
        procedure :: open => Hard_Spheres_Units_open
        procedure :: close => Hard_Spheres_Units_close
    
    end type Hard_Spheres_Units
    
    type, extends(Hard_Spheres_Units), public :: Dipolar_Hard_Spheres_Units
    
        integer :: rotate_delta

        integer :: snapIni_orientations
        integer :: snapFin_orientations
        integer :: snap_orientations
        
        integer :: potential_energy_real

        integer :: waveVectors
        integer :: structure_modulus
        
        integer :: totalMoment_modulus
    
    end type Dipolar_Hard_Spheres_Units
    
    type, public :: Between_Hard_Spheres_Units
    
        integer :: observables_thermalisation
        integer :: observables_equilibrium
        
        integer :: potential_energy_tabulation
        integer :: report
        
    contains
        
        procedure :: open => Between_Hard_Spheres_Units_open
        procedure :: close => Between_Hard_Spheres_Units_close
    
    end type Between_Hard_Spheres_Units
    
contains

    subroutine Hard_Spheres_Units_open(this, name)
    
        class(Hard_Spheres_Units), intent(out) :: this
        character(len=*), intent(in) :: name
        
        open(newunit=this%observables_thermalisation, recl=4096, &
             file=name//"_observables_thermalisation.out", status='new', action='write')
        open(newunit=this%observables_equilibrium, recl=4096, &
             file=name//"_observables_equilibrium.out", status='new', action='write')
        open(newunit=this%move_delta, recl=4096, file=name//"_move_delta.out", status='new', &
            action='write')
        
        open(newunit=this%snapIni_positions, recl=4096, file=name//"_snapIni_positions.out", &
             status='new', action='write')
        open(newunit=this%snapFin_positions, recl=4096, file=name//"_snapFin_positions.out", &
             status='new', action='write')
        open(newunit=this%snap_positions, recl=4096, file=name//"_snap_positions.shots", &
             status='new', action='write')
        
        open(newunit=this%report, recl=4096, file=name//"_report.txt", status='new', action='write')
        open(newunit=this%potential_energy, recl=4096, &
             file=name//"_potential_energy.tmp", status='new', action='write')
        
        select type (this)

            type is (Dipolar_Hard_Spheres_Units)
            
                open(newunit=this%rotate_delta, recl=4096, file=name//"_rotate_delta.out", &
                    status='new', action='write')
                
                open(newunit=this%snapIni_orientations, recl=4096, &
                     file=name//"_snapIni_orientations.out", status='new', action='write')
                open(newunit=this%snapFin_orientations, recl=4096, &
                     file=name//"_snapFin_orientations.out", status='new', action='write')
                open(newunit=this%snap_orientations, recl=4096, &
                     file=name//"_snap_orientations.shots", status='new', action='write')
                     
                open(newunit=this%potential_energy_real, recl=4096, &
                     file=name//"_potential_energy_real.tmp", status='new', action='write')

                open(newunit=this%structure_modulus, recl=4096, file=name//"_structure_modulus.out", &
                     status='new', action='write')
                open(newunit=this%waveVectors, recl=4096, file=name//"_waveVectors.tmp", &
                     status='new', action='write')
                     
                open(newunit=this%totalMoment_modulus, recl=4096, &
                     file=name//"_totalMoment_modulus.out", status='new', action='write')
                
                write(this%observables_equilibrium, *) "#", 4 ! 4 observables
                
            class default
            
                write(this%observables_equilibrium, *) "#", 3 ! 3 observables
                
        end select
        
    end subroutine Hard_Spheres_Units_open
    
    subroutine Between_Hard_Spheres_Units_open(this, name)
    
        class(Between_Hard_Spheres_Units), intent(out) :: this
        character(len=*), intent(in) :: name
    
        open(newunit=this%observables_thermalisation, recl=4096, &
             file=name//"_observables_thermalisation.out", status='new', action='write')
        open(newunit=this%observables_equilibrium, recl=4096, &
             file=name//"_observables_equilibrium.out", &
             status='new', action='write')
        write(this%observables_equilibrium, *) "#", 1 ! 1 observable: energy
        
        open(newunit=this%potential_energy_tabulation, recl=4096, &
             file=name//"_potential_energy_tabulation.tmp", status='new', action='write')
        open(newunit=this%report, recl=4096, file=name//"_report.txt", status='new', action='write')
        
    end subroutine Between_Hard_Spheres_Units_open
    
    subroutine Hard_Spheres_Units_close(this)
    
        class(Hard_Spheres_Units), intent(inout) :: this
        
        close(this%observables_thermalisation)
        close(this%observables_equilibrium)
        close(this%move_delta)

        close(this%snapIni_positions)
        close(this%snapFin_positions)
        close(this%snap_positions)
        
        close(this%report)
        close(this%potential_energy)
        
        select type (this)
            
            type is (Dipolar_Hard_Spheres_Units)
            
                close(this%rotate_delta)
        
                close(this%snapIni_orientations)
                close(this%snapFin_orientations)
                close(this%snap_orientations)
                
                close(this%potential_energy_real)

                close(this%waveVectors)
                close(this%structure_modulus)
                
                close(this%totalMoment_modulus)
                
        end select
    
    end subroutine Hard_Spheres_Units_close
    
    subroutine Between_Hard_Spheres_Units_close(this)
    
        class(Between_Hard_Spheres_Units), intent(inout) :: this
        
        close(this%report)
        close(this%potential_energy_tabulation)
        close(this%observables_thermalisation)
        close(this%observables_equilibrium)
        
    end subroutine Between_Hard_Spheres_Units_close

end module class_hard_spheres_units
