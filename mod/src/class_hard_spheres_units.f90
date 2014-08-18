!> \brief Description of the Hard_Spheres_Units class

module class_hard_spheres_units

use module_geometry, only: geometry

implicit none
private

    type, public :: Hard_Spheres_Units
    
        integer :: observables_thermalisation
        integer :: observables_equilibrium
        integer :: move_delta

        integer :: snap_initial_positions
        integer :: snap_final_positions
        integer :: snap_equilibrium_positions
        
        integer :: potential_energy
    
    contains
    
        procedure :: open => Hard_Spheres_Units_open
        procedure :: close => Hard_Spheres_Units_close
    
    end type Hard_Spheres_Units
    
    type, extends(Hard_Spheres_Units), public :: Dipolar_Hard_Spheres_Units
    
        integer :: rotate_delta

        integer :: snap_initial_orientations
        integer :: snap_final_orientations
        integer :: snap_equilibrium_orientations
        
        integer :: potential_energy_real

        integer :: wave_vectors
        integer :: ELC_wave_vectors
        integer :: structure_modulus
        integer :: ELC_structure_modulus
        
        integer :: total_moment_modulus
    
    end type Dipolar_Hard_Spheres_Units
    
    type, public :: Between_Hard_Spheres_Units
    
        integer :: observables_thermalisation
        integer :: observables_equilibrium
        
        integer :: potential_energy_tabulation
        
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
        
        open(newunit=this%snap_initial_positions, recl=4096, file=name//"_snap_initial_positions.out", &
             status='new', action='write')
        open(newunit=this%snap_final_positions, recl=4096, file=name//"_snap_final_positions.out", &
             status='new', action='write')
        open(newunit=this%snap_equilibrium_positions, recl=4096, &
             file=name//"_snap_equilibrium_positions.shots", status='new', action='write')
        
        open(newunit=this%potential_energy, recl=4096, &
             file=name//"_potential_energy.tmp", status='new', action='write')
        
        select type (this)

            type is (Dipolar_Hard_Spheres_Units)
            
                open(newunit=this%rotate_delta, recl=4096, file=name//"_rotate_delta.out", &
                    status='new', action='write')
                
                open(newunit=this%snap_initial_orientations, recl=4096, &
                     file=name//"_snap_initial_orientations.out", status='new', action='write')
                open(newunit=this%snap_final_orientations, recl=4096, &
                     file=name//"_snap_final_orientations.out", status='new', action='write')
                open(newunit=this%snap_equilibrium_orientations, recl=4096, &
                     file=name//"_snap_equilibrium_orientations.shots", status='new', action='write')
                     
                open(newunit=this%potential_energy_real, recl=4096, &
                     file=name//"_potential_energy_real.tmp", status='new', action='write')

                open(newunit=this%structure_modulus, recl=4096, file=name//"_structure_modulus.out", &
                     status='new', action='write')
                open(newunit=this%wave_vectors, recl=4096, file=name//"_wave_vectors.tmp", &
                     status='new', action='write')
                if (geometry%slab) then
                    open(newunit=this%ELC_structure_modulus, recl=4096, &
                         file=name//"_ELC_structure_modulus.out", status='new', action='write')
                    open(newunit=this%ELC_wave_vectors, recl=4096, file=name//"_ELC_wave_vectors.tmp", &
                        status='new', action='write')
                end if
                     
                open(newunit=this%total_moment_modulus, recl=4096, &
                     file=name//"_total_moment_modulus.out", status='new', action='write')
                
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
        
    end subroutine Between_Hard_Spheres_Units_open
    
    subroutine Hard_Spheres_Units_close(this)
    
        class(Hard_Spheres_Units), intent(inout) :: this
        
        close(this%observables_thermalisation)
        close(this%observables_equilibrium)
        close(this%move_delta)

        close(this%snap_initial_positions)
        close(this%snap_final_positions)
        close(this%snap_equilibrium_positions)
        
        close(this%potential_energy)
        
        select type (this)
            
            type is (Dipolar_Hard_Spheres_Units)
            
                close(this%rotate_delta)
        
                close(this%snap_initial_orientations)
                close(this%snap_final_orientations)
                close(this%snap_equilibrium_orientations)
                
                close(this%potential_energy_real)
                
                close(this%structure_modulus)
                close(this%wave_vectors)
                if (geometry%slab) then
                    close(this%ELC_structure_modulus)
                    close(this%ELC_wave_vectors)
                end if
                
                close(this%total_moment_modulus)
                
        end select
    
    end subroutine Hard_Spheres_Units_close
    
    subroutine Between_Hard_Spheres_Units_close(this)
    
        class(Between_Hard_Spheres_Units), intent(inout) :: this
        
        close(this%potential_energy_tabulation)
        close(this%observables_thermalisation)
        close(this%observables_equilibrium)
        
    end subroutine Between_Hard_Spheres_Units_close

end module class_hard_spheres_units
