!> \brief Description of the Physical System class

module class_physical_system

use, intrinsic :: iso_fortran_env, only: output_unit
use data_precisions, only: DP
use json_module, only: json_file
use module_types_micro, only: Box_Dimensions, Monte_Carlo_Arguments
use module_physics_micro, only: NwaveVectors
use module_data, only: test_data_found
use module_types_macro, only: Hard_Spheres_Macro, Dipolar_Hard_Spheres_Macro
use class_hard_spheres, only: Hard_Spheres, Dipolar_Hard_Spheres
use class_mixing_potential_energy, only: Mixing_Potential_Energy
use class_hard_spheres_observables, only: Hard_Spheres_Observables, Dipolar_Hard_Spheres_Observables
use class_hard_spheres_units, only: Hard_Spheres_Units, Dipolar_Hard_Spheres_Units
use module_monte_carlo_arguments, only: read_arguments
use module_physics_macro, only: init_random_seed, set_initial_configuration, &
                                init_spheres, init_hard_potential_energy, init_cells, &
                                set_ewald, &
                                total_energy, &
                                final_spheres, &
                                init_mix, mix_final, &
                                adapt_move, adapt_rotation, test_consist
use module_algorithms, only: move, widom, switch, rotate
use module_write, only: open_units, write_data, mix_open_units, write_results, mix_write_results, &
                        write_spheres_density

implicit none

private

    type, public :: Physical_System
    
        private
        
        character(len=5) :: name
        
        ! Box
        type(Box_Dimensions) :: Box
        integer :: num_particles
        
        ! Monte-Carlo
        integer :: num_thermalisation_steps, period_adaptation, num_equilibrium_steps
        integer :: decorrelation_factor, num_changes, num_moves, num_switches, num_rotations
        
        ! Type 1: Dipolar spheres
        type(Dipolar_Hard_Spheres) :: type1_spheres
        type(Dipolar_Hard_Spheres_Macro) :: type1_macro
        type(Dipolar_Hard_Spheres_Observables) :: type1_observables
        type(Dipolar_Hard_Spheres_Units) :: type1_units
        
        ! Type 2: Hard spheres
        type(Hard_Spheres) :: type2_spheres
        type(Hard_Spheres_Macro) :: type2_macro
        type(Hard_Spheres_Observables) :: type2_observables
        type(Hard_Spheres_Units) :: type2_units
        
        ! Mixing potential_energy
        type(Mixing_Potential_Energy) :: mix
        real(DP) :: mix_potential_energy, mix_potential_energy_sum, mix_potential_energy_conf
        integer :: mix_report_unit
        integer :: mix_potential_energy_tab_unit
        integer :: mix_observables_thermalisation_unit, mix_observables_equilibrium_unit
        
        ! Observables: write to files
        real(DP) :: potential_energy, potential_energy_sum
        integer :: report_unit
        integer :: observables_thermalisation_unit, observables_equilibrium_unit
        
        ! Switch
        integer :: switch_num_hits, switch_num_rejections
        real(DP) :: switch_rejection_rate, switch_sum_rejection
        
        logical :: write_potential_energy, snap
        integer :: reset_i_step
        real(DP) :: time_start, time_end
    
    contains
    
        !> Construction & destruction of the class
        procedure :: construct => Physical_System_construct
        procedure, private :: set_box => Physical_System_set_box
        procedure, private :: set_monte_carlo_steps => Physical_System_set_monte_carlo_steps
        procedure, private :: set_monte_carlo_changes => &
                              Physical_System_set_monte_carlo_changes        
        procedure :: destroy => Physical_System_destroy
        
        !> Initialization & Finalisation
        procedure :: init => Physical_System_init
        procedure, private :: open_all_units => Physical_System_open_all_units
        procedure, private :: write_all_reports => Physical_System_write_all_reports
        procedure, private :: write_report => Physical_System_write_report
        procedure, private :: init_switch => Physical_System_init_switch              
        procedure :: final => Physical_System_final
        procedure, private :: write_all_results => Physical_System_write_all_results
        procedure, private :: write_results => Physical_System_write_results
        procedure, private :: close_units => Physical_System_close_units
        
        !> Accessors & Mutators
        procedure :: get_num_thermalisation_steps => Physical_System_get_num_thermalisation_steps
        procedure :: get_num_equilibrium_steps => Physical_System_get_num_equilibrium_steps
        procedure :: set_time_start => Physical_System_set_time_start
        procedure :: set_time_end => Physical_System_set_time_end
        
        !> Simulation
        procedure :: random_changes => Physical_System_random_changes
        procedure :: update_rejections => Physical_System_update_rejections
        procedure :: adapt_changes => Physical_System_adapt_changes
        procedure :: write_observables_thermalisation => &
                     Physical_System_write_observables_thermalisation
        procedure :: fix_changes => Physical_System_fix_changes
        procedure :: measure_chemical_potentials => Physical_System_measure_chemical_potentials
        procedure :: accumulate_observables => Physical_System_accumulate_observables
        procedure :: write_observables_equilibrium => Physical_System_write_observables_equilibrium
        procedure :: take_snapshots => Physical_System_take_snapshots
        procedure :: reinitialize_quantites => Physical_System_reinitialize_quantites
    
    end type Physical_System
    
contains

    ! Construction
    
    subroutine Physical_System_construct(this, json)
            
        class(Physical_System), intent(out) :: this
        type(json_file), intent(inout) :: json
        
        character(len=4096) :: data_name
        logical :: found
        
        this%name = "sys"
        write(output_unit, *) this%name, " class construction"
        
        call this%set_box(json)
        call this%set_monte_carlo_steps(json)
        
        data_name = "Potential Energy.write"
        call json%get(data_name, this%write_potential_energy, found)
        call test_data_found(data_name, found)
        
        call this%type1_spheres%construct(json)
        call this%type2_spheres%construct(json)
        call this%mix%construct(json, this%type1_spheres%get_diameter(), &
                                      this%type2_spheres%get_diameter())
        
        call this%set_monte_carlo_changes(json)
        
    end subroutine Physical_System_construct
    
    subroutine Physical_System_set_box(this, json)
    
        class(Physical_System), intent(inout) :: this    
        type(json_file), intent(inout) :: json    
        
        character(len=4096) :: data_name
        logical :: found
        real(DP), dimension(:), allocatable :: Box_size
        integer, dimension(:), allocatable :: Box_wave
        
        data_name = "Box.size"
        call json%get(data_name, Box_size, found)
        call test_data_found(data_name, found)
        if (size(Box_size) /= size (this%Box%size)) error stop "Box size dimension"
        this%Box%size(:) = Box_size(:)
        
        data_name = "Box.wave"
        call json%get(data_name, Box_wave, found)
        call test_data_found(data_name, found)
        if (size(Box_wave) /= size (this%Box%wave)) error stop "Box wave dimension"
        this%Box%wave(:) = Box_wave(:)
        
        data_name = "Box.temperature"
        call json%get(data_name, this%Box%temperature, found)
        call test_data_found(data_name, found)                
        
    end subroutine Physical_System_set_box
    
    subroutine Physical_System_set_monte_carlo_steps(this, json)
        class(Physical_System), intent(inout) :: this
        type(json_file), intent(inout) :: json 
        
        character(len=4096) :: data_name
        logical :: found
        
        data_name = "Monte Carlo.number of thermalisation steps"
        call json%get(data_name, this%num_thermalisation_steps, found)
        call test_data_found(data_name, found)
        
        data_name = "Monte Carlo.period of adaptation"
        call json%get(data_name, this%period_adaptation, found)
        call test_data_found(data_name, found)
        
        data_name = "Monte Carlo.number of equilibrium steps"
        call json%get(data_name, this%num_equilibrium_steps, found)
        call test_data_found(data_name, found)
        
    end subroutine Physical_System_set_monte_carlo_steps
    
    subroutine Physical_System_set_monte_carlo_changes(this, json)
    
        class(Physical_System), intent(inout) :: this
        type(json_file), intent(inout) :: json
        
        integer :: switch_factor
        real(DP) :: type1_move_delta, type1_move_rejection
        real(DP) :: type1_rotation_delta, type1_rotation_delta_max, type1_rotation_rejection
        real(DP) :: type2_move_delta, type2_move_rejection
        
        character(len=4096) :: data_name
        logical :: found
        
        this%num_particles = this%type1_spheres%get_num_particles() + &
                             this%type2_spheres%get_num_particles()
        data_name = "Monte Carlo.decorrelation factor"
        call json%get(data_name, this%decorrelation_factor, found)
        call test_data_found(data_name, found)
        
        this%num_moves = this%decorrelation_factor * this%num_particles
        
        data_name = "Monte Carlo.switch factor"
        call json%get(data_name, switch_factor, found)
        call test_data_found(data_name, found)
        
        this%num_switches = switch_factor * this%decorrelation_factor * &
                            this%type1_spheres%get_num_particles()
        this%num_rotations = this%decorrelation_factor * this%type1_spheres%get_num_particles()
        this%num_changes = this%num_moves + this%num_switches + this%num_rotations
        
        data_name = "Monte Carlo.Dipoles.move.initial delta"
        call json%get(data_name, type1_move_delta, found)
        call test_data_found(data_name, found)
        
        data_name = "Monte Carlo.Dipoles.move.wanted rejection"
        call json%get(data_name, type1_move_rejection, found)
        call test_data_found(data_name, found)
        
        call this%type1_macro%move%init(type1_move_delta, type1_move_rejection)
        
        data_name = "Monte Carlo.Dipoles.rotation.initial delta"
        call json%get(data_name, type1_rotation_delta, found)
        call test_data_found(data_name, found)
        
        data_name = "Monte Carlo.Dipoles.rotation.maximum delta"
        call json%get(data_name, type1_rotation_delta_max, found)
        call test_data_found(data_name, found)
        
        data_name = "Monte Carlo.Dipoles.rotation.wanted rejection"
        call json%get(data_name, type1_rotation_rejection, found)
        call test_data_found(data_name, found)
        
        call this%type1_macro%rotation%init(type1_rotation_delta, type1_rotation_delta_max, &
                                            type1_rotation_rejection)
                      
        data_name = "Monte Carlo.Hard Spheres.move.initial delta"
        call json%get(data_name, type2_move_delta, found)
        call test_data_found(data_name, found)
        
        data_name = "Monte Carlo.Hard Spheres.move.wanted rejection"
        call json%get(data_name, type2_move_rejection, found)
        call test_data_found(data_name, found)
        
        call this%type2_macro%move%init(type2_move_delta, type2_move_rejection)
        
    end subroutine Physical_System_set_monte_carlo_changes
    
    ! Initialization
    
    subroutine Physical_System_init(this, json, args)   
     
        class(Physical_System), intent(inout) :: this
        type(json_file), intent(inout) :: json
        type(Monte_Carlo_Arguments), intent(in) :: args
        
        real(DP) :: potential_energy_conf      
        
        character(len=4096) :: data_name
        logical :: found  
        
        data_name = "Distribution.take snapshot"
        call json%get(data_name, this%snap, found)
        call test_data_found(data_name, found)
        
        data_name = "Monte Carlo.period of reset"
        call json%get(data_name, this%reset_i_step, found)
        call test_data_found(data_name, found)
        this%reset_i_step = this%reset_i_step / this%decorrelation_factor
        
        write(output_unit, *) "Monte-Carlo Simulation: Canonical ensemble"
        
        call this%open_all_units()
        call this%init_switch()
        
        call init_random_seed(args%random, this%report_unit)
        call set_initial_configuration(this%Box%size, args%initial, this%type1_spheres, &
                                       this%type2_spheres, this%mix%get_diameter(), &
                                       this%report_unit)
                                       
        this%mix_potential_energy_sum = 0._DP        
        call init_mix(this%Box%size, this%mix, this%type1_spheres, this%type2_spheres, &
                      this%write_potential_energy, this%mix_potential_energy_tab_unit, this%mix_potential_energy)
        call this%mix%write_report(this%mix_report_unit)
        
        call init_spheres(this%Box, this%type1_spheres, this%type1_units)
        call init_hard_potential_energy(this%type1_macro%hard_potential_energy, "Dipoles", &
                                 this%type1_spheres%get_diameter(), json)
        call init_cells(this%Box%size, this%type1_spheres, this%type1_macro, this%type2_spheres, &
                        this%mix)
        call set_ewald(this%Box, this%type1_spheres, this%type1_macro, json, this%type1_units)
        this%type1_observables%potential_energy = total_energy(this%Box, this%type1_spheres, this%type1_macro)
                        
        call init_spheres(this%Box, this%type2_spheres, this%type2_units)
        call init_hard_potential_energy(this%type2_macro%hard_potential_energy, "Hard Spheres", &
                                 this%type2_spheres%get_diameter(), json)
        call init_cells(this%Box%size, this%type2_spheres, this%type2_macro, this%type1_spheres, &
                        this%mix)        
        this%type2_observables%potential_energy = total_energy(this%Box, this%type2_spheres, this%type2_macro)
                        
        call this%write_all_reports()
        if (this%write_potential_energy) then
            call this%type1_macro%hard_potential_energy%write(this%type1_units%potential_energy)
            call this%type1_macro%ewald_real%write(this%type1_units%potential_energy_real)
            call this%type2_macro%hard_potential_energy%write(this%type2_units%potential_energy)
        end if
        
        this%potential_energy_sum = 0._DP
        potential_energy_conf = this%type1_observables%potential_energy + this%type2_observables%potential_energy + &
                         this%mix_potential_energy
        write(output_unit, *) "Initial potential_energy energy =", potential_energy_conf
        write(this%observables_thermalisation_unit, *) 0, potential_energy_conf
    
    end subroutine Physical_System_init
    
    subroutine Physical_System_open_all_units(this)    
        class(Physical_System), intent(inout) :: this
        
        call open_units(this%report_unit, this%observables_thermalisation_unit, &
                        this%observables_equilibrium_unit)        
        call this%type1_units%open(this%type1_spheres%get_name())           
        call this%type2_units%open(this%type2_spheres%get_name())       
        call mix_open_units(this%mix_report_unit, this%mix_potential_energy_tab_unit, &
                            this%mix_observables_thermalisation_unit, &
                            this%mix_observables_equilibrium_unit)        
    
    end subroutine Physical_System_open_all_units
    
    subroutine Physical_System_write_all_reports(this)    
        class(Physical_System), intent(in) :: this
        
        call write_data(this%report_unit)
        call this%write_report(this%report_unit)
        
        call write_spheres_density(this%Box%size, this%type1_spheres, this%num_particles, &
                                   this%type1_units%report)
        call this%type1_spheres%write_report(this%type1_units%report)
        
        call write_spheres_density(this%Box%size, this%type2_spheres, this%num_particles, &
                                   this%type2_units%report)
        call this%type2_spheres%write_report(this%type2_units%report)
        
        call this%mix%write_report(this%mix_report_unit)
    
    end subroutine Physical_System_write_all_reports
    
    subroutine Physical_System_write_report(this, report_unit)    
        class(Physical_System), intent(in) :: this
        integer, intent(in) :: report_unit

        write(report_unit, *) "Data macro: "
        
        write(report_unit ,*) "    Box_size(:) = ", this%Box%size(:)
        write(report_unit ,*) "    Volume = ", product(this%Box%size)
        write(report_unit ,*) "    Box_wave(:) = ", this%Box%wave(:)
        write(report_unit ,*) "    NwaveVectors =", NwaveVectors(this%Box%wave)
        write(report_unit ,*) "    Temperature = ", this%Box%temperature
        write(report_unit ,*) "    num_particles = ", this%num_particles        
        
        write(report_unit, *) "    num_equilibrium_steps = ", this%num_equilibrium_steps
        write(report_unit, *) "    num_thermalisation_steps = ", this%num_thermalisation_steps
        write(report_unit, *) "    decorrelation_factor = ", this%decorrelation_factor
        write(report_unit, *) "    num_moves = ", this%num_moves
        write(report_unit, *) "    num_switches = ", this%num_switches
        write(report_unit, *) "    num_rotations = ", this%num_rotations
        
        write(report_unit, *) "    reset_i_step = ", this%reset_i_step
        write(report_unit, *) "    write_potential_energy = ", this%write_potential_energy
    
    end subroutine Physical_System_write_report
    
    pure subroutine Physical_System_init_switch(this)
        class(Physical_System), intent(inout) :: this
        
        this%switch_num_hits = 0
        this%switch_num_rejections = 0
        this%switch_rejection_rate = 0._DP
        this%switch_sum_rejection = 0._DP
        
    end subroutine Physical_System_init_switch
    
    ! Finalisation
    
    subroutine Physical_System_final(this, json)    
    
        class(Physical_System), intent(inout) :: this
        type(json_file), intent(inout) :: json
        
        real(DP) :: type1_energy, type2_energy
        
        call final_spheres(this%Box, this%type1_spheres, this%type1_units)
        call set_ewald(this%Box, this%type1_spheres, this%type1_macro, json, this%type1_units)
        type1_energy = total_energy(this%Box, this%type1_spheres, this%type1_macro)
        call test_consist(this%type1_observables%potential_energy, type1_energy, this%type1_units%report)
        
        call final_spheres(this%Box, this%type2_spheres, this%type2_units)
        type2_energy = total_energy(this%Box, this%type2_spheres, this%type2_macro)
        call test_consist(this%type2_observables%potential_energy, type2_energy, this%type2_units%report)
        
        call mix_final(this%Box%size, this%mix, this%type1_spheres, this%type2_spheres, &
                       this%mix_report_unit, this%mix_potential_energy, this%mix_potential_energy_conf)
        
        call this%write_all_results()
        call this%close_units()
    
    end subroutine Physical_System_final
    
    subroutine Physical_System_write_all_results(this)
        class(Physical_System), intent(inout) :: this
        
        call this%type1_observables%write_results(this%Box%temperature, this%num_equilibrium_steps, this%type1_units%report)
        call this%type2_observables%write_results(this%Box%temperature, this%num_equilibrium_steps, this%type2_units%report)
        call mix_write_results(this%num_equilibrium_steps, this%mix_potential_energy_sum, this%mix_report_unit)
        
        call this%write_results()
    
    end subroutine Physical_System_write_all_results
    
    subroutine Physical_System_write_results(this)    
        class(Physical_System), intent(inout) :: this
        
        real(DP) :: potential_energy, potential_energy_conf
        real(DP) :: duration
        
        potential_energy = this%type1_observables%potential_energy + this%type2_observables%potential_energy + &
                    this%mix_potential_energy
        potential_energy_conf = total_energy(this%Box, this%type1_spheres, this%type1_macro) + &
                    total_energy(this%Box, this%type2_spheres, this%type2_macro) + &
                    this%mix_potential_energy_conf
        call test_consist(potential_energy, potential_energy_conf, this%report_unit)
        this%potential_energy_sum = this%type1_observables%potential_energy_sum + &
                             this%type2_observables%potential_energy_sum + &
                             this%mix_potential_energy_sum
        duration = this%time_end - this%time_start
        call write_results(this%num_particles, this%num_equilibrium_steps, this%potential_energy_sum, &
                           this%switch_sum_rejection, duration, this%report_unit)
    
    end subroutine Physical_System_write_results
    
    subroutine Physical_System_close_units(this)    
        class(Physical_System), intent(inout) :: this
    
        call this%type1_units%close()
        call this%type2_units%close()
        
        close(this%mix_report_unit)
        close(this%mix_potential_energy_tab_unit)
        close(this%mix_observables_thermalisation_unit)
        close(this%mix_observables_equilibrium_unit)
        
        close(this%report_unit)
        close(this%observables_thermalisation_unit)
        close(this%observables_equilibrium_unit)
    
    end subroutine Physical_System_close_units
    
    ! Destruction
    
    subroutine Physical_System_destroy(this)    
        class(Physical_System), intent(inout) :: this
        
        call this%mix%destroy()
        
        call this%type2_macro%mix_cells%destroy()
        call this%type2_macro%same_cells%destroy()        
        call this%type2_spheres%destroy()
        
        call this%type1_macro%ewald_reci%destroy()
        call this%type1_macro%ewald_real%destroy()
        call this%type1_macro%mix_cells%destroy()
        call this%type1_macro%same_cells%destroy()        
        call this%type1_spheres%destroy()
        
        write(output_unit, *) this%name, " class destruction"
    
    end subroutine Physical_System_destroy
    
    ! Accessors
    
    pure function Physical_System_get_num_thermalisation_steps(this) &
         result(get_num_thermalisation_steps)
        class(Physical_System), intent(in) :: this
        integer :: get_num_thermalisation_steps
                
        get_num_thermalisation_steps = this%num_thermalisation_steps
        
    end function Physical_System_get_num_thermalisation_steps
    
    pure function Physical_System_get_num_equilibrium_steps(this) result(get_num_equilibrium_steps)
        class(Physical_System), intent(in) :: this
        integer :: get_num_equilibrium_steps
                
        get_num_equilibrium_steps = this%num_equilibrium_steps
        
    end function Physical_System_get_num_equilibrium_steps
    
    ! Mutators
    
    subroutine Physical_System_set_time_start(this)
        class(Physical_System), intent(inout) :: this  
              
        call cpu_time(this%time_start)
        
    end subroutine Physical_System_set_time_start
    
    subroutine Physical_System_set_time_end(this)
        class(Physical_System), intent(inout) :: this  
              
        call cpu_time(this%time_end)
        
    end subroutine Physical_System_set_time_end
    
    ! Random changes
    
    subroutine Physical_System_random_changes(this)    
        class(Physical_System), intent(inout) :: this
        
        integer :: iChange, iChangeRand, iColRand
        real(DP) :: rand
        
        MC_Change: do iChange = 1, this%num_changes
        
            ! Randomly choosing the change
            call random_number(rand)
            iChangeRand = int(rand*real(this%num_changes, DP)) + 1
            
            if (iChangeRand <= this%num_moves) then
                ! Randomly choosing the type
                call random_number(rand)
                iColRand = int(rand*real(this%num_particles, DP)) + 1
                if (iColRand <= this%type1_spheres%get_num_particles()) then                    
                    call move(this%Box, &
                              this%type1_spheres, this%type1_macro, this%type1_observables, &
                              this%type2_spheres, this%type2_macro%mix_cells, &
                              this%mix, this%mix_potential_energy)
                else
                    call move(this%Box, &
                              this%type2_spheres, this%type2_macro, this%type2_observables, &
                              this%type1_spheres, this%type1_macro%mix_cells, &
                              this%mix, this%mix_potential_energy)
                end if
            else if (iChangeRand <= this%num_moves + this%num_switches) then
                call switch(this%Box, &
                            this%type1_spheres, this%type1_macro, this%type1_observables, &
                            this%type2_spheres, this%type2_macro, this%type2_observables, &
                            this%mix, this%mix_potential_energy, &
                            this%switch_num_rejections)
                this%switch_num_hits = this%switch_num_hits + 1
            else
                call rotate(this%Box, &
                            this%type1_spheres, this%type1_macro, this%type1_observables)
            end if
            
        end do MC_Change
    
    end subroutine Physical_System_random_changes
    
    subroutine Physical_System_update_rejections(this)    
        class(Physical_System), intent(inout) :: this
    
        call this%type1_observables%update_rejections()
        call this%type2_observables%update_rejections()
        
        this%switch_rejection_rate = real(this%switch_num_rejections, DP)/real(this%switch_num_hits, DP)
        this%switch_num_rejections = 0
        this%switch_num_hits = 0
        
    end subroutine Physical_System_update_rejections
    
    subroutine Physical_System_adapt_changes(this, i_step)    
        class(Physical_System), intent(inout) :: this
        integer, intent(in) :: i_step
        
        if (mod(i_step, this%period_adaptation) /= 0) then ! Rejections accumulation
            this%type1_observables%move_rejection_adapt = this%type1_observables%move_rejection_adapt + &
                                                  this%type1_observables%move_rejection_rate
            this%type1_observables%rotate_rejection_adapt = this%type1_observables%rotate_rejection_adapt + &
                                                    this%type1_observables%rotate_rejection_rate
            this%type2_observables%move_rejection_adapt = this%type2_observables%move_rejection_adapt + &
                                                  this%type2_observables%move_rejection_rate
        else ! Average & adaptation
            call adapt_move(this%Box%size, &
                            this%type1_macro%move, &
                            this%period_adaptation, i_step, &
                            this%type1_observables, &
                            this%type1_units%move_delta)
            call adapt_rotation(this%type1_macro%rotation, &
                                this%period_adaptation, i_step, &
                                this%type1_observables, &
                                this%type1_units%rotate_delta)
            call adapt_move(this%Box%size, &
                            this%type2_macro%move, &
                            this%period_adaptation, i_step, &
                            this%type2_observables, &
                            this%type2_units%move_delta)
        end if
        
    end subroutine Physical_System_adapt_changes
    
    subroutine Physical_System_write_observables_thermalisation(this, i_step)    
        class(Physical_System), intent(inout) :: this
        integer, intent(in) :: i_step
    
        call this%type1_observables%write(i_step, this%type1_units%observables_thermalisation)
        call this%type2_observables%write(i_step, this%type2_units%observables_thermalisation)
        
        write(this%mix_observables_thermalisation_unit, *) i_step, this%mix_potential_energy
        write(this%observables_thermalisation_unit, *) i_step, this%type1_observables%potential_energy + this%type2_observables%potential_energy + this%mix_potential_energy
            
    end subroutine Physical_System_write_observables_thermalisation
    
    subroutine Physical_System_fix_changes(this)
    
        class(Physical_System), intent(inout) :: this
        
        character(len=5) :: type_name
        
        type_name = this%type1_spheres%get_name()
        call this%type1_macro%move%set_delta(type_name, this%Box%size, &
                                             this%type1_observables%move_rejection_average, &
                                             this%type1_units%report)
        call this%type1_macro%rotation%set_delta(type_name, &
                                                 this%type1_observables%rotate_rejection_average, &
                                                 this%type1_units%report)
        type_name = this%type2_spheres%get_name()
        call this%type2_macro%move%set_delta(type_name, this%Box%size, &
                                             this%type2_observables%move_rejection_average, &
                                             this%type2_units%report)
    
    end subroutine Physical_System_fix_changes
    
    subroutine Physical_System_measure_chemical_potentials(this)
    
        class(Physical_System), intent(inout) :: this
    
        call widom(this%Box, &
                   this%type1_spheres, this%type1_macro, this%type1_observables, &
                   this%type2_spheres, this%type2_macro%mix_cells, &
                   this%mix)
        call widom(this%Box, &
                   this%type2_spheres, this%type2_macro, this%type2_observables, &
                   this%type1_spheres, this%type1_macro%mix_cells, &
                   this%mix)
    
    end subroutine Physical_System_measure_chemical_potentials
    
    subroutine Physical_System_accumulate_observables(this)
    
        class(Physical_System), intent(inout) :: this
    
        call this%type1_observables%accumulate()
        call this%type2_observables%accumulate()
        
        this%mix_potential_energy_sum = this%mix_potential_energy_sum + this%mix_potential_energy
        this%switch_sum_rejection = this%switch_sum_rejection + this%switch_rejection_rate
            
    end subroutine Physical_System_accumulate_observables
    
    subroutine Physical_System_write_observables_equilibrium(this, i_step)
    
        class(Physical_System), intent(inout) :: this
        integer, intent(in) :: i_step
    
        call this%type1_observables%write(i_step, this%type1_units%observables_equilibrium)
        call this%type2_observables%write(i_step, this%type2_units%observables_equilibrium)
        
        write(this%mix_observables_equilibrium_unit, *) i_step, this%mix_potential_energy
        write(this%observables_equilibrium_unit, *) i_step, this%type1_observables%potential_energy + &
                                                            this%type2_observables%potential_energy + &
                                                            this%mix_potential_energy
            
    end subroutine Physical_System_write_observables_equilibrium
    
    subroutine Physical_System_take_snapshots(this, i_step)
    
        class(Physical_System), intent(inout) :: this
        integer, intent(in) :: i_step
        
        if (this%snap) then ! Snap shots of the configuration
            call this%type1_spheres%write_snap_positions(i_step, this%type1_units%snap_positions)
            call this%type1_spheres%write_snap_orientations(i_step, &
                                                            this%type1_units%snap_orientations)
            call this%type2_spheres%write_snap_positions(i_step, this%type2_units%snap_positions)
        end if
        
    end subroutine Physical_System_take_snapshots
    
    !> Reinitialize summed quantities to prevent them from drifting.
    subroutine Physical_System_reinitialize_quantites(this, i_step)
    
        class(Physical_System), intent(inout) :: this
        integer, intent(in) :: i_step
    
        if (modulo(i_step, this%reset_i_step) == 0) then
            call this%type1_macro%ewald_reci%reset_structure(this%Box, this%type1_spheres, i_step, &
                                                             this%type1_units%structure_modulus)
            call this%type1_macro%ewald_bound%reset_total_moment(this%type1_spheres, i_step, &
                 this%type1_units%totalMoment_modulus)
        end if
        
    end subroutine Physical_System_reinitialize_quantites

end module class_physical_system
