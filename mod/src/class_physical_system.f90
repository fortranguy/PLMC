!> \brief Description of the Physical System class

module class_physical_system

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file, json_initialize, json_destroy, &
                       json_value, json_value_create, to_object, json_value_add, &
                       json_print
use module_data, only: test_data_file_exists, test_data_found, test_empty_string
use module_geometry, only: set_geometry
use module_types_micro, only: Box_Parameters, Monte_Carlo_Arguments
use module_physics_micro, only: num_wave_vectors
use class_hard_spheres, only: Hard_Spheres, Dipolar_Hard_Spheres, Between_Hard_Spheres
use class_hard_spheres_potential, only: Between_Hard_Spheres_Potential_Energy
use module_types_macro, only: Hard_Spheres_Macro, Dipolar_Hard_Spheres_Macro
use class_discrete_observable, only: Discrete_Observables
use class_hard_spheres_observables, only: Hard_Spheres_Observables, Dipolar_Hard_Spheres_Observables, &
                                          Between_Hard_Spheres_Observables
use class_hard_spheres_units, only: Hard_Spheres_Units, Dipolar_Hard_Spheres_Units, &
                                    Between_Hard_Spheres_Units
use module_physics_macro, only: init_random_seed, set_initial_configuration, &
                                init_spheres, init_cells, set_ewald, total_energy, final_spheres, &
                                init_between_spheres_potential, final_between_spheres_potential, &
                                test_consistency
use module_algorithms, only: move, widom, switch, rotate
use module_write, only: write_results, between_spheres_write_results

implicit none

private

    character(len=*), parameter :: data_filename = "data.json"
    character(len=*), parameter :: report_filename = "report.json"

    type, public :: Physical_System
    
        private
        
        type(json_file) :: data_json
        
        character(len=:), allocatable :: name
        
        ! Box
        type(Box_Parameters) :: Box
        
        ! Monte-Carlo
        integer :: num_thermalisation_steps, period_adaptation, num_equilibrium_steps
        integer :: decorrelation_factor, num_changes, num_moves, num_switches, num_rotations
        
        ! Type 1: Dipolar spheres
        type(Dipolar_Hard_Spheres) :: type1_spheres
        type(Dipolar_Hard_Spheres_Macro) :: type1_macro
        type(Dipolar_Hard_Spheres_Observables) :: type1_observables
        type(Dipolar_Hard_Spheres_Units) :: type1_units
        type(json_value), pointer :: type1_report_json
        
        ! Type 2: Hard spheres
        type(Hard_Spheres) :: type2_spheres
        type(Hard_Spheres_Macro) :: type2_macro
        type(Hard_Spheres_Observables) :: type2_observables
        type(Hard_Spheres_Units) :: type2_units
        type(json_value), pointer :: type2_report_json
        
        ! Between Spheres potential_energy
        type(Between_Hard_Spheres) :: between_spheres
        type(Between_Hard_Spheres_Potential_Energy) :: between_spheres_potential
        type(Between_Hard_Spheres_Observables) :: between_spheres_observables
        type(Between_Hard_Spheres_Units) :: between_spheres_units
        type(json_value), pointer :: between_spheres_report_json
        
        ! Observables and files units
        real(DP) :: potential_energy, potential_energy_sum
        type(json_value), pointer :: report_json, system_json
        integer :: report_unit
        integer :: observables_thermalisation_unit, observables_equilibrium_unit
        
        type(Discrete_Observables) :: switch_observable
        
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
        procedure, private :: open_units => Physical_System_open_units
        procedure, private :: json_create_all_values => Physical_System_json_create_all_values
        procedure, private :: write_all_reports => Physical_System_write_all_reports
        procedure, private :: write_report => Physical_System_write_report
        procedure :: final => Physical_System_final
        procedure, private :: write_all_results => Physical_System_write_all_results
        procedure, private :: write_results => Physical_System_write_results
        procedure, private :: json_destroy_all_values => Physical_System_json_destroy_all_values
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
    
    subroutine Physical_System_construct(this, args)
            
        class(Physical_System), intent(out) :: this
        type(Monte_Carlo_Arguments), intent(in) :: args
        
        character(len=4096) :: data_name
        logical :: found
        character(len=:), allocatable :: this_name

        call test_data_file_exists(data_filename)

        call json_initialize()
        call this%data_json%load_file(filename = data_filename)
        
        call set_geometry(args%geometry)
        
        data_name = "Box.name"
        call this%data_json%get(data_name, this_name, found)
        call test_data_found(data_name, found)
        call test_empty_string(data_name, this_name)
        this%name = this_name
        if (allocated(this_name)) deallocate(this_name)
        write(output_unit, *) this%name, " class construction"
        
        call this%set_box()
        call this%set_monte_carlo_steps()
        
        data_name = "Potential Energy.write"
        call this%data_json%get(data_name, this%write_potential_energy, found)
        call test_data_found(data_name, found)
        
        call this%type1_spheres%construct(this%data_json)
        call this%type2_spheres%construct(this%data_json)
        call this%between_spheres%construct(this%data_json, &
                                            this%type1_spheres%get_diameter(), &
                                            this%type2_spheres%get_diameter())
        
        call this%set_monte_carlo_changes()
        
    end subroutine Physical_System_construct
    
    subroutine Physical_System_set_box(this)
    
        class(Physical_System), intent(inout) :: this
        
        character(len=4096) :: data_name
        logical :: found
        real(DP), dimension(:), allocatable :: Box_size
        integer, dimension(:), allocatable :: Box_wave
        
        data_name = "Box.size"
        call this%data_json%get(data_name, Box_size, found)
        call test_data_found(data_name, found)
        if (size(Box_size) /= size (this%Box%size)) error stop "Box size dimension"
        this%Box%size(:) = Box_size(:)
        if (allocated(Box_size)) deallocate(Box_size)
        
        data_name = "Box.wave"
        call this%data_json%get(data_name, Box_wave, found)
        call test_data_found(data_name, found)
        if (size(Box_wave) /= size (this%Box%wave)) error stop "Box wave dimension"
        this%Box%wave(:) = Box_wave(:)
        if (allocated(Box_wave)) deallocate(Box_wave)
        
        data_name = "Box.temperature"
        call this%data_json%get(data_name, this%Box%temperature, found)
        call test_data_found(data_name, found)
        
    end subroutine Physical_System_set_box
    
    subroutine Physical_System_set_monte_carlo_steps(this)
    
        class(Physical_System), intent(inout) :: this
        
        character(len=4096) :: data_name
        logical :: found
        
        data_name = "Monte Carlo.number of thermalisation steps"
        call this%data_json%get(data_name, this%num_thermalisation_steps, found)
        call test_data_found(data_name, found)
        
        data_name = "Monte Carlo.period of adaptation"
        call this%data_json%get(data_name, this%period_adaptation, found)
        call test_data_found(data_name, found)
        
        data_name = "Monte Carlo.number of equilibrium steps"
        call this%data_json%get(data_name, this%num_equilibrium_steps, found)
        call test_data_found(data_name, found)
        
    end subroutine Physical_System_set_monte_carlo_steps
    
    subroutine Physical_System_set_monte_carlo_changes(this)
    
        class(Physical_System), intent(inout) :: this
        
        integer :: switch_factor
        real(DP) :: type1_move_delta, type1_move_rejection
        real(DP) :: type1_rotation_delta, type1_rotation_delta_max, type1_rotation_rejection
        real(DP) :: type2_move_delta, type2_move_rejection
        
        character(len=4096) :: data_name
        logical :: found
        
        this%Box%num_particles = this%type1_spheres%get_num_particles() + &
                                 this%type2_spheres%get_num_particles()
        data_name = "Monte Carlo.decorrelation factor"
        call this%data_json%get(data_name, this%decorrelation_factor, found)
        call test_data_found(data_name, found)
        
        this%num_moves = this%decorrelation_factor * this%Box%num_particles
        
        data_name = "Monte Carlo.switch factor"
        call this%data_json%get(data_name, switch_factor, found)
        call test_data_found(data_name, found)
        
        this%num_switches = switch_factor * this%decorrelation_factor * &
                            this%type1_spheres%get_num_particles()
        this%num_rotations = this%decorrelation_factor * this%type1_spheres%get_num_particles()
        this%num_changes = this%num_moves + this%num_switches + this%num_rotations
        
        data_name = "Monte Carlo.Dipolar Hard Spheres.move.initial delta"
        call this%data_json%get(data_name, type1_move_delta, found)
        call test_data_found(data_name, found)
        
        data_name = "Monte Carlo.Dipolar Hard Spheres.move.wanted rejection"
        call this%data_json%get(data_name, type1_move_rejection, found)
        call test_data_found(data_name, found)
        
        call this%type1_macro%move%init(type1_move_delta, type1_move_rejection)
        
        data_name = "Monte Carlo.Dipolar Hard Spheres.rotation.initial delta"
        call this%data_json%get(data_name, type1_rotation_delta, found)
        call test_data_found(data_name, found)
        
        data_name = "Monte Carlo.Dipolar Hard Spheres.rotation.maximum delta"
        call this%data_json%get(data_name, type1_rotation_delta_max, found)
        call test_data_found(data_name, found)
        
        data_name = "Monte Carlo.Dipolar Hard Spheres.rotation.wanted rejection"
        call this%data_json%get(data_name, type1_rotation_rejection, found)
        call test_data_found(data_name, found)
        
        call this%type1_macro%rotation%init(type1_rotation_delta, type1_rotation_delta_max, &
                                            type1_rotation_rejection)
                      
        data_name = "Monte Carlo.Hard Spheres.move.initial delta"
        call this%data_json%get(data_name, type2_move_delta, found)
        call test_data_found(data_name, found)
        
        data_name = "Monte Carlo.Hard Spheres.move.wanted rejection"
        call this%data_json%get(data_name, type2_move_rejection, found)
        call test_data_found(data_name, found)
        
        call this%type2_macro%move%init(type2_move_delta, type2_move_rejection)
        
    end subroutine Physical_System_set_monte_carlo_changes
    
    ! Initialization
    
    subroutine Physical_System_init(this, args)
     
        class(Physical_System), intent(inout) :: this
        type(Monte_Carlo_Arguments), intent(in) :: args
        
        real(DP) :: potential_energy_conf
        
        character(len=4096) :: data_name
        logical :: found
        
        data_name = "Distribution.take snapshot"
        call this%data_json%get(data_name, this%snap, found)
        call test_data_found(data_name, found)
        
        data_name = "Monte Carlo.period of reset"
        call this%data_json%get(data_name, this%reset_i_step, found)
        call test_data_found(data_name, found)
        this%reset_i_step = this%reset_i_step / this%decorrelation_factor
        
        write(output_unit, *) "Monte-Carlo Simulation: Canonical ensemble "
        
        call this%open_all_units()
        call this%json_create_all_values()
        
        call init_random_seed(args%random, this%report_json)
        call set_initial_configuration(this%Box%size, args%initial, &
                                       this%type1_spheres, this%type2_spheres, &
                                       this%between_spheres%get_diameter(), &
                                       this%report_json)
                                       
        call this%between_spheres%test_overlap(this%Box%size, &
                                               this%type1_spheres, this%type2_spheres)
        call this%between_spheres_potential%construct(this%data_json, "Between Spheres", &
                                                      this%between_spheres%get_diameter())
        call init_between_spheres_potential(this%Box%size, this%between_spheres_potential, &
                                            this%type1_spheres, this%type2_spheres, &
                                            this%write_potential_energy, &
                                            this%between_spheres_observables%potential_energy, &
                                            this%between_spheres_units%potential_energy_tabulation)
        
        call init_spheres(this%Box, this%type1_spheres, this%type1_units)
        call this%type1_macro%hard_potential%construct(this%data_json, "Dipolar Hard Spheres", &
                                                       this%type1_spheres%get_diameter())
        call init_cells(this%Box%size, this%type1_spheres, this%type1_macro, this%type2_spheres, &
                        this%between_spheres_potential)
        call set_ewald(this%Box, this%type1_spheres, this%type1_macro, this%data_json, this%type1_units)
        this%type1_observables%potential_energy = total_energy(this%Box, this%type1_spheres, &
                                                               this%type1_macro)
                        
        call init_spheres(this%Box, this%type2_spheres, this%type2_units)
        call this%type2_macro%hard_potential%construct(this%data_json, "Hard Spheres", &
                                                       this%type2_spheres%get_diameter())
        call init_cells(this%Box%size, this%type2_spheres, this%type2_macro, this%type1_spheres, &
                        this%between_spheres_potential)
        this%type2_observables%potential_energy = total_energy(this%Box, this%type2_spheres, &
                                                               this%type2_macro)
                        
        call this%write_all_reports(args%geometry)
        if (this%write_potential_energy) then
            call this%type1_macro%hard_potential%write(this%type1_units%potential_energy)
            call this%type1_macro%ewald_real%write(this%type1_units%potential_energy_real)
            call this%type2_macro%hard_potential%write(this%type2_units%potential_energy)
        end if
        
        this%potential_energy_sum = 0._DP
        potential_energy_conf = this%type1_observables%potential_energy + &
                                this%type2_observables%potential_energy + &
                                this%between_spheres_observables%potential_energy
        write(this%observables_thermalisation_unit, *) 0, potential_energy_conf

        rewind(this%report_unit)
        call json_print(this%report_json, this%report_unit)
    
    end subroutine Physical_System_init
    
    subroutine Physical_System_open_all_units(this)
    
        class(Physical_System), intent(inout) :: this

        call this%open_units()
        call this%type1_units%open(this%type1_spheres%get_name())
        call this%type2_units%open(this%type2_spheres%get_name())
        call this%between_spheres_units%open(this%between_spheres%get_name())
    
    end subroutine Physical_System_open_all_units

    subroutine Physical_System_open_units(this)

        class(Physical_System), intent(inout) :: this

        open(newunit=this%report_unit, recl=4096, file=report_filename, &
             status='new', action='write')
        open(newunit=this%observables_thermalisation_unit, recl=4096, &
             file="observables_thermalisation.out", status='new', action='write')
        open(newunit=this%observables_equilibrium_unit, recl=4096, &
             file="observables_equilibrium.out", status='new', action='write')
        write(this%observables_equilibrium_unit, *) "#", 1 ! 1 observable: energy

    end subroutine Physical_System_open_units

    subroutine Physical_System_json_create_all_values(this)

        class(Physical_System), intent(inout) :: this

        call json_value_create(this%report_json)
        call to_object(this%report_json)

        call json_value_create(this%system_json)
        call to_object(this%system_json, "System")
        call json_value_add(this%report_json, this%system_json)

        call json_value_create(this%between_spheres_report_json)
        call to_object(this%between_spheres_report_json, "Between Spheres")
        call json_value_add(this%report_json, this%between_spheres_report_json)

        call json_value_create(this%type1_report_json)
        call to_object(this%type1_report_json, "Dipolar Hard Spheres")
        call json_value_add(this%report_json, this%type1_report_json)

        call json_value_create(this%type2_report_json)
        call to_object(this%type2_report_json, "Hard Spheres")
        call json_value_add(this%report_json, this%type2_report_json)
        
    end subroutine Physical_System_json_create_all_values
    
    subroutine Physical_System_write_all_reports(this, geometry)
    
        class(Physical_System), intent(in) :: this
        character(len=*), intent(in) :: geometry
        
        call this%write_report(geometry)
        
        call this%type1_spheres%write_report(this%Box, this%type1_report_json)
        call this%type1_macro%hard_potential%write_report(this%type1_report_json)
        
        call this%type2_spheres%write_report(this%Box, this%type2_report_json)
        call this%type2_macro%hard_potential%write_report(this%type2_report_json)
        
        call this%between_spheres_potential%write_report(this%between_spheres_report_json)
    
    end subroutine Physical_System_write_all_reports
    
    subroutine Physical_System_write_report(this, geometry)
    
        class(Physical_System), intent(in) :: this
        character(len=*), intent(in) :: geometry

        type(json_value), pointer :: box_json, changes_json

        call json_value_create(box_json)
        call to_object(box_json, "Box")
        call json_value_add(this%system_json, box_json)

        call json_value_add(box_json, "volume", product(this%Box%size))
        call json_value_add(box_json, "geometry", geometry)
        call json_value_add(box_json, "number of wave vectors", &
                                       num_wave_vectors(this%Box%wave))
        call json_value_add(box_json, "number of particles", this%Box%num_particles)

        nullify(box_json)

        call json_value_create(changes_json)
        call to_object(changes_json, "Changes")
        call json_value_add(this%system_json, changes_json)

        call json_value_add(changes_json, "number of moves", this%num_moves)
        call json_value_add(changes_json, "number of switches", this%num_switches)
        call json_value_add(changes_json, "number of rotations", this%num_rotations)

        nullify(changes_json)
        
    end subroutine Physical_System_write_report
    
    ! Finalisation
    
    subroutine Physical_System_final(this)
    
        class(Physical_System), intent(inout) :: this
        
        real(DP) :: type1_energy, type2_energy
        
        call final_spheres(this%Box, this%type1_spheres, this%type1_units)
        call set_ewald(this%Box, this%type1_spheres, this%type1_macro, this%data_json, this%type1_units)
        type1_energy = total_energy(this%Box, this%type1_spheres, this%type1_macro)
        call test_consistency(this%type1_observables%potential_energy, type1_energy, &
                              this%type1_report_json)
        
        call final_spheres(this%Box, this%type2_spheres, this%type2_units)
        type2_energy = total_energy(this%Box, this%type2_spheres, this%type2_macro)
        call test_consistency(this%type2_observables%potential_energy, type2_energy, &
                              this%type2_report_json)
        
        call this%between_spheres%test_overlap(this%Box%size, &
                                               this%type1_spheres, this%type2_spheres)
        call final_between_spheres_potential(this%Box%size, this%between_spheres_potential, &
                                             this%type1_spheres, this%type2_spheres, &
                                             this%between_spheres_observables%potential_energy, &
                                             this%between_spheres_report_json)
        
        call this%write_all_results()
        call this%json_destroy_all_values()
        call this%close_units()
    
    end subroutine Physical_System_final
    
    subroutine Physical_System_write_all_results(this)
        class(Physical_System), intent(inout) :: this
        
        call this%type1_observables%write_results(this%Box%temperature, this%num_equilibrium_steps, &
                                                  this%type1_report_json)
        call this%type2_observables%write_results(this%Box%temperature, this%num_equilibrium_steps, &
                                                  this%type2_report_json)
        call between_spheres_write_results(this%num_equilibrium_steps, &
                                           this%between_spheres_observables%potential_energy_sum, &
                                           this%between_spheres_report_json)
        
        call this%write_results()
    
    end subroutine Physical_System_write_all_results
    
    subroutine Physical_System_write_results(this)
        class(Physical_System), intent(inout) :: this
        
        real(DP) :: potential_energy, potential_energy_conf
        real(DP) :: duration
        
        potential_energy = this%type1_observables%potential_energy + &
                           this%type2_observables%potential_energy + &
                           this%between_spheres_observables%potential_energy
        potential_energy_conf = total_energy(this%Box, this%type1_spheres, this%type1_macro) + &
                                total_energy(this%Box, this%type2_spheres, this%type2_macro) + &
                                this%between_spheres_potential%total(this%Box%size, &
                                                                     this%type1_spheres, &
                                                                     this%type2_spheres)
        
        call test_consistency(potential_energy, potential_energy_conf, this%report_json)        
        this%potential_energy_sum = this%type1_observables%potential_energy_sum + &
                                    this%type2_observables%potential_energy_sum + &
                                    this%between_spheres_observables%potential_energy_sum
        duration = this%time_end - this%time_start
        call write_results(this%Box%num_particles, this%num_equilibrium_steps, &
                           this%potential_energy_sum, this%switch_observable%sum_rejection, duration, &
                           this%report_json)

        rewind(this%report_unit)
        call json_print(this%report_json, this%report_unit)
    
    end subroutine Physical_System_write_results

    subroutine Physical_System_json_destroy_all_values(this)

        class(Physical_System), intent(inout) :: this

        nullify(this%type2_report_json)
        nullify(this%type1_report_json)
        nullify(this%between_spheres_report_json)
        nullify(this%system_json)
        
        call json_destroy(this%report_json)
    
    end subroutine Physical_System_json_destroy_all_values
    
    subroutine Physical_System_close_units(this)
    
        class(Physical_System), intent(inout) :: this
    
        call this%type1_units%close()
        call this%type2_units%close()
        call this%between_spheres_units%close()
        
        close(this%report_unit)
        close(this%observables_thermalisation_unit)
        close(this%observables_equilibrium_unit)
    
    end subroutine Physical_System_close_units
    
    ! Destruction
    
    subroutine Physical_System_destroy(this)
     
        class(Physical_System), intent(inout) :: this
        
        write(output_unit, *) this%name, " class destruction"
        
        call this%between_spheres%destroy()
        
        call this%type2_macro%between_cells%destroy()
        call this%type2_macro%same_cells%destroy()
        call this%type2_spheres%destroy()
        
        call this%type1_macro%ewald_reci%destroy()
        call this%type1_macro%ewald_real%destroy()
        call this%type1_macro%between_cells%destroy()
        call this%type1_macro%same_cells%destroy()
        call this%type1_spheres%destroy()
        
        if (allocated(this%name)) deallocate(this%name)

        call this%data_json%destroy()
    
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
        
        integer :: i_change, i_change_rand, i_particule_rand
        real(DP) :: rand
        
        MC_Change: do i_change = 1, this%num_changes
        
            ! Randomly choosing the change
            call random_number(rand)
            i_change_rand = int(rand*real(this%num_changes, DP)) + 1
            
            if (i_change_rand <= this%num_moves) then
                ! Randomly choosing the type
                call random_number(rand)
                i_particule_rand = int(rand*real(this%Box%num_particles, DP)) + 1
                if (i_particule_rand <= this%type1_spheres%get_num_particles()) then
                    call move(this%Box, &
                              this%type1_spheres, this%type1_macro, this%type1_observables, &
                              this%type2_spheres, this%type2_macro%between_cells, &
                              this%between_spheres_potential, &
                              this%between_spheres_observables%potential_energy)
                else
                    call move(this%Box, &
                              this%type2_spheres, this%type2_macro, this%type2_observables, &
                              this%type1_spheres, this%type1_macro%between_cells, &
                              this%between_spheres_potential, &
                              this%between_spheres_observables%potential_energy)
                end if
            else if (i_change_rand <= this%num_moves + this%num_switches) then
                call switch(this%Box, &
                            this%type1_spheres, this%type1_macro, this%type1_observables, &
                            this%type2_spheres, this%type2_macro, this%type2_observables, &
                            this%between_spheres_potential, &
                            this%between_spheres_observables%potential_energy, &
                            this%switch_observable)
            else
                call rotate(this%Box, &
                            this%type1_spheres, this%type1_macro, this%type1_observables)
            end if
            
        end do MC_Change
    
    end subroutine Physical_System_random_changes
    
    subroutine Physical_System_update_rejections(this)
        class(Physical_System), intent(inout) :: this
    
        call this%type1_observables%move%update_rejection()
        call this%type1_observables%rotation%update_rejection()
        call this%type2_observables%move%update_rejection()
        call this%switch_observable%update_rejection()
        
    end subroutine Physical_System_update_rejections
    
    subroutine Physical_System_adapt_changes(this, i_step)
        class(Physical_System), intent(inout) :: this
        integer, intent(in) :: i_step
        
        if (mod(i_step, this%period_adaptation) /= 0) then
            call this%type1_observables%move%accumulate_rejection()
            call this%type1_observables%rotation%accumulate_rejection()
            call this%type2_observables%move%accumulate_rejection()
        else
            call this%type1_observables%move%average_rejection(this%period_adaptation)
            call this%type1_macro%move%adapt_delta(this%Box%size, &
                this%type1_observables%move%rejection_average)
            write(this%type1_units%move_delta, *) i_step, this%type1_macro%move%get_delta_scalar(), &
                                                  this%type1_observables%move%rejection_average
                            
            call this%type1_observables%rotation%average_rejection(this%period_adaptation)
            call this%type1_macro%rotation%adapt_delta(&
                 this%type1_observables%rotation%rejection_average)
            write(this%type1_units%rotate_delta, *) i_step, this%type1_macro%rotation%get_delta(), &
                                                    this%type1_observables%rotation%rejection_average
                                
            call this%type2_observables%move%average_rejection(this%period_adaptation)
            call this%type2_macro%move%adapt_delta(this%Box%size, &
                this%type2_observables%move%rejection_average)
            write(this%type2_units%move_delta, *) i_step, this%type2_macro%move%get_delta_scalar(), &
                                                  this%type2_observables%move%rejection_average
        end if
        
    end subroutine Physical_System_adapt_changes
    
    subroutine Physical_System_write_observables_thermalisation(this, i_step)
        class(Physical_System), intent(inout) :: this
        integer, intent(in) :: i_step
    
        call this%type1_observables%write(i_step, this%type1_units%observables_thermalisation)
        call this%type2_observables%write(i_step, this%type2_units%observables_thermalisation)
        call this%between_spheres_observables%write(i_step, &
            this%between_spheres_units%observables_thermalisation)
        
        write(this%observables_thermalisation_unit, *) i_step, &
            this%type1_observables%potential_energy + &
            this%type2_observables%potential_energy + &
            this%between_spheres_observables%potential_energy
            
    end subroutine Physical_System_write_observables_thermalisation
    
    subroutine Physical_System_fix_changes(this)
    
        class(Physical_System), intent(inout) :: this
        
        call this%type1_macro%move%set_delta(this%Box%size, &
                                             this%type1_spheres%get_name(), &
                                             this%type1_observables%move%rejection_average, &
                                             this%type1_report_json)
        call this%type1_macro%rotation%set_delta(this%type1_spheres%get_name(), &
                                                 this%type1_observables%rotation%rejection_average, &
                                                 this%type1_report_json)
        call this%type2_macro%move%set_delta(this%Box%size, &
                                             this%type2_spheres%get_name(), &
                                             this%type2_observables%move%rejection_average, &
                                             this%type2_report_json)
    
    end subroutine Physical_System_fix_changes
    
    subroutine Physical_System_measure_chemical_potentials(this)
    
        class(Physical_System), intent(inout) :: this
    
        call widom(this%Box, &
                   this%type1_spheres, this%type1_macro, this%type1_observables, &
                   this%type2_spheres, this%type2_macro%between_cells, &
                   this%between_spheres_potential)
        call widom(this%Box, &
                   this%type2_spheres, this%type2_macro, this%type2_observables, &
                   this%type1_spheres, this%type1_macro%between_cells, &
                   this%between_spheres_potential)
    
    end subroutine Physical_System_measure_chemical_potentials
    
    subroutine Physical_System_accumulate_observables(this)
    
        class(Physical_System), intent(inout) :: this
    
        call this%type1_observables%accumulate()
        call this%type2_observables%accumulate()
        call this%between_spheres_observables%accumulate()        
        this%switch_observable%sum_rejection = this%switch_observable%sum_rejection + &
                                               this%switch_observable%rejection_rate
            
    end subroutine Physical_System_accumulate_observables
    
    subroutine Physical_System_write_observables_equilibrium(this, i_step)
    
        class(Physical_System), intent(inout) :: this
        integer, intent(in) :: i_step
    
        call this%type1_observables%write(i_step, this%type1_units%observables_equilibrium)
        call this%type2_observables%write(i_step, this%type2_units%observables_equilibrium)
        call this%between_spheres_observables%write(i_step, &
                                                    this%between_spheres_units%observables_equilibrium)
        
        write(this%observables_equilibrium_unit, *) i_step, &
            this%type1_observables%potential_energy + &
            this%type2_observables%potential_energy + &
            this%between_spheres_observables%potential_energy
            
    end subroutine Physical_System_write_observables_equilibrium
    
    subroutine Physical_System_take_snapshots(this, i_step)
    
        class(Physical_System), intent(inout) :: this
        integer, intent(in) :: i_step
        
        if (this%snap) then ! Snap shots of the configuration
            call this%type1_spheres%write_snap_positions(i_step, &
                 this%type1_units%snap_equilibrium_positions)
            call this%type1_spheres%write_snap_orientations(i_step, &
                 this%type1_units%snap_equilibrium_orientations)
            call this%type2_spheres%write_snap_positions(i_step, &
                 this%type2_units%snap_equilibrium_positions)
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
                 this%type1_units%total_moment_modulus)
        end if
        
    end subroutine Physical_System_reinitialize_quantites

end module class_physical_system
