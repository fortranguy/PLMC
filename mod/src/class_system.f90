!> \brief Description of the Physical System class

module class_system

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use data_box, only: num_dimensions
use json_module, only: json_file, json_initialize, json_destroy, &
                       json_value, json_value_create, to_object, json_value_add, &
                       json_print
use module_data, only: data_filename, data_post_filename, report_filename, report_post_filename, &
                       test_file_exists, test_data_found, test_empty_string
use module_types_micro, only: Box_Parameters, System_Arguments, Argument_Configurations
use module_geometry, only: geometry, set_geometry
use module_physics_micro, only: num_wave_vectors
use class_external_field, only: External_Field
use class_hard_spheres, only: Hard_Spheres, Dipolar_Hard_Spheres, Between_Hard_Spheres
use class_hard_spheres_potential, only: Between_Hard_Spheres_Potential_Energy
use module_types_macro, only: Hard_Spheres_Macro, Dipolar_Hard_Spheres_Macro
use class_discrete_observable, only: Discrete_Observables
use class_hard_spheres_observables, only: Hard_Spheres_Monte_Carlo_Observables, &
                                          Dipolar_Hard_Spheres_Monte_Carlo_Observables, &
                                          Between_Hard_Spheres_Monte_Carlo_Observables, &
                                          Hard_Spheres_Post_Processing_Observables
use class_distribution_function, only: Distribution_Function
use class_hard_spheres_units, only: Hard_Spheres_Monte_Carlo_Units, &
                                    Dipolar_Hard_Spheres_Monte_Carlo_Units, &
                                    Between_Hard_Spheres_Monte_Carlo_Units, &
                                    Hard_Spheres_Post_Processing_Units
use module_physics_macro, only: init_random_seed, set_initial_configuration, &
                                init_spheres, init_cells, reset_cells, &
                                set_ewald, total_energy, final_spheres, &
                                init_between_spheres_potential, final_between_spheres_potential, &
                                test_consistency
use module_algorithms, only: move, widom, switch, rotate

implicit none

private

    type, public :: System
    
        private
        
        type(json_file) :: data_json
        
        character(len=:), allocatable :: name
        
        ! Box
        type(Box_Parameters) :: Box
        type(External_Field) :: ext_field
        
        ! Monte-Carlo
        integer :: num_thermalisation_steps, num_equilibrium_steps
        
        ! Type 1: Dipolar spheres
        type(Dipolar_Hard_Spheres) :: type1_spheres
        type(Dipolar_Hard_Spheres_Macro) :: type1_macro
        type(json_value), pointer :: type1_report_json
        
        ! Type 2: Hard spheres
        type(Hard_Spheres) :: type2_spheres
        type(Hard_Spheres_Macro) :: type2_macro
        type(json_value), pointer :: type2_report_json
        
        ! Between Spheres potential_energy
        type(Between_Hard_Spheres) :: between_spheres
        type(Between_Hard_Spheres_Potential_Energy) :: between_spheres_potential
        
        ! Observables and files units
        type(json_value), pointer :: report_json
        integer :: report_unit

        integer :: reset_i_step
        logical :: snap
        real(DP) :: time_start, time_end
    
    contains
    
        !> Construction & destruction of the class
        procedure :: construct => System_construct
        procedure, private :: set_box => System_set_box
        procedure, private :: set_monte_carlo_steps => System_set_monte_carlo_steps
        procedure :: destroy => System_destroy
        
        procedure, private :: json_create_all_values => System_json_create_all_values
        procedure, private :: json_destroy_all_values => System_json_destroy_all_values
        procedure, private :: json_write_results => System_json_write_results
        procedure, private :: set_potentials => System_set_potentials
        
        !> Accessors & Mutators
        procedure :: get_num_thermalisation_steps => System_get_num_thermalisation_steps
        procedure :: get_num_equilibrium_steps => System_get_num_equilibrium_steps
        procedure :: set_time_start => System_set_time_start
        procedure :: set_time_end => System_set_time_end
    
    end type System
    
    type, extends(System), public :: System_Monte_Carlo
    
        private

        integer :: period_adaptation
        integer :: decorrelation_factor, num_changes, num_moves, num_switches, num_rotations
        
        type(Dipolar_Hard_Spheres_Monte_Carlo_Units) :: type1_units
        type(Dipolar_Hard_Spheres_Monte_Carlo_Observables) :: type1_observables
        
        type(Hard_Spheres_Monte_Carlo_Units) :: type2_units
        type(Hard_Spheres_Monte_Carlo_Observables) :: type2_observables
        
        type(Between_Hard_Spheres_Monte_Carlo_Observables) :: between_spheres_observables
        type(Between_Hard_Spheres_Monte_Carlo_Units) :: between_spheres_units
        type(json_value), pointer :: between_spheres_report_json
        
        real(DP) :: potential_energy, potential_energy_sum
        type(json_value), pointer :: system_json
        integer :: observables_thermalisation_unit, observables_equilibrium_unit
        type(Discrete_Observables) :: switch_observable
        logical :: write_potential_energy
        
    contains
    
        procedure, private :: set_changes => System_Monte_Carlo_set_changes
        procedure :: init => System_Monte_Carlo_init
        procedure, private :: open_all_units => System_Monte_Carlo_open_all_units
        procedure, private :: open_units => System_Monte_Carlo_open_units
        procedure, private :: write_all_reports => System_Monte_Carlo_write_all_reports
        procedure, private :: write_report => System_Monte_Carlo_write_report
        procedure :: final => System_Monte_Carlo_final
        procedure, private :: write_all_results => System_Monte_Carlo_write_all_results
        procedure, private :: write_results => System_Monte_Carlo_write_results
        procedure, private :: json_destroy_all_values => System_Monte_Carlo_json_destroy_all_values
        procedure, private :: close_units => System_Monte_Carlo_close_units
                
        !> Simulation
        procedure :: random_changes => System_Monte_Carlo_random_changes
        procedure :: update_rejections => System_Monte_Carlo_update_rejections
        procedure :: adapt_changes => System_Monte_Carlo_adapt_changes
        procedure :: write_observables_thermalisation => &
                     System_Monte_Carlo_write_observables_thermalisation
        procedure :: fix_changes => System_Monte_Carlo_fix_changes
        procedure :: accumulate_observables => System_Monte_Carlo_accumulate_observables
        procedure :: write_observables_equilibrium => System_Monte_Carlo_write_observables_equilibrium
        procedure :: take_snapshots => System_Monte_Carlo_take_snapshots
        procedure :: reset_potentials => System_Monte_Carlo_reset_potentials
    
    end type System_Monte_Carlo
    
    type, extends(System), public :: System_Post_Processing
    
        private

        type(json_file) :: data_post_json

        integer :: num_steps
        type(json_file) :: data_report_json
        
        type(Hard_Spheres_Post_Processing_Observables) :: type1_observables, type2_observables
        type(Hard_Spheres_Post_Processing_Units) :: type1_units, type2_units
        type(Distribution_Function) :: type1_field_distribution
        integer :: type1_positions_unit, type1_orientations_unit, type2_positions_unit
        
        logical :: first_set
        
    contains
    
        !> Initialization & Finalisation
        procedure :: init => System_Post_Processing_init
        procedure :: set_first => System_Post_Processing_set_first
        procedure, private :: open_all_units_in => System_Post_Processing_open_all_units_in
        procedure, private :: open_all_units => System_Post_Processing_open_all_units
        procedure :: set_coordinates => System_Post_Processing_set_coordinates
        procedure :: final => System_Post_Processing_final
        procedure, private :: write_all_results => System_Post_Processing_write_all_results
        procedure, private :: close_units_in => System_Post_Processing_close_units_in
        procedure, private :: close_units => System_Post_Processing_close_units
        
        procedure :: get_first_set => System_Post_Processing_get_first_set
        
        procedure :: measure_chemical_potentials => &
                     System_Post_Processing_measure_chemical_potentials
        procedure :: accumulate_observables => System_Post_Processing_accumulate_observables
        procedure :: write_observables => System_Post_Processing_write_observables
        procedure :: reset_potentials => System_Post_Processing_reset_potentials
        
    end type System_Post_Processing
    
contains

    ! Construction & destruction of the class
    
    subroutine System_construct(this, args)
            
        class(System), intent(out) :: this
        type(System_Arguments), intent(in) :: args
        
        character(len=4096) :: data_name
        logical :: found
        character(len=:), allocatable :: this_name
        
        character(len=:), allocatable :: Box_geometry
        real(DP) :: Box_height

        call json_initialize()
        
        call test_file_exists(data_filename)
        call this%data_json%load_file(filename = data_filename)

        select type (this)
            type is (System_Monte_Carlo)
                call set_geometry(args%geometry)
            type is (System_Post_Processing)
            
                call test_file_exists(report_filename)
                call this%data_report_json%load_file(report_filename)
                
                data_name = "System.Box.geometry"
                call this%data_report_json%get(data_name, Box_geometry, found)
                call test_data_found(data_name, found)
                call set_geometry(Box_geometry)
                if (allocated(Box_geometry)) deallocate(Box_geometry)

                call test_file_exists(data_post_filename)
                call this%data_post_json%load_file(filename = data_post_filename)
                
        end select
        
        data_name = "Box.name"
        call this%data_json%get(data_name, this_name, found)
        call test_data_found(data_name, found)
        call test_empty_string(data_name, this_name)
        this%name = this_name
        if (allocated(this_name)) deallocate(this_name)
        write(output_unit, *) this%name, " class construction"
        
        call this%set_box()
        call this%set_monte_carlo_steps()
        
        call this%type1_spheres%construct(this%Box, this%data_json, "Dipolar Hard Spheres")
        call this%type2_spheres%construct(this%Box, this%data_json, "Hard Spheres")
        call this%between_spheres%construct(this%data_json, &
                                            this%type1_spheres%get_diameter(), &
                                            this%type2_spheres%get_diameter())

        this%Box%num_particles = this%type1_spheres%get_num_particles() + &
                                 this%type2_spheres%get_num_particles()
                                            
        select type (this)
            type is (System_Monte_Carlo)
            
                data_name = "Potential Energy.write"
                call this%data_json%get(data_name, this%write_potential_energy, found)
                call test_data_found(data_name, found)
                
                call this%set_changes()

            type is (System_Post_Processing)

                call this%type1_spheres%set_test_num_particles(this%data_post_json, &
                                                               "Dipolar Hard Spheres")
                call this%type2_spheres%set_test_num_particles(this%data_post_json, &
                                                               "Hard Spheres")
                if (geometry%bulk) then
                    Box_height = this%Box%size(3)
                else if (geometry%slab) then
                    Box_height = this%Box%height
                end if
                call this%type1_field_distribution%construct(this%data_post_json, Box_height, &
                                                             num_dimensions)
                
        end select
        
    end subroutine System_construct
    
    subroutine System_set_box(this)
    
        class(System), intent(inout) :: this
        
        character(len=4096) :: data_name
        logical :: found
        real(DP), dimension(:), allocatable :: Box_size
        integer, dimension(:), allocatable :: Box_wave
        real(DP), dimension(:), allocatable :: external_field

        real(DP) :: z_ratio
        
        data_name = "Box.size"
        call this%data_json%get(data_name, Box_size, found)
        call test_data_found(data_name, found)
        
        data_name = "Box.wave"
        call this%data_json%get(data_name, Box_wave, found)
        call test_data_found(data_name, found)
        
        if (geometry%bulk) then
            if (size(Box_size) /= num_dimensions) error stop "Box size dimension"
            this%Box%size(:) = Box_size(:)
            if (size(Box_wave) /= num_dimensions) error stop "Box wave dimension"
            this%Box%wave(:) = Box_wave(:)
        else if (geometry%slab) then
            if (size(Box_size) /= 2) error stop "Box size dimension"
            this%Box%size(1:2) = Box_size(:)
            if (size(Box_wave) /= 2) error stop "Box wave dimension"
            this%Box%wave(1:2) = Box_wave(:)
        end if
        
        if (allocated(Box_size)) deallocate(Box_size)
        if (allocated(Box_wave)) deallocate(Box_wave)        

        if (geometry%slab) then
        
            data_name = "Box.height"
            call this%data_json%get(data_name, this%Box%height, found)
            call test_data_found(data_name, found)
            
            data_name = "Box.ratio in z"
            call this%data_json%get(data_name, z_ratio, found)
            call test_data_found(data_name, found)
            
            this%Box%size(num_dimensions) = z_ratio * this%Box%size(1)
            this%Box%wave(num_dimensions) = ceiling(z_ratio * real(this%Box%wave(1), DP))
            
        end if
        
        data_name = "Box.temperature"
        call this%data_json%get(data_name, this%Box%temperature, found)
        call test_data_found(data_name, found)

        data_name = "Box.external field"
        call this%data_json%get(data_name, external_field, found)
        call test_data_found(data_name, found)
        if (size(external_field) /= num_dimensions) error stop "Box external field dimension"
        call this%ext_field%set(external_field)
        if (allocated(external_field)) deallocate(external_field)
        
    end subroutine System_set_box
    
    subroutine System_set_monte_carlo_steps(this)
    
        class(System), intent(inout) :: this
        
        character(len=4096) :: data_name
        logical :: found
        
        data_name = "Monte Carlo.number of thermalisation steps"
        call this%data_json%get(data_name, this%num_thermalisation_steps, found)
        call test_data_found(data_name, found)
        
        data_name = "Monte Carlo.number of equilibrium steps"
        call this%data_json%get(data_name, this%num_equilibrium_steps, found)
        call test_data_found(data_name, found)
        
    end subroutine System_set_monte_carlo_steps
    
    subroutine System_Monte_Carlo_set_changes(this)
    
        class(System_Monte_Carlo), intent(inout) :: this
        
        integer :: switch_factor
        real(DP) :: type1_move_delta, type1_move_rejection
        real(DP) :: type1_rotation_delta, type1_rotation_delta_max, type1_rotation_rejection
        real(DP) :: type2_move_delta, type2_move_rejection
        
        character(len=4096) :: data_name
        logical :: found
        
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
        
        data_name = "Monte Carlo.period of adaptation"
        call this%data_json%get(data_name, this%period_adaptation, found)
        call test_data_found(data_name, found)
        
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
        
    end subroutine System_Monte_Carlo_set_changes

    subroutine System_destroy(this)

        class(System), intent(inout) :: this

        write(output_unit, *) this%name, " class destruction"

        call this%between_spheres%destroy()

        call this%type2_macro%between_cells%destroy()
        call this%type2_macro%same_cells%destroy()
        call this%type2_spheres%destroy()

        if (geometry%slab) call this%type1_macro%elc%destroy()
        call this%type1_macro%ewald_reci%destroy()
        call this%type1_macro%ewald_real%destroy()
        call this%type1_macro%between_cells%destroy()
        call this%type1_macro%same_cells%destroy()
        call this%type1_spheres%destroy()

        if (allocated(this%name)) deallocate(this%name)

        call this%data_json%destroy()

        select type(this)
            type is (System_Post_Processing)
                call this%type1_field_distribution%destroy()
                call this%data_post_json%destroy()
        end select

    end subroutine System_destroy
    
    ! Initialization
    
    subroutine System_Monte_Carlo_init(this, args)
     
        class(System_Monte_Carlo), intent(inout) :: this
        type(System_Arguments), intent(in) :: args
        
        real(DP) :: potential_energy_conf
        
        character(len=4096) :: data_name
        logical :: found
        
        call this%open_all_units()
        call this%json_create_all_values()
        
        data_name = "Distribution.take snapshot"
        call this%data_json%get(data_name, this%snap, found)
        call test_data_found(data_name, found)
        
        data_name = "Monte Carlo.period of reset"
        call this%data_json%get(data_name, this%reset_i_step, found)
        call test_data_found(data_name, found)
        this%reset_i_step = this%reset_i_step / this%decorrelation_factor
        
        write(output_unit, *) "Monte-Carlo Simulation: Canonical ensemble "
        
        call init_random_seed(args%random, this%report_json)
        call set_initial_configuration(this%Box, args%conf, &
                                       this%type1_spheres, this%type2_spheres, &
                                       this%between_spheres%get_diameter(), &
                                       this%report_json)
        call this%set_potentials()
                                       
        call this%between_spheres%test_overlap(this%Box%size, &
                                               this%type1_spheres, this%type2_spheres)
        
        call init_between_spheres_potential(this%Box%size, this%between_spheres_potential, &
                                            this%type1_spheres, this%type2_spheres, &
                                            this%write_potential_energy, &
                                            this%between_spheres_observables%potential_energy, &
                                            this%between_spheres_units%potential_energy_tabulation)
        
        call init_spheres(this%Box, this%type1_spheres, this%type1_units)
        call set_ewald(this%Box, this%type1_spheres, this%type1_macro, this%data_json, &
                       this%type1_units)
        this%type1_observables%potential_energy = total_energy(this%Box, this%type1_spheres, &
                                                               this%type1_macro, this%ext_field)
                        
        call init_spheres(this%Box, this%type2_spheres, this%type2_units)
        this%type2_observables%potential_energy = total_energy(this%Box, this%type2_spheres, &
                                                               this%type2_macro, this%ext_field)
                        
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
    
    end subroutine System_Monte_Carlo_init
    
    subroutine System_Monte_Carlo_open_all_units(this)
    
        class(System_Monte_Carlo), intent(inout) :: this

        call this%open_units()
        call this%type1_units%open(this%type1_spheres%get_name())
        call this%type2_units%open(this%type2_spheres%get_name())
        call this%between_spheres_units%open(this%between_spheres%get_name())
    
    end subroutine System_Monte_Carlo_open_all_units

    subroutine System_Monte_Carlo_open_units(this)

        class(System_Monte_Carlo), intent(inout) :: this

        open(newunit=this%report_unit, recl=4096, file=report_filename, &
             status='new', action='write')
        open(newunit=this%observables_thermalisation_unit, recl=4096, &
             file="observables_thermalisation.out", status='new', action='write')
        open(newunit=this%observables_equilibrium_unit, recl=4096, &
             file="observables_equilibrium.out", status='new', action='write')
        write(this%observables_equilibrium_unit, *) "#", 1 ! 1 observable: energy

    end subroutine System_Monte_Carlo_open_units

    subroutine System_json_create_all_values(this)

        class(System), intent(inout) :: this

        call json_value_create(this%report_json)
        call to_object(this%report_json)
        
        select type (this)
            type is (System_Monte_Carlo)

                call json_value_create(this%system_json)
                call to_object(this%system_json, "System")
                call json_value_add(this%report_json, this%system_json)

                call json_value_create(this%between_spheres_report_json)
                call to_object(this%between_spheres_report_json, "Between Spheres")
                call json_value_add(this%report_json, this%between_spheres_report_json)
                
        end select

        call json_value_create(this%type1_report_json)
        call to_object(this%type1_report_json, "Dipolar Hard Spheres")
        call json_value_add(this%report_json, this%type1_report_json)

        call json_value_create(this%type2_report_json)
        call to_object(this%type2_report_json, "Hard Spheres")
        call json_value_add(this%report_json, this%type2_report_json)
        
    end subroutine System_json_create_all_values
    
    subroutine System_set_potentials(this)
    
        class(System), intent(inout) :: this
        
        call this%between_spheres_potential%construct(this%data_json, "Between Spheres", &
                                                      this%between_spheres%get_diameter())
                                                      
        call this%type1_macro%hard_potential%construct(this%data_json, "Dipolar Hard Spheres", &
                                                       this%type1_spheres%get_diameter())
        call init_cells(this%Box%size, this%type1_spheres, this%type1_macro, this%type2_spheres, &
                        this%between_spheres_potential)
        
        call this%type2_macro%hard_potential%construct(this%data_json, "Hard Spheres", &
                                                       this%type2_spheres%get_diameter())
        call init_cells(this%Box%size, this%type2_spheres, this%type2_macro, this%type1_spheres, &
                        this%between_spheres_potential)
        
    end subroutine System_set_potentials
    
    subroutine System_Monte_Carlo_write_all_reports(this, geometry)
    
        class(System_Monte_Carlo), intent(in) :: this
        character(len=*), intent(in) :: geometry
        
        call this%write_report(geometry)
        
        call this%type1_spheres%write_report(this%Box%num_particles, this%type1_report_json)
        call this%type1_macro%hard_potential%write_report(this%type1_report_json)
        
        call this%type2_spheres%write_report(this%Box%num_particles, this%type2_report_json)
        call this%type2_macro%hard_potential%write_report(this%type2_report_json)
        
        call this%between_spheres_potential%write_report(this%between_spheres_report_json)
    
    end subroutine System_Monte_Carlo_write_all_reports
    
    subroutine System_Monte_Carlo_write_report(this, geometry)
    
        class(System_Monte_Carlo), intent(in) :: this
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
        
    end subroutine System_Monte_Carlo_write_report
    
    subroutine System_Post_Processing_init(this, args)
    
        class(System_Post_Processing), intent(inout) :: this
        type(System_Arguments), intent(in) :: args
        
        call this%open_all_units_in(args%conf)
        call this%open_all_units()
        call this%json_create_all_values()

        this%num_steps = 0 
        this%reset_i_step = 1
        
        call this%type1_spheres%set_data(this%type1_positions_unit)
        call this%type1_spheres%set_data(this%type1_orientations_unit)
        call this%type2_spheres%set_data(this%type2_positions_unit)
        
        this%first_set = .true.
        
    end subroutine System_Post_Processing_init
    
    subroutine System_Post_Processing_set_first(this)
    
        class(System_Post_Processing), intent(inout) :: this
        
        call this%set_potentials()                        
        call set_ewald(this%Box, this%type1_spheres, this%type1_macro, this%data_json)
        
        this%first_set = .false.
        
    end subroutine System_Post_Processing_set_first
    
    subroutine System_Post_Processing_open_all_units_in(this, arg_conf)
    
        class(System_Post_Processing), intent(inout) :: this
        type(Argument_Configurations), intent(in) :: arg_conf
        
        character(len=4096) :: filename
        integer :: length
        
        filename = arg_conf%files(1)
        length = arg_conf%length(1)
        open(newunit=this%type1_positions_unit, recl=4096, file=filename(1:length), &
             status='old', action='read')
             
        filename = arg_conf%files(2)
        length = arg_conf%length(2)
        open(newunit=this%type1_orientations_unit, recl=4096, file=filename(1:length), &
             status='old', action='read')
             
        filename = arg_conf%files(3)
        length = arg_conf%length(3)
        open(newunit=this%type2_positions_unit, recl=4096, file=filename(1:length), &
             status='old', action='read')
        
    end subroutine System_Post_Processing_open_all_units_in
    
    subroutine System_Post_Processing_open_all_units(this)
    
        class(System_Post_Processing), intent(inout) :: this

        open(newunit=this%report_unit, recl=4096, file=report_post_filename, status='new', &
             action='write')
        call this%type1_units%open(this%type1_spheres%get_name())
        call this%type2_units%open(this%type2_spheres%get_name())
    
    end subroutine System_Post_Processing_open_all_units
    
    subroutine System_Post_Processing_set_coordinates(this, i_step, coordinates_set)
    
        class(System_Post_Processing), intent(inout) :: this
        integer, intent(in) :: i_step
        logical :: coordinates_set
        
        logical :: type1_positions_set, type1_orientations_set
        logical :: type2_positions_set
        
        call this%type1_spheres%set_all_positions(this%Box%size, i_step, this%type1_positions_unit, &
                                                  type1_positions_set)
        call this%type1_spheres%set_all_orientations(i_step, this%type1_orientations_unit, &
                                                     type1_orientations_set)
        call this%type2_spheres%set_all_positions(this%Box%size, i_step, this%type2_positions_unit, &
                                                  type2_positions_set)
                                                  
        coordinates_set = type1_positions_set .and. type1_orientations_set .and. type2_positions_set
    
    end subroutine System_Post_Processing_set_coordinates
    
    ! Finalisation
    
    subroutine System_Monte_Carlo_final(this)
    
        class(System_Monte_Carlo), intent(inout) :: this
        
        real(DP) :: type1_energy, type2_energy
        
        call final_spheres(this%Box, this%type1_spheres, this%type1_units)
        call set_ewald(this%Box, this%type1_spheres, this%type1_macro, this%data_json, this%type1_units)
        type1_energy = total_energy(this%Box, this%type1_spheres, this%type1_macro, this%ext_field)
        call test_consistency(this%type1_observables%potential_energy, type1_energy, &
                              this%type1_report_json)
        
        call final_spheres(this%Box, this%type2_spheres, this%type2_units)
        type2_energy = total_energy(this%Box, this%type2_spheres, this%type2_macro, this%ext_field)
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
    
    end subroutine System_Monte_Carlo_final
    
    subroutine System_Monte_Carlo_write_all_results(this)
    
        class(System_Monte_Carlo), intent(inout) :: this
        
        call this%type1_observables%write_results(this%num_equilibrium_steps, this%type1_report_json)
        call this%type2_observables%write_results(this%num_equilibrium_steps, this%type2_report_json)
        call this%between_spheres_observables%write_results(this%num_equilibrium_steps, &
                                                            this%between_spheres_report_json)
        
        call this%write_results()
    
    end subroutine System_Monte_Carlo_write_all_results
    
    subroutine System_Monte_Carlo_write_results(this)
    
        class(System_Monte_Carlo), intent(inout) :: this
        
        real(DP) :: potential_energy, potential_energy_conf
        
        potential_energy = this%type1_observables%potential_energy + &
                           this%type2_observables%potential_energy + &
                           this%between_spheres_observables%potential_energy
        potential_energy_conf = total_energy(this%Box, this%type1_spheres, this%type1_macro, &
                                             this%ext_field) + &
                                total_energy(this%Box, this%type2_spheres, this%type2_macro, &
                                             this%ext_field) + &
                                this%between_spheres_potential%total(this%Box%size, &
                                                                     this%type1_spheres, &
                                                                     this%type2_spheres)
        
        call test_consistency(potential_energy, potential_energy_conf, this%report_json)
        this%potential_energy_sum = this%type1_observables%potential_energy_sum + &
                                    this%type2_observables%potential_energy_sum + &
                                    this%between_spheres_observables%potential_energy_sum
                                    
        call this%json_write_results()
    
    end subroutine System_Monte_Carlo_write_results

    subroutine System_json_write_results(this)

        class(System), intent(inout) :: this

        type(json_value), pointer :: results_json

        call json_value_create(results_json)
        call to_object(results_json, "Results")
        call json_value_add(this%report_json, results_json)

        select type (this)
            type is (System_Monte_Carlo)
                call json_value_add(results_json, "average energy", &
                                                this%potential_energy_sum / &
                                                real(this%num_equilibrium_steps, DP))
                if (this%Box%num_particles > 0) then
                    call json_value_add(results_json, "average energy per particule", &
                                                    this%potential_energy_sum / &
                                                    real(this%num_equilibrium_steps, DP) / &
                                                    real(this%Box%num_particles, DP)) ! DHS only?
                end if
                call json_value_add(results_json, "switch rejection rate", &
                                                this%switch_observable%sum_rejection / &
                                                real(this%num_equilibrium_steps, DP))
            type is (System_Post_Processing)
                call json_value_add(results_json, "number of steps", this%num_steps)
        end select
        
        call json_value_add(results_json, "duration (min)", &
                                          (this%time_end - this%time_start) / 60._DP)

        nullify(results_json)

        rewind(this%report_unit)
        call json_print(this%report_json, this%report_unit)

    end subroutine System_json_write_results

    subroutine System_json_destroy_all_values(this)

        class(System), intent(inout) :: this

        nullify(this%type2_report_json)
        nullify(this%type1_report_json)
        
        call json_destroy(this%report_json)
    
    end subroutine System_json_destroy_all_values

    subroutine System_Monte_Carlo_json_destroy_all_values(this)

        class(System_Monte_Carlo), intent(inout) :: this

        nullify(this%between_spheres_report_json)
        nullify(this%system_json)

        call this%System%json_destroy_all_values()

    end subroutine System_Monte_Carlo_json_destroy_all_values

    subroutine System_Monte_Carlo_close_units(this)

        class(System_Monte_Carlo), intent(inout) :: this

        call this%type1_units%close()
        call this%type2_units%close()
        call this%between_spheres_units%close()
        
        close(this%observables_thermalisation_unit)
        close(this%observables_equilibrium_unit)

        close(this%report_unit)
    
    end subroutine System_Monte_Carlo_close_units

    subroutine System_Post_Processing_final(this)

        class(System_Post_Processing), intent(inout) :: this

        call this%write_all_results()
        call this%json_destroy_all_values()
        call this%close_units()
        call this%close_units_in()

    end subroutine System_Post_Processing_final

    subroutine System_Post_Processing_write_all_results(this)

        class(System_Post_Processing), intent(inout) :: this
        
        call this%type1_observables%write_results(this%Box%temperature, this%num_steps, &
                                                  this%type1_spheres%get_widom_num_particles(), &
                                                  this%type1_report_json)
        call this%type1_field_distribution%write(this%num_steps, this%type1_units%local_field)
        call this%type2_observables%write_results(this%Box%temperature, this%num_steps, &
                                                  this%type2_spheres%get_widom_num_particles(), &
                                                  this%type2_report_json)

        call this%json_write_results()

    end subroutine System_Post_Processing_write_all_results

    subroutine System_Post_Processing_close_units_in(this)

        class(System_Post_Processing), intent(inout) :: this

        close(this%type2_positions_unit)
        close(this%type1_orientations_unit)
        close(this%type1_positions_unit)

    end subroutine System_Post_Processing_close_units_in

    subroutine System_Post_Processing_close_units(this)

        class(System_Post_Processing), intent(inout) :: this

        call this%type1_units%close()
        call this%type2_units%close()

        close(this%report_unit)

    end subroutine System_Post_Processing_close_units
    
    ! Accessors
    
    pure function System_get_num_thermalisation_steps(this) &
                  result(get_num_thermalisation_steps)
        
        class(System), intent(in) :: this
        integer :: get_num_thermalisation_steps
                
        get_num_thermalisation_steps = this%num_thermalisation_steps
        
    end function System_get_num_thermalisation_steps
    
    pure function System_get_num_equilibrium_steps(this) result(get_num_equilibrium_steps)
    
        class(System), intent(in) :: this
        integer :: get_num_equilibrium_steps
                
        get_num_equilibrium_steps = this%num_equilibrium_steps
        
    end function System_get_num_equilibrium_steps

    pure function System_Post_Processing_get_first_set(this) result(get_first_set)

        class(System_Post_Processing), intent(in) :: this
        logical :: get_first_set

        get_first_set = this%first_set

    end function System_Post_Processing_get_first_set
    
    ! Mutators
    
    subroutine System_set_time_start(this)
        class(System), intent(inout) :: this
              
        call cpu_time(this%time_start)
        
    end subroutine System_set_time_start
    
    subroutine System_set_time_end(this)
        class(System), intent(inout) :: this
              
        call cpu_time(this%time_end)
        
    end subroutine System_set_time_end
    
    subroutine System_Monte_Carlo_random_changes(this)
    
        class(System_Monte_Carlo), intent(inout) :: this
        
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
                call rotate(this%Box, this%ext_field, &
                            this%type1_spheres, this%type1_macro, this%type1_observables)
            end if
            
        end do MC_Change
    
    end subroutine System_Monte_Carlo_random_changes
    
    subroutine System_Monte_Carlo_update_rejections(this)
    
        class(System_Monte_Carlo), intent(inout) :: this
    
        call this%type1_observables%move%update_rejection()
        call this%type1_observables%rotation%update_rejection()
        call this%type2_observables%move%update_rejection()
        call this%switch_observable%update_rejection()
        
    end subroutine System_Monte_Carlo_update_rejections
    
    subroutine System_Monte_Carlo_adapt_changes(this, i_step)
    
        class(System_Monte_Carlo), intent(inout) :: this
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
        
    end subroutine System_Monte_Carlo_adapt_changes
    
    subroutine System_Monte_Carlo_write_observables_thermalisation(this, i_step)
    
        class(System_Monte_Carlo), intent(inout) :: this
        integer, intent(in) :: i_step
    
        call this%type1_observables%write(i_step, this%type1_units%observables_thermalisation)
        call this%type2_observables%write(i_step, this%type2_units%observables_thermalisation)
        call this%between_spheres_observables%write(i_step, &
            this%between_spheres_units%observables_thermalisation)
        
        write(this%observables_thermalisation_unit, *) i_step, &
            this%type1_observables%potential_energy + &
            this%type2_observables%potential_energy + &
            this%between_spheres_observables%potential_energy
            
    end subroutine System_Monte_Carlo_write_observables_thermalisation
    
    subroutine System_Monte_Carlo_fix_changes(this)
    
        class(System_Monte_Carlo), intent(inout) :: this
        
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
    
    end subroutine System_Monte_Carlo_fix_changes

    subroutine System_Post_Processing_measure_chemical_potentials(this)

        class(System_Post_Processing), intent(inout) :: this

        call widom(this%Box, this%ext_field, &
                   this%type1_spheres, this%type1_macro, this%type1_observables, &
                   this%type2_spheres, this%type2_macro%between_cells, &
                   this%between_spheres_potential)
        call widom(this%Box, this%ext_field, &
                   this%type2_spheres, this%type2_macro, this%type2_observables, &
                   this%type1_spheres, this%type1_macro%between_cells, &
                   this%between_spheres_potential)

    end subroutine System_Post_Processing_measure_chemical_potentials
    
    subroutine System_Monte_Carlo_accumulate_observables(this)
    
        class(System_Monte_Carlo), intent(inout) :: this
    
        call this%type1_observables%accumulate()
        call this%type2_observables%accumulate()
        call this%between_spheres_observables%accumulate()
        this%switch_observable%sum_rejection = this%switch_observable%sum_rejection + &
                                               this%switch_observable%rejection_rate
            
    end subroutine System_Monte_Carlo_accumulate_observables

    subroutine System_Post_Processing_accumulate_observables(this)

        class(System_Post_Processing), intent(inout) :: this

        this%num_steps = this%num_steps + 1
        call this%type1_observables%accumulate()
        call this%type2_observables%accumulate()

    end subroutine System_Post_Processing_accumulate_observables
    
    subroutine System_Monte_Carlo_write_observables_equilibrium(this, i_step)
    
        class(System_Monte_Carlo), intent(inout) :: this
        integer, intent(in) :: i_step
    
        call this%type1_observables%write(i_step, this%type1_units%observables_equilibrium)
        call this%type2_observables%write(i_step, this%type2_units%observables_equilibrium)
        call this%between_spheres_observables%write(i_step, &
                                                    this%between_spheres_units%observables_equilibrium)
        
        write(this%observables_equilibrium_unit, *) i_step, &
            this%type1_observables%potential_energy + &
            this%type2_observables%potential_energy + &
            this%between_spheres_observables%potential_energy
            
    end subroutine System_Monte_Carlo_write_observables_equilibrium

    subroutine System_Post_Processing_write_observables(this, i_step)

        class(System_Post_Processing), intent(inout) :: this
        integer, intent(in) :: i_step

        call this%type1_observables%write(i_step, this%type1_units%inv_activity)
        call this%type2_observables%write(i_step, this%type2_units%inv_activity)

    end subroutine System_Post_Processing_write_observables
    
    subroutine System_Monte_Carlo_take_snapshots(this, i_step)
    
        class(System_Monte_Carlo), intent(inout) :: this
        integer, intent(in) :: i_step
        
        if (this%snap) then ! Snap shots of the configuration
            call this%type1_spheres%write_all_positions(i_step, &
                 this%type1_units%snap_equilibrium_positions)
            call this%type1_spheres%write_all_orientations(i_step, &
                 this%type1_units%snap_equilibrium_orientations)
            call this%type2_spheres%write_all_positions(i_step, &
                 this%type2_units%snap_equilibrium_positions)
        end if
        
    end subroutine System_Monte_Carlo_take_snapshots
    
    subroutine System_Monte_Carlo_reset_potentials(this, i_step)
    
        class(System_Monte_Carlo), intent(inout) :: this
        integer, intent(in) :: i_step
    
        if (modulo(i_step, this%reset_i_step) == 0) then
            call this%type1_macro%ewald_reci%reset_structure(this%Box, this%type1_spheres, i_step, &
                                                             this%type1_units%structure_modulus)
            call this%type1_macro%ewald_bound%reset_total_moment(this%type1_spheres, i_step, &
                 this%type1_units%total_moment_modulus)
            if (geometry%slab) then
                call this%type1_macro%elc%reset_structure(this%Box, this%type1_spheres, i_step, &
                                                          this%type1_units%ELC_structure_modulus)
            end if
        end if
        
    end subroutine System_Monte_Carlo_reset_potentials

    subroutine System_Post_Processing_reset_potentials(this, i_step)

        class(System_Post_Processing), intent(inout) :: this
        integer, intent(in) :: i_step

        if (modulo(i_step, this%reset_i_step) == 0) then
            call this%type1_macro%ewald_reci%reset_structure(this%Box, this%type1_spheres, i_step)
            call this%type1_macro%ewald_bound%reset_total_moment(this%type1_spheres, i_step)
            if (geometry%slab) then
                call this%type1_macro%elc%reset_structure(this%Box, this%type1_spheres, i_step)
            end if
            call reset_cells(this%type1_spheres, this%type1_macro, this%type2_spheres)
            call reset_cells(this%type2_spheres, this%type2_macro, this%type1_spheres)
        end if

    end subroutine System_Post_Processing_reset_potentials

end module class_system
