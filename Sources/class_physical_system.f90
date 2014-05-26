!> \brief Description of the Physical System class

module class_physical_system

use, intrinsic :: iso_fortran_env, only: output_unit
use data_precisions, only: DP
use json_module, only: json_file
use module_types_micro, only: Box_Dimensions, Monte_Carlo_Arguments
use module_physics_micro, only: NwaveVectors
use module_data, only: test_data_found
use module_types_macro, only: Hard_Spheres_Macro, Dipolar_Hard_Spheres_Macro
use class_hard_spheres
use class_dipolar_hard_spheres
use class_mixing_potential
use class_observables
use class_units
use module_monte_carlo_arguments, only: read_arguments
use module_physics_macro, only: init_random_seed, set_initial_configuration, &
                                init_spheres, init_hard_potential, init_cells, &
                                set_ewald, &
                                total_energy, &
                                final_spheres, &
                                init_mix, mix_final, &
                                adapt_move, adapt_rotation, test_consist
use module_algorithms, only: move, widom, switch, rotate
use module_write, only: open_units, write_data, mix_open_units, write_results, mix_write_results

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
        type(MoreObservables) :: type1_obs
        type(MoreUnits) :: type1_units
        
        ! Type 2: Hard spheres
        type(Hard_Spheres) :: type2_spheres
        type(Hard_Spheres_Macro) :: type2_macro
        type(Observables) :: type2_obs
        type(Units) :: type2_units
        
        ! Mixing potential
        type(Mixing_Potential) :: mix
        real(DP) :: mix_Epot, mix_EpotSum, mix_Epot_conf
        integer :: mix_report_unit
        integer :: mix_Epot_tab_unit
        integer :: mix_obsThermal_unit, mix_obsEquilib_unit
        
        ! Observables: write to files
        real(DP) :: Epot, EpotSum
        integer :: report_unit
        integer :: obsThermal_unit, obsEquilib_unit
        
        ! Switch
        integer :: switch_Nhit, switch_Nreject
        real(DP) :: switch_reject, switch_rejectSum
        
        logical :: write_potential, snap
        integer :: reset_iStep
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
        
        data_name = "Potential.write"
        call json%get(data_name, this%write_potential, found)
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
        
        this%num_switches = switch_factor * this%decorrelation_factor * this%type1_spheres%get_num_particles()
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
        
        real(DP) :: Epot_conf      
        
        character(len=4096) :: data_name
        logical :: found  
        
        data_name = "Distribution.take snapshot"
        call json%get(data_name, this%snap, found)
        call test_data_found(data_name, found)
        
        data_name = "Monte Carlo.period of reset"
        call json%get(data_name, this%reset_iStep, found)
        call test_data_found(data_name, found)
        this%reset_iStep = this%reset_iStep / this%decorrelation_factor
        
        write(output_unit, *) "Monte-Carlo Simulation: Canonical ensemble"
        
        call this%open_all_units()
        call this%init_switch()
        
        call init_random_seed(args%random, this%report_unit)
        call set_initial_configuration(this%Box%size, args%initial, this%type1_spheres, &
                                       this%type2_spheres, this%mix%get_diameter(), &
                                       this%report_unit)
                                       
        this%mix_EpotSum = 0._DP        
        call init_mix(this%Box%size, this%mix, this%type1_spheres, this%type2_spheres, &
                      this%write_potential, this%mix_Epot_tab_unit, this%mix_Epot)
        call this%mix%write_report(this%mix_report_unit)
        
        call init_spheres(this%Box, this%type1_spheres, this%type1_units)
        call init_hard_potential(this%type1_macro%hard_potential, "Dipoles", &
                                 this%type1_spheres%get_diameter(), json)
        call init_cells(this%Box%size, this%type1_spheres, this%type1_macro, this%type2_spheres, &
                        this%mix)
        call set_ewald(this%Box, this%type1_spheres, this%type1_macro, json, this%type1_units)
        this%type1_obs%Epot = total_energy(this%Box, this%type1_spheres, this%type1_macro)
                        
        call init_spheres(this%Box, this%type2_spheres, this%type2_units)
        call init_hard_potential(this%type2_macro%hard_potential, "Hard Spheres", &
                                 this%type2_spheres%get_diameter(), json)
        call init_cells(this%Box%size, this%type2_spheres, this%type2_macro, this%type1_spheres, &
                        this%mix)        
        this%type2_obs%Epot = total_energy(this%Box, this%type2_spheres, this%type2_macro)
                        
        call this%write_all_reports()
        if (this%write_potential) then
            call this%type1_macro%hard_potential%write(this%type1_units%Epot)
            call this%type1_macro%ewald_real%write(this%type1_units%Epot_real)
            call this%type2_macro%hard_potential%write(this%type2_units%Epot)
        end if
        
        this%EpotSum = 0._DP
        Epot_conf = this%type1_obs%Epot + this%type2_obs%Epot + this%mix_Epot
        write(output_unit, *) "Initial potential energy =", Epot_conf
        write(this%obsThermal_unit, *) 0, Epot_conf
    
    end subroutine Physical_System_init
    
    subroutine Physical_System_open_all_units(this)    
        class(Physical_System), intent(inout) :: this
        
        call open_units(this%report_unit, this%obsThermal_unit, this%obsEquilib_unit)        
        call this%type1_units%open(this%type1_spheres%get_name())           
        call this%type2_units%open(this%type2_spheres%get_name())       
        call mix_open_units(this%mix_report_unit, this%mix_Epot_tab_unit, this%mix_obsThermal_unit, &
                            this%mix_obsEquilib_unit)        
    
    end subroutine Physical_System_open_all_units
    
    subroutine Physical_System_write_all_reports(this)    
        class(Physical_System), intent(in) :: this
        
        call write_data(this%report_unit)
        call this%write_report(this%report_unit)
        
        call this%type1_spheres%write_density(this%Box%size, this%num_particles, &
                                              this%type1_units%report)
        call this%type1_spheres%write_report(this%type1_units%report)
        
        call this%type2_spheres%write_density(this%Box%size, this%num_particles, &
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
        
        write(report_unit, *) "    reset_iStep = ", this%reset_iStep
        write(report_unit, *) "    write_potential = ", this%write_potential
    
    end subroutine Physical_System_write_report
    
    pure subroutine Physical_System_init_switch(this)
        class(Physical_System), intent(inout) :: this
        
        this%switch_Nhit = 0
        this%switch_Nreject = 0
        this%switch_reject = 0._DP
        this%switch_rejectSum = 0._DP
        
    end subroutine Physical_System_init_switch
    
    ! Finalisation
    
    subroutine Physical_System_final(this, json)    
    
        class(Physical_System), intent(inout) :: this
        type(json_file), intent(inout) :: json
        
        real(DP) :: type1_energy, type2_energy
        
        call final_spheres(this%Box, this%type1_spheres, this%type1_units)
        call set_ewald(this%Box, this%type1_spheres, this%type1_macro, json, this%type1_units)
        type1_energy = total_energy(this%Box, this%type1_spheres, this%type1_macro)
        call test_consist(this%type1_obs%Epot, type1_energy, this%type1_units%report)
        
        call final_spheres(this%Box, this%type2_spheres, this%type2_units)
        type2_energy = total_energy(this%Box, this%type2_spheres, this%type2_macro)
        call test_consist(this%type2_obs%Epot, type2_energy, this%type2_units%report)
        
        call mix_final(this%Box%size, this%mix, this%type1_spheres, this%type2_spheres, &
                       this%mix_report_unit, this%mix_Epot, this%mix_Epot_conf)
        
        call this%write_all_results()
        call this%close_units()
    
    end subroutine Physical_System_final
    
    subroutine Physical_System_write_all_results(this)
        class(Physical_System), intent(inout) :: this
        
        call this%type1_obs%write_results(this%Box%temperature, this%num_equilibrium_steps, this%type1_units%report)
        call this%type2_obs%write_results(this%Box%temperature, this%num_equilibrium_steps, this%type2_units%report)
        call mix_write_results(this%num_equilibrium_steps, this%mix_EpotSum, this%mix_report_unit)
        
        call this%write_results()
    
    end subroutine Physical_System_write_all_results
    
    subroutine Physical_System_write_results(this)    
        class(Physical_System), intent(inout) :: this
        
        real(DP) :: Epot, Epot_conf
        real(DP) :: duration
        
        Epot = this%type1_obs%Epot + this%type2_obs%Epot + this%mix_Epot
        Epot_conf = total_energy(this%Box, this%type1_spheres, this%type1_macro) + &
                    total_energy(this%Box, this%type2_spheres, this%type2_macro) + &
                    this%mix_Epot_conf
        call test_consist(Epot, Epot_conf, this%report_unit)
        this%EpotSum = this%type1_obs%EpotSum + this%type2_obs%EpotSum + this%mix_EpotSum
        duration = this%time_end - this%time_start
        call write_results(this%num_particles, this%num_equilibrium_steps, this%EpotSum, this%switch_rejectSum,&
                           duration, this%report_unit)
    
    end subroutine Physical_System_write_results
    
    subroutine Physical_System_close_units(this)    
        class(Physical_System), intent(inout) :: this
    
        call this%type1_units%close()
        call this%type2_units%close()
        
        close(this%mix_report_unit)
        close(this%mix_Epot_tab_unit)
        close(this%mix_obsThermal_unit)
        close(this%mix_obsEquilib_unit)
        
        close(this%report_unit)
        close(this%obsThermal_unit)
        close(this%obsEquilib_unit)
    
    end subroutine Physical_System_close_units
    
    ! Destruction
    
    subroutine Physical_System_destroy(this)    
        class(Physical_System), intent(inout) :: this
        
        write(output_unit, *) this%name, " class destruction"
        
        call this%type1_spheres%destroy()
        call this%type2_spheres%destroy()
        call this%mix%destroy()
    
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
                              this%type1_spheres, this%type1_macro, this%type1_obs, &
                              this%type2_spheres, this%type2_macro%mix_cells, &
                              this%mix, this%mix_Epot)
                else
                    call move(this%Box, &
                              this%type2_spheres, this%type2_macro, this%type2_obs, &
                              this%type1_spheres, this%type1_macro%mix_cells, &
                              this%mix, this%mix_Epot)
                end if
            else if (iChangeRand <= this%num_moves + this%num_switches) then
                call switch(this%Box, &
                            this%type1_spheres, this%type1_macro, this%type1_obs, &
                            this%type2_spheres, this%type2_macro, this%type2_obs, &
                            this%mix, this%mix_Epot, &
                            this%switch_Nreject)
                this%switch_Nhit = this%switch_Nhit + 1
            else
                call rotate(this%Box, &
                            this%type1_spheres, this%type1_macro, this%type1_obs)
            end if
            
        end do MC_Change
    
    end subroutine Physical_System_random_changes
    
    subroutine Physical_System_update_rejections(this)    
        class(Physical_System), intent(inout) :: this
    
        call this%type1_obs%update_rejections()
        call this%type2_obs%update_rejections()
        
        this%switch_reject = real(this%switch_Nreject, DP)/real(this%switch_Nhit, DP)
        this%switch_Nreject = 0
        this%switch_Nhit = 0
        
    end subroutine Physical_System_update_rejections
    
    subroutine Physical_System_adapt_changes(this, iStep)    
        class(Physical_System), intent(inout) :: this
        integer, intent(in) :: iStep
        
        if (mod(iStep, this%period_adaptation) /= 0) then ! Rejections accumulation
            this%type1_obs%move_rejectAdapt = this%type1_obs%move_rejectAdapt + &
                                              this%type1_obs%move_reject
            this%type1_obs%rotate_rejectAdapt = this%type1_obs%rotate_rejectAdapt + &
                                                this%type1_obs%rotate_reject
            this%type2_obs%move_rejectAdapt = this%type2_obs%move_rejectAdapt + &
                                              this%type2_obs%move_reject
        else ! Average & adaptation
            call adapt_move(this%Box%size, &
                            this%type1_macro%move, &
                            this%period_adaptation, iStep, &
                            this%type1_obs, &
                            this%type1_units%move_delta)
            call adapt_rotation(this%type1_macro%rotation, &
                                this%period_adaptation, iStep, &
                                this%type1_obs, &
                                this%type1_units%rotate_delta)
            call adapt_move(this%Box%size, &
                            this%type2_macro%move, &
                            this%period_adaptation, iStep, &
                            this%type2_obs, &
                            this%type2_units%move_delta)
        end if
        
    end subroutine Physical_System_adapt_changes
    
    subroutine Physical_System_write_observables_thermalisation(this, iStep)    
        class(Physical_System), intent(inout) :: this
        integer, intent(in) :: iStep
    
        call this%type1_obs%write(iStep, this%type1_units%obsThermal)
        call this%type2_obs%write(iStep, this%type2_units%obsThermal)
        
        write(this%mix_obsThermal_unit, *) iStep, this%mix_Epot
        write(this%obsThermal_unit, *) iStep, this%type1_obs%Epot + this%type2_obs%Epot + this%mix_Epot
            
    end subroutine Physical_System_write_observables_thermalisation
    
    subroutine Physical_System_fix_changes(this)
    
        class(Physical_System), intent(inout) :: this
        
        character(len=5) :: type_name
        
        type_name = this%type1_spheres%get_name()
        call this%type1_macro%move%set_delta(type_name, this%Box%size, this%type1_obs%move_rejectAvg, &
                                       this%type1_units%report)
        call this%type1_macro%rotation%set_delta(type_name, this%type1_obs%rotate_rejectAvg, &
                                           this%type1_units%report)
        type_name = this%type2_spheres%get_name()
        call this%type2_macro%move%set_delta(type_name, this%Box%size, this%type2_obs%move_rejectAvg, &
                                       this%type2_units%report)
    
    end subroutine Physical_System_fix_changes
    
    subroutine Physical_System_measure_chemical_potentials(this)
    
        class(Physical_System), intent(inout) :: this
    
        call widom(this%Box, &
                   this%type1_spheres, this%type1_macro, this%type1_obs, &
                   this%type2_spheres, this%type2_macro%mix_cells, &
                   this%mix)
        call widom(this%Box, &
                   this%type2_spheres, this%type2_macro, this%type2_obs, &
                   this%type1_spheres, this%type1_macro%mix_cells, &
                   this%mix)
    
    end subroutine Physical_System_measure_chemical_potentials
    
    subroutine Physical_System_accumulate_observables(this)
    
        class(Physical_System), intent(inout) :: this
    
        call this%type1_obs%accumulate()
        call this%type2_obs%accumulate()
        
        this%mix_EpotSum = this%mix_EpotSum + this%mix_Epot
        this%switch_rejectSum = this%switch_rejectSum + this%switch_reject
            
    end subroutine Physical_System_accumulate_observables
    
    subroutine Physical_System_write_observables_equilibrium(this, iStep)
    
        class(Physical_System), intent(inout) :: this
        integer, intent(in) :: iStep
    
        call this%type1_obs%write(iStep, this%type1_units%obsEquilib)
        call this%type2_obs%write(iStep, this%type2_units%obsEquilib)
        
        write(this%mix_obsEquilib_unit, *) iStep, this%mix_Epot
        write(this%obsEquilib_unit, *) iStep, this%type1_obs%Epot + this%type2_obs%Epot + this%mix_Epot
            
    end subroutine Physical_System_write_observables_equilibrium
    
    subroutine Physical_System_take_snapshots(this, iStep)
    
        class(Physical_System), intent(inout) :: this
        integer, intent(in) :: iStep
        
        if (this%snap) then ! Snap shots of the configuration
            call this%type1_spheres%write_snap_positions(iStep, this%type1_units%snap_positions)
            call this%type1_spheres%write_snap_orientations(iStep, &
                                                            this%type1_units%snap_orientations)
            call this%type2_spheres%write_snap_positions(iStep, this%type2_units%snap_positions)
        end if
        
    end subroutine Physical_System_take_snapshots
    
    !> Reinitialize summed quantities to prevent them from drifting.
    subroutine Physical_System_reinitialize_quantites(this, iStep)
    
        class(Physical_System), intent(inout) :: this
        integer, intent(in) :: iStep
    
        if (modulo(iStep, this%reset_iStep) == 0) then
            call this%type1_macro%ewald_reci%reset_structure(this%Box, this%type1_spheres, iStep, &
                                                             this%type1_units%structure_modulus)
            call this%type1_macro%ewald_bound%reset_total_moment(this%type1_spheres, iStep, &
                 this%type1_units%totalMoment_modulus)
        end if
        
    end subroutine Physical_System_reinitialize_quantites

end module class_physical_system
