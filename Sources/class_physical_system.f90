!> \brief Description of the Physical System class

module class_physical_system

use, intrinsic :: iso_fortran_env, only: output_unit
use data_precisions, only: DP
use data_box, only: Ndim, Box_size, Box_wave
use data_monte_carlo, only: Temperature, decorrelFactor, switch_factor, Nthermal, Nadapt, Nstep, &
                            reset_iStep, &
                            dipol_move_delta, dipol_move_rejectFix, &
                            dipol_rotate_delta, dipol_rotate_deltaMax, dipol_rotate_rejectFix, &
                            hard_move_delta, hard_move_rejectFix
use data_potential, only: write_potential
use data_distribution, only: snap
use module_types, only: Box_Dimensions, Monte_Carlo_Arguments
use module_physics_micro, only: NwaveVectors
use class_hard_spheres
use class_small_move
use class_small_rotation
use class_dipolar_spheres
use class_mixing_potential
use class_observables
use class_units
use module_monte_carlo_arguments, only: read_arguments
use module_physics_macro, only: init_randomSeed, set_initialConfiguration, init, final, mix_init, &
                                mix_final, adapt_move, adapt_rotation, test_consist
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
        real(DP) :: Temperature
        integer :: Nthermal, Nadapt, Nstep !< Markov-chain Monte-Carlo
        integer :: decorrelFactor, Nchange, Nmove, Nswitch, Nrotate !< changes
        
        ! Type 1: Dipolar spheres
        type(Dipolar_Spheres) :: type1_spheres
        type(Small_Move) :: type1_move
        type(Small_Rotation) :: type1_rotation
        type(MoreObservables) :: type1_obs !< e.g. energy, inverse of activity (-> chemical potential)
        type(MoreUnits) :: type1_units !< files units
        
        ! Type 2: Hard spheres
        type(Hard_Spheres) :: type2_spheres
        type(Small_Move) :: type2_move
        type(Observables) :: type2_obs
        type(Units) :: type2_units
        
        ! Mixing potential
        type(Mixing_Potential) :: mix !< short-range potential
        real(DP) :: mix_Epot, mix_EpotSum, mix_Epot_conf !< potential energy
        integer :: mix_report_unit
        integer :: mix_Epot_tab_unit
        integer :: mix_obsThermal_unit, mix_obsEquilib_unit
        
        ! Observables: write to files
        real(DP) :: Epot, EpotSum
        integer :: report_unit !< data & results file
        integer :: obsThermal_unit, obsEquilib_unit !< observables files: thermalisation & equilibrium
        
        ! Switch
        integer :: switch_Nhit, switch_Nreject
        real(DP) :: switch_reject, switch_rejectSum
        
        logical :: write_potential, snap
        integer :: reset_iStep
        real(DP) :: time_start, time_end
    
    contains
    
        !> Construction & destruction of the class
        procedure, private :: set_box => Physical_System_set_box
        procedure, private :: set_monteCarlo => Physical_System_set_monteCarlo
        procedure, private :: set_changes => Physical_System_set_changes
        procedure :: construct => Physical_System_construct
        procedure :: destroy => Physical_System_destroy
        
        !> Initialization & Finalisation
        procedure, private :: open_units => Physical_System_open_units
        procedure, private :: init_switch => Physical_System_init_switch
        procedure, private :: init_observables => Physical_System_init_observables
        procedure, private :: init_spheres => Physical_System_init_spheres
        procedure :: write_report => Physical_System_write_report
        procedure :: init => Physical_System_init
        procedure, private :: final_spheres => Physical_System_final_spheres
        procedure, private :: write_results => Physical_System_write_results
        procedure, private :: close_units => Physical_System_close_units
        procedure :: final => Physical_System_final
        
        !> Accessors & Mutators
        procedure :: get_Nthermal => Physical_System_get_Nthermal
        procedure :: get_Nstep => Physical_System_get_Nstep
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
    
    pure subroutine Physical_System_set_box(this)
        class(Physical_System), intent(inout) :: this
        
        this%Box%size(:) = Box_size(:)
        this%Box%wave(:) = Box_wave(:)
    end subroutine Physical_System_set_box
    
    pure subroutine Physical_System_set_monteCarlo(this)
        class(Physical_System), intent(inout) :: this
        
        this%Temperature = Temperature
        this%Nthermal = Nthermal
        this%Nadapt = Nadapt
        this%Nstep = Nstep
    end subroutine Physical_System_set_monteCarlo
    
    pure subroutine Physical_System_set_changes(this)
        class(Physical_System), intent(inout) :: this
        
        this%num_particles = this%type1_spheres%get_num_particles() + &
                             this%type2_spheres%get_num_particles()
        this%decorrelFactor = decorrelFactor
        this%Nmove = this%decorrelFactor * this%num_particles
        this%Nswitch = switch_factor * this%decorrelFactor * this%type1_spheres%get_num_particles()
        this%Nrotate = this%decorrelFactor * this%type1_spheres%get_num_particles()
        this%Nchange = this%Nmove + this%Nswitch + this%Nrotate
        
        call this%type1_move%init(dipol_move_delta, dipol_move_rejectFix) ! ugly
        call this%type1_rotation%init(dipol_rotate_delta, dipol_rotate_deltaMax, &
                                      dipol_rotate_rejectFix) ! ugly
        call this%type2_move%init(hard_move_delta, hard_move_rejectFix) ! ugly
    end subroutine Physical_System_set_changes

    subroutine Physical_System_construct(this)
    
        class(Physical_System), intent(out) :: this
        
        this%name = "sys"
        write(output_unit, *) this%name, " class construction"
        
        call this%set_box()
        call this%set_monteCarlo()
        this%write_potential = write_potential
        
        call this%type1_spheres%construct()
        call this%type2_spheres%construct()
        call this%mix%construct(this%type1_spheres%get_diameter(), this%type2_spheres%get_diameter())
        
        call this%set_changes()
    
    end subroutine Physical_System_construct
    
    ! Initialization
    
    subroutine Physical_System_open_units(this)
    
        class(Physical_System), intent(inout) :: this
        
        call open_units(this%report_unit, this%obsThermal_unit, this%obsEquilib_unit)
        call write_data(this%report_unit)
        call this%write_report(this%report_unit)
        
        call this%type1_units%open(this%type1_spheres%get_name())
        call this%type1_spheres%write_density(this%Box%size, this%num_particles, this%type1_units%report)
        
        call this%type2_units%open(this%type2_spheres%get_name())
        call this%type2_spheres%write_density(this%Box%size, this%num_particles, this%type2_units%report)
        
        call mix_open_units(this%mix_report_unit, this%mix_Epot_tab_unit, this%mix_obsThermal_unit, &
                            this%mix_obsEquilib_unit)
        call this%mix%write_report(this%mix_report_unit)
    
    end subroutine Physical_System_open_units
    
    pure subroutine Physical_System_init_switch(this)
        class(Physical_System), intent(inout) :: this
        this%switch_Nhit = 0
        this%switch_Nreject = 0
        this%switch_reject = 0._DP
        this%switch_rejectSum = 0._DP
    end subroutine Physical_System_init_switch
    
    subroutine Physical_System_init_observables(this)
    
        class(Physical_System), intent(inout) :: this
        
        call this%type1_obs%init()
        call this%type2_obs%init()
        this%mix_EpotSum = 0._DP
        this%EpotSum = 0._DP
        
    end subroutine Physical_System_init_observables
    
    subroutine Physical_System_init_spheres(this)
    
        class(Physical_System), intent(inout) :: this
    
        call mix_init(this%Box%size, this%mix, this%type1_spheres, this%type2_spheres, &
                      this%write_potential, this%mix_Epot_tab_unit, this%mix_Epot)
        call this%mix%write_report(this%mix_report_unit)
        call init(this%Box, this%type1_spheres, this%type2_spheres, this%mix, &
                  this%write_potential, this%type1_units, this%type1_obs%Epot)
        call init(this%Box, this%type2_spheres, this%type1_spheres, this%mix, &
                  this%write_potential, this%type2_units, this%type2_obs%Epot)
        
    end subroutine Physical_System_init_spheres

    subroutine Physical_System_init(this, args)
    
        class(Physical_System), intent(inout) :: this
        type(Monte_Carlo_Arguments), intent(in) :: args
        
        real(DP) :: Epot_conf
        
        this%snap = snap
        this%reset_iStep = reset_iStep
        
        write(output_unit, *) "Monte-Carlo Simulation: Canonical ensemble"
        
        call this%open_units()
        call this%init_switch()
        
        call init_randomSeed(args%random, this%report_unit)
        call set_initialConfiguration(this%Box%size, args%initial, this%type1_spheres, &
                                      this%type2_spheres, this%mix%get_min_distance(), this%report_unit)
        call this%init_observables()
        call this%init_spheres()
        
        Epot_conf = this%type1_obs%Epot + this%type2_obs%Epot + this%mix_Epot
        write(output_unit, *) "Initial potential energy =", Epot_conf
        write(this%obsThermal_unit, *) 0, Epot_conf
    
    end subroutine Physical_System_init
    
    subroutine Physical_System_write_report(this, report_unit)
    
        class(Physical_System), intent(in) :: this
        integer, intent(in) :: report_unit

        write(report_unit, *) "Data macro: "
        
        write(report_unit ,*) "    Box_size(:) = ", this%Box%size(:)
        write(report_unit ,*) "    Volume = ", product(this%Box%size)
        write(report_unit ,*) "    Box_wave(:) = ", this%Box%wave(:)
        write(report_unit ,*) "    NwaveVectors =", NwaveVectors(this%Box%wave)
        write(report_unit ,*) "    num_particles = ", this%num_particles
        write(report_unit ,*) "    Temperature = ", this%Temperature
        
        write(report_unit, *) "    Nstep = ", this%Nstep
        write(report_unit, *) "    Nthermal = ", this%Nthermal
        write(report_unit, *) "    decorrelFactor = ", this%decorrelFactor
        write(report_unit, *) "    Nmove = ", this%Nmove
        write(report_unit, *) "    Nswitch = ", this%Nswitch
        write(report_unit, *) "    Nrotate = ", this%Nrotate
        
        write(report_unit, *) "    reset_iStep = ", this%reset_iStep
        write(report_unit, *) "    write_potential = ", this%write_potential
    
    end subroutine Physical_System_write_report
    
    ! Finalisation
    
    subroutine Physical_System_final_spheres(this)
    
        class(Physical_System), intent(inout) :: this
        
        call final(this%Box, this%type1_spheres, this%type1_units, this%type1_obs)
        call final(this%Box, this%type2_spheres, this%type2_units, this%type2_obs)
        call mix_final(this%Box%size, this%mix, this%type1_spheres, this%type2_spheres, &
                       this%mix_report_unit, this%mix_Epot, this%mix_Epot_conf)
        call mix_write_results(this%Nstep, this%mix_EpotSum, this%mix_report_unit)
        
    end subroutine Physical_System_final_spheres
    
    subroutine Physical_System_write_results(this)
    
        class(Physical_System), intent(inout) :: this
        
        real(DP) :: Epot, Epot_conf
        real(DP) :: duration
        
        Epot = this%type1_obs%Epot + this%type2_obs%Epot + this%mix_Epot
        Epot_conf = this%type1_spheres%Epot_conf(this%Box) + &
                    this%type2_spheres%Epot_conf(this%Box) + this%mix_Epot_conf
        call test_consist(Epot, Epot_conf, this%report_unit)
        this%EpotSum = this%type1_obs%EpotSum + this%type2_obs%EpotSum + this%mix_EpotSum
        duration = this%time_end - this%time_start
        call write_results(this%num_particles, this%Nstep, this%EpotSum, this%switch_rejectSum, duration,&
                           this%report_unit)
    
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
    
    subroutine Physical_System_final(this)
    
        class(Physical_System), intent(inout) :: this
        
        call this%final_spheres()
        call this%write_results()
        call this%close_units()
    
    end subroutine Physical_System_final
    
    ! Destruction
    
    subroutine Physical_System_destroy(this)
    
        class(Physical_System), intent(inout) :: this
        
        write(output_unit, *) this%name, " class destruction"
        
        call this%type1_spheres%destroy()
        call this%type2_spheres%destroy()
        call this%mix%destroy()
    
    end subroutine Physical_System_destroy
    
    ! Accessors
    
    pure function Physical_System_get_Nthermal(this) result(get_Nthermal)
        class(Physical_System), intent(in) :: this
        integer :: get_Nthermal
        get_Nthermal = this%Nthermal
    end function Physical_System_get_Nthermal
    
    pure function Physical_System_get_Nstep(this) result(get_Nstep)
        class(Physical_System), intent(in) :: this
        integer :: get_Nstep
        get_Nstep = this%Nstep
    end function Physical_System_get_Nstep
    
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
        
        MC_Change: do iChange = 1, this%Nchange
        
            ! Randomly choosing the change
            call random_number(rand)
            iChangeRand = int(rand*real(this%Nchange, DP)) + 1
            
            if (iChangeRand <= this%Nmove) then
                ! Randomly choosing the type
                call random_number(rand)
                iColRand = int(rand*real(this%num_particles, DP)) + 1
                if (iColRand <= this%type1_spheres%get_num_particles()) then
                    call move(this%Box, this%type1_spheres, this%type1_move, this%type1_obs, &
                              this%type2_spheres, this%mix, this%mix_Epot)
                else
                    call move(this%Box, this%type2_spheres, this%type2_move, this%type2_obs, &
                              this%type1_spheres, this%mix, this%mix_Epot)
                end if
            else if (iChangeRand <= this%Nmove + this%Nswitch) then
                call switch(this%Box, this%type1_spheres, this%type1_obs, this%type2_spheres, &
                            this%type2_obs, this%mix, this%mix_Epot, this%switch_Nreject)
                this%switch_Nhit = this%switch_Nhit + 1
            else
                call rotate(this%Box, this%type1_spheres, this%type1_rotation, this%type1_obs)
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
        
        if (mod(iStep, this%Nadapt) /= 0) then ! Rejections accumulation
            this%type1_obs%move_rejectAdapt = this%type1_obs%move_rejectAdapt + &
                                              this%type1_obs%move_reject
            this%type1_obs%rotate_rejectAdapt = this%type1_obs%rotate_rejectAdapt + &
                                                this%type1_obs%rotate_reject
            this%type2_obs%move_rejectAdapt = this%type2_obs%move_rejectAdapt + &
                                              this%type2_obs%move_reject
        else ! Average & adaptation
            call adapt_move(this%Box%size, this%type1_move, this%Nadapt, iStep, this%type1_obs, &
                            this%type1_units%move_delta)
            call adapt_rotation(this%type1_rotation, this%Nadapt, iStep, this%type1_obs, &
                                this%type1_units%rotate_delta)
            call adapt_move(this%Box%size, this%type2_move, this%Nadapt, iStep, this%type2_obs, &
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
        call this%type1_move%set_delta(type_name, this%Box%size, this%type1_obs%move_rejectAvg, &
                                       this%type1_units%report)
        call this%type1_rotation%set_delta(type_name, this%type1_obs%rotate_rejectAvg, &
                                           this%type1_units%report)
        type_name = this%type2_spheres%get_name()
        call this%type2_move%set_delta(type_name, this%Box%size, this%type2_obs%move_rejectAvg, &
                                       this%type2_units%report)
    
    end subroutine Physical_System_fix_changes
    
    subroutine Physical_System_measure_chemical_potentials(this)
    
        class(Physical_System), intent(inout) :: this
    
        call widom(this%Box, this%type1_spheres, this%type1_obs, this%type2_spheres, this%mix)
        call widom(this%Box, this%type2_spheres, this%type2_obs, this%type1_spheres, this%mix)
    
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
            call this%type1_spheres%snap_positions(iStep, this%type1_units%snap_positions)
            call this%type1_spheres%snap_orientations(iStep, this%type1_units%snap_orientations)
            call this%type2_spheres%snap_positions(iStep, this%type2_units%snap_positions)
        end if
        
    end subroutine Physical_System_take_snapshots
    
    !> Reinitialize summed quantities to prevent them from drifting.
    subroutine Physical_System_reinitialize_quantites(this, iStep)
    
        class(Physical_System), intent(inout) :: this
        integer, intent(in) :: iStep
    
        if (modulo(iStep, this%reset_iStep) == 0) then
            call this%type1_spheres%reset_Epot_reci_structure(this%Box, iStep, &
                                                              this%type1_units%structure_modulus)
            call this%type1_spheres%reset_totalMoment(iStep, this%type1_units%totalMoment_modulus)
        end if
        
    end subroutine Physical_System_reinitialize_quantites

end module class_physical_system