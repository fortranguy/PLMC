!> \brief Description of the Physical System class

module class_physicalSystem

use, intrinsic :: iso_fortran_env, only: output_unit
use data_precisions, only: DP
use data_box, only: Ndim, Lsize, Kmax
use data_monteCarlo, only: Temperature, decorrelFactor, switch_factor, Nthermal, Nadapt, Nstep, &
                           reset_iStep
use data_distribution, only: snap
use module_types, only: argument_seed, argument_initial
use class_hardSpheres
use class_dipolarSpheres
use class_mixingPotential
use class_observables
use class_units
use module_monteCarlo_arguments, only: read_arguments
use module_physics_macro, only: init_randomSeed, set_initialCondition
use module_algorithms, only: move, widom, switch, rotate
use module_tools, only: open_units, mix_open_units, print_report, init, final, mix_init, mix_final, &
                        test_consist, print_results, adapt_move, adapt_rotate

implicit none

private

    type, public :: PhysicalSystem
    
        private
        
        character(len=5) :: name
        
        ! Box
        real(DP), dimension(Ndim) :: Lsize !< box size
        integer, dimension(Ndim) :: Kmax !< number of wave vectors
        integer :: Ncol !< number of particles
        
        ! Monte-Carlo
        real(DP) :: Temperature
        integer :: Nthermal, Nadapt, Nstep !< Markov-chain Monte-Carlo
        integer :: Nchange, Nmove, Nswitch, Nrotate !< changes
        
        ! Type 1: Dipolar spheres
        type(DipolarSpheres) :: type1_spheres !< physical properties and Monte-Carlo subroutines
        type(MoreObservables) :: type1_obs !< e.g. energy, inverse of activity (-> chemical potential)
        type(MoreUnits) :: type1_units !< files units
        
        ! Type 2: Hard spheres
        type(HardSpheres) :: type2_spheres
        type(Observables) :: type2_obs
        type(Units) :: type2_units
        
        ! Mixing potential
        type(MixingPotential) :: mix !< short-range potential
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
        
        logical :: snap
        integer :: reset_iStep
        real(DP) :: time_start, time_end
    
    contains
    
        !> Construction & destruction of the class
        procedure, private :: set_box => PhysicalSystem_set_box
        procedure, private :: set_monteCarlo => PhysicalSystem_set_monteCarlo
        procedure :: construct => PhysicalSystem_construct
        procedure :: destroy => PhysicalSystem_destroy
        
        !> Initialization & Finalisation     
        procedure, private :: construct_spheres => PhysicalSystem_construct_spheres
        procedure, private :: set_changes => PhysicalSystem_set_changes    
        procedure, private :: open_units => PhysicalSystem_open_units
        procedure, private :: init_switch => PhysicalSystem_init_switch
        procedure, private :: init_spheres => PhysicalSystem_init_spheres
        procedure, private :: init_observables => PhysicalSystem_init_observables
        procedure :: init => PhysicalSystem_init
        procedure, private :: final_spheres => PhysicalSystem_final_spheres
        procedure, private :: print_results => PhysicalSystem_print_results
        procedure, private :: close_units => PhysicalSystem_close_units
        procedure :: final => PhysicalSystem_final
        
        !> Accessors & Mutators
        procedure :: get_Nthermal => PhysicalSystem_get_Nthermal
        procedure :: get_Nstep => PhysicalSystem_get_Nstep
        procedure :: set_time_start => PhysicalSystem_set_time_start
        procedure :: set_time_end => PhysicalSystem_set_time_end
        
        !> Simulation
        procedure :: random_changes => PhysicalSystem_random_changes
        procedure :: update_rejections => PhysicalSystem_update_rejections
        procedure :: adapt_changes => PhysicalSystem_adapt_changes
        procedure :: write_observables_thermalisation => PhysicalSystem_write_observables_thermalisation
        procedure :: fix_changes => PhysicalSystem_fix_changes
        procedure :: measure_chemical_potentials => PhysicalSystem_measure_chemical_potentials
        procedure :: accumulate_observables => PhysicalSystem_accumulate_observables
        procedure :: write_observables_equilibrium => PhysicalSystem_write_observables_equilibrium
        procedure :: take_snapshots => PhysicalSystem_take_snapshots
        procedure :: reinitialize_quantites => PhysicalSystem_reinitialize_quantites
    
    end type PhysicalSystem
    
contains

    ! Construction
    
    pure subroutine PhysicalSystem_set_box(this)    
        class(PhysicalSystem), intent(inout) :: this        
        this%Lsize(:) = Lsize(:)
        this%Kmax(:) = Kmax(:)        
    end subroutine PhysicalSystem_set_box
    
    pure subroutine PhysicalSystem_set_monteCarlo(this)    
        class(PhysicalSystem), intent(inout) :: this        
        this%Temperature = Temperature
        this%Nthermal = Nthermal
        this%Nadapt = Nadapt
        this%Nstep = Nstep    
    end subroutine PhysicalSystem_set_monteCarlo
    
    pure subroutine PhysicalSystem_set_changes(this)    
        class(PhysicalSystem), intent(inout) :: this    
        this%Ncol = this%type1_spheres%get_Ncol() + this%type2_spheres%get_Ncol()
        this%Nmove = decorrelFactor * this%Ncol
        this%Nswitch = switch_factor * decorrelFactor * this%type1_spheres%get_Ncol()
        this%Nrotate = decorrelFactor * this%type1_spheres%get_Ncol()    
        this%Nchange = this%Nmove + this%Nswitch + this%Nrotate
    end subroutine PhysicalSystem_set_changes

    subroutine PhysicalSystem_construct(this)
    
        class(PhysicalSystem), intent(out) :: this
        
        this%name = "sys"
        write(output_unit, *) this%name, " class construction"
        
        call this%set_box()
        call this%set_monteCarlo()
    
    end subroutine PhysicalSystem_construct
    
    ! Initialization
    
    subroutine PhysicalSystem_construct_spheres(this)
    
        class(PhysicalSystem), intent(out) :: this
    
        call this%type1_spheres%construct()
        call this%type2_spheres%construct()
        call this%mix%construct(this%type1_spheres%get_sigma(), this%type2_spheres%get_sigma())
    
    end subroutine PhysicalSystem_construct_spheres
    
    subroutine PhysicalSystem_open_units(this)
    
        class(PhysicalSystem), intent(inout) :: this
        
        call open_units(this%report_unit, this%obsThermal_unit, this%obsEquilib_unit)        
        call print_report(this%Ncol, this%Nmove, this%Nswitch, this%Nrotate, this%reset_iStep, &
                          this%report_unit)        
        
        call this%type1_units%open(this%type1_spheres%get_name())
        call this%type1_spheres%print_density(this%Ncol, this%type1_units%report)        
        
        call this%type2_units%open(this%type2_spheres%get_name())
        call this%type2_spheres%print_density(this%Ncol, this%type2_units%report)    
        
        call mix_open_units(this%mix_report_unit, this%mix_Epot_tab_unit, this%mix_obsThermal_unit, &
                            this%mix_obsEquilib_unit)
        call this%mix%print_report(this%mix_report_unit)
    
    end subroutine PhysicalSystem_open_units
    
    pure subroutine PhysicalSystem_init_switch(this)    
        class(PhysicalSystem), intent(inout) :: this    
        this%switch_Nhit = 0
        this%switch_Nreject = 0
        this%switch_reject = 0._DP
        this%switch_rejectSum = 0._DP        
    end subroutine PhysicalSystem_init_switch
    
    subroutine PhysicalSystem_init_spheres(this)
    
        class(PhysicalSystem), intent(out) :: this
    
        call mix_init(this%mix, this%type1_spheres, this%type2_spheres, this%mix_Epot_tab_unit, &
                      this%mix_Epot)
        call this%mix%print_report(this%mix_report_unit)
        call init(this%type1_spheres, this%type2_spheres, this%mix, this%type1_units, &
                  this%type1_obs%Epot)
        call init(this%type2_spheres, this%type1_spheres, this%mix, this%type2_units, &
                  this%type2_obs%Epot)
        
    end subroutine PhysicalSystem_init_spheres
    
    subroutine PhysicalSystem_init_observables(this)
    
        class(PhysicalSystem), intent(inout) :: this
        
        real(DP) :: Epot_conf
        
        call this%type1_obs%init()
        call this%type2_obs%init()
        this%mix_EpotSum = 0._DP
        this%EpotSum = 0._DP
        
        Epot_conf = this%type1_obs%Epot + this%type2_obs%Epot + this%mix_Epot
        write(output_unit, *) "Initial potential energy =", Epot_conf
        write(this%obsThermal_unit, *) 0, Epot_conf
        
    end subroutine PhysicalSystem_init_observables    

    subroutine PhysicalSystem_init(this, arg_seed, arg_init)
    
        class(PhysicalSystem), intent(out) :: this
        type(argument_seed), intent(in) :: arg_seed
        type(argument_initial), intent(in) :: arg_init
        
        call this%construct_spheres()
        call this%set_changes()
        
        write(output_unit, *) "Monte-Carlo Simulation: Canonical ensemble"        
        
        call this%open_units()        
        call this%init_switch()        
        
        call init_randomSeed(arg_seed, this%report_unit)
        call set_initialCondition(arg_init, this%type1_spheres, this%type2_spheres, &
                                  this%mix%get_sigma(), this%report_unit)
        call this%init_spheres()
        call this%init_observables()
        
        this%snap = snap
        this%reset_iStep = reset_iStep
    
    end subroutine PhysicalSystem_init
    
    ! Finalisation
    
    subroutine PhysicalSystem_final_spheres(this)
    
        class(PhysicalSystem), intent(inout) :: this
        
        call final(this%type1_spheres, this%type1_units, this%type1_obs)
        call final(this%type2_spheres, this%type2_units, this%type2_obs)
        call mix_final(this%mix, this%type1_spheres, this%type2_spheres, this%mix_report_unit, &
                       this%mix_Epot, this%mix_EpotSum, this%mix_Epot_conf)
        
    end subroutine PhysicalSystem_final_spheres
    
    subroutine PhysicalSystem_print_results(this)
    
        class(PhysicalSystem), intent(inout) :: this
        
        real(DP) :: Epot, Epot_conf
        real(DP) :: duration
        
        Epot = this%type1_obs%Epot + this%type2_obs%Epot + this%mix_Epot
        Epot_conf = this%type1_spheres%Epot_conf() + this%type2_spheres%Epot_conf() + &
                    this%mix_Epot_conf
        call test_consist(Epot, Epot_conf, this%report_unit)
        this%EpotSum = this%type1_obs%EpotSum + this%type2_obs%EpotSum + this%mix_EpotSum
        duration = this%time_end - this%time_start
        call print_results(this%Ncol, this%EpotSum, this%switch_rejectSum, duration, this%report_unit)
    
    end subroutine PhysicalSystem_print_results
    
    subroutine PhysicalSystem_close_units(this)
    
        class(PhysicalSystem), intent(inout) :: this
    
        call this%type1_units%close()
        call this%type2_units%close()
        
        close(this%mix_report_unit)
        close(this%mix_Epot_tab_unit)
        close(this%mix_obsThermal_unit)
        close(this%mix_obsEquilib_unit)
        
        close(this%report_unit)
        close(this%obsThermal_unit)
        close(this%obsEquilib_unit)
    
    end subroutine PhysicalSystem_close_units
    
    subroutine PhysicalSystem_final(this)
    
        class(PhysicalSystem), intent(inout) :: this
        
        call this%final_spheres()
        call this%print_results()
        call this%close_units()
    
    end subroutine PhysicalSystem_final
    
    ! Destruction
    
    subroutine PhysicalSystem_destroy(this)
    
        class(PhysicalSystem), intent(inout) :: this
        
        write(output_unit, *) this%name, " class destruction"
        
        call this%type1_spheres%destroy()
        call this%type2_spheres%destroy()
        call this%mix%destroy()
    
    end subroutine PhysicalSystem_destroy
    
    ! Accessors
    
    pure function PhysicalSystem_get_Nthermal(this) result(get_Nthermal)
        class(PhysicalSystem), intent(in) :: this
        integer :: get_Nthermal
        get_Nthermal = this%Nthermal
    end function PhysicalSystem_get_Nthermal
    
    pure function PhysicalSystem_get_Nstep(this) result(get_Nstep)
        class(PhysicalSystem), intent(in) :: this
        integer :: get_Nstep
        get_Nstep = this%Nstep
    end function PhysicalSystem_get_Nstep
    
    ! Mutators
    
    subroutine PhysicalSystem_set_time_start(this)
        class(PhysicalSystem), intent(inout) :: this
        call cpu_time(this%time_start)
    end subroutine PhysicalSystem_set_time_start
    
    subroutine PhysicalSystem_set_time_end(this)
        class(PhysicalSystem), intent(inout) :: this
        call cpu_time(this%time_end)
    end subroutine PhysicalSystem_set_time_end        
    
    ! Random changes
    
    subroutine PhysicalSystem_random_changes(this)
    
        class(PhysicalSystem), intent(inout) :: this
        
        integer :: iChange, iChangeRand, iColRand
        real(DP) :: rand
        
        MC_Change: do iChange = 1, this%Nchange
        
            ! Randomly choosing the change
            call random_number(rand)
            iChangeRand = int(rand*real(this%Nchange, DP)) + 1
            
            if (iChangeRand <= this%Nmove) then            
                ! Randomly choosing the type
                call random_number(rand)
                iColRand = int(rand*real(this%Ncol, DP)) + 1                
                if (iColRand <= this%type1_spheres%get_Ncol()) then
                    call move(this%type1_spheres, this%type1_obs, this%type2_spheres, this%mix, &
                              this%mix_Epot)
                else
                    call move(this%type2_spheres, this%type2_obs, this%type1_spheres, this%mix, &
                              this%mix_Epot)
                end if                
            else if (iChangeRand <= this%Nmove + this%Nswitch) then            
                call switch(this%type1_spheres, this%type1_obs, this%type2_spheres, this%type2_obs, &
                            this%mix, this%mix_Epot, this%switch_Nreject)
                this%switch_Nhit = this%switch_Nhit + 1                
            else     
                call rotate(this%type1_spheres, this%type1_obs)                
            end if
            
        end do MC_Change
    
    end subroutine PhysicalSystem_random_changes
    
    subroutine PhysicalSystem_update_rejections(this)
    
        class(PhysicalSystem), intent(inout) :: this
    
        call this%type1_obs%update_rejections()
        call this%type2_obs%update_rejections()
        
        this%switch_reject = real(this%switch_Nreject, DP)/real(this%switch_Nhit, DP)
        this%switch_Nreject = 0
        this%switch_Nhit = 0 
        
    end subroutine PhysicalSystem_update_rejections
    
    subroutine PhysicalSystem_adapt_changes(this, iStep)
    
        class(PhysicalSystem), intent(inout) :: this
        integer, intent(in) :: iStep
        
        if (mod(iStep, this%Nadapt) /= 0) then ! Rejections accumulation
            this%type1_obs%move_rejectAdapt = this%type1_obs%move_rejectAdapt + &
                                              this%type1_obs%move_reject
            this%type1_obs%rotate_rejectAdapt = this%type1_obs%rotate_rejectAdapt + &
                                                this%type1_obs%rotate_reject
            this%type2_obs%move_rejectAdapt = this%type2_obs%move_rejectAdapt + &
                                              this%type2_obs%move_reject
        else ! Average & adaptation
            call adapt_move(this%type1_spheres, iStep, this%type1_obs, this%type1_units%move_delta)
            call adapt_rotate(this%type1_spheres, iStep, this%type1_obs, this%type1_units%rotate_delta)
            call adapt_move(this%type2_spheres, iStep, this%type2_obs, this%type2_units%move_delta)
        end if
        
    end subroutine PhysicalSystem_adapt_changes
    
    subroutine PhysicalSystem_write_observables_thermalisation(this, iStep)
    
        class(PhysicalSystem), intent(inout) :: this
        integer, intent(in) :: iStep
    
        call this%type1_obs%write(iStep, this%type1_units%obsThermal)
        call this%type2_obs%write(iStep, this%type2_units%obsThermal)
        
        write(this%mix_obsThermal_unit, *) iStep, this%mix_Epot
        write(this%obsThermal_unit, *) iStep, this%type1_obs%Epot + this%type2_obs%Epot + this%mix_Epot
            
    end subroutine PhysicalSystem_write_observables_thermalisation
    
    subroutine PhysicalSystem_fix_changes(this)
    
        class(PhysicalSystem), intent(inout) :: this
        
        call this%type1_spheres%set_move_delta(this%type1_obs%move_rejectAvg, this%type1_units%report)
        call this%type1_spheres%set_rotate_delta(this%type1_obs%rotate_rejectAvg, &
                                                 this%type1_units%report)
        call this%type2_spheres%set_move_delta(this%type2_obs%move_rejectAvg, this%type2_units%report)
    
    end subroutine PhysicalSystem_fix_changes
    
    subroutine PhysicalSystem_measure_chemical_potentials(this)
    
        class(PhysicalSystem), intent(inout) :: this
    
        call widom(this%type1_spheres, this%type1_obs, this%type2_spheres, this%mix)
        call widom(this%type2_spheres, this%type2_obs, this%type1_spheres, this%mix)
    
    end subroutine PhysicalSystem_measure_chemical_potentials
    
    subroutine PhysicalSystem_accumulate_observables(this)
    
        class(PhysicalSystem), intent(inout) :: this
    
        call this%type1_obs%accumulate()
        call this%type2_obs%accumulate()
        
        this%mix_EpotSum = this%mix_EpotSum + this%mix_Epot
        this%switch_rejectSum = this%switch_rejectSum + this%switch_reject
            
    end subroutine PhysicalSystem_accumulate_observables
    
    subroutine PhysicalSystem_write_observables_equilibrium(this, iStep)
    
        class(PhysicalSystem), intent(inout) :: this
        integer, intent(in) :: iStep
    
        call this%type1_obs%write(iStep, this%type1_units%obsEquilib)
        call this%type2_obs%write(iStep, this%type2_units%obsEquilib)
        
        write(this%mix_obsEquilib_unit, *) iStep, this%mix_Epot
        write(this%obsEquilib_unit, *) iStep, this%type1_obs%Epot + this%type2_obs%Epot + this%mix_Epot
            
    end subroutine PhysicalSystem_write_observables_equilibrium
    
    subroutine PhysicalSystem_take_snapshots(this, iStep)
    
        class(PhysicalSystem), intent(inout) :: this
        integer, intent(in) :: iStep
        
        if (this%snap) then ! Snap shots of the configuration
            call this%type1_spheres%snap_positions(iStep, this%type1_units%snap_positions)
            call this%type1_spheres%snap_orientations(iStep, this%type1_units%snap_orientations)
            call this%type2_spheres%snap_positions(iStep, this%type2_units%snap_positions)
        end if
        
    end subroutine PhysicalSystem_take_snapshots
    
    !> Reinitialize summed quantities to prevent them from drifting.
    subroutine PhysicalSystem_reinitialize_quantites(this, iStep)
    
        class(PhysicalSystem), intent(inout) :: this
        integer, intent(in) :: iStep
    
        if (modulo(iStep, this%reset_iStep) == 0) then
            call this%type1_spheres%reset_Epot_reci_structure(iStep, this%type1_units%structure_modulus)
            call this%type1_spheres%reset_totalMoment(iStep, this%type1_units%totalMoment_modulus)
        end if
        
    end subroutine PhysicalSystem_reinitialize_quantites

end module class_physicalSystem
