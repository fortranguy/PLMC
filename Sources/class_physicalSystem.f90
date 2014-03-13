!> \brief Description of the Physical System class

module class_physicalSystem

use, intrinsic :: iso_fortran_env, only: output_unit
use data_precisions, only: DP
use data_box, only: Ndim, Lsize, Kmax
use data_monteCarlo, only: Temperature, decorrelFactor, switch_factor, Nthermal, Nadapt, Nstep
use module_types, only: argument_seed, argument_initial
use class_hardSpheres
use class_dipolarSpheres
use class_mixingPotential
use class_observables
use class_units
use module_monteCarlo_arguments, only: read_arguments
use module_physics_macro, only: init_randomSeed, set_initialCondition
use module_tools, only: open_units, mix_open_units, print_report, init, final, mix_init, mix_final, &
                        test_consist, print_results

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
        integer :: Nmove, Nswitch, Nrotate !< changes
        
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
           
        real(DP) :: time_start, time_end
        
        ! Switch
        integer :: switch_Nhit, switch_Nreject
        real(DP) :: switch_reject, switch_rejectSum
    
    contains
    
        !> Construction and destruction of the class
        procedure :: set_box => PhysicalSystem_set_box
        procedure :: set_monteCarlo => PhysicalSystem_set_monteCarlo
        procedure :: set_changes => PhysicalSystem_set_changes
        procedure :: construct => PhysicalSystem_construct
        procedure :: destroy => PhysicalSystem_destroy
        
        !> Initialization                
        procedure :: open_units => PhysicalSystem_open_units
        procedure :: init_switch => PhysicalSystem_init_switch
        procedure :: init_spheres => PhysicalSystem_init_spheres
        procedure :: init_observables => PhysicalSystem_init_observables
        procedure :: init => PhysicalSystem_init
        
        !> Finalisation
        procedure :: final_spheres => PhysicalSystem_final_spheres
        procedure :: print_results => PhysicalSystem_print_results
        procedure :: final => PhysicalSystem_final
    
    end type PhysicalSystem
    
contains
    
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
    end subroutine PhysicalSystem_set_changes

    subroutine PhysicalSystem_construct(this)
    
        class(PhysicalSystem), intent(inout) :: this
        
        this%name = "sys"
        write(output_unit, *) this%name, " class construction"
        
        call this%set_box()
        call this%set_monteCarlo()
        
        call this%type1_spheres%construct()
        call this%type2_spheres%construct()
        call this%mix%construct(this%type1_spheres%get_sigma(), this%type2_spheres%get_sigma())
        
        call this%set_changes()
    
    end subroutine PhysicalSystem_construct
    
    subroutine PhysicalSystem_open_units(this)
    
        class(PhysicalSystem), intent(inout) :: this
        
        call open_units(this%report_unit, this%obsThermal_unit, this%obsEquilib_unit)        
        call print_report(this%Ncol, this%Nmove, this%Nswitch, this%Nrotate, this%report_unit)        
        
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
        
        write(output_unit, *) "Monte-Carlo Simulation: Canonical ensemble"        
        
        call this%open_units()        
        call this%init_switch()        
        
        call init_randomSeed(arg_seed, this%report_unit)
        call set_initialCondition(arg_init, this%type1_spheres, this%type2_spheres, &
                                  this%mix%get_sigma(), this%report_unit)
        call this%init_spheres()
        call this%init_observables()
    
    end subroutine PhysicalSystem_init
    
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
    
    subroutine PhysicalSystem_final(this)
    
        class(PhysicalSystem), intent(inout) :: this
        
        call this%final_spheres()
        call this%print_results()
    
    end subroutine PhysicalSystem_final    
    
    subroutine PhysicalSystem_destroy(this)
    
        class(PhysicalSystem), intent(inout) :: this
        
        write(output_unit, *) this%name, " class destruction"
    
    end subroutine PhysicalSystem_destroy

end module class_physicalSystem
