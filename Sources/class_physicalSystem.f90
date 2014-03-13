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
use module_tools, only: open_units, mix_open_units, print_report, mix_init, init

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
        real(DP) :: mix_Epot, mix_EpotSum !< potential energy
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
        procedure :: construct => PhysicalSystem_construct
        procedure :: destroy => PhysicalSystem_destroy
        
        !> Initialization
        procedure :: init_box => PhysicalSystem_init_box
        procedure :: init_monteCarlo => PhysicalSystem_init_monteCarlo
        procedure :: init_switch => PhysicalSystem_init_switch
        procedure :: init_changes => PhysicalSystem_init_changes
        procedure :: open_units => PhysicalSystem_open_units
        procedure :: init_spheres => PhysicalSystem_init_spheres
        procedure :: init_observables => PhysicalSystem_init_observables
        procedure :: init => PhysicalSystem_init
    
    end type PhysicalSystem
    
contains

    subroutine PhysicalSystem_construct(this)
    
        class(PhysicalSystem), intent(inout) :: this
        
        this%name = "sys"
        write(output_unit, *) this%name, " class construction"
        
        call this%type1_spheres%construct()
        call this%type2_spheres%construct()
        call this%mix%construct(this%type1_spheres%get_sigma(), this%type2_spheres%get_sigma())
    
    end subroutine PhysicalSystem_construct

    pure subroutine PhysicalSystem_init_box(this)
    
        class(PhysicalSystem), intent(inout) :: this
        
        this%Lsize(:) = Lsize(:)
        this%Kmax(:) = Kmax(:)
        
    end subroutine PhysicalSystem_init_box
    
    pure subroutine PhysicalSystem_init_monteCarlo(this)
    
        class(PhysicalSystem), intent(inout) :: this
        
        this%Temperature = Temperature
        this%Nthermal = Nthermal
        this%Nadapt = Nadapt
        this%Nstep = Nstep
    
    end subroutine PhysicalSystem_init_monteCarlo
    
    pure subroutine PhysicalSystem_init_switch(this)
    
        class(PhysicalSystem), intent(inout) :: this
    
        this%switch_Nhit = 0
        this%switch_Nreject = 0
        this%switch_reject = 0._DP
        this%switch_rejectSum = 0._DP
        
    end subroutine PhysicalSystem_init_switch
    
    pure subroutine PhysicalSystem_init_changes(this)
    
        class(PhysicalSystem), intent(inout) :: this
    
        this%Ncol = this%type1_spheres%get_Ncol() + this%type2_spheres%get_Ncol()
        this%Nmove = decorrelFactor * this%Ncol
        this%Nswitch = switch_factor * decorrelFactor * this%type1_spheres%get_Ncol()
        this%Nrotate = decorrelFactor * this%type1_spheres%get_Ncol()
    
    end subroutine PhysicalSystem_init_changes
    
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
        
        call this%init_box()
        call this%init_monteCarlo()
        call this%init_changes()
        
        write(output_unit, *) "Monte-Carlo Simulation: Canonical ensemble"        
        
        call this%open_units()        
        call this%init_switch()        
        
        call init_randomSeed(arg_seed, this%report_unit)
        call set_initialCondition(arg_init, this%type1_spheres, this%type2_spheres, &
                                  this%mix%get_sigma(), this%report_unit)
        call this%init_spheres()
        call this%init_observables()
    
    end subroutine PhysicalSystem_init
    
    subroutine PhysicalSystem_final(this)
    
        class(PhysicalSystem), intent(inout) :: this
    
    end subroutine PhysicalSystem_final    
    
    subroutine PhysicalSystem_destroy(this)
    
        class(PhysicalSystem), intent(inout) :: this
        
        write(output_unit, *) this%name, " class destruction"
    
    end subroutine PhysicalSystem_destroy

end module class_physicalSystem
