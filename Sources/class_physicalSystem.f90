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

implicit none

private

    type, public :: PhysicalSystem
    
        private
        
        character(len=5) :: name
        
        ! Box
        real(DP), dimension(Ndim) :: Lsize !< Box size
        integer, dimension(Ndim) :: Kmax !< Number of wave vectors
        integer :: Ncol !< Number of particles
        
        ! Monte-Carlo
        real(DP) :: Temperature
        integer :: Nthermal, Nadapt, Nstep !< Markov-chain Monte-Carlo
        integer :: Nmove, Nswitch, Nrotate !< Changes
        
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
        integer :: mix_report_unit
        integer :: mix_Epot_tab_unit
        integer :: mix_obsThermal_unit, mix_obsEquilib_unit
        
        ! Observables: write to files
        integer :: report_unit !< data & results file
        integer :: obsThermal_unit, obsEquilib_unit !< observables files: thermalisation & equilibrium
           
        real(DP) :: time_start, time_end
        
        ! Switch
        integer :: switch_Nhit, switch_Nreject
        real(DP) :: switch_reject, switch_rejectSum
    
    contains
    
        !> Construction and destruction of the class
        procedure :: init_box => PhysicalSystem_init_box
        procedure :: init_monteCarlo => PhysicalSystem_init_monteCarlo
        procedure :: init_spheres => PhysicalSystem_init_spheres
        procedure :: init_switch => PhysicalSystem_init_switch
        procedure :: init_changes => PhysicalSystem_init_changes
        procedure :: construct => PhysicalSystem_construct
        procedure :: destroy => PhysicalSystem_destroy
    
    end type PhysicalSystem
    
contains

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
    
    subroutine PhysicalSystem_init_spheres(this)
    
        class(PhysicalSystem), intent(inout) :: this
        
        call this%type1_spheres%construct()
        call this%type2_spheres%construct()
        call this%mix%construct(this%type1_spheres%get_sigma(), this%type2_spheres%get_sigma())
    
    end subroutine PhysicalSystem_init_spheres
    
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

    subroutine PhysicalSystem_construct(this, arg_seed, arg_init)
    
        class(PhysicalSystem), intent(out) :: this
        type(argument_seed), intent(in) :: arg_seed
        type(argument_initial), intent(in) :: arg_init
        
        this%name = "sys"
        write(output_unit, *) this%name, " class construction"
        
        call this%init_box()
        call this%init_monteCarlo()
        call this%init_spheres()
        call this%init_changes()
        call this%init_switch()
    
    end subroutine PhysicalSystem_construct
    
    subroutine PhysicalSystem_destroy(this)
    
        class(PhysicalSystem), intent(inout) :: this
        
        write(output_unit, *) this%name, " class destruction"
    
    end subroutine PhysicalSystem_destroy

end module class_physicalSystem
