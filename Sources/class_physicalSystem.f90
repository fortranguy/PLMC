!> \brief Description of the Physical System class

module class_physicalSystem

use data_precisions, only: DP
use data_box, only: Ndim
use module_types, only: argument_seed, argument_initial
use class_hardSpheres
use class_dipolarSpheres
use class_mixingPotential
use class_observables
use class_units

implicit none

private

    type, public :: PhysicalSystem
    
        private
        
        character(len=5) :: name
        
        real(DP), dimension(Ndim) :: Lsize !< Box size
        integer, dimension(Ndim) :: Kmax !< Number of wave vectors
        integer :: Ncol !< Number of particles
        
        real(DP) :: Temperature
        integer :: Nthermal, Nadapt, Nstep !< Markov-chain Monte-Carlo
        integer :: Nmove, Nswitch, Nrotate !< Changes
        integer :: switch_Nhit, switch_Nreject
        real(DP) :: switch_reject, switch_rejectSum
        
        real(DP) :: Epot_step, Epot_stepSum, Epot_conf  !< potential energy
        integer :: report_unit !< data & results file
        integer :: obsThermal_unit, obsEquilib_unit !< observables files: thermalisation & equilibrium
        
        type(DipolarSpheres) :: type1_spheres !< physical properties and Monte-Carlo subroutines
        type(MoreObservables) :: type1_obs !< e.g. energy, inverse of activity (-> chemical potential)
        type(MoreUnits) :: type1_units !< files units
        
        type(HardSpheres) :: type2_spheres
        type(Observables) :: type2_obs
        type(Units) :: type2_units
        
        type(MixingPotential) :: mix !< short-range potential
        real(DP) :: mix_Epot_step, mix_Epot_stepSum, mix_Epot_conf
        integer :: mix_report_unit
        integer :: mix_Epot_tab_unit
        integer :: mix_obsThermal_unit, mix_obsEquilib_unit
        
        type(argument_seed) :: arg_seed
        type(argument_initial) :: arg_init
        
        real(DP) :: time_start, time_end
    
    contains
    
    end type PhysicalSystem
    
contains

end module class_physicalSystem
