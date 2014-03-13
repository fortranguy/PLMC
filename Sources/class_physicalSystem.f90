!> \brief Description of the Physical System class

module class_physicalSystem

use data_precisions, only: DP
use data_box, only: Ndim
use class_hardSpheres
use class_dipolarSpheres
use class_mixingPotential
use class_observables
use class_units

implicit none

private

    type, public :: PhysicalSystem
    
        private
        
        !> Box
        real(DP), dimension(Ndim) :: Lsize
        integer, dimension(Ndim) :: Kmax
        
        !> Number of particles
        integer :: Ncol
        
        !> Markov-chain Monte-Carlo
        real(DP) :: Temperature
        integer :: Nthermal, Nadapt, Nstep
        
        !> Changes
        integer :: Nmove, Nswitch, Nrotate
        
        !> Switch
        integer :: switch_Nhit, switch_Nreject
        real(DP) :: switch_reject, switch_rejectSum
        
        !> Observables
        real(DP) :: Epot_step, Epot_stepSum, Epot_conf  !< potential energy
        integer :: report_unit !< data & results file
        integer :: obsThermal_unit, obsEquilib_unit !< observables files: thermalisation & equilibrium
        
        ! Type 1: Dipolar spheres: Ewald summation
        type(DipolarSpheres) :: type1_spheres !< physical properties and Monte-Carlo subroutines
        type(MoreObservables) :: type1_obs !< e.g. energy, inverse of activity (-> chemical potential)
        type(MoreUnits) :: type1_units !< files units
        
        ! Type 2: Hard spheres
        type(HardSpheres) :: type2_spheres
        type(Observables) :: type2_obs
        type(Units) :: type2_units
        
        ! Mixing potential between 2 types
        type(MixingPotential) :: mix !< short-range potential
        real(DP) :: mix_Epot_step, mix_Epot_stepSum, mix_Epot_conf
        integer :: mix_report_unit
        integer :: mix_Epot_tab_unit
        integer :: mix_obsThermal_unit, mix_obsEquilib_unit
        
        !> Time
        real(DP) :: time_start, time_end
    
    contains
    
    end type PhysicalSystem
    
contains

end module class_physicalSystem
