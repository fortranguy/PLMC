!> \brief Description of the Mixing Potential class

module class_mixingPotential

use data_constants
use data_particles
use data_potentiel

implicit none

private

	type, public :: MixingPotential
	
		private
		
		real(DP) :: rMin !< minimum distance between two particles
		real(DP) :: rCut !< short-range cut
		real(DP) :: dr !< discretisation step
        integer :: iMin !< minimum index of tabulation : minimum distance
        integer :: iCut !< maximum index of tabulation : until potential cut
        real(DP) :: epsilon !< factor in Yukawa
        real(DP) :: alpha !< coefficient in Yukawa
        real(DP), dimension(:), allocatable :: ePot_tab !< tabulation
	
	contains
	
		procedure :: ePot_init => Interacting_ePot_init
        procedure :: ePot => Interacting_ePot
        procedure :: ePot_neigh => Interacting_ePot_neigh
        procedure :: ePot_total => Interacting_ePot_total
        procedure :: consistTest => Interacting_consistTest
	
	end type
	
contains

end module class_mixingPotential
