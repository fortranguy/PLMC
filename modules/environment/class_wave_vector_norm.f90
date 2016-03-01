module class_wave_vector_norm

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use class_periodic_box, only: Abstract_Periodic_Box
use class_reciprocal_lattice, only: Abstract_Reciprocal_Lattice

implicit none

private

    type, abstract, public :: Abstract_Wave_Vector_Norm

    end type Abstract_Wave_Vector_Norm

end module class_wave_vector_norm
