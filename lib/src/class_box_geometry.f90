module class_box_geometry

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions

implicit none

private

    type, abstract, public :: Abstract_Box_Geometry
    private
        real(DP), dimension(num_dimensions) :: size
        integer, dimension(num_dimensions) :: wave
    contains
        procedure :: get_size => Abstract_Box_Geometry_get_size
        procedure :: get_wave => Abstract_Box_Geometry_get_wave
    end type Abstract_Box_Geometry
    
contains

end module class_box_geometry
