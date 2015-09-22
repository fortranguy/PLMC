module types_field_parameters

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions

implicit none

private

    type, abstract, public :: Abstract_Field_Parameters

    end type Abstract_Field_Parameters

    type, extends(Abstract_Field_Parameters), public :: Null_Field_Parameters

    end type Null_Field_Parameters

    type, extends(Abstract_Field_Parameters), public :: Constant_Field_Parameters
        real(DP) :: vector(num_dimensions)
    end type Constant_Field_Parameters

end module types_field_parameters
