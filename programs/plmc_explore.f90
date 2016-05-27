program plmc_explore

use, intrinsic :: iso_fortran_env, only: output_unit
use types_physical_model_wrapper, only: Physical_Model_Wrapper
use types_exploring_observables_wrapper, only: Exploring_Observables_Wrapper
use types_exploring_io, only: Exploring_IO_Wrapper
use procedures_plmc_factory, only: plmc_create, plmc_destroy, plmc_set
use procedures_plmc_reset, only: plmc_reset
use procedures_plmc_write, only: plmc_write
use procedures_plmc_help, only: plmc_catch_exploring_help

implicit none

    type(Physical_Model_Wrapper) :: physical_model
    type(Exploring_Observables_Wrapper) :: observables
    type(Exploring_IO_Wrapper) :: io

    call plmc_catch_exploring_help()

    call plmc_create(io%generating_data, io%exploring_data)
    call plmc_create(physical_model, io%generating_data)

    call plmc_create(observables, physical_model%mixture%components)
    call plmc_destroy(io%generating_data, io%exploring_data)

    call plmc_destroy(observables)
    call plmc_destroy(physical_model)

end program plmc_explore
