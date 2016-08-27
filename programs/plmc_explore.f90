program plmc_explore

use, intrinsic :: iso_fortran_env, only: output_unit
use types_physical_model_wrapper, only: Physical_Model_Wrapper
use types_markov_chain_explorer_wrapper, only: Markov_Chain_Explorer_Wrapper
use types_exploring_observables_wrapper, only: Exploring_Observables_Wrapper
use types_exploring_io, only: Exploring_IO_Wrapper
use procedures_markov_chain_explorer_factory, markov_chain_explorer_create => create, &
    markov_chain_explorer_destroy => destroy
use procedures_plmc_factory, only: plmc_create, plmc_destroy, plmc_set
use procedures_plmc_reset, only: plmc_reset
use procedures_plmc_visit, only: plmc_visit_set, plmc_visit
use procedures_plmc_write, only: plmc_write
use procedures_plmc_help, only: plmc_catch_exploring_help

implicit none

    type(Physical_Model_Wrapper) :: physical_model
    type(Markov_Chain_Explorer_Wrapper) :: markov_chain_explorer
    type(Exploring_Observables_Wrapper) :: observables
    type(Exploring_IO_Wrapper) :: io

    integer :: num_snaps, i_snap
    logical :: visit_energies

    call plmc_catch_exploring_help()

    call plmc_create(io%generating_data, io%exploring_data)
    call plmc_create(physical_model, io%generating_data, io%exploring_data)
    call plmc_set(io%generating_data)
    call plmc_set(num_snaps, io%generating_data)
    call plmc_visit_set(visit_energies, io%exploring_data, "Check.")
    call markov_chain_explorer_create(markov_chain_explorer, physical_model, visit_energies, io%&
        exploring_data)
        !! @bug [[markov_chain_explorer_create]] must be in [[plmc_create]]: ifort bug?
    call plmc_create(observables, physical_model%mixture%components)
    call plmc_create(io%readers, io%writers, physical_model, markov_chain_explorer, visit_energies)
    call plmc_destroy(io%generating_data, io%exploring_data)
    call plmc_create(io%json, io%report_data)
    call plmc_write(io%json, io%report_data)
    call plmc_destroy(io%json, io%report_data)

    write(output_unit, *) "Iterations start."
    do i_snap = 1, num_snaps
        call plmc_set(io%readers, i_snap)
        call plmc_reset(physical_model)
        call plmc_visit(observables%energies, physical_model, visit_energies)
        call markov_chain_explorer%maximum_box_compression_explorer%try(observables)
        call markov_chain_explorer%volume_change_method%try(observables)
        call markov_chain_explorer%particle_insertion_method%try(observables)
        call plmc_set(observables)
        call plmc_write(io%writers, observables, i_snap)
    end do
    write(output_unit, *) "Iterations end."

    call plmc_destroy(io%readers, io%writers)
    call plmc_destroy(observables)
    call markov_chain_explorer_destroy(markov_chain_explorer)
    call plmc_destroy(physical_model)

end program plmc_explore
