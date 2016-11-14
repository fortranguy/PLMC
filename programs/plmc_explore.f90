!> @note While [[plmc_generate]] accepts an arbitrary number of boxes,
!> [[plmc_explore]] overrides it to 1. Hence each box can be explored separately.
!>
!> @bug [[procedures_markov_chain_explorer_factory:create]] must be in
!> [[procedures_plmc_factory:plmc_create]] but it doesn't work. Is it an ifort bug?
program plmc_explore

use, intrinsic :: iso_fortran_env, only: output_unit
use procedures_property_inquirers, only: logical_from_json
use types_physical_model_wrapper, only: Physical_Model_Wrapper
use types_markov_chain_explorer_wrapper, only: Markov_Chain_Explorer_Wrapper
use types_exploring_observables_wrapper, only: Exploring_Observables_Wrapper
use types_exploring_io, only: Exploring_IO_Wrapper
use procedures_markov_chain_explorer_factory, markov_chain_explorer_create => create, &
    markov_chain_explorer_destroy => destroy
use procedures_plmc_factory, only: plmc_create, plmc_destroy, plmc_set
use procedures_plmc_resetter, only: plmc_reset
use procedures_plmc_visitor, only: plmc_visit
use procedures_plmc_writer, only: plmc_write
use procedures_plmc_help, only: plmc_catch_exploring_help

implicit none

    type(Physical_Model_Wrapper) :: physical_model
    type(Markov_Chain_Explorer_Wrapper) :: markov_chain_explorer
    type(Exploring_Observables_Wrapper) :: observables
    type(Exploring_IO_Wrapper) :: io

    integer :: i_box, num_snaps, i_snap
    logical :: visit_energies

    call plmc_catch_exploring_help()

    call plmc_create(io%generating_data, io%exploring_data)
    call plmc_create(physical_model, io%generating_data, io%exploring_data, unique_box=.true.)
    call plmc_set(io%generating_data)
    call plmc_set(num_snaps, io%generating_data)
    visit_energies = logical_from_json(io%exploring_data, "Check.visit energies")
    call markov_chain_explorer_create(markov_chain_explorer, physical_model, visit_energies, io%&
        exploring_data)
    call plmc_create(observables, physical_model%mixture%components)
    call plmc_create(io%readers, io%writers, physical_model, markov_chain_explorer, visit_energies,&
        io%generating_data)
    call plmc_destroy(io%generating_data, io%exploring_data)
    call plmc_create(io%json, io%report_data)
    call plmc_write(io%json, io%report_data)
    call plmc_destroy(io%json, io%report_data)

    write(output_unit, *) "Iterations start."
    do i_snap = 1, num_snaps
        call plmc_set(io%readers, i_snap)
        call plmc_reset(physical_model)
        if (visit_energies) then
            call plmc_visit(observables%energies, physical_model, use_cells=.true.)
        end if
        do i_box = 1, size(markov_chain_explorer%changed_boxes_size_ratio)
            call markov_chain_explorer%maximum_boxes_compression_explorer(i_box)%reset()
            call markov_chain_explorer%maximum_boxes_compression_explorer(i_box)%&
                try(observables%maximum_boxes_compression_delta(i_box))
            call markov_chain_explorer%changed_boxes_size_ratio(i_box)%&
                set(observables%maximum_boxes_compression_delta(i_box))
        end do
        call markov_chain_explorer%volume_change_method%try(observables)
        call markov_chain_explorer%particle_insertion_method%try(observables)
        do i_box = 1, size(markov_chain_explorer%dipolar_neighbourhoods_visitors)
            call markov_chain_explorer%dipolar_neighbourhoods_visitors(i_box)%&
                try(observables%adjacency_matrices(:, :, i_box))
        end do
        call plmc_set(observables)
        call plmc_write(io%writers, observables, i_snap)
    end do
    write(output_unit, *) "Iterations end."

    call plmc_destroy(io%readers, io%writers)
    call plmc_destroy(observables)
    call markov_chain_explorer_destroy(markov_chain_explorer)
    call plmc_destroy(physical_model)

end program plmc_explore
