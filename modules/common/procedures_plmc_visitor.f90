module procedures_plmc_visitor

use json_module, only: json_file
use procedures_errors, only: error_exit
use procedures_short_interactions_visitor, only: short_interactions_visit => visit, &
    short_interactions_visit_cells => visit_cells
use procedures_dipolar_interactions_visitor, only: dipolar_interactions_visit => visit
use types_physical_model_wrapper, only: Physical_Model_Wrapper
use types_observables_energies, only: Concrete_Observables_Energies
use procedures_property_inquirers, only: logical_from_json

implicit none

private
public :: plmc_visit_set, plmc_visit

interface plmc_visit_set
    module procedure :: set_visit_energies
end interface plmc_visit_set

interface plmc_visit
    module procedure :: visit_generating, visit_exploring
end interface plmc_visit

contains

    subroutine set_visit_energies(visit_energies, exploring_data, prefix)
        logical, intent(out) :: visit_energies
        type(json_file), intent(inout) :: exploring_data
        character(len=*), intent(in) :: prefix

        visit_energies = logical_from_json(exploring_data, prefix//"visit energies")
    end subroutine set_visit_energies

    subroutine visit_generating(energies, physical_model)
        type(Concrete_Observables_Energies), intent(inout) :: energies
        type(Physical_Model_Wrapper), intent(in) :: physical_model

        logical :: overlap

        call short_interactions_visit(overlap, energies%walls_energies, physical_model%&
            mixture%components, physical_model%short_interactions)
        if (overlap) call error_exit("procedures_plmc_visitor: visit_generating: "//&
            "short_interactions_visit: walls: overlap.")
        call short_interactions_visit(overlap, energies%short_energies, physical_model%&
            mixture%components, physical_model%short_interactions)
        if (overlap) call error_exit("procedures_plmc_visitor: visit_generating: "//&
            "short_interactions_visit: short: overlap.")
        call dipolar_interactions_visit(energies%field_energies, physical_model%environment%&
            external_field, physical_model%mixture%components)
        call dipolar_interactions_visit(energies%dipolar_energies, energies%dipolar_shared_energy, &
            physical_model%mixture%components, physical_model%dipolar_interactions_dynamic)
    end subroutine visit_generating

    subroutine visit_exploring(energies, physical_model, visit_energies)
        type(Concrete_Observables_Energies), intent(inout) :: energies
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        logical, intent(in) :: visit_energies

        logical :: overlap

        if (.not.visit_energies) return
        call short_interactions_visit(overlap, energies%walls_energies, physical_model%&
            mixture%components, physical_model%short_interactions)
        if (overlap) call error_exit("procedures_plmc_visitor: visit_exploring: "//&
            "short_interactions_visit: walls: overlap.")
        call short_interactions_visit_cells(overlap, energies%short_energies, physical_model%&
            mixture%components,physical_model%short_interactions)
        if (overlap) call error_exit("procedures_plmc_visitor: visit_exploring: "//&
            "short_interactions_visit_cells: overlap.")
        call dipolar_interactions_visit(energies%field_energies, physical_model%environment%&
            external_field, physical_model%mixture%components)
        call dipolar_interactions_visit(energies%dipolar_energies, energies%dipolar_shared_energy, &
            physical_model%mixture%components, physical_model%dipolar_interactions_dynamic)
    end subroutine visit_exploring

end module procedures_plmc_visitor
