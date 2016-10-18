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
    module procedure :: visit
end interface plmc_visit

contains

    subroutine set_visit_energies(visit_energies, exploring_data, prefix)
        logical, intent(out) :: visit_energies
        type(json_file), intent(inout) :: exploring_data
        character(len=*), intent(in) :: prefix

        visit_energies = logical_from_json(exploring_data, prefix//"visit energies")
    end subroutine set_visit_energies

    subroutine visit(energies, physical_model, visit_energies)
        type(Concrete_Observables_Energies), intent(inout) :: energies(:)
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        logical, optional, intent(in) :: visit_energies

        logical :: overlap
        integer :: i_box

        if (present(visit_energies)) then
            if (.not.visit_energies) return
        end if

        do i_box = 1, size(energies)
            call short_interactions_visit(overlap, energies(i_box)%walls_energies, physical_model%&
                mixture%components(:, i_box), physical_model%short_interactions%&
                walls_visitors(i_box), physical_model%short_interactions%wall_pairs)
            if (overlap) call error_exit("procedures_plmc_visitor: visit: "//&
                "short_interactions_visit: walls: overlap.")
            call short_interactions_visit(overlap, energies(i_box)%short_energies, physical_model%&
                mixture%components(:, i_box), physical_model%short_interactions%components_visitors(i_box), &
                physical_model%short_interactions%components_pairs)
            if (overlap) call error_exit("procedures_plmc_visitor: visit: "//&
                "short_interactions_visit: short: overlap.")
            call dipolar_interactions_visit(energies(i_box)%field_energies, physical_model%environment%&
                external_fields(i_box), physical_model%mixture%components(:, i_box))
            call dipolar_interactions_visit(energies(i_box)%dipolar_energies, energies(i_box)%dipolar_shared_energy, &
                physical_model%mixture%components(:, i_box), physical_model%dipolar_interactions_dynamic(i_box))
            end do
    end subroutine visit

end module procedures_plmc_visitor
