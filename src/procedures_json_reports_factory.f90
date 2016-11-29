module procedures_json_reports_factory

use, intrinsic :: iso_fortran_env, only: error_unit
use json_module, only: json_core
use data_output_objects, only: random_number_generator_object, generating_algorithms_weights_object
use types_json_report, only: Generating_JSON_Report, Exploring_JSON_Report

implicit none

private
public :: create, destroy

interface create
    module procedure :: create_generating, create_exploring
end interface create

interface destroy
    module procedure :: destroy_generating, destroy_exploring
end interface destroy

contains

    subroutine create_generating(json, report)
        type(json_core), intent(inout) :: json
        type(Generating_JSON_Report), intent(inout) :: report

        call json%create_object(report%random_seed, random_number_generator_object)
        call json%add(report%root, report%random_seed)

        call json%create_object(report%algorithms_weight, &
            generating_algorithms_weights_object)
        call json%add(report%root, report%algorithms_weight)
    end subroutine create_generating

    subroutine destroy_generating(json, report)
        type(json_core), intent(inout) :: json
        type(Generating_JSON_Report), intent(inout) :: report

        report%algorithms_weight => null()
        report%random_seed => null()
        call json%destroy(report%root)
        if (json%failed()) call json%print_error_message(error_unit)
    end subroutine destroy_generating

    subroutine create_exploring(json, report)
        type(json_core), intent(inout) :: json
        type(Exploring_JSON_Report), intent(inout) :: report

        call json%create_object(report%random_seed, random_number_generator_object)
        call json%add(report%root, report%random_seed)
    end subroutine create_exploring

    subroutine destroy_exploring(json, report)
        type(json_core), intent(inout) :: json
        type(Exploring_JSON_Report), intent(inout) :: report

        report%random_seed => null()
        call json%destroy(report%root)
        if (json%failed()) call json%print_error_message(error_unit)
    end subroutine destroy_exploring

end module procedures_json_reports_factory
