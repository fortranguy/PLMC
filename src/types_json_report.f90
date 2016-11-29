module types_json_report

use json_module, only: json_value

implicit none

private

    type, public :: Generating_JSON_Report
        type(json_value), pointer :: root => null()
        type(json_value), pointer :: random_seed => null()
        type(json_value), pointer :: algorithms_weight => null()
    end type Generating_JSON_Report

    type, public :: Exploring_JSON_Report
        type(json_value), pointer :: root => null()
        type(json_value), pointer :: random_seed => null()
    end type Exploring_JSON_Report

end module types_json_report
