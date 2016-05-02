module classes_plmc_factory

implicit none

private

    type, abstract, public :: Abstract_PLMC_Factory
    contains
        !procedure :: create_physical
    end type Abstract_PLMC_Factory

end module classes_plmc_factory
