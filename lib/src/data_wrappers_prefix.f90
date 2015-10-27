module data_wrappers_prefix

implicit none

private
public :: environment_prefix, mixture_prefix, changes_prefix, &
    short_interactions_prefix, long_interactions_prefix, writers_prefix

    character(len=*), parameter :: environment_prefix = "Environment."
    character(len=*), parameter :: mixture_prefix = "Mixture."
    character(len=*), parameter :: changes_prefix = "Changes."
    character(len=*), parameter :: short_interactions_prefix = "Short Interactions."
    character(len=*), parameter :: long_interactions_prefix = "Long Interactions."
    character(len=*), parameter :: writers_prefix = "Output."

end module data_wrappers_prefix
