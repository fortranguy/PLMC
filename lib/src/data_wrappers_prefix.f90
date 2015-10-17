module data_wrappers_prefix

implicit none

private
public :: environment_prefix, mixture_prefix, changes_prefix, &
    short_potentials_prefix, ewalds_prefix

    character(len=*), parameter :: environment_prefix = "Environment."
    character(len=*), parameter :: mixture_prefix = "Mixture."
    character(len=*), parameter :: changes_prefix = "Changes."
    character(len=*), parameter :: short_potentials_prefix = "Short Potentials."
    character(len=*), parameter :: ewalds_prefix = "Ewalds."

end module data_wrappers_prefix
