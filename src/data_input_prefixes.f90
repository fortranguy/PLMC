module data_input_prefixes

implicit none

private
public :: environment_prefix, mixture_prefix, short_interactions_prefix, &
    dipolar_interactions_prefix, random_number_generator_prefix, changes_prefix, readers_prefix, &
    writers_prefix, particle_insertion_prefix, volume_change_prefix, density_prefix, radial_prefix,&
    dipoles_graph_prefix

    character(len=*), parameter :: environment_prefix = "Environment."
    character(len=*), parameter :: mixture_prefix = "Mixture."
    character(len=*), parameter :: short_interactions_prefix = "Short Interactions."
    character(len=*), parameter :: dipolar_interactions_prefix = "Dipolar Interactions."
    character(len=*), parameter :: random_number_generator_prefix = "Random Number Generator."
    character(len=*), parameter :: changes_prefix = "Changes."
    character(len=*), parameter :: readers_prefix = "Input."
    character(len=*), parameter :: writers_prefix = "Output."
    character(len=*), parameter :: particle_insertion_prefix = "Particle Insertion."
    character(len=*), parameter :: volume_change_prefix = "Volume Change."
    character(len=*), parameter :: density_prefix = "Density."
    character(len=*), parameter :: radial_prefix = "Radial."
    character(len=*), parameter :: dipoles_graph_prefix = "Dipoles Graph."

end module data_input_prefixes
