project: Physics of Liquids Monte-Carlo
summary: A simple Monte-Carlo program for colloids written in Modern Fortran
favicon: ./plmc_logo.png
src_dir: ./programs
    ./modules
extra_mods: iso_fortran_env:https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fFORTRAN_005fENV.html
    json_module: http://jacobwilliams.github.io/json-fortran/
output_dir: ./doc
page_dir: ./modules
predocmark: >
search: true
graph: true
coloured_edges: true
display: public
    protected
    private
sort: permission
version: Green Mile
year: 2016
author: Salomon CHUNG
email: salomon.chung@u-pe.fr
website: http://www.plmc.u-pec.fr/
exclude: procedures_moved_component_coordinates_factory.f90
exclude: procedures_logical_factory.f90
exclude: procedures_plmc_factory.f90
exclude: procedures_des_real_factory.f90
exclude: procedures_des_reci_factory.f90
exclude: procedures_des_self_factory.f90
exclude: procedures_dlc_factory.f90
exclude: procedures_boxes_factory.f90
exclude: procedures_fields_factory.f90
exclude: procedures_walls_factory.f90
exclude: procedures_composition_factory.f90
exclude: procedures_coordinates_factory.f90
exclude: procedures_hard_core_factory.f90
exclude: procedures_mixture_factory.f90
exclude: procedures_exploring_observables_factory.f90
exclude: procedures_generating_observables_factory.f90
exclude: procedures_observables_energies_factory.f90
exclude: procedures_random_coordinates_factory.f90
exclude: procedures_component_coordinates_reader_factory.f90
exclude: procedures_readers_factory.f90
exclude: procedures_cells_factory.f90
exclude: procedures_pairs_factory.f90
exclude: procedures_component_coordinates_writer_factory.f90
exclude: procedures_changes_success_writer_factory.f90
exclude: procedures_line_writer_factory.f90
exclude: procedures_triangle_writer_factory.f90
