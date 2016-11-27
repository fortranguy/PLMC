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
exclude: procedures_cells_factory.f90
exclude: procedures_pairs_factory.f90
