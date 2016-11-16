module procedures_exploration_inquirers

use json_module, only: json_file
use procedures_property_inquirers, only: logical_from_json
use classes_volume_change_method, only: Abstract_Volume_Change_Method, Concrete_Volume_Change_Method
use classes_particle_insertion_method, only: Abstract_Particle_Insertion_Method, &
    Concrete_Particle_Insertion_Method
use classes_dipolar_neighbourhoods_visitor, only: Abstract_Dipolar_Neighbourhoods_Visitor, &
    Concrete_Dipolar_Neighbourhoods_Visitor

implicit none

private
public :: measure_pressure, measure_chemical_potentials, make_dipoles_graph

interface measure_pressure
    module procedure :: measure_pressure_from_json
    module procedure :: measure_pressure_from_method
end interface measure_pressure

interface measure_chemical_potentials
    module procedure :: measure_chemical_potentials_from_json
    module procedure :: measure_chemical_potentials_from_method
end interface measure_chemical_potentials

interface make_dipoles_graph
    module procedure :: make_dipoles_graph_from_json
    module procedure :: make_dipoles_graph_from_visitor
end interface make_dipoles_graph

contains

    logical function measure_pressure_from_json(exploring_data, prefix) result(measure_pressure)
        type(json_file), intent(inout) :: exploring_data
        character(len=*), intent(in) :: prefix

        measure_pressure = logical_from_json(exploring_data, prefix//"measure pressure")
    end function measure_pressure_from_json

    pure logical function measure_pressure_from_method(volume_change_method) &
        result(measure_pressure)
        class(Abstract_Volume_Change_Method), intent(in) :: volume_change_method

        select type (volume_change_method)
            type is (Concrete_Volume_Change_Method)
                measure_pressure = .true.
            class default
                measure_pressure = .false.
        end select
    end function measure_pressure_from_method

    logical function measure_chemical_potentials_from_json(exploring_data, prefix)&
        result(measure_chemical_potentials)
        type(json_file), intent(inout) :: exploring_data
        character(len=*), intent(in) :: prefix

        measure_chemical_potentials = logical_from_json(exploring_data, &
            prefix//"measure chemical potentials")
    end function measure_chemical_potentials_from_json

    pure logical function measure_chemical_potentials_from_method(particle_insertion_method) &
        result(measure_chemical_potentials)
        class(Abstract_Particle_Insertion_Method), intent(in) :: particle_insertion_method

        select type (particle_insertion_method)
            type is (Concrete_Particle_Insertion_Method)
                measure_chemical_potentials = .true.
            class default
                measure_chemical_potentials = .false.
        end select
    end function measure_chemical_potentials_from_method

    logical function make_dipoles_graph_from_json(exploring_data, prefix) result(make_dipoles_graph)
        type(json_file), intent(inout) :: exploring_data
        character(len=*), intent(in) :: prefix

        make_dipoles_graph = logical_from_json(exploring_data, prefix//"make graph")
    end function make_dipoles_graph_from_json

    elemental logical function make_dipoles_graph_from_visitor(dipolar_neighbourhoods_visitor) &
        result(make_dipoles_graph)
        class(Abstract_Dipolar_Neighbourhoods_Visitor), intent(in) :: dipolar_neighbourhoods_visitor

        select type (dipolar_neighbourhoods_visitor)
            type is (Concrete_Dipolar_Neighbourhoods_Visitor)
                make_dipoles_graph = .true.
            class default
                make_dipoles_graph = .false.
        end select
    end function make_dipoles_graph_from_visitor

end module procedures_exploration_inquirers
