module procedures_exploration_inquirers

use json_module, only: json_file
use procedures_property_inquirers, only: logical_from_json
use classes_maximum_box_compression_explorer, only: Abstract_Maximum_Box_Compression_Explorer, &
    Concrete_Maximum_Box_Compression_Explorer
use classes_volume_change_method, only: Abstract_Volume_Change_Method, Concrete_Volume_Change_Method
use classes_particle_insertion_method, only: Abstract_Particle_Insertion_Method, &
    Concrete_Particle_Insertion_Method

implicit none

private
public :: measure_maximum_compression, measure_pressure, measure_chemical_potentials

interface measure_maximum_compression
    module procedure :: measure_maximum_compression_from_json
    module procedure :: measure_maximum_compression_from_explorer
end interface measure_maximum_compression

interface measure_pressure
    module procedure :: measure_pressure_from_json
    module procedure :: measure_pressure_from_method
end interface measure_pressure

interface measure_chemical_potentials
    module procedure :: measure_chemical_potentials_from_json
    module procedure :: measure_chemical_potentials_from_method
end interface measure_chemical_potentials

contains

    logical function measure_maximum_compression_from_json(exploring_data, prefix) &
        result(measure_maximum_compression)
        type(json_file), intent(inout) :: exploring_data
        character(len=*), intent(in) :: prefix

        measure_maximum_compression = logical_from_json(exploring_data, &
            prefix//"measure maximum compression")
    end function measure_maximum_compression_from_json

    pure logical function measure_maximum_compression_from_explorer(&
        maximum_box_compression_explorer) result(measure_maximum_compression)
        class(Abstract_Maximum_Box_Compression_Explorer), intent(in) :: &
            maximum_box_compression_explorer

        select type (maximum_box_compression_explorer)
            type is (Concrete_Maximum_Box_Compression_Explorer)
                measure_maximum_compression = .true.
            class default
                measure_maximum_compression = .false.
        end select
    end function measure_maximum_compression_from_explorer

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

end module procedures_exploration_inquirers
