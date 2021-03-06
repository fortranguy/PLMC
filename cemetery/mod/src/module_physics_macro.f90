!> Subroutines for Physics / macro: after

module module_physics_macro

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit, error_unit, iostat_end
use data_precisions, only: real_zero, input_tiny, consistency_tiny, field_consistency_tiny
use data_box, only: num_dimensions
use json_module, only: json_file, json_value, json_value_create, to_object, json_value_add
use module_data, only: test_data_found
use module_types_micro, only: Box_Parameters, Argument_Random, Argument_Configurations
use module_geometry, only: geometry
use module_physics_micro, only: PBC_distance, random_surface
use class_external_field, only: External_Field
use class_hard_spheres, only: Hard_Spheres, Dipolar_Hard_Spheres
use class_small_move, only: Small_Move
use class_small_rotation, only: Small_Rotation
use class_hard_spheres_potential, only: Between_Hard_Spheres_Potential_Energy
use module_types_macro, only: Hard_Spheres_Macro, Dipolar_Hard_Spheres_Macro
use class_hard_spheres_observables, only: Hard_Spheres_Monte_Carlo_Observables, &
                                          Dipolar_Hard_Spheres_Monte_Carlo_Observables
use class_hard_spheres_units, only: Hard_Spheres_Monte_Carlo_Units, &
                                    Dipolar_Hard_Spheres_Monte_Carlo_Units

implicit none
private
public init_random_seed, set_initial_configuration, &
       init_spheres, init_cells, reset_cells, &
       set_ewald, &
       total_energy, check_field_energy, &
       final_spheres, init_between_spheres_potential, final_between_spheres_potential, &
       test_consistency, test_particles_inside

contains

    !> Random number generator: seed
    
    subroutine init_random_seed(arg_rand, report_json)
    
        type(Argument_Random) :: arg_rand
        type(json_value), pointer, intent(inout) :: report_json
    
        integer :: i, n, clock
        integer, dimension(:), allocatable :: seed
        type(json_value), pointer :: random_json

        call random_seed(size = n)
        allocate(seed(n))

        call json_value_create(random_json)
        call to_object(random_json, "RNG seed")
        call json_value_add(report_json, random_json)

        select case (arg_rand%choice)
        
            case ('v')

                call system_clock(count=clock)
                seed(:) = clock + 37 * [ (i - 1, i = 1, n) ]
                call random_seed(put = seed)
                ! not sufficent ? cf.
                ! http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html
                call json_value_add(random_json, "type", "variable")

            case ('p')
                call random_seed(put = arg_rand%seed)
                call json_value_add(random_json, "type", "put")
                
            case ('f')
                call json_value_add(random_json, "type", "fixed")

            case default
                error stop "Error: init_random_seed"

        end select

        call random_seed(get = seed)
        call json_value_add(random_json, "size", n)
        call json_value_add(random_json, "value", seed)

        deallocate(seed)
        nullify(random_json)
        
    end subroutine init_random_seed
    
    !> Initial configuration
    
    subroutine set_initial_configuration(Box, arg_conf, dipolar, spherical, &
                                         between_spheres_min_distance, &
                                         report_json)
        
        type(Box_Parameters), intent(in) :: Box
        type(Argument_Configurations), intent(in) :: arg_conf
        class(Dipolar_Hard_Spheres), intent(inout) :: dipolar
        class(Hard_Spheres), intent(inout) :: spherical
        real(DP), intent(in) :: between_spheres_min_distance
        type(json_value), pointer, intent(inout) :: report_json
        
        select case (arg_conf%choice)
        
            case ('r')
                call random_depositions(Box, dipolar, spherical, between_spheres_min_distance)
                call random_orientations(dipolar, dipolar%get_num_particles())
                call json_value_add(report_json, "initial configuration", &
                                                 "random depositions + random orientations")
                
            case ('f')
                call json_value_add(report_json, "initial configuration", &
                                                 "old configuration")
                call old_configuration(arg_conf%files(1), arg_conf%length(1), dipolar, &
                                       norm2(Box%size), "positions")
                call old_configuration(arg_conf%files(2), arg_conf%length(2), dipolar, &
                                       1._DP, "orientations")
                call old_configuration(arg_conf%files(3), arg_conf%length(3), spherical, &
                                       norm2(Box%size), "positions")
            
            case default
                error stop "Error: in setting new configuration"
                
        end select
        
    end subroutine set_initial_configuration

    !> Random depositions configuration
    
    subroutine random_depositions(Box, spheres1, spheres2, between_spheres_min_distance)

        type(Box_Parameters), intent(in) :: Box
        class(Hard_Spheres), intent(inout) :: spheres1, spheres2
        real(DP), intent(in) :: between_spheres_min_distance

        integer :: i_particle, i_particle_test
        real(DP), dimension(num_dimensions) :: random_vector, position, position_test
        real(DP) :: rTest, position_z
        
        ! Type 1
        do i_particle = 1, spheres1%get_num_particles()
        
7101        continue
            call random_number(random_vector)
            if (geometry%bulk) then
                call spheres1%set_position(i_particle, random_vector*Box%size)
            else if (geometry%slab) then
                call spheres1%set_position_2d(i_particle, random_vector(1:2)*Box%size(1:2))
                position_z = (Box%height - spheres1%get_diameter()) * random_vector(3) + &
                             spheres1%get_diameter()/2._DP
                call spheres1%set_position_z(i_particle, position_z)
            end if
            
            position(:) = spheres1%get_position(i_particle)
            do i_particle_test = 1, i_particle-1
                position_test(:) = spheres1%get_position(i_particle_test)
                rTest = PBC_distance(Box%size, position, position_test)
                if (rTest < spheres1%get_diameter()) then
                    goto 7101
                end if
            end do
        
        end do
        
        ! Type 2
        do i_particle = 1, spheres2%get_num_particles()
        
7102        continue
            call random_number(random_vector)
            if (geometry%bulk) then
                call spheres2%set_position(i_particle, random_vector*Box%size)
            else if (geometry%slab) then
                call spheres2%set_position_2d(i_particle, random_vector(1:2)*Box%size(1:2))
                position_z = (Box%height - spheres2%get_diameter()) * random_vector(3) + &
                             spheres2%get_diameter()/2._DP
                call spheres2%set_position_z(i_particle, position_z)
            end if
            
            position(:) = spheres2%get_position(i_particle)
            do i_particle_test = 1, spheres1%get_num_particles()
                position_test(:) = spheres1%get_position(i_particle_test)
                rTest = PBC_distance(Box%size, position, position_test)
                if (rTest < between_spheres_min_distance) then
                    goto 7102
                end if
            end do
            
            do i_particle_test = 1, i_particle-1
                position_test(:) = spheres2%get_position(i_particle_test)
                rTest = PBC_distance(Box%size, position, position_test)
                if (rTest < spheres2%get_diameter()) then
                    goto 7102
                end if
            end do
        
        end do
    
    end subroutine random_depositions
    
    !> Uniform (gaussian) orientations
    
    subroutine random_orientations(spheres, num_particles)
    
        class(Dipolar_Hard_Spheres), intent(inout) :: spheres
        integer, intent(in) :: num_particles
        
        integer :: i_particle
        
        do i_particle = 1, num_particles
            call spheres%set_orientation(i_particle, random_surface())
        end do
    
    end subroutine random_orientations
    
    !> From an old configuration
    
    subroutine old_configuration(file, length, spheres, norm_max, coordinates_name)
    
        character(len=*), intent(in) :: file
        integer, intent(in) :: length
        class(Hard_Spheres), intent(inout) :: spheres
        real(DP), intent(in) :: norm_max
        character(len=*), intent(in) :: coordinates_name
        
        integer :: file_unit, readStat

        integer :: i_particle
        real(DP), dimension(num_dimensions) :: coordinate
        real(DP) :: coordinate_norm
        
        write(output_unit, *) spheres%get_name()//"."//coordinates_name, " <- ", file(1:length)
        open(newunit=file_unit, recl=4096, file=file(1:length), status='old', action='read')
        
        i_particle = 0
        do
            read(file_unit, fmt=*, iostat=readStat) coordinate(:)
            if (readStat == iostat_end) exit
            i_particle = i_particle + 1
        end do
        
        if (i_particle == spheres%get_num_particles()) then
            rewind(file_unit)
            do i_particle = 1, spheres%get_num_particles()
            
                read(file_unit, *) coordinate(:)
                select case(coordinates_name)
                    case("positions")
                        call spheres%set_position(i_particle, coordinate)
                        coordinate_norm = norm2(spheres%get_position(i_particle))
                    case("orientations")
                        select type(spheres)
                            type is (Dipolar_Hard_Spheres)
                                call spheres%set_orientation(i_particle, coordinate)
                                coordinate_norm = norm2(spheres%get_orientation(i_particle))
                        end select
                    case default
                        write(error_unit) "Error: unknown coordinate name"
                        error stop
                end select
                
                if (coordinate_norm > norm_max+input_tiny) then
                    write(error_unit, *) "Norm error in file: ", file(1:length)
                    write(error_unit, *) "Index ", i_particle
                    write(error_unit, *) "Norm =", coordinate_norm
                    error stop
                end if
            end do
        else
            write(error_unit, *) "Error reading: ", file(1:length)
            write(error_unit, *) "i_particle", i_particle, " /= ", "num_particles", &
                                 spheres%get_num_particles()
            error stop
        end if
        
        close(file_unit)
        
    end subroutine old_configuration
    
    !> Spheres initialisations
    
    subroutine init_spheres(Box, this_spheres, this_units)
        
        type(Box_Parameters), intent(in) :: Box
        class(Hard_Spheres), intent(inout) :: this_spheres
        class(Hard_Spheres_Monte_Carlo_Units), intent(in) :: this_units

        call check_spheres_in_box(Box, this_spheres)
        call this_spheres%test_overlap(Box%size)
        call this_spheres%write_data(this_units%snap_equilibrium_positions)
        call this_spheres%write_all_positions(0, this_units%snap_initial_positions, &
                                               double_precision=.true.)
        
        select type (this_spheres)
            type is (Dipolar_Hard_Spheres)
                select type (this_units)
                    type is (Dipolar_Hard_Spheres_Monte_Carlo_Units)
                        call this_spheres%write_data(this_units%snap_equilibrium_orientations)
                        call this_spheres%write_all_orientations(0, &
                             this_units%snap_initial_orientations, double_precision=.true.)
                end select
        end select
    
    end subroutine init_spheres

    subroutine check_spheres_in_box(Box, this_spheres)

        type(Box_Parameters), intent(in) :: Box
        class(Hard_Spheres), intent(in) :: this_spheres

        logical :: lower_bounds_ok, upper_bounds_ok

        integer :: i_particle

        do i_particle = 1, this_spheres%get_num_particles()

            if (geometry%bulk) then
                lower_bounds_ok = all(0._DP < this_spheres%get_position(i_particle))
                upper_bounds_ok = all(this_spheres%get_position(i_particle) < Box%size)
            end if

            if (geometry%slab) then
                lower_bounds_ok = all(0._DP < this_spheres%get_position_2d(i_particle)) .and. &
                                (this_spheres%get_diameter()/2._DP < &
                                this_spheres%get_position_z(i_particle))
                upper_bounds_ok = all(this_spheres%get_position_2d(i_particle) < Box%size(1:2)) .and. &
                                (this_spheres%get_position_z(i_particle) < &
                                Box%height - this_spheres%get_diameter()/2._DP)
            end if

            if (.not.lower_bounds_ok .or. .not.upper_bounds_ok) then
                write(error_unit, *) this_spheres%get_name(), " number ", i_particle, &
                                     " is not in the box."
                error stop
            end if

        end do

    end subroutine check_spheres_in_box
    
    subroutine init_cells(Box_size, this_spheres, this_macro, other_spheres, between_spheres_potential)
      
        real(DP), dimension(:), intent(in) :: Box_size
        class(Hard_Spheres), intent(in) :: this_spheres, other_spheres
        class(Hard_Spheres_Macro), intent(inout) :: this_macro
        class(Between_Hard_Spheres_Potential_Energy), intent(in) :: between_spheres_potential
        
        real(DP), dimension(num_dimensions) :: same_cell_size, between_cell_size
        
        same_cell_size(:) = this_macro%hard_potential%get_cutoff()
        call this_macro%same_cells%construct(Box_size, same_cell_size, &
                                             this_macro%hard_potential%get_cutoff())
        call this_macro%same_cells%all_particles_to_cells(this_spheres)
        
        between_cell_size = between_spheres_potential%get_cutoff()
        call this_macro%between_cells%construct(Box_size, between_cell_size, &
                                                between_spheres_potential%get_cutoff())
        call this_macro%between_cells%all_particles_to_cells(other_spheres)

    end subroutine init_cells

    subroutine reset_cells(this_spheres, this_macro, other_spheres)

        class(Hard_Spheres), intent(in) :: this_spheres, other_spheres
        class(Hard_Spheres_Macro), intent(inout) :: this_macro

        call this_macro%same_cells%dealloc_nodes()
        call this_macro%same_cells%alloc_nodes()
        call this_macro%same_cells%all_particles_to_cells(this_spheres)

        call this_macro%between_cells%dealloc_nodes()
        call this_macro%between_cells%alloc_nodes()
        call this_macro%between_cells%all_particles_to_cells(other_spheres)

    end subroutine reset_cells
    
    subroutine set_ewald(Box, this_spheres, this_macro, data_json, this_units)
    
        type(Box_Parameters), intent(in) :: Box
        class(Dipolar_Hard_Spheres), intent(in) :: this_spheres
        class(Dipolar_Hard_Spheres_Macro), intent(inout) :: this_macro
        type(json_file), intent(inout) :: data_json
        class(Dipolar_Hard_Spheres_Monte_Carlo_Units), intent(in), optional :: this_units
        
        character(len=4096) :: data_name
        logical :: found
        
        real(DP) :: alpha_factor, alpha
        real(DP) :: min_distance
        
        data_name = "Potential Energy.Dipolar Hard Spheres.Ewald summation.alpha factor"
        call data_json%get(data_name, alpha_factor, found)
        call test_data_found(data_name, found)
        
        alpha = alpha_factor / Box%size(1)
        min_distance = this_macro%hard_potential%get_min_distance()
        call this_macro%ewald_real%construct(Box%size, alpha, min_distance, data_json)
        call this_macro%ewald_reci%construct(Box, alpha, this_spheres)
        if (present(this_units)) then
            call this_macro%ewald_reci%count_wave_vectors(Box%wave, this_units%wave_vectors)
        end if
        call this_macro%ewald_self%set_alpha(alpha)
        call this_macro%ewald_bound%set_total_moment(this_spheres)
        if (geometry%slab) then
            call this_macro%dlc%construct(Box, this_spheres)
            if (present(this_units)) then
                call this_macro%dlc%count_wave_vectors(Box%wave, this_units%ELC_wave_vectors)
            end if
        end if
    
    end subroutine set_ewald
    
    pure function total_energy(Box, this_spheres, this_macro, ext_field)
    
        type(Box_Parameters), intent(in) :: Box
        class(Hard_Spheres), intent(in) :: this_spheres
        class(Hard_Spheres_Macro), intent(in) :: this_macro
        class(External_Field), intent(in) :: ext_field
        real(DP) :: total_energy
        
        total_energy = this_macro%hard_potential%total(Box%size, this_spheres)
        
        select type (this_spheres)
            type is (Dipolar_Hard_Spheres)
                select type (this_macro)
                    type is (Dipolar_Hard_Spheres_Macro)
                        total_energy = total_energy + &
                                       this_macro%ewald_real%total_energy(Box%size, this_spheres) + &
                                       this_macro%ewald_reci%total_energy(Box) - &
                                       this_macro%ewald_self%total_energy(this_spheres) + &
                                       this_macro%ewald_bound%total_energy(Box%size) + &
                                       ext_field%total_energy(this_macro%ewald_bound%get_total_moment())
                        if (geometry%slab) then
                            total_energy = total_energy - this_macro%dlc%total_energy(Box)
                        end if
                end select
        end select
        
    end function total_energy
    
    pure function total_energy_field(Box, this_spheres, this_macro, ext_field)
    
        type(Box_Parameters), intent(in) :: Box
        class(Dipolar_Hard_Spheres), intent(in) :: this_spheres
        class(Dipolar_Hard_Spheres_Macro), intent(in) :: this_macro
        class(External_Field), intent(in) :: ext_field
        real(DP) :: total_energy_field
        
        logical, parameter :: using_field = .true.
        
        total_energy_field = this_macro%hard_potential%total(Box%size, this_spheres)

        total_energy_field = total_energy_field + &
            this_macro%ewald_real%total_energy(Box%size, this_spheres, using_field) + &
            this_macro%ewald_reci%total_energy(Box, using_field, this_spheres) - &
            this_macro%ewald_self%total_energy(this_spheres, using_field) + &
            this_macro%ewald_bound%total_energy(Box%size, using_field) + &
            ext_field%total_energy(this_macro%ewald_bound%get_total_moment())
        if (geometry%slab) then
            total_energy_field = total_energy_field - &
                this_macro%dlc%total_energy(Box, using_field, this_spheres)
        end if

    end function total_energy_field
    
    subroutine check_field_energy(Box, this_spheres, this_macro, ext_field)
    
        type(Box_Parameters), intent(in) :: Box
        class(Dipolar_Hard_Spheres), intent(in) :: this_spheres
        class(Dipolar_Hard_Spheres_Macro), intent(in) :: this_macro
        class(External_Field), intent(in) :: ext_field
        
        real(DP) :: energy, energy_field
        
        energy = total_energy(Box, this_spheres, this_macro, ext_field)
        energy_field = total_energy_field(Box, this_spheres, this_macro, ext_field)
        
        if (abs(energy - energy_field) > field_consistency_tiny) then
            write(error_unit, *) this_spheres%get_name(), &
                                 ": Energy with and without field don't match: "
            error stop
        end if
        
    end subroutine
    
    !> Spheres finalizations
    
    subroutine final_spheres(Box, this_spheres, this_units)
    
        type(Box_Parameters), intent(in) :: Box
        class(Hard_Spheres), intent(inout) :: this_spheres
        class(Hard_Spheres_Monte_Carlo_Units), intent(in) :: this_units

        call check_spheres_in_box(Box, this_spheres)
        call this_spheres%test_overlap(Box%size)
        call this_spheres%write_all_positions(0, this_units%snap_final_positions, &
                                               double_precision=.true.)
        
        select type (this_spheres)
            type is (Dipolar_Hard_Spheres)
                select type (this_units)
                    type is (Dipolar_Hard_Spheres_Monte_Carlo_Units)
                        call this_spheres%write_all_orientations(0, &
                             this_units%snap_final_orientations, double_precision=.true.)
                end select
        end select
    
    end subroutine final_spheres
    
    !> Mix initialisation
    
    subroutine init_between_spheres_potential(Box_size, between_spheres_potential, spheres1, spheres2, &
                                              write_potential_energy, &
                                              potential_energy, potential_energy_unit)
    
        real(DP), dimension(:), intent(in) :: Box_size
        class(Between_Hard_Spheres_Potential_Energy), intent(inout) :: between_spheres_potential
        class(Hard_Spheres), intent(in) :: spheres1, spheres2
        logical, intent(in) :: write_potential_energy
        real(DP), intent(out) :: potential_energy
        integer, intent(in) :: potential_energy_unit
        
        if (write_potential_energy) then
            call between_spheres_potential%write(potential_energy_unit)
        end if
        potential_energy = between_spheres_potential%total(Box_size, spheres1, spheres2)
    
    end subroutine init_between_spheres_potential
    
    !> Mix finalization
    
    subroutine final_between_spheres_potential(Box_size, between_spheres_potential, &
                                               spheres1, spheres2, potential_energy, &
                                               report_json)
    
        real(DP), dimension(:), intent(in) :: Box_size
        class(Between_Hard_Spheres_Potential_Energy), intent(inout) :: between_spheres_potential
        class(Hard_Spheres), intent(in) :: spheres1, spheres2
        real(DP), intent(in) :: potential_energy
        type(json_value), pointer, intent(in) :: report_json
        
        real(DP) :: potential_energy_conf
        
        potential_energy_conf = between_spheres_potential%total(Box_size, spheres1, spheres2)
        call test_consistency(potential_energy, potential_energy_conf, report_json)
    
    end subroutine final_between_spheres_potential
    
    !> Consistency test
    
    subroutine test_consistency(potential_energy, potential_energy_conf, report_json)
    
        real(DP), intent(in) :: potential_energy, potential_energy_conf
        type(json_value), pointer, intent(in) :: report_json
        
        real(DP) :: difference
        type(json_value), pointer :: consistency_json

        call json_value_create(consistency_json)
        call to_object(consistency_json, "Consistency")
        call json_value_add(report_json, consistency_json)
        
        if (abs(potential_energy_conf) < real_zero) then
            difference = abs(potential_energy_conf-potential_energy)
            call json_value_add(consistency_json, "absolute difference", difference)
        else
            difference = abs((potential_energy_conf-potential_energy)/potential_energy_conf)
            call json_value_add(consistency_json, "relative difference", difference)
        end if
        
        if (difference > consistency_tiny) then ! not sufficient for HS ?
            call json_value_add(consistency_json, "pass", .false.)
        else
            call json_value_add(consistency_json, "pass", .true.)
        end if

        nullify(consistency_json)
    
    end subroutine test_consistency
    
    subroutine test_particles_inside(Box_lower_bounds, Box_upper_bounds, &
                                     num_particles, all_positions, particles_inside, num_particles_step)
                                     
        real(DP), dimension(:), intent(in) :: Box_lower_bounds, Box_upper_bounds
        integer, intent(in) :: num_particles
        real(DP), dimension(:, :), intent(in) :: all_positions
        logical, dimension(:), intent(out) :: particles_inside
        integer, intent(out) :: num_particles_step
        
        integer :: i_particle
        logical, dimension(num_dimensions) :: position_inside
        
        num_particles_step = 0
        do i_particle = 1, num_particles
            
            position_inside(:) = (Box_lower_bounds(:) < all_positions(:, i_particle)) .and. &
                                 (all_positions(:, i_particle) < Box_upper_bounds(:))                                 
            
            if (all(position_inside)) then
                particles_inside(i_particle) = .true.
                num_particles_step = num_particles_step + 1
            else
                particles_inside(i_particle) = .false.
            end if
            
        end do
    
    end subroutine test_particles_inside

end module module_physics_macro
