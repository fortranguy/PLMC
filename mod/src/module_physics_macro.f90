!> Subroutines for Physics / macro: after

module module_physics_macro

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit, error_unit, iostat_end
use data_precisions, only: real_zero, io_tiny, consist_tiny
use data_box, only: num_dimensions
use json_module, only: json_file
use module_data, only: test_data_found
use module_types_micro, only: Box_Parameters, Argument_Random, Argument_Initial
use module_physics_micro, only: PBC_distance, random_surface
use class_hard_spheres, only: Hard_Spheres, Dipolar_Hard_Spheres
use class_small_move, only: Small_Move
use class_small_rotation, only: Small_Rotation
use class_hard_spheres_potential, only: Between_Hard_Spheres_Potential_Energy
use module_types_macro, only: Hard_Spheres_Macro, Dipolar_Hard_Spheres_Macro
use class_hard_spheres_observables, only: Hard_Spheres_Observables, Dipolar_Hard_Spheres_Observables
use class_hard_spheres_units, only: Hard_Spheres_Units, Dipolar_Hard_Spheres_Units

implicit none
private
public init_random_seed, set_initial_configuration, &
       init_spheres, init_cells, &
       set_ewald, &
       total_energy, &
       final_spheres, init_between_spheres_potential, final_between_spheres_potential, &
       test_consist

contains

    !> Random number generator: seed
    
    subroutine init_random_seed(arg_rand, report_unit)
    
        type(Argument_Random) :: arg_rand
        integer, intent(in) :: report_unit
    
        integer :: i, n, clock
        integer, dimension(:), allocatable :: seed

        call random_seed(size = n)
        allocate(seed(n))

        select case (arg_rand%choice)
        
            case ('v')

                call system_clock(count=clock)
                seed(:) = clock + 37 * [ (i - 1, i = 1, n) ]
                call random_seed(put = seed)
                ! not sufficent ? cf.
                ! http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html
                write(report_unit, *) "Random number generator: variable"

            case ('p')
                call random_seed(put = arg_rand%seed)
                write(report_unit, *) "Random number generator: put"
                
            case ('f')
                write(report_unit, *) "Random number generator: fix"

            case default
                error stop "Error: init_random_seed"

        end select

        call random_seed(get = seed)
        write(report_unit ,*) "    size = ", n
        write(report_unit ,*) "    seed = ", seed(:)

        deallocate(seed)
        
    end subroutine init_random_seed
    
    !> Initial configuration
    
    subroutine set_initial_configuration(Box_size, arg_init, dipolar, spherical, &
                                         between_spheres_min_distance, &
                                         report_unit)
        
        real(DP), dimension(:), intent(in) :: Box_size
        type(Argument_Initial), intent(in) :: arg_init
        class(Dipolar_Hard_Spheres), intent(inout) :: dipolar
        class(Hard_Spheres), intent(inout) :: spherical
        real(DP), intent(in) :: between_spheres_min_distance
        integer, intent(in) :: report_unit
        
        write(report_unit, *) "Initial configuration: "
        
        select case (arg_init%choice)
        
            case ('r')
                call random_depositions(Box_size, dipolar, spherical, between_spheres_min_distance)
                call random_orientations(dipolar, dipolar%get_num_particles())
                write(output_unit, *) "Random depositions + random orientations"
                write(report_unit, *) "    Random depositions + random orientations"
                
            case ('f')
                write(output_unit, *) "Old configuration"
                write(report_unit, *) "    Old configuration"
                call old_configuration(arg_init%files(1), arg_init%length(1), dipolar, &
                                       norm2(Box_size), "positions")
                call old_configuration(arg_init%files(2), arg_init%length(2), dipolar, &
                                       1._DP, "orientations")
                call old_configuration(arg_init%files(3), arg_init%length(3), spherical, &
                                       norm2(Box_size), "positions")
            
            case default
                error stop "Error: in setting new configuration"
                
        end select
        
    end subroutine set_initial_configuration

    !> Random depositions configuration
    
    subroutine random_depositions(Box_size, spheres1, spheres2, between_spheres_min_distance)

        real(DP), dimension(:), intent(in) :: Box_size
        class(Hard_Spheres), intent(inout) :: spheres1, spheres2
        real(DP), intent(in) :: between_spheres_min_distance

        integer :: i_particle, i_particle_test
        real(DP), dimension(num_dimensions) :: random_position, position, position_test
        real(DP) :: rTest
        
        ! Type 1
        do i_particle = 1, spheres1%get_num_particles()
        
7101        continue
            call random_number(random_position)
            call spheres1%set_position(i_particle, random_position*Box_size)
            
            position(:) = spheres1%get_position(i_particle)
            do i_particle_test = 1, i_particle-1
                position_test(:) = spheres1%get_position(i_particle_test)
                rTest = PBC_distance(Box_size, position, position_test)
                if (rTest < spheres1%get_diameter()) then
                    goto 7101
                end if
            end do
        
        end do
        
        ! Type 2
        do i_particle = 1, spheres2%get_num_particles()
        
7102        continue
            call random_number(random_position)
            call spheres2%set_position(i_particle, random_position*Box_size)
            
            position(:) = spheres2%get_position(i_particle)
            do i_particle_test = 1, spheres1%get_num_particles()
                position_test(:) = spheres1%get_position(i_particle_test)
                rTest = PBC_distance(Box_size, position, position_test)
                if (rTest < between_spheres_min_distance) then
                    goto 7102
                end if
            end do
            
            do i_particle_test = 1, i_particle-1
                position_test(:) = spheres2%get_position(i_particle_test)
                rTest = PBC_distance(Box_size, position, position_test)
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
        
        write(output_unit, *) spheres%get_name()//coordinates_name, " <- ", file(1:length)
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
                
                if (coordinate_norm > norm_max+io_tiny) then
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
        class(Hard_Spheres_Units), intent(in) :: this_units
        
        call this_spheres%test_overlap(Box%size)
        call this_spheres%write_snap_data(this_units%snap_equilibrium_positions)
        call this_spheres%write_snap_positions(0, this_units%snap_initial_positions, &
                                               double_precision=.true.)
        
        select type (this_spheres)
            type is (Dipolar_Hard_Spheres)
                select type (this_units)
                    type is (Dipolar_Hard_Spheres_Units)
                        call this_spheres%write_snap_data(this_units%snap_equilibrium_orientations)
                        call this_spheres%write_snap_orientations(0, &
                             this_units%snap_initial_orientations, double_precision=.true.)
                end select
        end select
    
    end subroutine init_spheres
    
    subroutine init_cells(Box_size, this_spheres, this_macro, other_spheres, between_spheres_potential)
      
        real(DP), dimension(:), intent(in) :: Box_size
        class(Hard_Spheres), intent(in) :: this_spheres, other_spheres
        class(Hard_Spheres_Macro), intent(inout) :: this_macro
        class(Between_Hard_Spheres_Potential_Energy), intent(in) :: between_spheres_potential
        
        real(DP), dimension(num_dimensions) :: same_cell_size, between_cell_size
        
        same_cell_size(:) = this_macro%hard_potential%get_range_cut()
        call this_macro%same_cells%construct(Box_size, same_cell_size, &
                                             this_macro%hard_potential%get_range_cut())
        call this_macro%same_cells%all_cols_to_cells(this_spheres%get_num_particles(), this_spheres)
        
        between_cell_size = between_spheres_potential%get_range_cut()
        call this_macro%between_cells%construct(Box_size, between_cell_size, &
                                                between_spheres_potential%get_range_cut())
        call this_macro%between_cells%all_cols_to_cells(other_spheres%get_num_particles(), &
                                                        other_spheres)

    end subroutine init_cells
    
    subroutine set_ewald(Box, this_spheres, this_macro, json, this_units)
    
        type(Box_Parameters), intent(in) :: Box
        class(Dipolar_Hard_Spheres), intent(in) :: this_spheres
        class(Dipolar_Hard_Spheres_Macro), intent(inout) :: this_macro
        type(json_file), intent(inout) :: json
        class(Dipolar_Hard_Spheres_Units), intent(in) :: this_units
        
        character(len=4096) :: data_name
        logical :: found
        
        real(DP) :: alpha_factor, alpha
        real(DP) :: min_distance
        
        data_name = "Potential Energy.Dipolar Hard Spheres.Ewald summation.alpha factor"
        call json%get(data_name, alpha_factor, found)
        call test_data_found(data_name, found)
        
        alpha = alpha_factor / Box%size(1)
        min_distance = this_macro%hard_potential%get_min_distance()
        call this_macro%ewald_real%construct(Box%size, alpha, min_distance, json)
        call this_macro%ewald_reci%construct(Box, alpha, this_spheres)
        call this_macro%ewald_reci%count_wave_vectors(Box%wave, this_units%wave_vectors)
        call this_macro%ewald_self%set_alpha(alpha)
        call this_macro%ewald_bound%set_total_moment(this_spheres)
    
    end subroutine set_ewald
    
    pure function total_energy(Box, this_spheres, this_macro)
    
        type(Box_Parameters), intent(in) :: Box
        class(Hard_Spheres), intent(in) :: this_spheres
        class(Hard_Spheres_Macro), intent(in) :: this_macro
        real(DP) :: total_energy
        
        total_energy = this_macro%hard_potential%total(Box%size, this_spheres)
        
        select type (this_spheres)
            type is (Dipolar_Hard_Spheres)
                select type (this_macro)
                    type is (Dipolar_Hard_Spheres_Macro)
                        total_energy = total_energy + &
                                       this_macro%ewald_real%total(Box%size, this_spheres) + &
                                       this_macro%ewald_reci%total(Box) - &
                                       this_macro%ewald_self%total(this_spheres) + &
                                       this_macro%ewald_bound%total(Box%size)
                end select
        end select
        
    end function total_energy
    
    !> Spheres finalizations
    
    subroutine final_spheres(Box, this_spheres, this_units)
    
        type(Box_Parameters), intent(in) :: Box
        class(Hard_Spheres), intent(inout) :: this_spheres
        class(Hard_Spheres_Units), intent(in) :: this_units
        
        call this_spheres%test_overlap(Box%size)
        call this_spheres%write_snap_positions(0, this_units%snap_final_positions, &
                                               double_precision=.true.)
        
        select type (this_spheres)
            type is (Dipolar_Hard_Spheres)
                select type (this_units)
                    type is (Dipolar_Hard_Spheres_Units)
                        call this_spheres%write_snap_orientations(0, &
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
                                               between_spheres_report_unit)
    
        real(DP), dimension(:), intent(in) :: Box_size
        class(Between_Hard_Spheres_Potential_Energy), intent(inout) :: between_spheres_potential
        class(Hard_Spheres), intent(in) :: spheres1, spheres2
        real(DP), intent(in) :: potential_energy
        integer, intent(in) :: between_spheres_report_unit
        
        real(DP) :: potential_energy_conf
        
        potential_energy_conf = between_spheres_potential%total(Box_size, spheres1, spheres2)
        call test_consist(potential_energy, potential_energy_conf, between_spheres_report_unit)
    
    end subroutine final_between_spheres_potential
    
    !> Consistency test
    
    subroutine test_consist(potential_energy, potential_energy_conf, report_unit)
    
        real(DP), intent(in) :: potential_energy, potential_energy_conf
        integer, intent(in) :: report_unit
        
        real(DP) :: difference
        
        write(report_unit, *) "Consistency test: "
        write(report_unit, *) "    potential_energy = ", potential_energy
        write(report_unit, *) "    potential_energy_conf = ", potential_energy_conf
        
        if (abs(potential_energy_conf) < real_zero) then
            difference = abs(potential_energy_conf-potential_energy)
            write(report_unit, *) "    absolute difference = ", difference
        else
            difference = abs((potential_energy_conf-potential_energy)/potential_energy_conf)
            write(report_unit, *) "    relative difference = ", difference
        end if
        
        if (difference > consist_tiny) then ! not sufficient for HS ?
            write(report_unit, *) "    WARNING !"
        else
            write(report_unit, *) "    OK !"
        end if
    
    end subroutine test_consist

end module module_physics_macro
