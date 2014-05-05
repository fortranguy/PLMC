!> Subroutines for Physics / macro: after

module module_physics_macro

use, intrinsic :: iso_fortran_env, only: output_unit, error_unit, iostat_end
use data_precisions, only: DP, real_zero, io_tiny, consist_tiny
use data_box, only: Ndim
use module_types, only: Box_Dimensions, Argument_Random, Argument_Initial
use module_physics_micro, only: dist_PBC, random_surface
use class_neighbour_cells
use class_hard_spheres
use class_hard_spheres_potential
use class_small_move
use class_small_rotation
use class_dipolar_spheres
use class_mixing_potential
use class_observables
use class_units

implicit none
private
public init_random_seed, set_initial_configuration, &
       init_spheres, init_cells, final_spheres, init_mix, mix_final, &
       adapt_move, adapt_rotation, test_consist

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
    
    subroutine set_initial_configuration(Box_size, arg_init, dipolar, spherical, mix_min_distance, &
                                         report_unit)
        
        real(DP), dimension(:), intent(in) :: Box_size
        type(Argument_Initial), intent(in) :: arg_init
        class(Dipolar_Spheres), intent(inout) :: dipolar
        class(Hard_Spheres), intent(inout) :: spherical
        real(DP), intent(in) :: mix_min_distance
        integer, intent(in) :: report_unit
        
        write(report_unit, *) "Initial configuration: "
        
        select case (arg_init%choice)
        
            case ('r')
                call random_depositions(Box_size, dipolar, spherical, mix_min_distance)
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
    
    subroutine random_depositions(Box_size, spheres1, spheres2, mix_min_distance)

        real(DP), dimension(:), intent(in) :: Box_size
        class(Hard_Spheres), intent(inout) :: spheres1, spheres2
        real(DP), intent(in) :: mix_min_distance

        integer :: i_particle, i_particle_test
        real(DP), dimension(Ndim) :: xRand, position, position_test
        real(DP) :: rTest
        
        ! Type 1
        do i_particle = 1, spheres1%get_num_particles()
        
7101        continue
            call random_number(xRand)
            call spheres1%set_position(i_particle, xRand*Box_size)
            
            position(:) = spheres1%get_position(i_particle)
            do i_particle_test = 1, i_particle-1
                position_test(:) = spheres1%get_position(i_particle_test)
                rTest = dist_PBC(Box_size, position, position_test)
                if (rTest < spheres1%get_diameter()) then
                    goto 7101
                end if
            end do
        
        end do
        
        ! Type 2
        do i_particle = 1, spheres2%get_num_particles()
        
7102        continue
            call random_number(xRand)
            call spheres2%set_position(i_particle, xRand*Box_size)
            
            position(:) = spheres2%get_position(i_particle)
            do i_particle_test = 1, spheres1%get_num_particles()
                position_test(:) = spheres1%get_position(i_particle_test)
                rTest = dist_PBC(Box_size, position, position_test)
                if (rTest < mix_min_distance) then
                    goto 7102
                end if
            end do
            
            do i_particle_test = 1, i_particle-1
                position_test(:) = spheres2%get_position(i_particle_test)
                rTest = dist_PBC(Box_size, position, position_test)
                if (rTest < spheres2%get_diameter()) then
                    goto 7102
                end if
            end do
        
        end do
    
    end subroutine random_depositions
    
    !> Uniform (gaussian) orientations
    
    subroutine random_orientations(spheres, num_particles)
    
        class(Dipolar_Spheres), intent(inout) :: spheres
        integer, intent(in) :: num_particles
        
        integer :: i_particle
        
        do i_particle = 1, num_particles
            call spheres%set_orientation(i_particle, random_surface())
        end do
    
    end subroutine random_orientations
    
    !> From an old configuration
    
    subroutine old_configuration(file, length, spheres, normMax, coordinates_name)
    
        character(len=*), intent(in) :: file
        integer, intent(in) :: length
        class(Hard_Spheres), intent(inout) :: spheres
        real(DP), intent(in) :: normMax
        character(len=*), intent(in) :: coordinates_name
        
        integer :: file_unit, readStat

        integer :: i_particle
        real(DP), dimension(Ndim) :: coordinate
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
                            type is (Dipolar_Spheres)
                                call spheres%set_orientation(i_particle, coordinate)
                                coordinate_norm = norm2(spheres%get_orientation(i_particle))
                        end select
                    case default
                        write(error_unit) "Error: unknown coordinate name"
                        error stop
                end select
                
                if (coordinate_norm > normMax+io_tiny) then
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
    
    subroutine init_spheres(Box, this, this_potential, write_potential, this_units, &
                            this_Epot) ! separate ?
    
        type(Box_Dimensions), intent(in) :: Box
        class(Hard_Spheres), intent(inout) :: this
        class(Hard_Spheres_Potential), intent(inout) :: this_potential
        logical, intent(in) :: write_potential
        class(Units), intent(in) :: this_units
        real(DP), intent(out) :: this_Epot
        
        call this%test_overlap(Box%size)
        call this%write_snap_data(this_units%snap_positions)
        call this%write_snap_positions(0, this_units%snapIni_positions)
        call this_potential%construct(this%get_diameter())
        this_Epot = this_potential%conf()
        
        if (write_potential) then
            call this_potential%write(this_units%Epot)
        end if
        select type (this)
            type is (Dipolar_Spheres)
                this_Epot = this_Epot + this%Epot_conf(Box) ! temp
                select type (this_units)
                    type is (MoreUnits)
                        call this%write_snap_data(this_units%snap_orientations)
                        call this%write_snap_orientations(0, this_units%snapIni_orientations)
                        if (write_potential) then
                            call this%write_Epot_real(this_units%Epot_real)
                        end if
                        call this%Epot_reci_count_waveVectors(Box%wave, this_units%waveVectors)
                end select
        end select
    
    end subroutine init_spheres
    
    subroutine init_cells(Box_size, this_spheres, sameCells, other_spheres, mixCells, this_potential, &
                          mix)
    
        real(DP), dimension(:), intent(in) :: Box_size
        class(Hard_Spheres), intent(in) :: this_spheres, other_spheres
        class(Neighbour_Cells), intent(out) :: sameCells, mixCells
        class(Hard_Spheres_Potential), intent(in) :: this_potential
        class(Mixing_Potential), intent(in) :: mix
        
        real(DP), dimension(Ndim) :: same_cell_size
        
        same_cell_size(:) = this_potential%get_range_cut()
        call sameCells%construct(Box_size, same_cell_size, this_potential%get_range_cut())
        call sameCells%all_cols_to_cells(this_spheres%get_num_particles(), &
                                              this_spheres%get_all_positions())
        
        call mixCells%construct(Box_size, mix%get_cell_size(), mix%get_rCut())
        call mixCells%all_cols_to_cells(other_spheres%get_num_particles(), &
                                             other_spheres%get_all_positions())

    end subroutine init_cells
    
    !> Spheres finalizations
    
    subroutine final_spheres(Box, this, this_units, this_obs)
    
        type(Box_Dimensions), intent(in) :: Box
        class(Hard_Spheres), intent(inout) :: this
        class(Units), intent(in) :: this_units
        class(Observables), intent(in) :: this_obs
        
        call this%test_overlap(Box%size)
        select type (this)
            type is (Dipolar_Spheres)
                call this%set_Epot(Box) ! ugly !
                call test_consist(this_obs%Epot, this%Epot_conf(Box), this_units%report)
                ! for HS too ?
        end select
        call this%write_snap_positions(0, this_units%snapFin_positions)
        
        select type (this)
            type is (Dipolar_Spheres)
                select type (this_units)
                    type is (MoreUnits)
                        call this%write_snap_orientations(0, this_units%snapFin_orientations)
                end select
        end select
    
    end subroutine final_spheres
    
    !> Mix initialisation
    
    subroutine init_mix(Box_size, mix, spheres1, spheres2, write_potential, mix_Epot_unit, mix_Epot)
    
        real(DP), dimension(:), intent(in) :: Box_size
        class(Mixing_Potential), intent(inout) :: mix
        class(Hard_Spheres), intent(in) :: spheres1, spheres2
        logical, intent(in) :: write_potential
        integer, intent(in) :: mix_Epot_unit
        real(DP), intent(out) :: mix_Epot
    
        call mix%test_overlap(Box_size, spheres1, spheres2)
        call mix%set_Epot()
        if (write_potential) then
            call mix%write_Epot(mix_Epot_unit)
        end if
        call mix%set_cell_size()
        mix_Epot = mix%Epot_conf(Box_size, spheres1, spheres2)
    
    end subroutine init_mix
    
    !> Mix finalization
    
    subroutine mix_final(Box_size, mix, spheres1, spheres2, mix_report_unit, mix_Epot, mix_Epot_conf)
    
        real(DP), dimension(:), intent(in) :: Box_size
        class(Mixing_Potential), intent(inout) :: mix
        class(Hard_Spheres), intent(in) :: spheres1, spheres2
        integer, intent(in) :: mix_report_unit
        real(DP), intent(in) :: mix_Epot
        real(DP), intent(out) :: mix_Epot_conf
        
        call mix%test_overlap(Box_size, spheres1, spheres2)
        call mix%set_Epot()
        mix_Epot_conf = mix%Epot_conf(Box_size, spheres1, spheres2)
        call test_consist(mix_Epot, mix_Epot_conf, mix_report_unit)
    
    end subroutine mix_final
    
    !> Change: average & adaptation
    
    subroutine adapt_move(Box_size, this, Nadapt, iStep, obs, move_unit)
    
        real(DP), dimension(:), intent(in) :: Box_size
        class(Small_Move), intent(inout) :: this
        integer, intent(in) :: Nadapt, iStep
        class(Observables), intent(inout) :: obs
        integer, intent(in) :: move_unit
    
        obs%move_rejectAvg = obs%move_rejectAdapt / real(Nadapt-1, DP)
        obs%move_rejectAdapt = 0._DP
        call this%adapt_delta(Box_size, obs%move_rejectAvg)
        write(move_unit, *) iStep, this%get_delta_scalar(), obs%move_rejectAvg
    
    end subroutine adapt_move
    
    subroutine adapt_rotation(this, Nadapt, iStep, obs, rotate_unit)
    
        class(Small_Rotation), intent(inout) :: this
        integer, intent(in) :: Nadapt, iStep
        class(MoreObservables), intent(inout) :: obs
        integer, intent(in) :: rotate_unit
        
        obs%rotate_rejectAvg = obs%rotate_rejectAdapt / real(Nadapt-1, DP)
        obs%rotate_rejectAdapt = 0._DP
        call this%adapt_delta(obs%rotate_rejectAvg)
        write(rotate_unit, *) iStep, this%get_delta(), obs%rotate_rejectAvg
        
    end subroutine adapt_rotation
    
    !> Consistency test
    
    subroutine test_consist(Epot, Epot_conf, report_unit)
    
        real(DP), intent(in) :: Epot, Epot_conf
        integer, intent(in) :: report_unit
        
        real(DP) :: difference
        
        write(report_unit, *) "Consistency test: "
        write(report_unit, *) "    Epot = ", Epot
        write(report_unit, *) "    Epot_conf = ", Epot_conf
        
        if (abs(Epot_conf) < real_zero) then
            difference = abs(Epot_conf-Epot)
            write(report_unit, *) "    absolute difference = ", difference
        else
            difference = abs((Epot_conf-Epot)/Epot_conf)
            write(report_unit, *) "    relative difference = ", difference
        end if
        
        if (difference > consist_tiny) then ! not sufficient for HS ?
            write(report_unit, *) "    WARNING !"
        else
            write(report_unit, *) "    OK !"
        end if
    
    end subroutine test_consist

end module module_physics_macro
