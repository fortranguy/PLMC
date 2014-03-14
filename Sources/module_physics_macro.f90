!> Subroutines for Physics / macro: after

module module_physics_macro

use, intrinsic :: iso_fortran_env, only: output_unit, error_unit, iostat_end
use data_precisions, only: DP, real_zero, io_tiny, consist_tiny
use data_box, only: Ndim, Lsize
use data_monteCarlo, only: Nadapt
use data_potential, only: write_potential
use module_types, only: argument_random, argument_initial
use module_physics_micro, only: dist_PBC, random_surface
use class_hardSpheres
use class_dipolarSpheres
use class_mixingPotential
use class_observables
use class_units

implicit none
private
public init_randomSeed, set_initialCondition, init, final, mix_init, mix_final, &
       adapt_move, adapt_rotate, test_consist

contains

    !> Random number generator: seed
    
    subroutine init_randomSeed(arg_rand, report_unit)
    
        type(argument_random) :: arg_rand
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
                error stop "Error: init_randomSeed"

        end select

        call random_seed(get = seed)
        write(report_unit ,*) "    size = ", n
        write(report_unit ,*) "    seed = ", seed(:)

        deallocate(seed)
        
    end subroutine init_randomSeed

    !> Random depositions configuration
    
    subroutine randomDepositions(type1, type2, mix_sigma)

        class(HardSpheres), intent(inout) :: type1, type2
        real(DP), intent(in) :: mix_sigma

        integer :: iCol, iColTest
        real(DP), dimension(Ndim) :: xRand
        real(DP) :: rTest
        
        ! Type 1
        do iCol = 1, type1%get_Ncol()
        
7101        continue
            call random_number(xRand)
            type1%positions(:, iCol) = xRand*Lsize(:)
            
            do iColTest = 1, iCol-1
                rTest = dist_PBC(type1%positions(:, iColTest), type1%positions(:, iCol))
                if (rTest < type1%get_sigma()) then
                    goto 7101
                end if
            end do
        
        end do
        
        ! Type 2
        do iCol = 1, type2%get_Ncol()
        
7102        continue
            call random_number(xRand)
            type2%positions(:, iCol) = xRand*Lsize(:)
            
            do iColTest = 1, type1%get_Ncol()
                rTest = dist_PBC(type1%positions(:, iColTest), type2%positions(:, iCol))
                if (rTest < mix_sigma) then
                    goto 7102
                end if
            end do
            
            do iColTest = 1, iCol-1
                rTest = dist_PBC(type2%positions(:, iColTest), type2%positions(:, iCol))
                if (rTest < type2%get_sigma()) then
                    goto 7102
                end if
            end do
        
        end do
    
    end subroutine randomDepositions
    
    !> Uniform (gaussian) orientations
    
    subroutine randomOrientations(orientations, Ncol)
    
        real(DP), dimension(:, :), intent(out) :: orientations
        integer, intent(in) :: Ncol
        
        integer :: iCol
        
        do iCol = 1, Ncol
            orientations(:, iCol) = random_surface()
        end do
    
    end subroutine randomOrientations
    
    !> From an old configuration
    
    subroutine oldConfiguration(file, length, type_name, type_Ncol, type_coords, normMax)
    
        character(len=*), intent(in) :: file
        integer, intent(in) :: length
        character(len=*), intent(in) :: type_name
        integer, intent(in) :: type_Ncol
        real(DP), dimension(:, :), intent(out) :: type_coords
        real(DP), intent(in) :: normMax
        
        integer :: file_unit, readStat

        integer :: iCol
        real(DP), dimension(Ndim) :: vecDummy
        
        write(output_unit, *) type_name, " <- ", file(1:length)
        open(newunit=file_unit, recl=4096, file=file(1:length), status='old', action='read')
        
        iCol = 0
        do
            read(file_unit, fmt=*, iostat=readStat) vecDummy(:)
            if (readStat == iostat_end) exit
            iCol = iCol + 1
        end do
        
        if (iCol == type_Ncol) then
            rewind(file_unit)
            do iCol = 1, type_Ncol
                read(file_unit, *) type_coords(:, iCol)
                if (norm2(type_coords(:, iCol)) > normMax+io_tiny) then
                    write(error_unit, *) "Norm error in file: ", file(1:length)
                    write(error_unit, *) "Coordinates ", type_coords(:, iCol)
                    write(error_unit, *) "Norm =", norm2(type_coords(:, iCol))
                    error stop
                end if
            end do
        else
            write(error_unit, *) "Error reading: ", file(1:length)
            write(error_unit, *) "iCol", iCol, " /= ", "type_Ncol", type_Ncol
            error stop
        end if
        
        close(file_unit)
        
    end subroutine oldConfiguration
    
    !> Initial condition
    
    subroutine set_initialCondition(arg_init, dipolar, spherical, mix_sigma, report_unit)

        type(argument_initial), intent(in) :: arg_init
        class(DipolarSpheres), intent(inout) :: dipolar
        class(HardSpheres), intent(inout) :: spherical
        real(DP), intent(in) :: mix_sigma
        integer, intent(in) :: report_unit
        
        write(report_unit, *) "Initial condition: "
        
        select case (arg_init%choice)
        
            case ('r')
                call randomDepositions(dipolar, spherical, mix_sigma)
                call randomOrientations(dipolar%orientations, dipolar%get_Ncol())
                write(output_unit, *) "Random depositions + random orientations"
                write(report_unit, *) "    Random depositions + random orientations"
                
            case ('f')
                write(output_unit, *) "Old configuration"
                write(report_unit, *) "    Old configuration"
                call oldConfiguration(arg_init%files(1), arg_init%length(1), &
                                      dipolar%get_name()//"_positions", dipolar%get_Ncol(), &
                                      dipolar%positions, norm2(Lsize))
                call oldConfiguration(arg_init%files(2), arg_init%length(2), &
                                      dipolar%get_name()//"_orientations", dipolar%get_Ncol(), &
                                      dipolar%orientations, 1._DP)
                call oldConfiguration(arg_init%files(3), arg_init%length(3), &
                                      spherical%get_name()//"_positions", spherical%get_Ncol(), &
                                      spherical%positions, norm2(Lsize))
            
            case default
                error stop "Error: set_initialCondition"
                
        end select
        
    end subroutine set_initialCondition
    
    !> Spheres initialisations
    
    subroutine init(this, other, mix, this_units, this_Epot)
    
        class(HardSpheres), intent(inout) :: this
        class(HardSpheres), intent(in) :: other
        class(MixingPotential), intent(in) :: mix
        class(Units), intent(in) :: this_units
        real(DP), intent(inout) :: this_Epot
        
        call this%test_overlap()
        call this%snap_data(this_units%snap_positions)
        call this%snap_positions(0, this_units%snapIni_positions)
        call this%set_Epot()
        
        if (write_potential) then
            call this%write_Epot(this_units%Epot)
        end if
        select type (this)
            type is (DipolarSpheres)
                select type (this_units)
                    type is (MoreUnits)
                        call this%snap_data(this_units%snap_orientations)
                        call this%snap_orientations(0, this_units%snapIni_orientations)
                        if (write_potential) then
                            call this%write_Epot_real(this_units%Epot_real)
                        end if
                        call this%Epot_reci_count_waveVectors(this_units%waveVectors)
                end select
        end select
        this_Epot = this%Epot_conf()
        
        call this%construct_cells(other, mix%get_cell_size(), mix%get_rCut())
        call this%write_report(this_units%report)
    
    end subroutine init
    
    !> Spheres finalizations
    
    subroutine final(this, this_units, this_obs)
    
        class(HardSpheres), intent(inout) :: this
        class(Units), intent(in) :: this_units
        class(Observables), intent(in) :: this_obs
        
        call this%test_overlap()
        call this%set_Epot()
        call test_consist(this_obs%Epot, this%Epot_conf(), this_units%report)
        call this%snap_positions(0, this_units%snapFin_positions)
        call this_obs%write_results(this_units%report)
        
        select type (this)
            type is (DipolarSpheres)
                select type (this_units)
                    type is (MoreUnits)
                        call this%snap_orientations(0, this_units%snapFin_orientations)
                end select
        end select
    
    end subroutine final
    
    !> Mix initialisation
    
    subroutine mix_init(mix, type1, type2, mix_Epot_unit, mix_Epot)
    
        class(MixingPotential), intent(inout) :: mix
        class(HardSpheres), intent(in) :: type1, type2
        integer, intent(in) :: mix_Epot_unit
        real(DP), intent(out) :: mix_Epot
    
        call mix%test_overlap(type1, type2)
        call mix%set_Epot()
        if (write_potential) then
            call mix%write_Epot(mix_Epot_unit)
        end if
        call mix%set_cell_size()
        mix_Epot = mix%Epot_conf(type1, type2)
    
    end subroutine mix_init
    
    !> Mix finalization
    
    subroutine mix_final(mix, type1, type2, mix_report_unit, mix_Epot, mix_Epot_conf)
    
        class(MixingPotential), intent(inout) :: mix
        class(HardSpheres), intent(in) :: type1, type2
        integer, intent(in) :: mix_report_unit
        real(DP), intent(in) :: mix_Epot
        real(DP), intent(out) :: mix_Epot_conf
        
        call mix%test_overlap(type1, type2)
        call mix%set_Epot()
        mix_Epot_conf = mix%Epot_conf(type1, type2)
        call test_consist(mix_Epot, mix_Epot_conf, mix_report_unit)
    
    end subroutine mix_final
    
    !> Change: average & adaptation
    
    subroutine adapt_move(this, iStep, obs, move_unit)
    
        class(HardSpheres), intent(inout) :: this
        integer, intent(in) :: iStep
        class(Observables), intent(inout) :: obs
        integer, intent(in) :: move_unit
    
        obs%move_rejectAvg = obs%move_rejectAdapt / real(Nadapt-1, DP)
        obs%move_rejectAdapt = 0._DP
        call this%adapt_move_delta(obs%move_rejectAvg)
        write(move_unit, *) iStep, this%get_move_delta(), obs%move_rejectAvg
    
    end subroutine adapt_move
    
    subroutine adapt_rotate(this, iStep, obs, rotate_unit)
    
        class(DipolarSpheres), intent(inout) :: this
        integer, intent(in) :: iStep
        class(MoreObservables), intent(inout) :: obs
        integer, intent(in) :: rotate_unit
        
        obs%rotate_rejectAvg = obs%rotate_rejectAdapt / real(Nadapt-1, DP)
        obs%rotate_rejectAdapt = 0._DP
        call this%adapt_rotate_delta(obs%rotate_rejectAvg)
        write(rotate_unit, *) iStep, this%get_rotate_delta(), obs%rotate_rejectAvg
        
    end subroutine adapt_rotate
    
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
