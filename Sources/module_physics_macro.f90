!> Subroutines for Physics / macro: after

module module_physics_macro

use, intrinsic :: iso_fortran_env, only : output_unit, error_unit, iostat_end
use data_precisions, only: DP, io_tiny
use data_box, only: Ndim, Lsize
use module_types, only: argument_seed, argument_initial
use module_physics_micro, only: dist_PBC, random_surface
use class_hardSpheres
use class_dipolarSpheres
use class_mixingPotential

implicit none
private
public init_randomSeed, set_initialCondition

contains

    !> Random number generator : seed
    
    subroutine init_randomSeed(arg_seed, report_unit)
    
        type(argument_seed) :: arg_seed
        integer, intent(in) :: report_unit
    
        integer :: i, n, clock
        integer, dimension(:), allocatable :: seed

        call random_seed(size = n)
        allocate(seed(n))

        select case (arg_seed%choice)
        
            case ('v')

                call system_clock(count=clock)
                seed(:) = clock + 37 * [ (i - 1, i = 1, n) ]
                call random_seed(put = seed)
                ! not sufficent ? cf.
                ! http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html
                write(report_unit, *) "Random number generator: variable"

            case ('p')
                call random_seed(put = arg_seed%seed)
                write(report_unit, *) "Random number generator: put"
                
            case ('f')
                write(report_unit, *) "Random number generator: fix"

            case default
                error stop "Error : init_randomSeed"

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
        
        write(report_unit, *) "Initial condition :"
        
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
                error stop "Error : set_initialCondition"
                
        end select
        
    end subroutine set_initialCondition

end module module_physics_macro
