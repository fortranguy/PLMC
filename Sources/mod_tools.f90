!> \brief Subroutines as tools for the main program

module mod_tools

use data_constants
use data_particles
use data_mc
use mod_physics

implicit none

contains

    !> Random number generator : seed
    
    subroutine init_random_seed(unitReport)
    
        integer, intent(in) :: unitReport
    
        integer :: i, n, clock
        integer, dimension(:), allocatable :: seed

        call random_seed(size = n)
        allocate(seed(n))

        call system_clock(count=clock)

        seed(:) = clock + 37 * [ (i - 1, i = 1, n) ]
        call random_seed(put = seed)
        
        write(unitReport, *) "RNG :"
        write(unitReport ,*) "    n = ", n
        write(unitReport ,*) "    seed(:) = ", seed(:)

        deallocate(seed)
        
    end subroutine init_random_seed
    
    !> Initial condition
    
    subroutine initialCondition(type1_X, type1_rMin, type2_X, type2_rMin, 
        mix_rMin, unitReport)
    
        real(DP), dimension(:, :), intent(inout) :: type1_X, type2_X
        real(DP), intent(in) :: type1_rMin, type2_rMin, mix_rMin
        integer, intent(in) :: unitReport
        

        character(len=20) :: init
        integer :: longueur, statut
        
        call get_command_argument(1, init, longueur, statut)
        if (statut /= 0) stop "error get_command_argument"
        if (command_argument_count() > 1) stop "Too many arguments"
            
        write(unitReport, *) "Initial condition :"
        
        select case (init)
            case ("rand")
                call randomDeposition(sph_X, sph_rMin)
                write(unitReport, *) "    Random deposition"
            case default
                write(*, *) "Enter the initial condition : "
                write(*, *) "   'rand'."
                stop
        end select
        
    end subroutine initialCondition
    
    !> Random deposition configuration
    
    subroutine randomDeposition(type1_X, type1_rMin, type2_X, type2_rMin, 
        mix_rMin)
    
        real(DP), dimension(:, :), intent(inout) :: type1_X, type2_X
        real(DP), intent(in) :: type1_rMin, type2_rMin, mix_rMin
    
        integer :: iCol, NcolOK, test_nOK
        integer :: type1_Ncol, type2_Ncol
        real(DP), dimension(Dim) :: xRand, xTest
        real(DP) :: rTest
        
        write(*, *) "Random deposition"
        
        ! Type 1
        
        !   First
        call random_number(xRand)
        type1_X(:, 1) = xRand*Lsize(:)
        NcolOK = 1
        
        !   Others
        type1_Ncol = size(type1_X, 2)
        do while (NcolOK < type1_Ncol)
        
            call random_number(xRand)
            xTest(:) = xRand(:)*Lsize(:)
            
            test_nOK = 0
            do iCol = 1, NcolOK
                rTest = dist(type1_X(:, iCol), xTest(:))
                if (rTest >= type1_rMin) then
                    test_nOK = test_nOK + 1
                else
                    exit
                end if
            end do
            
            if (test_nOK == NcolOK) then
                NcolOK = NcolOK + 1
                type1_X(:, NcolOK) = xTest(:)
                write(*, *) "    Type 1 particle n°", NcolOK, "OK"
            end if
            
        end do
        
        ! Type 2
        
        !   First
        test_nOK = 0
        do while (test_nOK == type1_Ncol)
        
            call random_number(xRand)
            xTest = xRand*Lsize(:)
            
            do iCol = 1, type1_Ncol
                rTest = dist(type1_X(:, iCol), xTest(:))
                if (rTest >= mix_rMin) then
                    test_nOK = test_nOK + 1
                else
                    exit
                end if
            end do
        
        end do        
        type2_X(:, 1) = xTest(:)        
        NcolOK = 1
        
        !   Others
        type2_Ncol = size(type2_X, 2)
        do while (NcolOK < type2_Ncol)
        
            call random_number(xRand)
            xTest(:) = xRand(:)*Lsize(:)
            
            nOK = 0
            do iCol = 1, NcolOK
                rTest = dist(type1_X(:, iCol), xTest(:))
                if (rTest >= mix_rMin) then
                    test_nOK = test_nOK + 1
                else
                    exit
                end if
            end do
            
            do iCol = 1, NcolOK
                rTest = dist(type2_X(:, iCol), xTest(:))
                if (rTest >= mix_rMin) then
                    test_nOK = test_nOK + 1
                else
                    exit
                end if
            end do
            
            if (test_nOK == NcolOK) then
                NcolOK = NcolOK + 1
                type2_X(:, NcolOK) = xTest(:)
                write(*, *) "    Type 2 particle n°", NcolOK, "OK"
            end if
            
        end do
        
        do iCol = 1, type1_Ncol
            type1_X(:, iCol) = type1_X(:, iCol) + type1_rMin/2._DP
        end do
        
        do iCol = 1, type2_Ncol
            type2_X(:, iCol) = type2_X(:, iCol) + type2_rMin/2._DP
        end do
    
    end subroutine randomDeposition
    
    subroutine report(unitReport)
    
    	integer, intent(in) :: unitReport
    
    	write(unitReport, *) "Simulation MC_C :"
        write(unitReport ,*) "    Lsize(:) = ", Lsize(:)
        write(unitReport ,*) "    Vol = ", product(Lsize)
        write(unitReport, *) "    Nstep = ", Nstep
        write(unitReport, *) "    Ntherm = ", Ntherm
        write(unitReport, *) "    Nmove = ", Nmove
    
    end subroutine report
    
end module mod_tools
