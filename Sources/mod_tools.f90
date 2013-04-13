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
    
    subroutine initialCondition(inter_X, unitReport)
    
        real(DP), dimension(:, :), intent(inout) :: inter_X
        integer, intent(in) :: unitReport
        
        real(DP) :: compac, densite
        character(len=20) :: init
        integer :: longueur, statut
        
        call get_command_argument(1, init, longueur, statut)
        if (statut /= 0) stop "error get_command_argument"
        if (command_argument_count() > 1) stop "Too many arguments"
            
        write(unitReport, *) "Initial condition :"
        
        select case (init)
            case ("cube")
                call primitiveCubic(inter_X)
                write(unitReport, *) "    Primitive cubic"
            case ("rand")
                call randomDeposition(inter_X)
                write(unitReport, *) "    Random deposition"
            case default
                write(*, *) "Enter the initial condition : "
                write(*, *) "   'cube' or 'rand'."
                stop
        end select
        
        densite = real(inter_Ncol, DP) / product(Lsize)
        write(*, *) "    Density = ", densite
        write(unitReport, *) "    Density = ", densite
        
        compac = 4._DP/3._DP*PI*inter_radius**3 * densite
        write(*, *) "    Compacity = ", compac
        write(unitReport, *) "    Compacity = ", compac
        
    end subroutine initialCondition
    
    !> Primitive cubic configuration
    
    subroutine primitiveCubic(inter_X)
    
        real(DP), dimension(:, :), intent(inout) :: inter_X
    
        integer :: iDir
        integer :: i, j, k, iCol
        integer, dimension(Dim) :: nCols
        real(DP), dimension(Dim) :: ratio
        real(DP) :: oneThird = 1._DP/3._DP
        
        write(*, *) "Primitive cubic"
        
        ! Proportion according to the direction       
        nCols(1) = int( (inter_Ncol*Lsize(1)**2/Lsize(2)/Lsize(3))**oneThird )
        nCols(2) = int( (inter_Ncol*Lsize(2)**2/Lsize(3)/Lsize(1))**oneThird )
        nCols(3) = int( (inter_Ncol*Lsize(3)**2/Lsize(1)/Lsize(2))**oneThird )
        
        ! Check
        iDir = 1
        do while (product(nCols)<inter_Ncol)
            nCols(iDir) = nCols(iDir) + 1
            iDir = iDir + 1
        end do
        
        ratio(:) = Lsize(:)/real(nCols(:), DP) ! A vérifier
        do iDir = 1, Dim
            if ( inter_rMin*real(nCols(iDir), DP) > Lsize(iDir) ) then
                write(*, *) "    Error : too dense in the direction", iDir
                stop
            end if
        end do
        
        ! Filling
        do k = 1, nCols(3)
            do j = 1, nCols(2)
                do i = 1, nCols(1)            
                    iCol = i + nCols(1)*(j-1) + nCols(1)*nCols(2)*(k-1)
                    if (iCol <= inter_Ncol) then
                        inter_X(1, iCol) = ratio(1)*real(i, DP)
                        inter_X(2, iCol) = ratio(2)*real(j, DP)
                        inter_X(3, iCol) = ratio(3)*real(k, DP)
                    end if
                end do
            end do
        end do
    
        do iDir = 1, Dim
            inter_X(iDir, :) = inter_X(iDir, :) - 0.5_DP*ratio(iDir) 
            ! just inside
        end do
    
    end subroutine primitiveCubic
    
    !> Random deposition configuration
    
    subroutine randomDeposition(inter_X)
    
        real(DP), dimension(:, :), intent(inout) :: inter_X
    
        integer :: iCol, Ncols, nOK
        real(DP), dimension(Dim) :: xTest
        real(DP) :: rTest
    
        write(*, *) "Random deposition"
    
        call random_number(inter_X(:, 1))
        inter_X(:, 1) = inter_X(:, 1)*(Lsize(:)-2*inter_radius)
        Ncols = 1        
        
        do while (Ncols<inter_Ncol)
        
            call random_number(xTest)
            xTest(:) = xTest(:)*(Lsize(:)-2._DP*inter_radius)
            
            nOK = 0
            do iCol = 1, Ncols
                rTest = dist(inter_X(:, iCol), xTest(:))
                if (rTest >= inter_rMin) then
                    nOK = nOK + 1
                else
                    exit
                end if
            end do
            
            if (nOK == Ncols) then
                Ncols = Ncols + 1
                inter_X(:, Ncols) = xTest(:)
                write(*, *) "    Particule", Ncols, "OK"
            end if
            
        end do
        
        do iCol = 1, inter_Ncol
            inter_X(:, iCol) = inter_X(:, iCol) + inter_radius
        end do
    
    end subroutine randomDeposition
    
end module mod_tools
