!> \brief Description of the Sphere class

module class_spheres

use data_constants
use data_cell
use data_mc
use data_neighbours
use mod_physics
use class_neighbours

implicit none

private

    type, public :: Spheres
    
        !private

        ! Particles
        real(DP) :: radius !< radius of a particle
        real(DP) :: rMin !< minimum distance between two particles
        integer ::  Ncol !< number of a component particles
        real(DP), dimension(:, :), allocatable :: X !< position of a particle

        ! Monte-Carlo
        real(DP), dimension(Dim) :: dx !< displacement
        integer :: Nwidom

        ! Potential
        real(DP) :: rCut !< short-range cut    
        
        ! Neighbours (cell/grid scheme)   
        type(Neighbours) :: same !< same kind
        
    contains
        
        !> Take a snap shot of the configuration
        procedure :: snapShot => Spheres_snapShot
        
        !> Do an overlap test
        procedure :: overlapTest => Spheres_overlapTest
        
        !> Assign all particles to cells
        procedure :: cols_to_cells => Spheres_cols_to_cells
                
        !> Adapt the displacement dx during thermalisation
        procedure :: adaptDx => Spheres_adaptDx
        procedure :: getDx => Spheres_getDx
        
    end type Spheres
    
contains
    
    !> Configuration state
      
    subroutine Spheres_snapShot(this, unitSnap)
        
        class(Spheres), intent(in) :: this
        integer, intent(in) :: unitSnap
    
        integer :: iCol
        
        do iCol = 1, this%Ncol
            write(unitSnap, *) this%X(:, iCol)
        end do    

    end subroutine Spheres_snapShot
    
    !> Overlapt test
    
    subroutine Spheres_overlapTest(this)
    
        class(Spheres), intent(in) :: this
    
        integer :: jCol, iCol
        real(DP) :: r_ij
    
        do jCol = 1, this%Ncol
            do iCol = 1, this%Ncol
                if (iCol /= jCol) then
                    
                    r_ij = dist(this%X(:, iCol), this%X(:, jCol))
                    if (r_ij < this%rMin) then
                        write(*, *) "    Overlap !", iCol, jCol
                        write(*, * ) "    r_ij = ", r_ij
                        stop
                    end if
                    
                end if
            end do
        end do
        
        write(*, *) "    Overlap test : OK !"
    
    end subroutine Spheres_overlapTest
    
    !> Fill cells with colloids
    
    subroutine Spheres_cols_to_cells(this)
    
        class(Spheres), intent(inout) :: this
        
        call this%same%all_col_to_cell(this%Ncol, this%X)
    
    end subroutine Spheres_cols_to_cells
    
    !> Adaptation of dx during the thermalisation
    
    subroutine Spheres_adaptDx(this, iStep, rejectsRateSum, unitReport)
    
        class(Spheres), intent(inout) :: this 
        integer, intent(in) :: iStep, unitReport
        real(DP), intent(in) :: rejectsRateSum    
        
        integer, parameter :: multiple = 2**2
        real(DP) :: rejectsRate
        real(DP), parameter :: rejectsRateFix = 0.5_DP
        real(DP), parameter :: dx_eps = 0.05_DP, taux_eps = 0.05_DP
        real(DP), parameter :: more = 1._DP+dx_eps, less = 1._DP-dx_eps
        
        rejectsRate = 0._DP
        
        if (mod(iStep, multiple) == 0 .and. iStep>2) then
        
            rejectsRate = rejectsRateSum/real(iStep-1, DP)
        
            if (rejectsRate < rejectsRateFix - taux_eps) then            
                this%dx(:) = this%dx(:) * more
                this%dx(:) = modulo(this%dx(:), Lsize(:))
            else if (rejectsRate > rejectsRateFix + taux_eps) then
                this%dx(:) = this%dx(:) * less
                this%dx(:) = modulo(this%dx(:), Lsize(:))
            end if

        end if
        
        if (iStep == Ntherm) then
        
            if (rejectsRate == 0._DP) then
                write(*, *) "Error : dx adaptation problem"
                stop
            end if
            
            write(unitReport, *) "Displacement :"
            write(unitReport, *) "    dx(:) = ", this%dx(:)
            write(unitReport, *) "    rejection relative difference = ", &
            							! wrong translation ?
                abs(rejectsRate - rejectsRateFix)/rejectsRateFix
            
        end if
    
    end subroutine Spheres_adaptDx
    
    function Spheres_getDx(this)
        
        class(Spheres), intent(in) :: this
        
        real(DP) :: Spheres_getDx
        ! average dx of 3 vector components
        Spheres_getDx = sum(this%dx)/size(this%dx)
        
    end function Spheres_getDx    

end module class_spheres
