!> \brief Description of the Sphere class

module class_spheres

use iso_fortran_env
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
        character(len=5) :: name

        ! Particles
        real(DP) :: radius !< radius of a particle
        real(DP) :: rMin !< minimum distance between two particles
        integer ::  Ncol !< number of a component particles
        real(DP), dimension(:, :), allocatable, public :: X !< positions of all particles
        
        ! Snashot
        integer :: snap_factor

        ! Monte-Carlo
        real(DP), dimension(Dim) :: dx !< displacement
        real(DP), dimension(Dim) :: dxSave
        real(DP) :: rejFix
        integer :: Nadapt
        integer :: Nwidom

        ! Potential
        real(DP) :: rCut !< short-range cut
        
        ! Neighbours (cell/grid scheme)
        type(Neighbours) :: same !< same kind
        type(Neighbours) :: mix !< other kind
        
    contains
    
        !> Accessors
        procedure :: getName => Spheres_getName
        procedure :: getNcol => Spheres_getNcol
        procedure :: getRmin => Spheres_getRmin
        procedure :: getNadapt => Spheres_getNadapt
        
        procedure :: printInfo => Spheres_printInfo
        
        !> Take a snap shot of the configuration : positions
        procedure :: snapShot_X => Spheres_snapShot_X
        
        !> Do an overlap test
        procedure :: overlapTest => Spheres_overlapTest
        
        !> Assign all particles to cells
        procedure :: cols_to_cells => Spheres_cols_to_cells
                
        !> Adapt the displacement dx during thermalisation
        procedure :: adaptDx => Spheres_adaptDx
        procedure :: definiteDx => Spheres_definiteDx
        procedure :: getDx => Spheres_getDx
        
    end type Spheres
    
contains

    !> Accessor : name

    function Spheres_getName(this) result(getName)
    
        class(Spheres), intent(in) :: this        
        character(len=5) :: getName
        
        getName = this%name
    
    end function Spheres_getName

    !> Accessor : Ncol

    function Spheres_getNcol(this) result(getNcol)
    
        class(Spheres), intent(in) :: this        
        integer :: getNcol
        
        getNcol = this%Ncol
    
    end function Spheres_getNcol
    
    !> Accessor : rMin
    
    function Spheres_getRmin(this) result(getRmin)
    
        class(Spheres), intent(in) :: this        
        real(DP) :: getRmin
        
        getRmin = this%rMin
    
    end function Spheres_getRmin
    
    !> Accessor : Nadapt
    
    function Spheres_getNadapt(this) result(getNadapt)
    
        class(Spheres), intent(in) :: this        
        integer :: getNadapt
        
        getNadapt = this%Nadapt
        
    end function Spheres_getNadapt
    
    !> Print density and compacity
    
    subroutine Spheres_printInfo(this, report_unit)
    
        class(Spheres), intent(in) :: this
        integer, intent(in) :: report_unit
        
        real(DP) :: density, compac
        
        density = real(this%Ncol, DP) / product(Lsize)
        compac = 4._DP/3._DP*PI*this%radius**3 * density
        
        write(output_unit, *) this%name, " : ", "density = ", density, "compacity = ", compac
        
        write(report_unit, *) "    density = ", density
        write(report_unit, *) "    compacity = ", compac
    
    end subroutine Spheres_printInfo
    
    !> Configuration state : positions
      
    subroutine Spheres_snapShot_X(this, iStep, snap_unit)
        
        class(Spheres), intent(in) :: this
        integer, intent(in) :: iStep
        integer, intent(in) :: snap_unit
    
        integer :: iCol
        
        if (modulo(iStep, this%snap_factor) == 0) then
        
            do iCol = 1, this%Ncol
                write(snap_unit, *) this%X(:, iCol)
            end do
            
        end if            

    end subroutine Spheres_snapShot_X
    
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
                        write(error_unit, *) this%name, "    Overlap !", iCol, jCol
                        write(error_unit, *) "    r_ij = ", r_ij
                        stop
                    end if
                    
                end if
            end do
        end do
        
        write(output_unit, *) this%name, " :    Overlap test : OK !"
    
    end subroutine Spheres_overlapTest
    
    !> Fill cells with colloids
    
    subroutine Spheres_cols_to_cells(this, other_X)
    
        class(Spheres), intent(inout) :: this
        real(DP), dimension(:, :), intent(in) :: other_X
        
        integer :: other_Ncol 
        
        call this%same%cols_to_cells(this%Ncol, this%X)
        
        other_Ncol = size(other_X, 2)
        call this%mix%cols_to_cells(other_Ncol, other_X)
    
    end subroutine Spheres_cols_to_cells
    
    !> Adaptation of dx during the thermalisation
    
    subroutine Spheres_adaptDx(this, rej)
    
        class(Spheres), intent(inout) :: this
        real(DP), intent(in) :: rej
        
        real(DP), parameter :: eps_dx = 0.05_DP
        real(DP), parameter :: eps_rej = 0.1_DP * eps_dx
        real(DP), parameter :: more = 1._DP+eps_dx
        real(DP), parameter :: less = 1._DP-eps_dx
        
        real(DP) :: dx_normSqr, Lsize_normSqr
        
        Lsize_normSqr = dot_product(Lsize, Lsize)
        
        if (rej < this%rejFix - eps_rej) then
        
            this%dx(:) = this%dx(:) * more
            
            dx_normSqr = dot_product(this%dx, this%dx)
            if (dx_normSqr > Lsize_normSqr) then
                this%dx(:) = Lsize(:)
            end if
            
        else if (rej > this%rejFix + eps_rej) then
        
            this%dx(:) = this%dx(:) * less
            
        end if
    
    end subroutine Spheres_adaptDx
    
    subroutine Spheres_definiteDx(this, rej, report_unit)
    
        class(Spheres), intent(inout) :: this    
        real(DP), intent(in) :: rej
        integer, intent(in) :: report_unit
        
        real(DP) :: dx_normSqr, Lsize_normSqr
        
            if (rej == 0._DP) then
                write(error_unit, *) this%name, " :    Warning : dx adaptation problem."
                this%dx(:) = this%dxSave(:)
                write(error_unit, *) "default dx :", this%dx(:)
            end if
            
            dx_normSqr = dot_product(this%dx, this%dx)
            Lsize_normSqr = dot_product(Lsize, Lsize)
            if (dx_normSqr >= Lsize_normSqr) then
                write(error_unit, *) this%name, " :   Warning : dx too big."
                this%dx(:) = Lsize(:)
                write(error_unit, *) "big dx :", this%dx(:)
            end if
            
            write(output_unit, *) this%name, " :    Thermalisation : over"
            
            write(report_unit, *) "Displacement :"
            write(report_unit, *) "    dx(:) = ", this%dx(:)
            write(report_unit, *) "    rejection relative difference = ", &
                                       abs(rej-this%rejFix)/this%rejFix
    
    end subroutine Spheres_definiteDx
    
    function Spheres_getDx(this) result(getDx)
        
        class(Spheres), intent(in) :: this        
        real(DP) :: getDx
        
        ! average dx of 3 vector components
        getDx = sum(this%dx)/size(this%dx)
        
    end function Spheres_getDx    

end module class_spheres
