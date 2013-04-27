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
        real(DP), dimension(:, :), allocatable :: X !< position of a particle

        ! Monte-Carlo
        real(DP), dimension(Dim) :: dx !< displacement
        real(DP), dimension(Dim) :: dx_save
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
        
        procedure :: printInfo => Spheres_printInfo
        
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
    
    !> Print density and compacity
    
    subroutine Spheres_printInfo(this, report_unit)
    
        class(Spheres), intent(in) :: this
        integer, intent(in) :: report_unit
        
        real(DP) :: density, compac
        
        write(output_unit, *) this%name, " :"
        
        density = real(this%Ncol, DP) / product(Lsize)
        write(output_unit, *) "    density = ", density
        write(report_unit, *) "    density = ", density
        
        compac = 4._DP/3._DP*PI*this%radius**3 * density
        write(output_unit, *) "    compacity = ", compac
        write(report_unit, *) "    compacity = ", compac
    
    end subroutine Spheres_printInfo
    
    !> Configuration state
      
    subroutine Spheres_snapShot(this, snap_unit)
        
        class(Spheres), intent(in) :: this
        integer, intent(in) :: snap_unit
    
        integer :: iCol
        
        do iCol = 1, this%Ncol
            write(snap_unit, *) this%X(:, iCol)
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
                        write(output_unit, *) this%name, "    Overlap !", iCol, jCol
                        write(output_unit, *) "    r_ij = ", r_ij
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
    
    subroutine Spheres_adaptDx(this, iStep, rej, report_unit)
    
        class(Spheres), intent(inout) :: this 
        integer, intent(in) :: iStep
        real(DP), intent(in) :: rej
        integer, intent(in) :: report_unit
        
        real(DP), parameter :: rejFix = 0.5_DP
        real(DP), parameter :: dx_eps = 0.05_DP, taux_eps = 0.05_DP
        real(DP), parameter :: more = 1._DP+dx_eps, less = 1._DP-dx_eps
        
        if (mod(iStep, this%Nadapt) == 0 .and. iStep>2) then
        
            if (rej < rejFix - taux_eps) then            
                this%dx(:) = this%dx(:) * more
                this%dx(:) = modulo(this%dx(:), Lsize(:))
            else if (rej > rejFix + taux_eps) then
                this%dx(:) = this%dx(:) * less
                this%dx(:) = modulo(this%dx(:), Lsize(:))
            end if

        end if
        
        if (iStep == Ntherm) then
        
            if (rej == 0._DP) then
                write(output_unit, *) this%name, "    Warning : dx adaptation problem."
                this%dx(:) = this%dx_save(:)
                write(output_unit, *) "default dx :", this%dx(:)
                write(output_unit, *)
            end if
            
            write(output_unit, *) this%name, " :    Thermalisation : over"
            
            write(report_unit, *) "Displacement :"
            write(report_unit, *) "    dx(:) = ", this%dx(:)
            write(report_unit, *) "    rejection relative difference = ", abs(rej-rejFix)/rejFix
                                       ! wrong translation ?
                                       
        end if
    
    end subroutine Spheres_adaptDx
    
    function Spheres_getDx(this)
        
        class(Spheres), intent(in) :: this
        
        real(DP) :: Spheres_getDx
        ! average dx of 3 vector components
        Spheres_getDx = sum(this%dx)/size(this%dx)
        
    end function Spheres_getDx    

end module class_spheres
