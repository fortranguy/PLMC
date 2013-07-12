!> \brief Description of the Sphere class

module class_spheres

use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
use data_precisions, only : DP
use data_constants, only : PI
use data_cell, only : Ndim, Lsize, Volume
use mod_physics, only : dist_PBC
use class_neighbourCells

implicit none
private

    type, public :: Spheres
    
        ! private
        ! The attributes must be private according to the encapsulation principle.
        ! Nevertheless, it is public for inheritance.
        ! The class must not be instanciated in the main program.
    
        character(len=5) :: name

        ! Particles
        real(DP) :: radius !< radius of a particle
        real(DP) :: rMin !< minimum distance between two particles
        integer ::  Ncol !< number of a component particles
        real(DP), dimension(:, :), allocatable, public :: positions !< positions of all particles
                                                                    !< Warning : use carefully !
        
        ! Snashot
        integer :: snap_factor

        ! Monte-Carlo
        real(DP), dimension(Ndim) :: deltaX !< displacement
        real(DP), dimension(Ndim) :: deltaXsave
        real(DP) :: rejectFix
        integer :: Nadapt
        integer :: Nwidom

        ! Potential
        real(DP) :: rCut !< short-range cut
        
        ! Neighbours (cell/grid scheme)
        type(NeighbourCells) :: same !< same kind
        type(NeighbourCells) :: mix !< other kind
        
    contains
    
        !> Accessors
        procedure :: getName => Spheres_getName
        procedure :: getNcol => Spheres_getNcol
        procedure :: getRmin => Spheres_getRmin
        procedure :: getNadapt => Spheres_getNadapt
        
        procedure :: printDensity => Spheres_printDensity
        
        !> Take a snap shot of the configuration : positions
        procedure :: snapShot_positions => Spheres_snapShot_positions
        
        !> Do an overlap test
        procedure :: overlapTest => Spheres_overlapTest
        
        !> Assign all particles to cells
        procedure :: all_cols_to_cells => Spheres_all_cols_to_cells
                
        !> Adapt the displacement deltaX during thermalisation
        procedure :: adaptDeltaX => Spheres_adaptDeltaX
        procedure :: definiteDeltaX => Spheres_definiteDeltaX
        procedure :: getDeltaX => Spheres_getDeltaX
        
    end type Spheres
    
contains

    !> Accessor : name

    pure function Spheres_getName(this) result(getName)
    
        class(Spheres), intent(in) :: this        
        character(len=5) :: getName
        
        getName = this%name
    
    end function Spheres_getName

    !> Accessor : Ncol

    pure function Spheres_getNcol(this) result(getNcol)
    
        class(Spheres), intent(in) :: this        
        integer :: getNcol
        
        getNcol = this%Ncol
    
    end function Spheres_getNcol
    
    !> Accessor : rMin
    
    pure function Spheres_getRmin(this) result(getRmin)
    
        class(Spheres), intent(in) :: this        
        real(DP) :: getRmin
        
        getRmin = this%rMin
    
    end function Spheres_getRmin
    
    !> Accessor : Nadapt
    
    pure function Spheres_getNadapt(this) result(getNadapt)
    
        class(Spheres), intent(in) :: this        
        integer :: getNadapt
        
        getNadapt = this%Nadapt
        
    end function Spheres_getNadapt
    
    !> Print density and compacity
    
    subroutine Spheres_printDensity(this, report_unit)
    
        class(Spheres), intent(in) :: this
        integer, intent(in) :: report_unit
        
        real(DP) :: density, compacity
        
        density = real(this%Ncol, DP) / Volume
        compacity = 4._DP/3._DP*PI*this%radius**3 * density
        
        write(output_unit, *) this%name, " : ", "density = ", density, "compacity = ", compacity
        
        write(report_unit, *) "    density = ", density
        write(report_unit, *) "    compacity = ", compacity
    
    end subroutine Spheres_printDensity
    
    !> Configuration state : positions
      
    subroutine Spheres_snapShot_positions(this, iStep, snap_unit)
        
        class(Spheres), intent(in) :: this
        integer, intent(in) :: iStep
        integer, intent(in) :: snap_unit
    
        integer :: iCol
        
        if (modulo(iStep, this%snap_factor) == 0) then
        
            do iCol = 1, this%Ncol
                write(snap_unit, *) this%positions(:, iCol)
            end do
            
        end if            

    end subroutine Spheres_snapShot_positions
    
    !> Overlapt test
    
    subroutine Spheres_overlapTest(this)
    
        class(Spheres), intent(in) :: this
    
        integer :: jCol, iCol
        real(DP) :: r_ij
    
        do jCol = 1, this%Ncol
            do iCol = 1, this%Ncol
                if (iCol /= jCol) then
                    
                    r_ij = dist_PBC(this%positions(:, iCol), this%positions(:, jCol))
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
    
    subroutine Spheres_all_cols_to_cells(this, other)
    
        class(Spheres), intent(inout) :: this
        class(Spheres), intent(in) :: other
        
        call this%same%all_cols_to_cells(this%Ncol, this%positions)
        call this%mix%all_cols_to_cells(other%Ncol, other%positions)
    
    end subroutine Spheres_all_cols_to_cells
    
    !> Adaptation of deltaX during the thermalisation
    
    subroutine Spheres_adaptDeltaX(this, reject)
    
        class(Spheres), intent(inout) :: this
        real(DP), intent(in) :: reject
        
        real(DP), parameter :: eps_deltaX = 0.05_DP
        real(DP), parameter :: eps_reject = 0.1_DP * eps_deltaX
        real(DP), parameter :: more = 1._DP+eps_deltaX
        real(DP), parameter :: less = 1._DP-eps_deltaX
        
        if (reject < this%rejectFix - eps_reject) then
        
            this%deltaX(:) = this%deltaX(:) * more
            
            if (norm2(this%deltaX) > norm2(Lsize)) then
                this%deltaX(:) = Lsize(:)
            end if
            
        else if (reject > this%rejectFix + eps_reject) then
        
            this%deltaX(:) = this%deltaX(:) * less
            
        end if
    
    end subroutine Spheres_adaptDeltaX
    
    subroutine Spheres_definiteDeltaX(this, reject, report_unit)
    
        class(Spheres), intent(inout) :: this    
        real(DP), intent(in) :: reject
        integer, intent(in) :: report_unit

        if (reject == 0._DP) then
            write(error_unit, *) this%name, " :    Warning : deltaX adaptation problem."
            this%deltaX(:) = this%deltaXsave(:)
            write(error_unit, *) "default deltaX :", this%deltaX(:)
        end if

        if (norm2(this%deltaX) > norm2(Lsize)) then
            write(error_unit, *) this%name, " :   Warning : deltaX too big."
            this%deltaX(:) = Lsize(:)
            write(error_unit, *) "big deltaX :", this%deltaX(:)
        end if

        write(output_unit, *) this%name, " :    Thermalisation : over"

        write(report_unit, *) "Displacement :"
        write(report_unit, *) "    deltaX(:) = ", this%deltaX(:)
        write(report_unit, *) "    rejection relative difference = ", &
                                    abs(reject-this%rejectFix)/this%rejectFix
    
    end subroutine Spheres_definiteDeltaX
    
    pure function Spheres_getDeltaX(this) result(getDeltaX)
        
        class(Spheres), intent(in) :: this        
        real(DP) :: getDeltaX
        
        ! average deltaX of 3 vector components
        getDeltaX = sum(this%deltaX)/size(this%deltaX)
        
    end function Spheres_getDeltaX

end module class_spheres
