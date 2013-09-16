!> \brief Description of the Sphere class

module class_spheres

use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
use data_precisions, only : DP
use data_constants, only : PI
use data_box, only : Ndim, Lsize, Volume
use module_physics, only : dist_PBC
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
        real(DP) :: rMin !< minimum distance between two particles
        real(DP) :: radius !< radius of a particle
        integer ::  Ncol !< number of a component particles
        real(DP), dimension(:, :), allocatable, public :: positions !< positions of all particles
                                                                    !< Warning : use carefully !
        
        ! Snashot
        integer :: snap_factor

        ! Monte-Carlo
        real(DP), dimension(Ndim) :: move_delta !< displacement
        real(DP), dimension(Ndim) :: move_deltaSave
        real(DP) :: move_rejectFix
        integer :: move_Nadapt
        integer :: Nwidom

        ! Potential
        real(DP) :: rCut !< short-range cut
        
        ! Neighbours (cell/grid scheme)
        type(NeighbourCells) :: sameCells !< same kind
        type(NeighbourCells) :: mixCells !< other kind
        
    contains
    
        !> Accessors
        procedure :: getName => Spheres_getName
        procedure :: getNcol => Spheres_getNcol
        procedure :: getRmin => Spheres_getRmin
        procedure :: getRcut => Spheres_getRcut        
        procedure :: getMove_Nadapt => Spheres_getMove_Nadapt
        
        procedure :: printDensity => Spheres_printDensity
        
        !> Take a snap shot of the configuration : positions
        procedure :: snapShot_positions => Spheres_snapShot_positions
        
        !> Do an overlap test
        procedure :: overlapTest => Spheres_overlapTest
        
        !> Assign all particles to cells
        procedure :: all_cols_to_cells => Spheres_all_cols_to_cells
                
        !> Adapt the displacement move_delta during thermalisation
        procedure :: adaptMove_delta => Spheres_adaptMove_delta
        procedure :: definiteMove_delta => Spheres_definiteMove_delta
        procedure :: getMove_delta => Spheres_getMove_delta
        
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
    
    !> Accessor : rCut
    
    pure function Spheres_getRcut(this) result(getRcut)
    
        class(Spheres), intent(in) :: this        
        real(DP) :: getRcut
        
        getRcut = this%rCut
    
    end function Spheres_getRcut
    
    !> Accessor : move_Nadapt
    
    pure function Spheres_getMove_Nadapt(this) result(getMove_Nadapt)
    
        class(Spheres), intent(in) :: this        
        integer :: getMove_Nadapt
        
        getMove_Nadapt = this%move_Nadapt
        
    end function Spheres_getMove_Nadapt
    
    !> Print density and compacity
    
    subroutine Spheres_printDensity(this, report_unit)
    
        class(Spheres), intent(in) :: this
        integer, intent(in) :: report_unit
        
        real(DP) :: density, compacity
        
        density = real(this%Ncol + 1, DP) / Volume ! cheating ? cf. Widom
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
        
        call this%sameCells%all_cols_to_cells(this%Ncol, this%positions)
        call this%mixCells%all_cols_to_cells(other%Ncol, other%positions)
    
    end subroutine Spheres_all_cols_to_cells
    
    !> Adaptation of move_delta during the thermalisation
    
    subroutine Spheres_adaptMove_delta(this, reject)
    
        class(Spheres), intent(inout) :: this
        real(DP), intent(in) :: reject
        
        real(DP), parameter :: move_delta_eps = 0.05_DP
        real(DP), parameter :: move_reject_eps = 0.1_DP * move_delta_eps
        real(DP), parameter :: more = 1._DP+move_delta_eps
        real(DP), parameter :: less = 1._DP-move_delta_eps
        
        if (reject < this%move_rejectFix - move_reject_eps) then
        
            this%move_delta(:) = this%move_delta(:) * more
            
            if (norm2(this%move_delta) > norm2(Lsize)) then
                this%move_delta(:) = Lsize(:)
            end if
            
        else if (reject > this%move_rejectFix + move_reject_eps) then
        
            this%move_delta(:) = this%move_delta(:) * less
            
        end if
    
    end subroutine Spheres_adaptMove_delta
    
    subroutine Spheres_definiteMove_delta(this, reject, report_unit)
    
        class(Spheres), intent(inout) :: this    
        real(DP), intent(in) :: reject
        integer, intent(in) :: report_unit

        if (reject == 0._DP) then
            write(error_unit, *) this%name, " :    Warning : move_delta adaptation problem."
            this%move_delta(:) = this%move_deltaSave(:)
            write(error_unit, *) "default move_delta :", this%move_delta(:)
        end if

        if (norm2(this%move_delta) > norm2(Lsize)) then
            write(error_unit, *) this%name, " :   Warning : move_delta too big."
            this%move_delta(:) = Lsize(:)
            write(error_unit, *) "big move_delta :", this%move_delta(:)
        end if

        write(output_unit, *) this%name, " :    Thermalisation : over"

        write(report_unit, *) "Displacement :"
        write(report_unit, *) "    move_delta(:) = ", this%move_delta(:)
        write(report_unit, *) "    rejection relative difference = ", &
                                    abs(reject-this%move_rejectFix)/this%move_rejectFix
    
    end subroutine Spheres_definiteMove_delta
    
    pure function Spheres_getMove_delta(this) result(getMove_delta)
        
        class(Spheres), intent(in) :: this        
        real(DP) :: getMove_delta
        
        ! average move_delta of 3 vector components
        getMove_delta = sum(this%move_delta)/size(this%move_delta)
        
    end function Spheres_getMove_delta

end module class_spheres
