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
        procedure :: get_name => Spheres_get_name
        procedure :: get_Ncol => Spheres_get_Ncol
        procedure :: get_rMin => Spheres_get_rMin
        procedure :: get_rCut => Spheres_get_rCut        
        procedure :: get_move_Nadapt => Spheres_get_move_Nadapt
        !> Specifier
        procedure :: construct_mixCells => Spheres_construct_mixCells
        
        procedure :: print_density => Spheres_print_density
        
        !> Take a snap shot of the configuration : positions
        procedure :: snapShot_positions => Spheres_snapShot_positions
        
        !> Do an overlap test
        procedure :: test_overlap => Spheres_test_overlap
        
        !> Assign all particles to cells
        procedure :: all_cols_to_cells => Spheres_all_cols_to_cells
                
        !> Adapt the displacement move_delta during thermalisation
        procedure :: adapt_move_delta => Spheres_adapt_move_delta
        procedure :: set_move_delta => Spheres_set_move_delta
        procedure :: get_move_delta => Spheres_get_move_delta
        
    end type Spheres
    
contains

    !> Accessor : name

    pure function Spheres_get_name(this) result(get_name)
    
        class(Spheres), intent(in) :: this        
        character(len=5) :: get_name
        
        get_name = this%name
    
    end function Spheres_get_name

    !> Accessor : Ncol

    pure function Spheres_get_Ncol(this) result(get_Ncol)
    
        class(Spheres), intent(in) :: this        
        integer :: get_Ncol
        
        get_Ncol = this%Ncol
    
    end function Spheres_get_Ncol
    
    !> Accessor : rMin
    
    pure function Spheres_get_rMin(this) result(get_rMin)
    
        class(Spheres), intent(in) :: this        
        real(DP) :: get_rMin
        
        get_rMin = this%rMin
    
    end function Spheres_get_rMin
    
    !> Accessor : rCut
    
    pure function Spheres_get_rCut(this) result(get_rCut)
    
        class(Spheres), intent(in) :: this        
        real(DP) :: get_rCut
        
        get_rCut = this%rCut
    
    end function Spheres_get_rCut
    
    !> Accessor : move_Nadapt
    
    pure function Spheres_get_move_Nadapt(this) result(get_move_Nadapt)
    
        class(Spheres), intent(in) :: this
        integer :: get_move_Nadapt
        
        get_move_Nadapt = this%move_Nadapt
        
    end function Spheres_get_move_Nadapt
    
    !> Specifier : mixCells construction
    
    subroutine Spheres_construct_mixCells(this, mix_cell_size, mix_rCut)
    
        class(Spheres), intent(inout) :: this
        real(DP), dimension(:), intent(in) :: mix_cell_size
        real(DP), intent(in) :: mix_rCut
        
        write(output_unit, *) this%name, ": mixCells construction"
        
        call this%mixCells%construct(mix_cell_size, mix_rCut)
    
    end subroutine Spheres_construct_mixCells
    
    !> Print density and compacity
    
    subroutine Spheres_print_density(this, report_unit)
    
        class(Spheres), intent(in) :: this
        integer, intent(in) :: report_unit
        
        real(DP) :: density, compacity
        
        density = real(this%Ncol + 1, DP) / Volume ! cheating ? cf. Widom
        compacity = 4._DP/3._DP*PI*this%radius**3 * density
        
        write(output_unit, *) this%name, " : ", "density = ", density, "compacity = ", compacity
        
        write(report_unit, *) "    density = ", density
        write(report_unit, *) "    compacity = ", compacity
    
    end subroutine Spheres_print_density
    
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
    
    subroutine Spheres_test_overlap(this)
    
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
    
    end subroutine Spheres_test_overlap
    
    !> Fill cells with colloids
    
    subroutine Spheres_all_cols_to_cells(this, other)
    
        class(Spheres), intent(inout) :: this
        class(Spheres), intent(in) :: other
        
        call this%sameCells%all_cols_to_cells(this%Ncol, this%positions)
        call this%mixCells%all_cols_to_cells(other%Ncol, other%positions)
    
    end subroutine Spheres_all_cols_to_cells
    
    !> Adaptation of move_delta during the thermalisation
    
    subroutine Spheres_adapt_move_delta(this, reject)
    
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
    
    end subroutine Spheres_adapt_move_delta
    
    subroutine Spheres_set_move_delta(this, reject, report_unit)
    
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
    
    end subroutine Spheres_set_move_delta
    
    pure function Spheres_get_move_delta(this) result(get_move_delta)
        
        class(Spheres), intent(in) :: this        
        real(DP) :: get_move_delta
        
        ! average move_delta of 3 vector components
        get_move_delta = sum(this%move_delta)/size(this%move_delta)
        
    end function Spheres_get_move_delta

end module class_spheres
