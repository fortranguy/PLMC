!***********************************************************************
!* MODULE: Component class                                              *
!***********************************************************************

module class_component

use data_cell
use data_neighbours
use mod_neighbours
use mod_physique

implicit none

private

    type, public :: Component

        ! Particles

        real(DP) :: radius
        real(DP) :: rmin
        integer ::  Ncol

        ! Monte-Carlo
        
        real(DP), dimension(Dim) :: dx

        ! Potential domain

        real(DP) :: rcut
        real(DP) :: pas
        integer :: iMin
        integer :: Ntab
        
        ! Potential shape
        
        real(DP) :: epsilon
        real(DP) :: alpha
        real(DP), dimension(:), allocatable :: Vtab
        
    contains
    
        procedure :: ePot => component_ePot
        procedure :: ePotNeigh => component_ePotNeigh
        procedure :: mcMove => component_mcMove
        procedure :: widom => component_widom
        
    end type Component
    
contains

    ! Energie potentielle -----------------------------------------------------

    function component_ePot(this, r) result(ePot)
        
        class(Component), intent(in) :: this
        real(DP), intent(in) :: r
        
        integer :: i
        real(DP) :: r_i, ePot
       
        if (r < this%rcut) then
       
            i = int(r/this%pas)
            r_i = real(i, DP)*this%pas
            ePot = this%Vtab(i) + (r-r_i)/this%pas * &
                (this%Vtab(i+1)-this%Vtab(i))
           
        else
       
            ePot = 0.
           
        end if
        
    end function component_ePot
    
    ! -----------------------
    
    subroutine component_ePotNeigh(this, iCol, xCol, iCell, overlap, energ)
        
        class(Component), intent(in) :: this        
        integer, intent(in) :: iCol, iCell
        real(DP), dimension(Dim), intent(in) :: xCol
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energ
    
        integer :: iNeigh,  iCell_neigh
        real(DP) :: r
    
        type(Link), pointer :: courant => null(), suivant => null()
        
        overlap = .false.
        energ = 0._DP
    
        do iNeigh = 1, cell_neighs_nb
        
            iCell_neigh = cell_neighs(iNeigh, iCell)
            courant => cellsBegin(iCell_neigh)%particle%next            
            if (.not. associated(courant%next)) cycle
            
            do
            
                suivant => courant%next
            
                if (courant%iCol /= iCol) then
                
                    r = dist(xCol(:), X(:, courant%iCol))
                    if (r < this%rmin) then
                        overlap = .true.
                        return
                    end if
                    energ = energ + this%ePot(r)
       
                end if
                
                if (.not. associated(suivant%next)) exit
                
                courant => suivant
            
            end do            
            
        end do
    
    end subroutine component_ePotNeigh
    
    ! Déplacement d'une particule ---------------------------------------------
    
    subroutine component_mcMove(this, enTot, Nrejects)
    
    	class(Component), intent(in) :: this
        real(DP), intent(inout) :: enTot
        integer, intent(inout) :: Nrejects
        
        logical :: overlap
        integer :: iOld
        real(DP) :: rand
        real(DP), dimension(Dim) :: xNew
        integer :: iCellBefore, iCellAfter
        real(DP) :: eNew, eOld, dEn
        
        call random_number(rand)
        iOld = int(rand*Ncol) + 1
        
        call random_number(xNew)
        xNew(:) = X(:, iOld) + (xNew(:)-0.5_DP)*this%dx(:)
        xNew(:) = modulo(xNew(:), Lsize(:))
        iCellAfter = position_to_cell(xNew)
        call this%ePotNeigh(iOld, xNew, iCellAfter, overlap, eNew)
        
        if (.not. overlap) then
        
            iCellBefore = position_to_cell(X(:, iOld))
            call this%ePotNeigh(iOld, X(:, iOld), iCellBefore, overlap, eOld)
        	
            dEn = eNew - eOld
        
            call random_number(rand)
            if ( rand < exp(-dEn/Tstar) ) then
                X(:, iOld) = xNew(:)
                enTot = enTot + dEn
                
                if ( iCellBefore /= iCellAfter ) then                
                    call remove_cell_col(iOld, iCellBefore)
                    call add_cell_col(iOld, iCellAfter)
                end if
                
            else
                Nrejects = Nrejects + 1
            end if
            
        else
        
            Nrejects = Nrejects + 1
            
        end if
    
    end subroutine component_mcMove
    
    ! Méthode de Widom --------------------------------------------------------

    subroutine component_widom(this, nWidom, activExInv)
        
        class(Component), intent(in) :: this
        integer, intent(in) :: nWidom
        real(DP), intent(inOut) :: activExInv 
        
        integer :: iWid
        real(DP) :: widTestSum
        real(DP), dimension(Dim) :: xTest
        integer :: iCellTest
        logical :: overlap        
        real(DP) :: enTest
        
        widTestSum = 0._DP
        
        do iWid = 1, nWidom           
            
            call random_number(xTest)
            xTest(:) = Lsize(:) * xTest(:)    
            iCellTest = position_to_cell(xTest)
            call this%ePotNeigh(0, xTest, iCellTest, overlap, enTest) 
            
            if (.not. overlap) then
                widTestSum = widTestSum + exp(-enTest/Tstar)
            end if
            
        end do
        
        activExInv = widTestSum/real(nWidom, DP)
        
    end subroutine component_widom

end module class_component
