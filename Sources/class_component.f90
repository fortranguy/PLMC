!***********************************************************************
!* MODULE: Component class                                              *
!***********************************************************************

module class_component

use data_cell
use data_neighbours
use mod_physique

implicit none

private
public :: sph, sph_init

    type Link
        integer :: iCol
        type(Link), pointer :: next => null()
    end type Link
    
    type LinkedList
        type(Link), pointer :: particle => null()
    end type LinkedList

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
        
        ! Neighbours
        
        real(DP), dimension(Dim) :: cell_Lsize
        integer, dimension(Dim) :: cell_coordMax
        integer, dimension(:, :), allocatable :: cell_neighs
        type(LinkedList), allocatable, dimension(:) :: cells, cellsNext
        type(LinkedList), allocatable, dimension(:) :: cellsBegin
        
    contains

        procedure :: overlapTest => component_overlapTest
        
        procedure :: alloc_Cells => component_alloc_Cells
        procedure :: libere_chaine => component_libere_chaine
        procedure :: dealloc_Cells => component_dealloc_Cells
        procedure :: check_CellsSize => component_check_CellsSize
        procedure :: position_to_cell => component_position_to_cell
        procedure :: all_col_to_cell => component_all_col_to_cell
        procedure :: remove_cell_col => component_remove_cell_col
        procedure :: add_cell_col => component_add_cell_col
        procedure :: cell_coord_to_ind => component_cell_coord_to_ind
        procedure :: cell_neigh_coord_to_ind => component_cell_neigh_coord_to_ind
        procedure :: cell_period => component_cell_period
        procedure :: ini_cell_neighs => component_ini_cell_neighs
        
        procedure :: ePot => component_ePot
        procedure :: ePotNeigh => component_ePotNeigh
        procedure :: enTotCalc => component_enTotCalc
        
        procedure :: mcMove => component_mcMove
        procedure :: widom => component_widom
        
    end type Component
    
    type(Component), protected :: sph
    
contains

    subroutine sph_init()
    
        integer, dimension(:, :), allocatable :: sph_cell_neighs
        
        allocate(sph_cell_neighs(cell_neighs_nb, int(Lsize(1)/rcut) * &
            int(Lsize(2)/rcut) * int(Lsize(3)/rcut)))
    
        ! Component initialization
        
        call ePotIni()
        
        ! Construction
                
        sph =   Component(&        
                    radius = radius, &
                    rmin = rmin, &
                    Ncol = Ncol, &
                    dx = dx, &
                    rcut = rcut, &
                    pas = pas, &
                    iMin = iMin, &
                    Ntab = Ntab, &
                    epsilon = epsilon, &
                    alpha = alpha, &
                    Vtab = Vtab, &
                    cell_Lsize = [rcut, rcut, rcut], &
                    cell_coordMax = [int(Lsize(1)/rcut), &
                        int(Lsize(2)/rcut), int(Lsize(3)/rcut)], &
                    cell_neighs = sph_cell_neighs &
                )
        
    end subroutine sph_init
    
    ! Test d'overlap ----------------------------------------------------------
    
    subroutine component_overlapTest(this)
    
        class(Component), intent(in) :: this
    
        integer :: jCol, iCol
        real(DP) :: r_ij
    
        do jCol = 1, Ncol
            do iCol = 1, Ncol
                if (iCol /= jCol) then
                    
                    r_ij = dist(X(:, iCol), X(:, jCol))
                    if (r_ij < this%rmin) then
                        write(*, *) "    Overlap !", iCol, jCol
                        write(*, * ) "    r_ij = ", r_ij
                        stop
                    end if
                    
                end if
            end do
        end do
        
        write(*, *) "    Overlap test : OK !"
    
    end subroutine component_overlapTest
    
    ! Linked-list allocation
    
    subroutine component_alloc_Cells(this)
    
        class(Component), intent(inout) :: this
    
        integer :: iCell, nCells
        
        nCells = product(this%cell_coordMax)

        allocate(this%cellsBegin(nCells))
        allocate(this%cells(nCells))
        allocate(this%cellsNext(nCells))
        
        do iCell = 1, nCells

            allocate(this%cellsBegin(iCell)%particle)
            this%cells(iCell)%particle => this%cellsBegin(iCell)%particle
            this%cells(iCell)%particle%iCol = 0
            
            allocate(this%cellsNext(iCell)%particle)
            this%cellsNext(iCell)%particle%iCol = 0            
            this%cells(iCell)%particle%next => this%cellsNext(iCell)%particle
            this%cells(iCell)%particle => this%cellsNext(iCell)%particle     
    
        end do
        
    end subroutine component_alloc_Cells
    
    ! Libération
    
    recursive subroutine component_libere_chaine(this, courant)
    
        class(Component), intent(inout) :: this
        type(Link), pointer :: courant
        
        if (associated(courant%next)) then
            call this%libere_chaine(courant%next)
        end if
        deallocate(courant)
        
    end subroutine component_libere_chaine
    
    subroutine component_dealloc_Cells(this)
    
        class(Component), intent(inout) :: this
    
        integer :: iCell
        integer :: nCells
    
        nCells = this%cell_coordMax(1) * this%cell_coordMax(2) * &
            this%cell_coordMax(3)
        do iCell = 1, nCells
            if (associated(this%cellsBegin(iCell)%particle)) then
                call libere_chaine(this%cellsBegin(iCell)%particle)
            end if
        end do
    
    end subroutine component_dealloc_Cells
    
    ! Vérification de la taille des cellules (voisines)
    
    subroutine component_check_CellsSize(this)
    
        class(Component), intent(in) :: this
        
        integer :: iDir
        
        do iDir = 1, Dim
        
            if (this%cell_Lsize(iDir) < this%rcut) then
                write(*, *) "Cellule trop petite dans la direction", iDir, ":"
                write(*, *) this%cell_Lsize(iDir), "<", this%rcut
                stop
            end if
            
            if (this%cell_coordMax(iDir) < cell_neigh_coordMax(iDir)) then
                write(*, *) "Trop peu de cellules dans la direction", iDir, ":"
                write(*, *) this%cell_coordMax(iDir), "<",&
                    cell_neigh_coordMax(iDir)
                stop
            end if
            
        end do
        
    end subroutine component_check_CellsSize
    
    ! Assignation : particule -> cellule
    
    function component_position_to_cell(this, xCol) &
        result(position_to_cell)
    
        class(Component), intent(in) :: this
        real(DP), dimension(Dim), intent(in) :: xCol
        
        integer, dimension(Dim) :: cell_coord
        integer :: position_to_cell
    
        cell_coord(:) = int( xCol(:)/this%cell_Lsize(:) ) + 1
        position_to_cell = cell_coord(1) + this%cell_coordMax(1) * &
            (cell_coord(2)-1) + this%cell_coordMax(1) * &
            this%cell_coordMax(2) * (cell_coord(3)-1)
    
    end function component_position_to_cell
    
    subroutine component_all_col_to_cell(this)
    
        class(Component), intent(inout) :: this
    
        integer :: iCol
        integer :: iCell, nCells
        
        nCells = this%cell_coordMax(1) * this%cell_coordMax(2) * &
            this%cell_coordMax(3)
    
        do iCol = 1, this%Ncol
    
            iCell = this%position_to_cell(X(:,iCol))
            this%cells(iCell)%particle%iCol = iCol
            
            allocate(this%cellsNext(iCell)%particle)
            this%cellsNext(iCell)%particle%iCol = 0
            this%cells(iCell)%particle%next => &
                this%cellsNext(iCell)%particle
            this%cells(iCell)%particle => this%cellsNext(iCell)%particle
            
        end do
        
        do iCell = 1, nCells
            
            this%cells(iCell)%particle%next => null()
            
        end do
        
    end subroutine component_all_col_to_cell
    
    ! Mise à jour des TV
    
    subroutine component_remove_cell_col(this, iCol, iCellBefore)
    
        class(Component), intent(inout) :: this
    
        integer, intent(in) :: iCol, iCellBefore
        
        type(Link), pointer :: courant => null()
        type(Link), pointer :: suivant => null(), precedent => null()
    
        precedent => this%cellsBegin(iCellBefore)%particle
        courant => precedent%next
        
        do
        
            suivant => courant%next
        
            if ( courant%iCol == iCol ) then
            
                precedent%next => courant%next
                deallocate(courant)
                courant => suivant
                exit
                
            else
            
                precedent => courant
                
            end if
            
            courant => suivant
        
        end do
            
    end subroutine component_remove_cell_col   
    
    subroutine component_add_cell_col(this, iCol, iCellAfter)
    
        class(Component), intent(inout) :: this
    
        integer, intent(in) :: iCol, iCellAfter
    
        type(Link), pointer :: nouveau => null()
        type(Link), pointer :: suivant => null(), precedent => null()           
          
        
        precedent => this%cellsBegin(iCellAfter)%particle
        
        do
        
            suivant => precedent%next
            
            if (.not. associated(suivant%next)) then
            
                allocate(nouveau)
                nouveau%next => precedent%next
                precedent%next => nouveau
                nouveau%iCol = iCol
                exit
                
            end if
            
            precedent => suivant
            
        end do
            
    end  subroutine component_add_cell_col
    
! -----------------------------------------------------------------------------
! VOISINS :
! -----------------------------------------------------------------------------
    
    function component_cell_coord_to_ind(this, coord) result(cell_coord_to_ind)
        
        class(Component), intent(in) :: this    
        integer, dimension(Dim), intent(in) :: coord
        
        integer :: cell_coord_to_ind
        
        cell_coord_to_ind = coord(1) + this%cell_coordMax(1)*(coord(2)-1) + &
            this%cell_coordMax(1)*this%cell_coordMax(2)*(coord(3)-1)
    
    end function component_cell_coord_to_ind
    
    function component_cell_neigh_coord_to_ind(this, neigh_coord) &
        result(cell_neigh_coord_to_ind)
    
        class(Component), intent(in) :: this
        integer, dimension(Dim), intent(in) :: neigh_coord
        
        integer :: cell_neigh_coord_to_ind
        
        cell_neigh_coord_to_ind = neigh_coord(1) + &
            cell_neigh_coordMax(1) * (neigh_coord(2)-1) &
            + cell_neigh_coordMax(1) * cell_neigh_coordMax(2) * &
            (neigh_coord(3)-1)
    
    end function component_cell_neigh_coord_to_ind
    
    function component_cell_period(this, coord) result(cell_period)
    
        class(Component), intent(in) :: this    
        integer, dimension(Dim), intent(in) :: coord
        
        integer, dimension(Dim) :: cell_period
        
        cell_period(:) = coord(:)
        
        where (cell_period(:) < 1)
            cell_period(:) = cell_period(:) + this%cell_coordMax(:)
        end where
        
        where (cell_period(:) > this%cell_coordMax(:))
            cell_period(:) = cell_period(:) - this%cell_coordMax(:)
        end where
    
    end function component_cell_period
    
    subroutine component_ini_cell_neighs(this)
    
        class(Component), intent(inout) :: this 
    
        integer :: i, j, k, ind
        integer :: neigh_i, neigh_j, neigh_k, neigh_ind
        integer, dimension(Dim) :: coord, neigh_coord
        
        do i = 1, this%cell_coordMax(1)
        do j = 1, this%cell_coordMax(2)
        do k = 1, this%cell_coordMax(3)
            
            ind = this%cell_coord_to_ind([i, j, k])

            do neigh_i = 1, cell_neigh_coordMax(1)
            do neigh_j = 1, cell_neigh_coordMax(2)
            do neigh_k = 1, cell_neigh_coordMax(3)
            
                neigh_coord(:) = [neigh_i, neigh_j, neigh_k]                
                neigh_ind = this%cell_neigh_coord_to_ind(neigh_coord(:))          
                neigh_coord(:) = neigh_coord(:) - cell_neigh_coordMax(:) + 1
                    ! Par rapport au centre [i, j, k]
                
                coord(:) = [i, j, k] + neigh_coord(:)
                
                this%cell_neighs(neigh_ind, ind) = &
                    this%cell_coord_to_ind( this%cell_period(coord(:)) )
                    
            end do
            end do
            end do
        
        end do
        end do
        end do
            
    end subroutine component_ini_cell_neighs

    ! Energie potentielle -------------------------------------------------

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
        
            iCell_neigh = this%cell_neighs(iNeigh, iCell)
            courant => this%cellsBegin(iCell_neigh)%particle%next            
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
    
    ! Déplacement d'une particule -----------------------------------------
    
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
        iCellAfter = this%position_to_cell(xNew)
        call this%ePotNeigh(iOld, xNew, iCellAfter, overlap, eNew)
        
        if (.not. overlap) then
        
            iCellBefore = this%position_to_cell(X(:, iOld))
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
    
    ! Méthode de Widom ----------------------------------------------------

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
            iCellTest = this%position_to_cell(xTest)
            call this%ePotNeigh(0, xTest, iCellTest, overlap, enTest) 
            
            if (.not. overlap) then
                widTestSum = widTestSum + exp(-enTest/Tstar)
            end if
            
        end do
        
        activExInv = widTestSum/real(nWidom, DP)
        
    end subroutine component_widom

	! Energie potentielle totale
    
    function component_enTotCalc(this) result(enTotCalc)
    
    	class(Component), intent(in) :: this
    	
        integer :: iCol, jCol
        real(DP) :: r_ij
        real(DP) :: enTotCalc
    
        enTotCalc = 0._DP
        
        do jCol = 1, this%Ncol
            do iCol = 1, this%Ncol
                if (iCol /= jCol) then
                
                    r_ij = dist(X(:, iCol), X(:, jCol))
                    enTotCalc = enTotCalc + this%ePot(r_ij)
                    
                end if
            end do
        end do
        
        enTotCalc = 0.5_DP*enTotCalc
    
    end function component_enTotCalc

end module class_component
