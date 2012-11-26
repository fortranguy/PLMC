!***********************************************************************
!* MODULE: Constants                                                   *
!***********************************************************************
module data_constants

        integer, parameter :: DP = selected_real_kind(15, 307)
            ! double precision
        real(DP), parameter :: PI = acos(-1._DP)
        
end module data_constants
!***********************************************************************
!* MODULE : Cell                                                       *
!* PURPOSE : declaration of the cell parameters                        *   
!* COMMENT : size1,2=x, y; size3=z                                     *
!***********************************************************************
module data_cell

use data_constants
    
    integer, parameter :: Dim = 3
    real(DP), parameter :: Lsize1 = 12._DP
    real(DP), parameter :: Lsize2 = Lsize1
    real(DP), parameter :: Lsize3 = 12._DP
    real(DP), dimension(Dim), parameter :: Lsize = &
        [Lsize1, Lsize2, Lsize3]
    real(DP), dimension(Dim), parameter :: LsizeMi = 0.5_DP * Lsize
    real(DP), parameter::Tstar = 1._DP
    
end module data_cell
!***********************************************************************
!* MODULE: Neighbours                                                  *
!***********************************************************************
module data_neighbours

use data_cell

    ! Type
    type Particle
        integer :: iCol
        type(Particle), pointer :: next => null()
    end type Particle
    
    type ContainerParticle
        type(Particle), pointer :: particle => null()
    end type ContainerParticle
    
    ! Cellules
    integer, parameter :: cell_iMax = 3, cell_jMax = cell_iMax, cell_kMax = 3
    integer, dimension(dim), parameter :: cell_coordMax = &
        [cell_iMax, cell_jMax, cell_kMax]
    real, parameter :: cell_Lsize1 = Lsize1/real(cell_iMax), &
        cell_Lsize2 = Lsize2/real(cell_jMax), &
        cell_Lsize3 = Lsize3/real(cell_kMax)
    real, dimension(Dim), parameter :: cell_Lsize = &
        [cell_Lsize1, cell_Lsize2, cell_Lsize3]
    type(ContainerParticle), allocatable, dimension(:) :: cells, cellsNext
    type(ContainerParticle), allocatable, dimension(:), protected :: cellsBegin
                                            !Attention : ne pas déplacer.
    ! Voisins
    integer, parameter :: cell_neigh_iMax = 3, &
        cell_neigh_jMax = cell_neigh_iMax, &
        cell_neigh_kMax = cell_neigh_iMax ! évident
    integer, dimension(dim), parameter :: cell_neigh_coordMax = &
        [cell_neigh_iMax, cell_neigh_jMax, cell_neigh_kMax]
    integer, parameter :: cell_neighs_nb = &
        cell_neigh_iMax*cell_neigh_jMax*cell_neigh_kMax ! inclus soi-mêm
    integer, dimension(cell_iMax*cell_jMax*cell_kMax, cell_neighs_nb) :: &
        cell_neighs
    
contains
    
    ! Allocation
    
    subroutine alloc_Cells()
    
        integer :: iCell, nCells
        
        nCells = product(cell_coordMax)

        allocate(cellsBegin(nCells))
        allocate(cells(nCells))
        allocate(cellsNext(nCells))
        
        do iCell = 1, nCells

            allocate(cellsBegin(iCell)%particle)
            cells(iCell)%particle => cellsBegin(iCell)%particle
            cells(iCell)%particle%iCol = 0
            
            allocate(cellsNext(iCell)%particle)
            cellsNext(iCell)%particle%iCol = 0            
            cells(iCell)%particle%next => cellsNext(iCell)%particle
            cells(iCell)%particle => cellsNext(iCell)%particle     
    
        end do
        
    end subroutine alloc_Cells
    
    ! Libération
    
    recursive subroutine libere_chaine(courant)
    
        type(Particle), pointer :: courant
        
        if (associated(courant%next)) then
            call libere_chaine(courant%next)
        end if
        deallocate(courant)
        
    end subroutine libere_chaine
    
    subroutine dealloc_Cells()
    
        integer :: iCell
        integer :: nCells = cell_coordMax(1)*cell_coordMax(2)*cell_coordMax(3)
    
        do iCell = 1, nCells
            if (associated(cellsBegin(iCell)%particle)) then
                call libere_chaine(cellsBegin(iCell)%particle)
            end if
        end do
    
    end subroutine dealloc_Cells
    
end module data_neighbours
!***********************************************************************
!* MODULE : Particles                                                  *
!* PURPOSE : declaration of the particles parameters                   *
!* COMMENT : 1=big particles ; 2=small particles                       *
!***********************************************************************
module data_particles
    use data_cell
    use data_constants
    real(DP), parameter :: rayon1 = .5_DP
    real(DP), parameter :: rmin = 1._DP
    integer, parameter ::  Ncol1 = 270 ! Vs Ncolmax
    integer, parameter :: Ncolmax = 5000 
    real(DP), dimension(Dim, Ncolmax) :: X
end module data_particles
!***********************************************************************
    
!***********************************************************************
!* MODULE : MC                                                         *
!* PURPOSE : declaration of the MC parameters                          *
!* COMMENT : Field-induced layer formation in dipolar nanofilms        *
!***********************************************************************
module data_mc
    use data_constants
    use data_particles
    integer, parameter :: Nstep = 10**5
    integer, parameter :: Ntherm = 40
    integer, parameter :: Nmove = 4*Ncol1 ! new
    real(DP), dimension(Dim) :: dx = 3._DP ! new, à tester ?
end module data_mc
!***********************************************************************

!***********************************************************************
!* MODULE : Potentiel                                                  *
!* COMMENT : HS + Yukawa + cut                                         *
!***********************************************************************
module data_potentiel

use data_cell
use data_particles
use data_mc
use data_constants
    real(DP), parameter :: rcut11 = 4._DP ! new
    real(DP), parameter :: pas11 = 5.E-5_DP ! new
    real(DP), parameter :: surpas11 = 1._DP/pas11 ! new
    integer, parameter :: Ntab11 = int(rcut11*surpas11) ! new
    integer, parameter :: iMin = int( rmin/rcut11*real(Ntab11, DP) ) ! new
    real(DP), dimension(iMin:Ntab11), protected :: Vtab11
    real(DP), parameter :: epsilon11 = 1._DP, alpha11 = 5._DP ! new
    
contains
    
    subroutine ePotIni()

        integer :: i
        real(DP) :: r_i
       
	    ! cut
        do i = iMin, Ntab11       
            r_i = real(i, DP)*pas11
            Vtab11(i) = epsilon11*exp(-alpha11*(r_i-rmin))/r_i           
        end do
            ! problème pour V(rcut) ?
        
        ! shift        
        Vtab11(:) = Vtab11(:) - epsilon11*exp(-alpha11*(rcut11-rmin))/rcut11

    end subroutine ePotIni
        
end module data_potentiel
!***********************************************************************
