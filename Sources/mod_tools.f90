module mod_tools

use data_constants
use data_particles
use data_mc
use data_potentiel
use data_neighbours
use class_component

implicit none

    contains
    
    function sph_constructor
    
        type(Component) :: sph_constructor
    
        ! Component initialization
        
        call ePotIni()
        
        ! Construction                

        sph_constructor%radius = radius
        sph_constructor%rmin = rmin
        sph_constructor%Ncol = Ncol
        sph_constructor%dx = dx
        sph_constructor%rcut = rcut
        sph_constructor%pas = pas
        sph_constructor%iMin = iMin
        sph_constructor%Ntab = Ntab
        sph_constructor%epsilon = epsilon
        sph_constructor%alpha = alpha        
        allocate(sph_constructor%Vtab(iMin:Ntab))
        sph_constructor%Vtab(:) = Vtab(:)
        sph_constructor%cell_Lsize(:) = [rcut, rcut, rcut]
        sph_constructor%cell_coordMax(:) = int(Lsize(:)/rcut)
        allocate(sph_constructor%cell_neighs(cell_neighs_nb, &
            product( int(Lsize(:)/rcut) )))
    
    end function sph_constructor

    ! Générateurs de nombres aléatoires : graine ------------------------------
    
    subroutine init_random_seed(unitRapport)
    
        integer, intent(in) :: unitRapport
    
        integer :: i, n, clock
        integer, dimension(:), allocatable :: seed

        call random_seed(size = n)
        allocate(seed(n))

        call system_clock(count=clock)

        seed(:) = clock + 37 * [ (i - 1, i = 1, n) ]
        call random_seed(put = seed)
        
        write(unitRapport, *) "Graine :"
        write(unitRapport ,*) "    n = ", n
        write(unitRapport ,*) "    seed(:) = ", seed(:)

        deallocate(seed)
        
    end subroutine init_random_seed
    
    ! Etat de la configuration ------------------------------------------------
      
    subroutine snapShot(unitSnap)
        
        integer, intent(in) :: unitSnap
    
        integer :: iCol
        
        do iCol = 1, Ncol
            write(unitSnap, *) X(:, iCol)
        end do    

    end subroutine
    
    ! Rapport -----------------------------------------------------------------
    
    subroutine rapport(comp, nWidom, unitRapport)
    
        type(Component), intent(in) :: comp
        integer, intent(in) :: nWidom
        integer, intent(in) :: unitRapport    
        
        write(unitRapport, *) "Simulation MC_C :"
        write(unitRapport ,*) "    Lsize(:) = ", Lsize(:)
        write(unitRapport ,*) "    Vol = ", product(Lsize)
        write(unitRapport ,*) "    Ncol = ", Ncol
        write(unitRapport ,*) "    nWidom = ", nWidom
        write(unitRapport, *) "    Nstep = ", Nstep
        write(unitRapport, *) "    Ntherm = ", Ntherm
        write(unitRapport, *) "    Nmove = ", Nmove
        write(unitRapport, *) "    epsilon = ", epsilon
        write(unitRapport, *) "    alpha = ", alpha
        write(unitRapport, *) "    rcut = ", rcut
        write(unitRapport, *) "    pas = ", pas
        write(unitRapport, *) "    cell_coordMax(:) = ", comp%cell_coordMax(:)
        write(unitRapport, *) "    cell_Lsize(:) = ", comp%cell_Lsize(:)
        
    end subroutine rapport
    
    ! Résultats ---------------------------------------------------------------
        
    subroutine mcResults(enTotSum, activExInvSum, tauxRejectsSum, duree,&
    	unitRapport)

        real(DP), intent(in) :: enTotSum, activExInvSum     
        real(DP), intent(in) :: tauxRejectsSum
        real(DP), intent(in) :: duree
        integer, intent(in) :: unitRapport
        
        real(DP) :: realNstep = real(Nstep, DP)
        real(DP) :: potChiId, potChiEx
    
        write(unitRapport, *) "Résultats :"
        write(unitRapport, *) "    Energie moyenne = ", &
            enTotSum/realNstep
        write(unitRapport, *) "    Energie moyenne par particule = ", &
            enTotSum/realNstep/real(Ncol, DP)
        potChiId = -Tstar*log( product(Lsize)/real(Ncol+1,DP) )
        write(unitRapport, *) "    Potentiel chimique idéal = ", potChiId
        potChiEx = -Tstar*log( activExInvSum/realNstep )
        write(unitRapport, *) "    Potentiel chimique (excès) moyen = ", &
            potChiEx           
        write(unitRapport, *) "    potChi.moy = ", potChiId + potChiEx
        write(unitRapport, *) "    Taux rejets = ", &
            tauxRejectsSum/real(Nstep+Ntherm, DP)
        write(unitRapport, *) "    Durée =", duree/60._DP, "min"        
    
    end subroutine mcResults
    
end module mod_tools
    
