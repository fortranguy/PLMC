module mod_tools

use data_constants
use data_particles
use data_mc
use data_potentiel
use data_neighbours

implicit none

    contains

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
        
        do iCol = 1, Ncol1
            write(unitSnap, *) X(:, iCol)
        end do    

    end subroutine
    
    ! Rapport -----------------------------------------------------------------
    
    subroutine rapport(nWidom, unitRapport)
    
        integer, intent(in) :: nWidom
        integer, intent(in) :: unitRapport    
        
        write(unitRapport, *) "Simulation MC_C :"
        write(unitRapport ,*) "    Lsize(:) = ", Lsize(:)
        write(unitRapport ,*) "    Vol = ", product(Lsize)
        write(unitRapport ,*) "    Ncol1 = ", Ncol1
        write(unitRapport ,*) "    nWidom = ", nWidom
        write(unitRapport, *) "    Nstep = ", Nstep
        write(unitRapport, *) "    Ntherm = ", Ntherm
        write(unitRapport, *) "    Nmove = ", Nmove
        write(unitRapport, *) "    epsilon = ", epsilon
        write(unitRapport, *) "    alpha = ", alpha
        write(unitRapport, *) "    rcut = ", rcut
        write(unitRapport, *) "    pas = ", pas
        write(unitRapport, *) "    cell_coordMax(:) = ", cell_coordMax(:)
        write(unitRapport, *) "    cell_Lsize(:) = ", cell_Lsize(:)
        
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
            enTotSum/realNstep/real(Ncol1, DP)
        potChiId = -Tstar*log( product(Lsize)/real(Ncol1+1,DP) )
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
    
