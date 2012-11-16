module mod_tools

use data_particles
use data_mc
use data_potentiel
use data_constants
    
    contains

    ! Générateurs de nombres aléatoires : graine ------------------------------
    
    subroutine init_random_seed()
    
        integer :: i, n, clock
        integer, dimension(:), allocatable :: seed

        call random_seed(size = n)
        allocate(seed(n))

        call system_clock(count=clock)

        seed = clock + 37 * [ (i - 1, i = 1, n) ]
        call random_seed(put = seed)

        deallocate(seed)
        
    end subroutine init_random_seed
    
    ! Etat de la configuration ------------------------------------------------
      
    subroutine snapShot(unitSnap)
        
        integer, intent(in) :: unitSnap
    
        integer :: iPart
        
        do iPart = 1, Ncol1
            write(unitSnap, *) X(:, iPart)
        end do    

    end subroutine
    
    ! Résultats ---------------------------------------------------------------
        
    subroutine mcResults(enTotSum, activExInvSum, tauxRejects, unitRapport)

        real(DP), intent(in) :: enTotSum, activExInvSum     
        real(DP), intent(in) :: tauxRejects
        integer, intent(in) :: unitRapport
        
        real(DP) :: realNstep = real(Nstep, DP)
    
        write(unitRapport, *) "Résultats :"
        write(unitRapport, *) "    Energie moyenne = ", &
            enTotSum/realNstep
        write(unitRapport, *) "    Energie moyenne par particule = ", &
            enTotSum/realNstep/real(Ncol1, DP)
		write(unitRapport, *) "    Potentiel chimique idéal = ", &
            -Tstar*log( product(Lsize)/real(Ncol1+1,DP) )
        write(unitRapport, *) "    Potentiel chimique (excès) moyen = ", &
            -Tstar*log( activExInvSum/realNstep )
        write(unitRapport, *) "    Taux rejets = ", &
            tauxRejects/real(Nstep+Ntherm, DP)
    
    end subroutine mcResults
    
    ! Rapport -----------------------------------------------------------------
    
    subroutine rapport(nWidom, Lratio, unitRapport)
    
        integer, intent(in) :: nWidom
        real(DP), intent(in) :: Lratio
        integer, intent(in) :: unitRapport    
        
        write(unitRapport, *) "Simulation MC_C :"
        write(unitRapport ,*) "    Lsize(:) = ", Lsize(:)
        write(unitRapport ,*) "    Vol = ", product(Lsize)
        write(unitRapport ,*) "    Ncol1 = ", Ncol1
        write(unitRapport ,*) "    nWidom = ", nWidom
        write(unitRapport ,*) "    Lratio = ", Lratio
        write(unitRapport, *) "    Nstep = ", Nstep
        write(unitRapport, *) "    Ntherm = ", Ntherm
        write(unitRapport, *) "    Nmove = ", Nmove
        write(unitRapport, *) "    dx(:) = ", dx(:)
        write(unitRapport, *) "    epsilon11 = ", epsilon11
        write(unitRapport, *) "    alpha11 = ", alpha11
        write(unitRapport, *) "    rcut11 = ", rcut11
        write(unitRapport, *) "    surpas11 = ", surpas11
        
    end subroutine rapport
    
end module mod_tools
    
