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
    
    ! Rapport -----------------------------------------------------------------
    
    subroutine rapport(nWidom, Lratio, unitRapport)
    
        integer, intent(in) :: nWidom
        real(DP), intent(in) :: Lratio
        integer, intent(in) :: unitRapport    
        
        write(unitRapport, *) "Simulation MC_C :"
        write(unitRapport ,*) "    Lsize(:) = ", Lsize(:)
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
        
        write(*, *) "Simulation MC_C :"
        write(* ,*) "    Lsize(:) = ", Lsize(:)
        write(* ,*) "    Ncol1 = ", Ncol1
        write(* ,*) "    nWidom = ", nWidom
        write(* ,*) "    Lratio = ", Lratio
        write(*, *) "    Nstep = ", Nstep
        write(*, *) "    Ntherm = ", Ntherm
        write(*, *) "    Nmove = ", Nmove
        write(*, *) "    dx(:) = ", dx(:)
        write(*, *) "    epsilon11 = ", epsilon11
        write(*, *) "    alpha11 = ", alpha11
        write(*, *) "    rcut11 = ", rcut11
        write(*, *) "    surpas11 = ", surpas11
        
    end subroutine rapport
    
end module mod_tools
    
