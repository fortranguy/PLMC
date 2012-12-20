module mod_dist

use data_constants
use data_distrib

implicit none

contains

    function sphereVol(iDist)
    
        integer, intent(in) :: iDist    
        real(DP) :: sphereVol
        
        sphereVol = 4._DP/3._DP * PI * ( real(iDist, DP)*deltaDist )**3
        
    end function sphereVol
    
end module mod_dist

program distribution

use data_constants
use data_distrib
use data_particles
use mod_dist
use mod_physique

implicit none

	real(DP), parameter :: densite = real(Ncol1, DP) / (Lsize1*Lsize2*Lsize3)
	integer, dimension(:), allocatable :: distrib
	integer, parameter :: unitSnapEnCours = 10, unitDistrib = 11, &
		unitEnerg = 12

	integer :: iStep
	integer :: iCol, jCol
	real(DP), dimension(Dim) :: DeltaX
	real(DP) :: r_ij
	integer :: iDist, iDistMin, iDistMax
	real(DP) :: r
	real(DP) :: numerat, denomin
	real(DP), dimension(:), allocatable :: fct_dist
	real(DP) :: energSum
	
	if (.not.snap) stop "Snap désactivé."
	
	call initDistriParams()
	allocate(distrib(Ndist))
	allocate(fct_dist(Ndist))
	
	distrib(:) = 0
	
	
	open(unit=unitSnapEnCours, recl=4096, file="snap.shot", status='old', &
		action='read')
	
	do iStep = 1, Nstep
	
		! Lecture :
			
		do iCol = 1, Ncol1
		    read(unitSnapEnCours, *) X(:, iCol)
	    end do
    
		! Traitement
	
		do iCol = 1, Ncol1
			do jCol = iCol + 1, Ncol1
		
				DeltaX(:) = X(:, jCol) - X(:, iCol)
		        call pbc_dif(DeltaX)
		        r_ij = sqrt(dot_product(DeltaX, DeltaX))
		        
		        iDist =  int( r_ij/deltaDist )
		        distrib(iDist) = distrib(iDist) + 1
		
			end do		
		end do
	
	end do
	
	close(unitSnapEnCours)
	
	! Ecriture
	
	iDistMin = 0
	iDistMax = 0
	
	open(unit=unitDistrib, file="fct_distrib.out", action="write")	
		do iDist = 1, Ndist
		
			r = (real(iDist, DP) + 0.5_DP) * deltaDist
			numerat = real(distrib(iDist), DP) / real(Nstep, DP)
			denomin = real(Ncol1, DP) * (sphereVol(iDist+1) - sphereVol(iDist))
			fct_dist(iDist) = 2._DP * numerat / denomin / densite
			write(unitDistrib, *) r, fct_dist(iDist)
			
			if (r>=rmin .and. r<=rcut11) then
				if (iDistMin == 0) then
					iDistMin = iDist
				end if
				iDistMax = iDist
			end if
			
		end do
	close(unitDistrib)
	
	! Energie par particule

	call ePotIni()
	
	energSum = 0._DP
	
	do iDist = iDistMin, iDistMax
		r = (real(iDist, DP) + 0.5_DP) * deltaDist
		energSum = energSum + ePot(r) * fct_dist(iDist) * 4._DP*PI*r**2	
	end do

	open(unit=unitEnerg, file="epp_dist.out", action="write")
		write(unitEnerg, *) "epp_dist =", &
			densite/2._DP * energSum * deltaDist
	close(unitEnerg)
	
	deallocate(fct_dist)
	deallocate(distrib)

end program distribution
