program distribution

use data_constants
use data_distrib
use data_particles
use mod_physique

implicit none

	real(DP), parameter :: densite = real(Ncol1, DP) / (Lsize1*Lsize2*Lsize3)
	integer, dimension(:, :), allocatable :: distrib
	integer, parameter :: unitSnapEnCours = 10, unitdistrib = 11

	integer :: iStep
	character(len=20) :: fileSnap, iSnap
	integer :: iCol, jCol
	real(DP), dimension(Dim) :: DeltaX
	real(DP) :: r_ij
	integer :: iDist
	real(DP) :: r
	real(DP) :: numerat, denomin
	real(DP) :: energSum 
	
	if (.not.snap) stop "Snap désactivé."
	
	call initDistriParams()
	allocate(distrib(Ndist, Ncol1))
	
	distrib(:, :) = 0
	
	do iStep = Ntherm-1, Nstep + Ntherm
	
		! Lecture :
	
		write(iSnap, "(i)") iStep
		fileSnap = trim("snap"//adjustl(iSnap))//".out"
		open(unit=unitSnapEnCours, recl=4096, file=fileSnap, &
			status='old', action='read')
			
			do iCol = 1, Ncol1
			    read(unitSnapEnCours, *) X(:, iCol)
		    end do
		    
		close(unitSnapEnCours)
    
		! Traitement
	
		do iCol = 1, Ncol1
		do jCol = 1, Ncol1
	
			if (jCol/=iCol) then
			
				DeltaX(:) = X(:, jCol) - X(:, iCol)
	            call pbc_dif(DeltaX)
	            r_ij = sqrt(dot_product(DeltaX, DeltaX))
	            
	            iDist =  int( r_ij/deltaDist )
	            distrib(iDist, iCol) = distrib(iDist, iCol) + 1
			
			end if
		
		end do		
		end do
	
	end do
	
	! Ecriture et Energie par particule
	
	energSum = 0._DP
	call ePotIni()
	
	open(unit=unitdistrib, file="fct_distrib.out", action="write")	
		do iDist = 1, Ndist
		
			r = (real(iDist, DP) + 0.5_DP) * deltaDist
			numerat = real(sum(distrib(iDist, :)), DP)
			numerat = numerat / real(Nstep, DP)
			denomin = real(Ncol1, DP) * (sphereVol(iDist+1) - sphereVol(iDist))
			denomin = denomin * densite
			write(unitdistrib, *) r, numerat / denomin
			
			if (r>=rmin .and. r<=rcut11) then
				energSum = energSum + ePot(r) * numerat/denomin * 4._DP*PI*r**2
			end if
			
		end do
	close(unitdistrib)
	
	write(*, *) "Energie par particule =", densite/2._DP * energSum * deltaDist
	
	deallocate(distrib)

end program distribution
