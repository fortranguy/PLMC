!***********************************************************************
!* MODULE: Constants                                                   *
!***********************************************************************
module data_constants

    integer, parameter :: DP = selected_real_kind(15, 307)
        ! double precision
    real(DP), parameter :: PI = acos(-1._DP)
        
end module data_constants
!***********************************************************************

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
    real(DP), parameter :: Lsize3 = Lsize1
    real(DP), dimension(Dim), parameter :: Lsize = &
        [Lsize1, Lsize2, Lsize3]
    real(DP), dimension(Dim), parameter :: LsizeMi = 0.5_DP * Lsize
    
end module data_cell
!***********************************************************************

!***********************************************************************
!* MODULE : Particles                                                  *
!* PURPOSE : declaration of the particles parameters                   *
!* COMMENT : 1=big particles ; 2=small particles                       *
!***********************************************************************
module data_particles

use data_constants
use data_cell
    
    real(DP), parameter :: rayon1 = .5_DP
    real(DP), parameter :: rmin = 1._DP
    integer, parameter ::  Ncol1 = 270
    integer, parameter :: Ncolmax = 2*Ncol1 
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
use data_cell

implicit none

    real(DP), parameter :: Tstar = 1._DP
    integer, parameter :: Nstep = 2**10
    integer, parameter :: Ntherm = 2**8
    integer, parameter :: Nmove = 2**2 * Ncol1 ! new
    real(DP), dimension(Dim), protected :: dx = 2._DP ! new, à modifier.
    
contains

    subroutine adapt_dx(iStep, tauxRejectsSum, unitRapport)
    
        integer, intent(in) :: iStep, unitRapport
        real(DP), intent(in) :: tauxRejectsSum    
        
        integer, parameter :: multiple = 2**2
        real(DP) :: tauxRejects
        real(DP), parameter :: tauxRejectsFix = 0.5_DP
        real(DP), parameter :: dx_eps = 0.05_DP, taux_eps = 0.005_DP
        real(DP), parameter :: more = 1._DP+dx_eps, less = 1._DP-dx_eps
        
        tauxRejects = 0._DP
        
        if (mod(iStep, multiple) == 0 .and. iStep>2) then
        
            tauxRejects = tauxRejectsSum/real(iStep-1, DP)
        
            if (tauxRejects < tauxRejectsFix - taux_eps) then            
                dx(:) = dx(:) * more
                dx(:) = modulo(dx(:), Lsize(:))
            else if (tauxRejects > tauxRejectsFix + taux_eps) then
                dx(:) = dx(:) * less
                dx(:) = modulo(dx(:), Lsize(:))
            end if

        end if
        
        if (iStep == Ntherm) then
        
            if (tauxRejects == 0._DP) then
                write(*, *) "Problème adaptation dx."
                stop
            end if
            
            write(unitRapport, *) "Déplacement :"
            write(unitRapport, *) "    dx(:) = ", dx(:)
            write(unitRapport, *) "    écart relatif rejet = ", &
                abs(tauxRejects - tauxRejectsFix)/tauxRejectsFix
            
        end if
    
    end subroutine adapt_dx

end module data_mc
!***********************************************************************

!***********************************************************************
!* MODULE : Potentiel                                                  *
!* COMMENT : HS + Yukawa + cut                                         *
!***********************************************************************
module data_potentiel

use data_constants
use data_particles

implicit none

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
        
        ! shift        
        Vtab11(:) = Vtab11(:) - epsilon11*exp(-alpha11*(rcut11-rmin))/rcut11

    end subroutine ePotIni
        
end module data_potentiel
!***********************************************************************

!***********************************************************************
!* MODULE: Distribution function                                       *
!***********************************************************************
module data_distrib

use data_constants
use data_cell

implicit none

	logical, parameter :: snap = .false.
	real(DP), parameter :: deltaDist = 0.01_DP
	real(DP), protected :: rMax
	integer, protected :: Ndist
	
contains
	
	subroutine initDistriParams()
	
	implicit none
	
		rMax = sqrt(dot_product(LsizeMi, LsizeMi))
		Ndist = int(rMax/deltaDist)
	
	end subroutine initDistriParams

end module data_distrib
!***********************************************************************
