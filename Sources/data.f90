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
!* COMMENT : size1=x, y; size2=z                                       *
!***********************************************************************
    module data_cell
    use data_constants
    integer, parameter :: Dim = 3
    real(DP), parameter :: Lsize1 = 16._DP
    real(DP), parameter :: Lsize2 = 16._DP
    real(DP), dimension(Dim), parameter :: Lsize = &
        [Lsize1, Lsize1, Lsize2]
    real(DP), dimension(Dim), parameter :: LsizeMi = 0.5_DP * Lsize
    real(DP), parameter::Tstar = 1._DP
    end module data_cell
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
    integer, parameter ::  Ncol1 = 80 ! Vs Ncolmax
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
    integer, parameter :: Nstep = 1000 ! 10000
    integer, parameter :: Ntherm = 1000 ! 1
    integer, parameter :: Nmove = 1000 ! 2**14 ! new
    real(DP), dimension(Dim) :: dx = 2._DP ! new, à tester ?
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
    ! comprendre : surpas ?
    integer, parameter :: Ntab11 = int(rcut11*surpas11) ! new
    real(DP), dimension(Ntab11) :: Vtab11 ! new 
    real(DP), parameter :: epsilon11=1._DP, alpha11=5._DP ! new
    integer, parameter :: iMin = int( rmin/rcut11*real(Ntab11, DP) ) + 1 ! new
    end module data_potentiel
!***********************************************************************
