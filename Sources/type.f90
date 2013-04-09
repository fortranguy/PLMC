type, public :: HardSpheres

    ! Particles

    real(DP), private :: radius !< radius of a particle
    real(DP), private :: rmin !< minimum distance between two particles
    integer, private ::  Ncol !< number of a component particles
    real(DP), dimension(:, :), allocatable :: X !< position of a particle

    ! Monte-Carlo
    
    real(DP), dimension(Dim), private :: dx !< displacement

    ! Potential

    real(DP), private :: rcut !< short-range cut
    
    
    ! Neighbours (cell/grid scheme)
    
    type(Neighbours), private :: same !< same kind
    
contains
    
    !> Destructor of the class
    procedure :: destructor => HardSpheres_destructor
    !> Print a report of the component in a file
    procedure :: report => HardSpheres_report
    !> Take a snap shot of the configuration
    procedure :: snapShot => HardSpheres_snapShot
    !> Do an overlap test
    procedure :: overlapTest => HardSpheres_overlapTest
    !> Assign all particles to cells
    procedure :: cols_to_cells => HardSpheres_cols_to_cells
    
    !> Adapt the displacement dx during thermalisation
    procedure :: adapt_dx => HardSpheres_adapt_dx
    procedure :: get_dx => HardSpheres_get_dx
    
    procedure :: ePotNeigh => HardSpheres_ePotNeigh
    
    procedure :: mcMove => HardSpheres_mcMove
    procedure :: widom => HardSpheres_widom
    
end type HardSpheres

type, public, extends(HardSphere) :: SphericalCharges

    ! Potential

    real(DP), private :: pas !< discretisation step
    integer, private :: iMin !< minimum index of tabulation
    integer, private :: Ntab !< maximum index of tabulation
    real(DP), private :: epsilon !< factor in Yukawa
    real(DP), private :: alpha !< coefficient in Yukawa
    real(DP), dimension(:), allocatable, private :: Vtab !< tabulation
    
contains
    
    !> Tabulate the potential
    procedure :: ePotIni => SphericalCharges_ePotIni
    procedure :: ePot => SphericalCharges_ePot
    procedure :: ePotNeigh => SphericalCharges_ePotNeigh
    procedure :: enTotCalc => SphericalCharges_enTotCalc
    
end type SphericalCharges
