type, public :: Spheres

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
    
    !> Take a snap shot of the configuration
    procedure :: snapShot => Spheres_snapShot
    !> Do an overlap test
    procedure :: overlapTest => Spheres_overlapTest
    !> Assign all particles to cells
    procedure :: cols_to_cells => Spheres_cols_to_cells
    
    !> Adapt the displacement dx during thermalisation
    procedure :: adapt_dx => Spheres_adapt_dx
    procedure :: get_dx => Spheres_get_dx
    
end type Spheres

type, extends(Spheres), public :: InteractingSpheres

    ! Potential

    real(DP), private :: pas !< discretisation step
    integer, private :: iMin !< minimum index of tabulation
    integer, private :: Ntab !< maximum index of tabulation
    real(DP), private :: epsilon !< factor in Yukawa
    real(DP), private :: alpha !< coefficient in Yukawa
    real(DP), dimension(:), allocatable, private :: Vtab !< tabulation
    
contains

    !> Destructor of the class
    procedure :: destructor => InteractingSpheres_destructor
    !> Print a report of the component in a file
    procedure :: report => InteractingSpheres_report
    
    !> Tabulate the potential
    procedure :: ePotIni => InteractingSpheres_ePotIni
    procedure :: ePot => InteractingSpheres_ePot
    procedure :: ePotNeigh => InteractingSpheres_ePotNeigh
    procedure :: enTotCalc => InteractingSpheres_enTotCalc

    procedure :: mcMove => InteractingSpheres_mcMove
    procedure :: widom => InteractingSpheres_widom
    
end type InteractingSpheres
