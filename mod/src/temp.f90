type1%positions(1:2, iCol) = Lsize(1:2) * xRand(1:2)
type1%positions(3, iCol) = (Height - type1%get_sigma()) * xRand(3) + type1%get_sigma()/2._DP

real(DP), parameter :: Height = Lsize1 ! u_length
real(DP), parameter :: Lsize_ratio = 1.5_DP
real(DP), parameter :: Lsize3 = Lsize_ratio * Lsize1 ! u_length

integer, parameter :: Kmax3 = ceiling(Lsize_ratio * real(Kmax1, DP)) ! 1/u_length

integer, parameter :: NnearCell_layer = NnearCell_dim(1) * NnearCell_dim(2)
integer, parameter :: NnearCell = NnearCell_layer * NnearCell_dim(3)
