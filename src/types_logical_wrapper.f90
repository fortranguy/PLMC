module types_logical_wrapper

implicit none

private

    type, public :: Logical_Line
        logical, allocatable :: line(:)
    end type Logical_Line

    type, public :: Logical_Triangle
        type(Logical_Line), allocatable :: triangle(:)
    end type Logical_Triangle

    type, public :: Logical_Rectangle
        logical, allocatable :: rectangle(:, :)
    end type Logical_Rectangle

end module types_logical_wrapper
