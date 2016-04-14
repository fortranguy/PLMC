module class_box_size_reader

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions, max_line_length, max_word_length
use procedures_checks, only: check_file_exists
use class_periodic_box, only: Abstract_Periodic_Box

implicit none

private

    type, abstract, public :: Abstract_Box_Size_Reader
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: read => Abstract_read
    end type Abstract_Box_Size_Reader

    type, extends(Abstract_Box_Size_Reader), public :: Concrete_Box_Size_Reader

    end type Concrete_Box_Size_Reader

contains

    subroutine Abstract_construct(this, periodic_box)
        class(Abstract_Box_Size_Reader), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box

        this%periodic_box => periodic_box
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Box_Size_Reader), intent(inout) :: this

        this%periodic_box => null()
    end subroutine Abstract_destroy

    subroutine Abstract_read(this, filename)
        class(Abstract_Box_Size_Reader), intent(in) :: this
        character(len=*), intent(in) :: filename

        real(DP) :: size(num_dimensions)
        character(len=max_word_length) :: comment_caracter
        integer :: file_unit

        call check_file_exists(filename)
        open(newunit=file_unit, recl=max_line_length, file=filename, status="old", action="read")
        read(file_unit, *) comment_caracter ! header
        read(file_unit, *) size
        close(file_unit)

        call this%periodic_box%set(size)
    end subroutine Abstract_read

end module class_box_size_reader
