module classes_hetero_couples

use procedures_checks, only: check_positive

implicit none

private

    type, abstract, public :: Abstract_Hetero_Couples
    private
        integer :: num_couples = 0
        integer, allocatable :: couples(:, :)
    contains
        procedure(Abstract_construct), deferred :: construct
        procedure :: destroy => Abstract_destroy
        procedure :: get_num => Abstract_get_num
        procedure :: get => Abstract_get
    end type Abstract_Hetero_Couples

    abstract interface

        subroutine Abstract_construct(this, num_partners)
        import :: Abstract_Hetero_Couples
            class(Abstract_Hetero_Couples), intent(out) :: this
            integer, intent(in) :: num_partners
        end subroutine Abstract_construct

    end interface

    type, extends(Abstract_Hetero_Couples), public :: Half_Hetero_Couples
    contains
        procedure :: construct => Half_construct
    end type Half_Hetero_Couples

    type, extends(Abstract_Hetero_Couples), public :: Full_Hetero_Couples
    contains
        procedure :: construct => Full_construct
    end type Full_Hetero_Couples

    type, extends(Abstract_Hetero_Couples), public :: Null_Hetero_Couples
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: get_num => Null_get_num
        procedure :: get => Null_get
    end type Null_Hetero_Couples

contains

!implementation Abstract_Hetero_Couples

    subroutine Abstract_destroy(this)
        class(Abstract_Hetero_Couples), intent(inout) :: this

        if (allocated(this%couples)) deallocate(this%couples)
    end subroutine Abstract_destroy

    pure integer function Abstract_get_num(this) result(num_couples)
        class(Abstract_Hetero_Couples), intent(in) :: this

        num_couples = this%num_couples
    end function Abstract_get_num

    pure function Abstract_get(this, i_couple) result(couple)
        class(Abstract_Hetero_Couples), intent(in) :: this
        integer, intent(in) :: i_couple
        integer :: couple(2)

        couple = this%couples(:, i_couple)
    end function Abstract_get

!implementation Abstract_Hetero_Couples

!implementation Half_Hetero_Couples

    subroutine Half_construct(this, num_partners)
        class(Half_Hetero_Couples), intent(out) :: this
        integer, intent(in) :: num_partners

        integer :: i_couple, i_partner, j_partner

        call check_positive("Half_Hetero_Couples: construct", "num_partners", num_partners)
        this%num_couples = num_partners * (num_partners - 1) / 2
        allocate(this%couples(2, this%num_couples))
        this%couples = 0
        i_couple = 0
        do j_partner = 2, num_partners
            do i_partner = 1, j_partner - 1
                i_couple = i_couple + 1
                this%couples(:, i_couple) = [j_partner, i_partner]
            end do
        end do
    end subroutine Half_construct

!end implementation Half_Hetero_Couples

!implementation Full_Hetero_Couples

    subroutine Full_construct(this, num_partners)
        class(Full_Hetero_Couples), intent(out) :: this
        integer, intent(in) :: num_partners

        integer :: i_couple, i_partner, j_partner
        call check_positive("Half_Hetero_Couples: construct", "num_partners", num_partners)
        this%num_couples = num_partners * (num_partners - 1)
        allocate(this%couples(2, this%num_couples))
        this%couples = 0
        i_couple = 0
        do j_partner = 1, num_partners
            do i_partner = 1, num_partners
                if (i_partner /= j_partner) then
                    i_couple = i_couple + 1
                    this%couples(:, i_couple) = [j_partner, i_partner]
                end if
            end do
        end do
    end subroutine Full_construct

!end implementation Full_Hetero_Couples

!implementation Null_Hetero_Couples

    subroutine Null_construct(this, num_partners)
        class(Null_Hetero_Couples), intent(out) :: this
        integer, intent(in) :: num_partners
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Hetero_Couples), intent(inout) :: this
    end subroutine Null_destroy

    pure integer function Null_get_num(this) result(num_couples)
        class(Null_Hetero_Couples), intent(in) :: this
        num_couples = 0
    end function Null_get_num

    pure function Null_get(this, i_couple) result(couple)
        class(Null_Hetero_Couples), intent(in) :: this
        integer, intent(in) :: i_couple
        integer :: couple(2)
        couple = 0
    end function Null_get

!implementation Null_Hetero_Couples

end module classes_hetero_couples
