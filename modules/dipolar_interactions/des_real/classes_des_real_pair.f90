module classes_des_real_pair

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions, PI
use procedures_errors, only: error_exit
use procedures_checks, only: check_positive, check_potential_domain
use classes_box_size_memento, only: Abstract_Box_Size_Memento
use procedures_box_size, only: box_size_edge => edge
use classes_permittivity, only: Abstract_Permittivity
use procedures_dipolar_interactions_micro, only: des_real_B, des_real_C
use types_potential_domain, only: Concrete_Potential_Domain
use types_potential_domain_selector, only: Concrete_Potential_Domain_Selector
use classes_des_convergence_parameter, only: Abstract_DES_Convergence_Parameter

implicit none

private

    type, abstract, public :: Abstract_DES_Real_Pair
    private
        class(Abstract_Box_Size_Memento), pointer :: box_size_memento => null()
        type(Concrete_Potential_Domain) :: domain
        real(DP) :: coulomb = 0._DP
        real(DP) :: alpha_x_box_edge = 0._DP
    contains
        procedure(Abstract_construct), deferred :: construct
        procedure(Abstract_destroy), deferred :: destroy
        procedure :: target => Abstract_target
        procedure(Abstract_reset), deferred :: reset
        procedure :: meet => Abstract_meet
        procedure(Abstract_expression), private, deferred :: expression
    end type Abstract_DES_Real_Pair

    abstract interface

        subroutine Abstract_construct(this, box_size_memento, permittivity, alpha, domain)
        import :: Abstract_Box_Size_Memento, Abstract_Permittivity, Concrete_Potential_Domain, &
            Abstract_DES_Convergence_Parameter, Abstract_DES_Real_Pair
            class(Abstract_DES_Real_Pair), intent(out) :: this
            class(Abstract_Box_Size_Memento), intent(in) :: box_size_memento
            class(Abstract_Permittivity), intent(in) :: permittivity
            class(Abstract_DES_Convergence_Parameter), intent(in) :: alpha
            type(Concrete_Potential_Domain), intent(in) :: domain
        end subroutine Abstract_construct

        subroutine Abstract_destroy(this)
        import :: Abstract_DES_Real_Pair
            class(Abstract_DES_Real_Pair), intent(inout) :: this
        end subroutine Abstract_destroy

        subroutine Abstract_reset(this)
        import :: Abstract_DES_Real_Pair
            class(Abstract_DES_Real_Pair), intent(inout) :: this
        end subroutine Abstract_reset

        !> \[ r \mapsto \frac{1}{4\pi\epsilon} [B(r), C(r)] \]
        pure function Abstract_expression(this, distance) result(expression)
        import :: DP, Abstract_DES_Real_Pair
            real(DP), dimension(2) :: expression
            class(Abstract_DES_Real_Pair), intent(in) :: this
            real(DP), intent(in) :: distance !! \( r \)
        end function Abstract_expression

    end interface

    type, extends(Abstract_DES_Real_Pair), public :: Tabulated_DES_Real_Pair
    private
        real(DP), dimension(:, :), allocatable :: tabulation
    contains
        procedure :: construct => Tabulated_construct
        procedure :: destroy => Tabulated_destroy
        procedure :: reset => Tabulated_reset
        procedure, private :: set_domain => Tabulated_set_domain
        procedure, private :: create_tabulation => Tabulated_create_tabulation
        procedure, private :: expression => Tabulated_expression
    end type Tabulated_DES_Real_Pair

    type, extends(Abstract_DES_Real_Pair), public :: Raw_DES_Real_Pair
    private
        real(DP) :: alpha = 0._DP
        real(DP) :: expression_domain_max(2) = 0._DP
    contains
        procedure :: construct => Raw_construct
        procedure :: destroy => Raw_destroy
        procedure :: reset => Raw_reset
        procedure, private :: set_domain => Raw_set_domain
        procedure, private :: expression => Raw_expression
    end type Raw_DES_Real_Pair

    type, extends(Abstract_DES_Real_Pair), public :: Null_DES_Real_Pair
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: target => Null_target
        procedure :: reset => Null_reset
        procedure :: meet => Null_meet
        procedure, private :: expression => Null_expression
    end type Null_DES_Real_Pair

contains

!implementation Abstract_DES_Real_Pair

    subroutine Abstract_target(this, box_size_memento)
        class(Abstract_DES_Real_Pair), intent(inout) :: this
        class(Abstract_Box_Size_Memento), target, intent(in) :: box_size_memento

        this%box_size_memento => box_size_memento
    end subroutine Abstract_target

    !> \[
    !>      u(\vec{r}_{ij}, \vec{\mu}_i, \vec{\mu}_j) = \frac{1}{4\pi \epsilon}
    !>          \left[ (\vec{\mu}_i \cdot \vec{\mu}_j) B_\alpha(r_{ij}) -
    !>              (\vec{\mu}_i \cdot \vec{r}_{ij}) (\vec{\mu}_j \cdot \vec{r}_{ij})
    !>                  C_\alpha(r_{ij}) \right]
    !> \]
    !> cf. [[procedures_dipolar_interactions_micro:des_real_B]] &
    !> [[procedures_dipolar_interactions_micro:des_real_C]]
    pure real(DP) function Abstract_meet(this, vector_ij, moment_i, moment_j) result(energy)
        class(Abstract_DES_Real_Pair), intent(in) :: this
        real(DP), dimension(:), intent(in) :: vector_ij
        real(DP), dimension(:), intent(in) :: moment_i, moment_j

        real(DP) :: expression(2)

        expression = this%expression(norm2(vector_ij))
        energy = dot_product(moment_i, moment_j) * expression(1) - &
            dot_product(moment_i, vector_ij) * dot_product(moment_j, vector_ij) * expression(2)
    end function Abstract_meet

!end implementation Abstract_DES_Real_Pair

!implementation Tabulated_DES_Real_Pair

    !> @note create_tabulation() is delayed in reset()
    subroutine Tabulated_construct(this, box_size_memento, permittivity, alpha, domain)
        class(Tabulated_DES_Real_Pair), intent(out) :: this
        class(Abstract_Box_Size_Memento), intent(in) :: box_size_memento
        class(Abstract_Permittivity), intent(in) :: permittivity
        class(Abstract_DES_Convergence_Parameter), intent(in) :: alpha
        type(Concrete_Potential_Domain), intent(in) :: domain

        call this%target(box_size_memento)
        this%coulomb = 1._DP / (4._DP*PI * permittivity%get())
        this%alpha_x_box_edge = alpha%get_times_box_edge()
        call this%set_domain(domain)
    end subroutine Tabulated_construct

    !> @note Since box_size_memento may not be initialised, max_over_box_edge won't be checked now.
    subroutine Tabulated_set_domain(this, domain)
        class(Tabulated_DES_Real_Pair), intent(inout) :: this
        type(Concrete_Potential_Domain), intent(in) :: domain

        type(Concrete_Potential_Domain_Selector) :: selector

        this%domain%min = domain%min
        this%domain%max_over_box_edge = domain%max_over_box_edge
        this%domain%max = this%domain%max_over_box_edge * box_size_edge(this%box_size_memento%get())
        this%domain%delta = domain%delta
        selector%check_max = .false.
        selector%check_max_over_box_edge = .true.
        selector%check_delta = .true.
        call check_potential_domain("Tabulated_DES_Real_Pair: set_domain", this%domain, selector)
    end subroutine Tabulated_set_domain

    pure subroutine Tabulated_create_tabulation(this)
        class(Tabulated_DES_Real_Pair), intent(inout) :: this

        real(DP) :: alpha, distance_i
        integer :: i_min, i_max, i_distance

        i_min = int(this%domain%min / this%domain%delta)
        i_max = int(this%domain%max / this%domain%delta) + 1
        allocate(this%tabulation(i_min:i_max, 2))
        alpha = this%alpha_x_box_edge / box_size_edge(this%box_size_memento%get())
        do i_distance = i_min, i_max
            distance_i = real(i_distance, DP) * this%domain%delta
            this%tabulation(i_distance, 1) = des_real_B(alpha, distance_i)
            this%tabulation(i_distance, 2) = des_real_C(alpha, distance_i)
        end do
        this%tabulation(:, 1) = this%coulomb * (this%tabulation(:, 1) - this%tabulation(i_max, 1))
        this%tabulation(:, 2) = this%coulomb * (this%tabulation(:, 2) - this%tabulation(i_max, 2))
    end subroutine Tabulated_create_tabulation

    subroutine Tabulated_destroy(this)
        class(Tabulated_DES_Real_Pair), intent(inout) :: this

        if (allocated(this%tabulation)) deallocate(this%tabulation)
        this%box_size_memento => null()
    end subroutine Tabulated_destroy

    subroutine Tabulated_reset(this)
        class(Tabulated_DES_Real_Pair), intent(inout) :: this

        this%domain%max = this%domain%max_over_box_edge * box_size_edge(this%box_size_memento%get())
        if (this%domain%min > this%domain%max) then
            call error_exit("Tabulated_DES_Real_Pair: reset: domain: min > max.")
        end if
        if (allocated(this%tabulation)) deallocate(this%tabulation)
        call this%create_tabulation()
    end subroutine Tabulated_reset

    !> Linear interpolation
    pure function Tabulated_expression(this, distance) result(expression)
        real(DP) :: expression(2)
        class(Tabulated_DES_Real_Pair), intent(in) :: this
        real(DP), intent(in) :: distance

        integer :: i_distance
        real(DP) :: distance_i

        if (distance < this%domain%max) then
            i_distance = int(distance/this%domain%delta)
            distance_i = real(i_distance, DP) * this%domain%delta
            expression = this%tabulation(i_distance, :) + &
                (distance - distance_i) * &
                (this%tabulation(i_distance + 1, :) - this%tabulation(i_distance, :)) / &
                this%domain%delta
        else
            expression = 0._DP
        end if
    end function Tabulated_expression

!end implementation Tabulated_DES_Real_Pair

!implementation Raw_DES_Real_Pair

    !> @note this%alpha and this%expression_domain_max will be set later, cf.
    !> [[classes_des_real_pair:Raw_reset]].
    subroutine Raw_construct(this, box_size_memento, permittivity, alpha, domain)
        class(Raw_DES_Real_Pair), intent(out) :: this
        class(Abstract_Box_Size_Memento), intent(in) :: box_size_memento
        class(Abstract_Permittivity), intent(in) :: permittivity
        class(Abstract_DES_Convergence_Parameter), intent(in) :: alpha
        type(Concrete_Potential_Domain), intent(in) :: domain

        call this%target(box_size_memento)
        this%coulomb = 1._DP / (4._DP * PI * permittivity%get())
        this%alpha_x_box_edge = alpha%get_times_box_edge()
        call this%set_domain(domain)
    end subroutine Raw_construct

    !> note cf. [[classes_des_real_pair:Tabulated_set_domain]]
    subroutine Raw_set_domain(this, domain)
        class(Raw_DES_Real_Pair), intent(inout) :: this
        type(Concrete_Potential_Domain), intent(in) :: domain

        type(Concrete_Potential_Domain_Selector) :: selector

        this%domain%min = domain%min
        this%domain%max_over_box_edge = domain%max_over_box_edge
        this%domain%max = this%domain%max_over_box_edge * box_size_edge(this%box_size_memento%get())
        this%domain%delta = 0._DP
        selector%check_max = .false.
        selector%check_max_over_box_edge = .true.
        selector%check_delta = .false.
        call check_potential_domain("Raw_DES_Real_Pair: set_domain", this%domain, selector)
    end subroutine Raw_set_domain

    subroutine Raw_destroy(this)
        class(Raw_DES_Real_Pair), intent(inout) :: this

        this%box_size_memento => null()
    end subroutine Raw_destroy

    subroutine Raw_reset(this)
        class(Raw_DES_Real_Pair), intent(inout) :: this

        this%domain%max = this%domain%max_over_box_edge * box_size_edge(this%box_size_memento%get())
        if (this%domain%min > this%domain%max) then
            call error_exit("Raw_DES_Real_Pair: reset: domain: min > max.")
        end if
        this%alpha = this%alpha_x_box_edge / box_size_edge(this%box_size_memento%get())
        this%expression_domain_max =  [des_real_B(this%alpha, this%domain%max), &
            des_real_C(this%alpha, this%domain%max)]
    end subroutine Raw_reset

    pure function Raw_expression(this, distance) result(expression)
        real(DP) :: expression(2)
        class(Raw_DES_Real_Pair), intent(in) :: this
        real(DP), intent(in) :: distance

        if (distance < this%domain%max) then
            expression = [des_real_B(this%alpha, distance), des_real_C(this%alpha, distance)]
            expression = this%coulomb * (expression - this%expression_domain_max)
        else
            expression = 0._DP
        end if
    end function Raw_expression

!end implementation Raw_DES_Real_Pair

!implementation Null_DES_Real_Pair

    subroutine Null_construct(this, box_size_memento, permittivity, alpha, domain)
        class(Null_DES_Real_Pair), intent(out) :: this
        class(Abstract_Box_Size_Memento), intent(in) :: box_size_memento
        class(Abstract_Permittivity), intent(in) :: permittivity
        class(Abstract_DES_Convergence_Parameter), intent(in) :: alpha
        type(Concrete_Potential_Domain), intent(in) :: domain
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_DES_Real_Pair), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_target(this, box_size_memento)
        class(Null_DES_Real_Pair), intent(inout) :: this
        class(Abstract_Box_Size_Memento), target, intent(in) :: box_size_memento
    end subroutine Null_target

    subroutine Null_reset(this)
        class(Null_DES_Real_Pair), intent(inout) :: this
    end subroutine Null_reset

    pure function Null_expression(this, distance) result(expression)
        real(DP) :: expression(2)
        class(Null_DES_Real_Pair), intent(in) :: this
        real(DP), intent(in) :: distance
        expression = 0._DP
    end function Null_expression

    pure real(DP) function Null_meet(this, vector_ij, moment_i, moment_j) result(energy)
        class(Null_DES_Real_Pair), intent(in) :: this
        real(DP), dimension(:), intent(in) :: vector_ij
        real(DP), dimension(:), intent(in) :: moment_i, moment_j
        energy = 0._DP
    end function Null_meet

!end implementation Null_DES_Real_Pair

end module classes_des_real_pair
