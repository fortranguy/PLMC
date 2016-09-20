module classes_des_real_pair

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions, PI
use procedures_errors, only: error_exit
use procedures_checks, only: check_positive, check_potential_domain
use classes_box_volume_memento, only: Abstract_Box_Volume_Memento
use classes_permittivity, only: Abstract_Permittivity
use procedures_dipolar_interactions_micro, only: des_real_B, des_real_C
use types_potential_domain, only: Dipolar_Potential_Domain
use classes_des_convergence_parameter, only: Abstract_DES_Convergence_Parameter

implicit none

private

    type, abstract, public :: Abstract_DES_Real_Pair
    private
        class(Abstract_Box_Volume_Memento), pointer :: box_volume_memento => null()
        type(Dipolar_Potential_Domain) :: domain
        real(DP) :: domain_max = 0._DP
        real(DP) :: coulomb = 0._DP
        class(Abstract_DES_Convergence_Parameter), pointer :: alpha => null()
    contains
        procedure(Abstract_construct), deferred :: construct
        procedure(Abstract_destroy), deferred :: destroy
        procedure :: target => Abstract_target
        procedure(Abstract_reset), deferred :: reset
        procedure :: meet => Abstract_meet
        procedure(Abstract_expression), private, deferred :: expression
    end type Abstract_DES_Real_Pair

    abstract interface

        subroutine Abstract_construct(this, box_volume_memento, permittivity, alpha, domain)
        import :: Abstract_Box_Volume_Memento, Abstract_Permittivity, Dipolar_Potential_Domain, &
            Abstract_DES_Convergence_Parameter, Abstract_DES_Real_Pair
            class(Abstract_DES_Real_Pair), intent(out) :: this
            class(Abstract_Box_Volume_Memento), intent(in) :: box_volume_memento
            class(Abstract_Permittivity), intent(in) :: permittivity
            class(Abstract_DES_Convergence_Parameter), intent(in) :: alpha
            type(Dipolar_Potential_Domain), intent(in) :: domain
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

    subroutine Abstract_target(this, alpha)
        class(Abstract_DES_Real_Pair), intent(inout) :: this
        class(Abstract_DES_Convergence_Parameter), target, intent(in) :: alpha

        this%alpha => alpha
    end subroutine Abstract_target

    !> \[
    !>      u(\vec{r}_{ij}, \vec{\mu}_i, \vec{\mu}_j) = \frac{1}{4\pi \epsilon}
    !>          \left[ (\vec{\mu}_i \cdot \vec{\mu}_j) B_\alpha(r_{ij}) -
    !>              (\vec{\mu}_i \cdot \vec{r}_{ij}) (\vec{\mu}_j \cdot \vec{r}_{ij})
    !>                  C_\alpha(r_{ij})
    !>          \right]
    !> \]
    !> cf. [[procedures_dipolar_interactions_micro:des_real_B]] &
    !> [[procedures_dipolar_interactions_micro:des_real_C]]
    pure real(DP) function Abstract_meet(this, vector_ij, moment_i, moment_j) result(energy)
        class(Abstract_DES_Real_Pair), intent(in) :: this
        real(DP), dimension(:), intent(in) :: vector_ij
        real(DP), dimension(:), intent(in) :: moment_i, moment_j

        real(DP), dimension(2) :: coefficient

        coefficient(1) = dot_product(moment_i, moment_j)
        coefficient(2) =-dot_product(moment_i, vector_ij) * dot_product(moment_j, vector_ij)
        energy = dot_product(coefficient, this%expression(norm2(vector_ij)))
    end function Abstract_meet

!end implementation Abstract_DES_Real_Pair

!implementation Tabulated_DES_Real_Pair

    subroutine Tabulated_construct(this, box_volume_memento, permittivity, alpha, domain)
        class(Tabulated_DES_Real_Pair), intent(out) :: this
        class(Abstract_Box_Volume_Memento), intent(in) :: box_volume_memento
        class(Abstract_Permittivity), intent(in) :: permittivity
        class(Abstract_DES_Convergence_Parameter), intent(in) :: alpha
        type(Dipolar_Potential_Domain), intent(in) :: domain

        this%coulomb = 1._DP / (4._DP*PI * permittivity%get())
        call this%target(alpha)
        call this%set_domain(domain)
        call this%create_tabulation()
    end subroutine Tabulated_construct

    subroutine Tabulated_set_domain(this, domain)
        class(Tabulated_DES_Real_Pair), intent(inout) :: this
        type(Dipolar_Potential_Domain), intent(in) :: domain

        call check_potential_domain("Tabulated_DES_Real_Pair: set_domain", domain, this%alpha%&
            get_box_edge())
        this%domain%min = domain%min
        this%domain%max_over_box = domain%max_over_box
        this%domain_max = this%domain%max_over_box * this%alpha%get_box_edge()
        this%domain%delta = domain%delta
    end subroutine Tabulated_set_domain

    pure subroutine Tabulated_create_tabulation(this)
        class(Tabulated_DES_Real_Pair), intent(inout) :: this

        real(DP) :: alpha, distance_i
        integer :: i_min, i_max, i_distance

        alpha = this%alpha%get()
        i_min = int(this%domain%min/this%domain%delta)
        i_max = int(this%domain_max/this%domain%delta) + 1
        allocate(this%tabulation(i_min:i_max, 2))
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
        this%alpha => null()
        this%box_volume_memento => null()
    end subroutine Tabulated_destroy

    subroutine Tabulated_reset(this)
        class(Tabulated_DES_Real_Pair), intent(inout) :: this

        if (allocated(this%tabulation)) deallocate(this%tabulation)
        this%domain_max = this%domain%max_over_box * this%alpha%get_box_edge()
        call this%create_tabulation()
    end subroutine Tabulated_reset

    !> Linear interpolation
    pure function Tabulated_expression(this, distance) result(expression)
        real(DP) :: expression(2)
        class(Tabulated_DES_Real_Pair), intent(in) :: this
        real(DP), intent(in) :: distance

        integer :: i_distance
        real(DP) :: distance_i

        if (distance < this%domain_max) then
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

    subroutine Raw_construct(this, box_volume_memento, permittivity, alpha, domain)
        class(Raw_DES_Real_Pair), intent(out) :: this
        class(Abstract_Box_Volume_Memento), intent(in) :: box_volume_memento
        class(Abstract_Permittivity), intent(in) :: permittivity
        class(Abstract_DES_Convergence_Parameter), intent(in) :: alpha
        type(Dipolar_Potential_Domain), intent(in) :: domain

        this%coulomb = 1._DP / (4._DP * PI * permittivity%get())
        call this%target(alpha)
        call this%set_domain(domain)
    end subroutine Raw_construct

    subroutine Raw_set_domain(this, domain)
        class(Raw_DES_Real_Pair), intent(inout) :: this
        type(Dipolar_Potential_Domain), intent(in) :: domain

        call check_positive("Raw_DES_Real_Pair", "domain%min", domain%min)
        call check_positive("Raw_DES_Real_Pair", "domain%max_over_box", domain%max_over_box)
        this%domain%min = domain%min
        this%domain%max_over_box = domain%max_over_box
        this%domain_max = this%domain%max_over_box * this%alpha%get_box_edge()
        if (this%domain%min > this%domain_max) then
            call error_exit("Raw_DES_Real_Pair: set_domain: min > max.")
        end if
        this%domain%delta = 0._DP
    end subroutine Raw_set_domain

    subroutine Raw_destroy(this)
        class(Raw_DES_Real_Pair), intent(inout) :: this

        this%alpha => null()
        this%box_volume_memento => null()
    end subroutine Raw_destroy

    subroutine Raw_reset(this)
        class(Raw_DES_Real_Pair), intent(inout) :: this

        this%domain_max = this%domain%max_over_box * this%alpha%get_box_edge()
    end subroutine Raw_reset

    pure function Raw_expression(this, distance) result(expression)
        real(DP) :: expression(2)
        class(Raw_DES_Real_Pair), intent(in) :: this
        real(DP), intent(in) :: distance

        real(DP) :: expression_domain_max(2)
        real(DP) :: alpha

        if (distance < this%domain_max) then
            alpha = this%alpha%get()
            expression = [des_real_B(alpha, distance), des_real_C(alpha, distance)]
            expression_domain_max =  [des_real_B(alpha, this%domain_max), &
                des_real_C(alpha, this%domain_max)]
            expression = this%coulomb * (expression - expression_domain_max)
        else
            expression = 0._DP
        end if
    end function Raw_expression

!end implementation Raw_DES_Real_Pair

!implementation Null_DES_Real_Pair

    subroutine Null_construct(this, box_volume_memento, permittivity, alpha, domain)
        class(Null_DES_Real_Pair), intent(out) :: this
        class(Abstract_Box_Volume_Memento), intent(in) :: box_volume_memento
        class(Abstract_Permittivity), intent(in) :: permittivity
        class(Abstract_DES_Convergence_Parameter), intent(in) :: alpha
        type(Dipolar_Potential_Domain), intent(in) :: domain
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_DES_Real_Pair), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_target(this, alpha)
        class(Null_DES_Real_Pair), intent(inout) :: this
        class(Abstract_DES_Convergence_Parameter), target, intent(in) :: alpha
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
