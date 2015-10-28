module procedures_parallelepiped_domain_macro

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use class_parallelepiped_domain, only: Abstract_Parallelepiped_Domain

implicit none

private
public :: domains_overlap

contains

    pure logical function domains_overlap(domain_1, domain_2)
        class(Abstract_Parallelepiped_Domain), intent(in) :: domain_1, domain_2

        real(DP), dimension(num_dimensions) :: domain_1_vertex, domain_2_vertex
        integer :: i_vertex, j_vertex, k_vertex

        domains_overlap = .false.
        do k_vertex = -1, 1, 2
            do j_vertex = -1, 1, 2
                do i_vertex = -1, 1, 2
                    domain_1_vertex = domain_1%get_origin() + &
                        real([i_vertex, j_vertex, k_vertex], DP)/2._DP * domain_1%get_size()
                    domain_2_vertex = domain_2%get_origin() + &
                        real([i_vertex, j_vertex, k_vertex], DP)/2._DP * domain_2%get_size()
                    domains_overlap = domain_2%is_inside(domain_1_vertex) .or. &
                        domain_1%is_inside(domain_2_vertex)
                    if (domains_overlap) return
                end do
            end do
        end do
    end function domains_overlap

    function domain_random_position(domain) result(random_position)
        class(Abstract_Parallelepiped_Domain), intent(in) :: domain
        real(DP), dimension(num_dimensions) :: random_position

        real(DP) :: rand_3d(num_dimensions)

        call random_number(rand_3d)
        random_position = domain%get_origin() + (rand_3d - 0.5_DP) * domain%get_size()
    end function domain_random_position

end module procedures_parallelepiped_domain_macro
