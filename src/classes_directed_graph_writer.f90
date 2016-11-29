module classes_directed_graph_writer

use data_strings, only: max_line_length
use types_logical_wrapper, only: Logical_Rectangle
use classes_number_to_string, only: Concrete_Number_to_String
use types_string_wrapper, only: String_Wrapper
use procedures_string_factory, only: string_destroy => destroy
use procedures_checks, only: check_string_not_empty

implicit none

private

    type, abstract, public :: Abstract_Directed_Graph_Writer
    contains
        procedure(Abstract_destroy), deferred :: destroy
        procedure(Abstract_write), deferred :: write
    end type Abstract_Directed_Graph_Writer

    abstract interface

        subroutine Abstract_destroy(this)
        import :: Abstract_Directed_Graph_Writer
            class(Abstract_Directed_Graph_Writer), intent(inout) :: this
        end subroutine Abstract_destroy

        subroutine Abstract_write(this, i_step, i_box, adjacency_matrices)
        import :: Logical_Rectangle, Abstract_Directed_Graph_Writer
            class(Abstract_Directed_Graph_Writer), intent(in) :: this
            integer, intent(in) :: i_step, i_box
            type(Logical_Rectangle), intent(in) ::adjacency_matrices(:, :)
        end subroutine Abstract_write

    end interface

    type, extends(Abstract_Directed_Graph_Writer), public :: Dot_Directed_Graph_Writer
    private
        type(String_Wrapper), allocatable :: paths(:)
        character(len=:), allocatable :: basename, graph_name
    contains
        procedure :: construct => Dot_construct
        procedure :: destroy => Dot_destroy
        procedure :: write => Dot_write
    end type Dot_Directed_Graph_Writer

    type, extends(Abstract_Directed_Graph_Writer), public :: Null_Directed_Graph_Writer
    contains
        procedure :: destroy => Null_destroy
        procedure :: write => Null_write
    end type Null_Directed_Graph_Writer

contains

!implementation Dot_Directed_Graph_Writer

    subroutine Dot_construct(this, paths, basename, graph_name)
        class(Dot_Directed_Graph_Writer), intent(out) :: this
        type(String_Wrapper), intent(in) :: paths(:)
        character(len=*), intent(in) :: basename, graph_name

        allocate(this%paths, source=paths)
        call check_string_not_empty("Dot_Directed_Graph_Writer: construct: basename", basename)
        this%basename = basename
        call check_string_not_empty("Dot_Directed_Graph_Writer: construct: graph_name", graph_name)
        this%graph_name = graph_name
    end subroutine Dot_construct

    subroutine Dot_destroy(this)
        class(Dot_Directed_Graph_Writer), intent(inout) :: this

        if (allocated(this%graph_name)) deallocate(this%graph_name)
        if (allocated(this%basename)) deallocate(this%basename)
        call string_destroy(this%paths)
    end subroutine Dot_destroy

    subroutine Dot_write(this, i_step, i_box, adjacency_matrices)
        class(Dot_Directed_Graph_Writer), intent(in) :: this
        integer, intent(in) :: i_step, i_box
        type(Logical_Rectangle), intent(in) ::adjacency_matrices(:, :)

        integer :: graph_unit
        type(Concrete_Number_to_String) :: string
        integer :: i_component, j_component, i_particle, j_particle

        open(newunit=graph_unit, recl=max_line_length, file=this%paths(i_box)%string//this%&
            basename//"_"//string%get(i_step)//".dot", action="write")
        write(graph_unit, *) "digraph ", this%graph_name, " {"
        do j_component = 1,  size(adjacency_matrices, 2)
        do i_component = 1, size(adjacency_matrices, 1)
            do j_particle = 1, size(adjacency_matrices(i_component, j_component)%rectangle, 2)
            do i_particle = 1, size(adjacency_matrices(i_component, j_component)%rectangle, 1)
                if (adjacency_matrices(i_component, j_component)%rectangle(i_particle, j_particle))&
                    then
                    write(graph_unit, *) "    ", "c"//string%get(i_component)//"_"//string%&
                        get(i_particle)//" -> c"//string%get(j_component)//"_"//string%&
                        get(j_particle)
                end if
            end do
            end do
        end do
        end do
        write(graph_unit, *) "}"
        close(graph_unit)
    end subroutine Dot_write

!end implementation Dot_Directed_Graph_Writer

!implementation Null_Directed_Graph_Writer

    subroutine Null_destroy(this)
        class(Null_Directed_Graph_Writer), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_write(this, i_step, i_box, adjacency_matrices)
        class(Null_Directed_Graph_Writer), intent(in) :: this
        integer, intent(in) :: i_step, i_box
        type(Logical_Rectangle), intent(in) ::adjacency_matrices(:, :)
    end subroutine Null_write

!end implementation Null_Directed_Graph_Writer

end module classes_directed_graph_writer
