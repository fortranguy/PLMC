module procedures_directed_graph_writer_factory

use procedures_errors, only: error_exit
use types_string_wrapper, only: String_Wrapper
use classes_directed_graph_writer, only: Abstract_Directed_Graph_Writer, Dot_Directed_Graph_Writer,&
    Null_Directed_Graph_Writer

implicit none

private
public :: create, destroy

contains

    subroutine create(writer, paths, basename, graph_name, write_graph)
        class(Abstract_Directed_Graph_Writer), allocatable, intent(out) :: writer
        type(String_Wrapper), intent(in) :: paths(:)
        character(len=*), intent(in) :: basename, graph_name
        logical, intent(in) :: write_graph

        if (write_graph) then
            allocate(Dot_Directed_Graph_Writer :: writer)
        else
            allocate(Null_Directed_Graph_Writer :: writer)
        end if

        select type (writer)
            type is (Dot_Directed_Graph_Writer)
                call writer%construct(paths, basename, graph_name)
            type is (Null_Directed_Graph_Writer)
            class default
                call error_exit("procedures_directed_graph_writer_factory: create: writer: "//&
                    "type unknown.")
        end select
    end subroutine create

    subroutine destroy(writer)
        class(Abstract_Directed_Graph_Writer), allocatable, intent(inout) :: writer

        if (allocated(writer)) then
            call writer%destroy()
            deallocate(writer)
        end if
    end subroutine destroy

end module procedures_directed_graph_writer_factory
