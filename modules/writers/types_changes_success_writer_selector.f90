module types_changes_success_writer_selector

implicit none

private

    type, public :: Changes_Success_Writer_Selector
        logical :: write_translations = .false.
        logical :: write_rotations = .false.
        logical :: write_exchanges = .false.
    end type Changes_Success_Writer_Selector

end module types_changes_success_writer_selector
