MODULE version_data

  IMPLICIT NONE
  INTEGER, PARAMETER :: i4  = SELECTED_INT_KIND(9)  ! 4-byte 2^31 ~ 10^9

  CHARACTER(LEN=*), PARAMETER :: c_code_name = 'Lare3D'
  INTEGER(i4) :: c_version, c_revision, c_minor_rev
  INTEGER(i4), PARAMETER :: c_code_io_version = 1
  CHARACTER(LEN=*), PARAMETER :: c_commit_id = _COMMIT
  CHARACTER(LEN=*), PARAMETER :: c_compile_machine = _MACHINE
  CHARACTER(LEN=*), PARAMETER :: c_compile_flags = 'unknown'
  INTEGER(i4), PARAMETER :: c_compile_date = _DATE
  CHARACTER(LEN=16) :: version_string
  CHARACTER(LEN=70) :: ascii_header

END MODULE version_data
