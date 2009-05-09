MODULE output_arb

  USE shared_data
  USE iocommon
  USE output

  IMPLICIT NONE

CONTAINS

  ! This subroutine is used to wrap a block containing program specific data
  ! Which there is no general way of allowing other programs to read
  ! It permits the use of a single string to idenitify the program that wrote it

  SUBROUTINE cfd_write_arb_block(name, class, generator_desc, &
      data_length, writer)

    CHARACTER(LEN = *), INTENT(IN) :: name, class, generator_desc
    INTEGER(8), INTENT(IN) :: data_length
    INTERFACE
      SUBROUTINE writer(filehandle, current_displacement)
        USE shared_data
        INTEGER, INTENT(IN) :: filehandle
        INTEGER(KIND = MPI_OFFSET_KIND), INTENT(IN) :: current_displacement
      END SUBROUTINE writer
    END INTERFACE

    INTEGER(KIND = 8) :: mdlength, blocklength
    INTEGER(KIND = MPI_OFFSET_KIND) :: initial_displacement

    ! Outputs general block header as described in cfd_write_block_header
    ! and then a single string

    mdlength = 1 * max_string_len
    blocklength = mdlength + data_length

    CALL cfd_write_block_header(name, class, TYPE_ARB_DB, blocklength, &
        mdlength, default_rank)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_CHARACTER, MPI_CHARACTER, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank == default_rank) CALL cfd_safe_write_string(generator_desc)

    current_displacement = current_displacement + max_string_len

    CALL writer(cfd_filehandle, current_displacement)
    current_displacement = current_displacement + data_length

  END SUBROUTINE cfd_write_arb_block

END MODULE output_arb
