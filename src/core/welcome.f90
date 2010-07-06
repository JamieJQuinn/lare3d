MODULE welcome

  USE shared_data

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: welcome_message

CONTAINS

  SUBROUTINE welcome_message

    INTEGER, PARAMETER :: LOGOX = 34, LOGOY = 11
    INTEGER, DIMENSION(LOGOX, LOGOY) :: LOGO
    CHARACTER(LOGOX*2+1) :: LOGOSTRING
    CHARACTER, DIMENSION(5) :: LOGOELS
    CHARACTER, DIMENSION(4) :: clrstr = (/ ' ', '[', '2', 'J' /)
    INTEGER :: ix, iy

    IF (rank .NE. 0) RETURN

    clrstr(1) = CHAR(27)
    WRITE(*, '(1x, 4a1)') clrstr

    LOGOELS = (/ ' ', '@', " ", " ", " " /)

    PRINT *, ""
    PRINT *, ""
    LOGO(:,  1) = (/ 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, &
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 /)
    LOGO(:,  2) = (/ 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4 /)
    LOGO(:,  3) = (/ 3, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, &
        1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 4 /)
    LOGO(:,  4) = (/ 3, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, &
        0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 4 /)
    LOGO(:,  5) = (/ 3, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, &
        0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 4 /)
    LOGO(:,  6) = (/ 3, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, &
        1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 4 /)
    LOGO(:,  7) = (/ 3, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, &
        0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 4 /)
    LOGO(:,  8) = (/ 3, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, &
        0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 4 /)
    LOGO(:,  9) = (/ 3, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, &
        1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 4 /)
    LOGO(:, 10) = (/ 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4 /)
    LOGO(:, 11) = (/ 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, &
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 /)

    LOGOSTRING = " "

    DO iy = 1, LOGOY*2+1
      DO ix = 1, LOGOX
        LOGOSTRING(ix*2-1:ix*2-1) = LOGOELS(LOGO(ix, MAX(iy/2, 1)) + 1)
        LOGOSTRING(ix*2:ix*2) = LOGOELS(LOGO(ix, MAX(iy/2, 1)) + 1)
      END DO
      WRITE(*, *), LOGOSTRING
    END DO

    WRITE(*, *) ""
    WRITE(*, '("Welcome to Lare3D Version ", I1, ".", I1)'), &
        version, revision
    WRITE(*, *) ""

  END SUBROUTINE welcome_message

END MODULE welcome
