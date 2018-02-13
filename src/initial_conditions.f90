MODULE initial_conditions

  USE shared_data
  USE neutral

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: set_initial_conditions

CONTAINS

  !****************************************************************************
  ! This function sets up the initial condition for the code
  ! The variables which must be set are:
  !   rho - density
  !   v{x,y,z} - Velocities in x, y, z
  !   b{x,y,z} - Magnetic fields in x, y, z
  !   energy - Specific internal energy
  !
  ! You may also need the neutral fraction. This can be calculated by a
  ! function call to get_neutral(temperature, rho, z). This routine is in
  ! core/neutral.f90 and requires the local temperature and mass density.
  ! For example to set xi_n to the neutral fraction use:
  !   xi_n = get_neutral(temperature, rho, z)
  !****************************************************************************

  SUBROUTINE set_initial_conditions

    ! Straight field

    INTEGER:: ix, iy, iz
    REAL(num) :: T0

    vx = 0.0_num
    vy = 0.0_num
    vz = 0.0_num

    bx = 0.0_num
    by = 0.0_num
    bz = 0.0_num

    DO iz = -1, nz+2
      DO iy = -1, ny+2
        DO ix = -1, nx+2
          bx(ix,iy,iz) = xb(ix)
          by(ix,iy,iz) = yb(iy)
          bz(ix,iy,iz) = -2.0_num*zb(iz)
        END DO
      END DO
    END DO

    T0 = (B0**2)*mf*mh_si/(kb_si*mu0_si*RHO0)

    energy = 1.0_num/(gamma-1.0_num)*1.0e6_num/T0
    rho = 1.0_num
    grav = 0.0_num

  END SUBROUTINE set_initial_conditions

END MODULE initial_conditions
