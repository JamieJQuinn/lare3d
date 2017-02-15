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

    vx = 0.0_num
    vy = 0.0_num
    vz = 0.0_num
    bx = 0.0_num
    by = 0.0_num
    bz = 0.0_num
    rho = 1.0_num
    energy = 1.0_num

    DO ix = -2, nx + 2
      bx(ix,:,:) = xb(ix)
    END DO
    DO iy = -2, ny + 2
      by(:,iy,:) = yc(iy)
    END DO
    DO iz = -2, nz + 2
      bz(:,:,iz) = -2.0_num*zc(iz)
    END DO

  END SUBROUTINE set_initial_conditions

END MODULE initial_conditions
