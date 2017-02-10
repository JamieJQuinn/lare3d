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

    ! This is about the most complicated example for initial conditions
    ! used here as it covers including gravity and neutrals.
    ! The normalisation assumed is that from the defauls control.f90

    INTEGER :: loop
    INTEGER :: ix, iy, iz
    REAL(num) :: x, y, z

    vx = 0.0_num
    vy = 0.0_num
    vz = 0.0_num

    DO iz = -1, nz + 2
      DO iy = -1, ny + 2
        DO ix = -1, nx + 2
          x = x_min + ix/nx*(x_max - xmin)
          y = y_min + iy/ny*(y_max - ymin)
          z = z_min + iz/nz*(z_max - zmin)
          rho(ix,iy,iz) = RHO0
          energy(ix,iy,iz) = 0.0_num
          bx(ix,iy,iz) = x
          by(ix,iy,iz) = y
          bz(ix,iy,iz) = -2.0_num*z
        END DO
      END DO
    END DO

  END SUBROUTINE set_initial_conditions

END MODULE initial_conditions
