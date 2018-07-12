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

    INTEGER :: iz, izm, iz1
    REAL(num) :: z_photosphere = 0._num,&
      z_transition      = 10._num,&
      z_corona          = 20._num
    REAL(num), DIMENSION(:), ALLOCATABLE :: zc_global, dzb_global, dzc_global
    REAL(num), DIMENSION(:), ALLOCATABLE :: temp_ref, rho_ref
    ! Tube variables
    REAL(num) :: r2, buoyant_factor, lambda = 10.0_num, p, b1
    REAL(num) :: r0 = 2.5_num, z0 = -12._num, x0 = 0._num, alpha = 0.5_num, bFlux = 5._num
    ! Ambient dipole variables
    REAL(num) :: zd = -100._num, bd = 7.5e3_num

    ALLOCATE( zc_global(-1:nz_global+1))
    ALLOCATE(dzb_global(-1:nz_global+1))
    ALLOCATE(dzc_global(-1:nz_global+1))
    ALLOCATE(  temp_ref(-1:nz_global+2))
    ALLOCATE(   rho_ref(-1:nz_global+2))

    vx = 0.0_num
    vy = 0.0_num
    vz = 0.0_num
    bx = 0.0_num
    by = 0.0_num
    bz = 0.0_num

    grav = 1._num

    ! Fill in zc_global with the positions central to the zb_global points
    DO iz = -1, nz_global + 1
      zc_global(iz) = 0.5_num * (zb_global(iz-1) + zb_global(iz))
    END DO

    ! Fill in dzb_global and dzc_global
    DO iz = -1, nz_global
      dzb_global(iz) = zb_global(iz) - zb_global(iz-1)
      dzc_global(iz) = zc_global(iz+1) - zc_global(iz)
    END DO

    temp_ref = 1.0_num
    ! Describe temperature profile
    DO iz = -1, nz_global+1
      IF (zc_global(iz) .LE. z_photosphere) THEN
        temp_ref(iz) = 1.0_num - zc_global(iz)*(gamma-1.0_num)/gamma
      ELSEIF (zc_global(iz) .LE. z_transition) THEN
        temp_ref(iz) = 1.0_num
      ELSEIF (zc_global(iz) .LE. z_corona) THEN
        temp_ref(iz) = (150.0_num)**((zc_global(iz)-z_transition)/(z_corona-z_transition))
      ELSE
        temp_ref(iz) = 150.0_num
      ENDIF
    ENDDO

    temp_ref(nz_global+1:nz_global+2) = temp_ref(nz_global)

    ! hydrostatic profile + ideal gas law
    rho_ref = 1._num
    DO iz = nz_global, 0, -1
      IF (zc_global(iz) < z_photosphere) THEN
        izm = iz - 1
        rho_ref(izm) = rho_ref(iz ) &
          * (2._num*temp_ref(iz ) + dzb_global(iz )) &
          / (2._num*temp_ref(izm) - dzb_global(izm))
      END IF
    END DO

    DO iz = 0, nz_global
      IF (zc_global(iz) >= z_photosphere) THEN
        izm = iz - 1
        rho_ref(iz ) = rho_ref(izm) &
          * (2._num*temp_ref(izm) - dzb_global(izm)) &
          / (2._num*temp_ref(iz ) + dzb_global(iz ))
      END IF
    END DO

    rho_ref(nz_global+1:nz_global+2) = rho_ref(nz_global)
    rho_ref(-1) = rho_ref(0)

    ! Set z index local to MPI process
    iz1 = n_global_min(3) - 1

    ! Fill in energy and density
    DO iz = -1, nz + 2
      rho(:,:,iz) = rho_ref(iz1)
      energy(:,:,iz) = temp_ref(iz1)/(gamma-1.0_num)

      iz1 = iz1 + 1
    END DO

    ! Insert ambient field, symmetric in y
    DO ix = -1, nx+1
      DO iy = -1, ny+1
        DO iz = -1, nz+1
          r2 = SQRT((xb(ix))**2 + (zc(iz)-zd)**2)
          b1 = -bd*r2**(-3)*(1._num - 3._num*(zc(iz) - zd)**2*r2**(-2))
          bx(ix,iy,iz) = b1 + bx(ix,iy,iz)

          r2 = SQRT((xc(ix))**2 + (zb(iz)-zd)**2)
          b1 = -3._num*bd*(zb(iz)-zd)*xc(ix)*r2**(-5)
          bz(ix,iy,iz) = b1 + bz(ix,iy,iz)
        END DO
      END DO
    END DO

    ! Insert tube with axis pointing in y-dir
    DO ix = -1, nx+1
      DO iz = -1, nz+1
        DO iy = -1, ny+1
          r2 = (xc(ix)-x0)**2 + (zc(iz)-z0)**2
          b1 = bFlux*EXP(-r2/r0**2)
          by(ix,iy,iz) = b1 + by(ix,iy,iz)

          r2 = (xb(ix)-x0)**2 + (zc(iz)-z0)**2
          b1 = bFlux*EXP(-r2/r0**2)
          bx(ix,iy,iz) = -b1*alpha*(zc(iz)-z0) + bx(ix,iy,iz)

          r2 = (xc(ix)-x0)**2 + (zb(iz)-z0)**2
          b1 = bFlux*EXP(-r2/r0**2)
          bz(ix,iy,iz) = b1*alpha*(xc(ix)-x0) + bz(ix,iy,iz)

          r2 = (xc(ix)-x0)**2 + (zc(iz)-z0)**2
          b1 = bFlux**2 * EXP(-2.0_num*r2/r0**2) * &    ! Not field but pexc.
               (1.0_num + r2*alpha**2 - alpha**2*r0**2*0.5_num) * 0.5_num
          p = energy(ix,iy,iz)*rho(ix,iy,iz)*(gamma-1.0_num)
          buoyant_factor = EXP(-((yc(iy)/lambda)**2))
          rho(ix,iy,iz) = rho(ix,iy,iz) - rho(ix,iy,iz) * &
               b1/p*buoyant_factor
          p = p - b1
          energy(ix,iy,iz) = p/(rho(ix,iy,iz)*(gamma-1.0_num))
        END DO
      END DO
    END DO

    DEALLOCATE(zc_global, dzb_global, dzc_global)
    DEALLOCATE(temp_ref, rho_ref)
  END SUBROUTINE set_initial_conditions

END MODULE initial_conditions
