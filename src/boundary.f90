!******************************************************************************
! This module contains the boundary conditions for the entire code
! Any new boundary conditions should be added here
!******************************************************************************

MODULE boundary

  USE shared_data
  USE mpiboundary

  IMPLICIT NONE

CONTAINS

  !****************************************************************************
  ! Set up any necessary variables for the chosen boundary conditions
  !****************************************************************************

  SUBROUTINE set_boundary_conditions

    ! Must be called twice
    LOGICAL, SAVE :: first_call = .TRUE.

    IF (first_call) THEN
      any_open = .FALSE.
      IF (xbc_min == BC_OPEN .OR. xbc_max == BC_OPEN &
          .OR. ybc_min == BC_OPEN .OR. ybc_max == BC_OPEN &
          .OR. zbc_min == BC_OPEN .OR. zbc_max == BC_OPEN) any_open = .TRUE.
      first_call = .FALSE.
    ELSE
      IF (xbc_min == BC_OPEN) THEN
        bx(-2,:,:) = bx(-1,:,:)
        by(-1,:,:) = by( 0,:,:)
        bz(-1,:,:) = bz( 0,:,:)
      END IF
      IF (xbc_max == BC_OPEN) THEN
        bx(nx+2,:,:) = bx(nx+1,:,:)
        by(nx+2,:,:) = by(nx+1,:,:)
        bz(nx+2,:,:) = bz(nx+1,:,:)
      END IF
      IF (ybc_min == BC_OPEN) THEN
        bx(:,-1,:) = bx(:, 0,:)
        by(:,-2,:) = by(:,-1,:)
        bz(:,-1,:) = bz(:, 0,:)
      END IF
      IF (ybc_max == BC_OPEN) THEN
        bx(:,ny+2,:) = bx(:,ny+1,:)
        by(:,ny+2,:) = by(:,ny+1,:)
        bz(:,ny+2,:) = bz(:,ny+1,:)
      END IF
      IF (zbc_min == BC_OPEN) THEN
        bx(:,:,-1) = bx(:,:, 0)
        by(:,:,-1) = by(:,:, 0)
        bz(:,:,-2) = bz(:,:,-1)
      END IF
      IF (zbc_max == BC_OPEN) THEN
        bx(:,:,nz+2) = bx(:,:,nz+1)
        by(:,:,nz+2) = by(:,:,nz+1)
        bz(:,:,nz+2) = bz(:,:,nz+1)
      END IF
    END IF

  END SUBROUTINE set_boundary_conditions



  !****************************************************************************
  ! Call all of the boundaries needed by the core Lagrangian solver
  !****************************************************************************

  SUBROUTINE boundary_conditions

    CALL bfield_bcs
    CALL energy_bcs
    CALL density_bcs
    CALL velocity_bcs
    CALL damp_boundaries

  END SUBROUTINE boundary_conditions



  !****************************************************************************
  ! Boundary conditions for magnetic field through plane
  !****************************************************************************

  SUBROUTINE bfield_bcs

    CALL bfield_mpi

    IF (proc_x_min == MPI_PROC_NULL .AND. xbc_min == BC_OTHER) THEN
      bx(-1,:,:) = 0._num
      bx(-2,:,:) = 0._num
      by( 0,:,:) = 0._num
      by(-1,:,:) = 0._num
      bz( 0,:,:) = 0._num
      bz(-1,:,:) = 0._num
    END IF

    IF (proc_x_max == MPI_PROC_NULL .AND. xbc_max == BC_OTHER) THEN
      bx(nx+1,:,:) = 0._num
      bx(nx+2,:,:) = 0._num
      by(nx+1,:,:) = 0._num
      by(nx+2,:,:) = 0._num
      bz(nx+1,:,:) = 0._num
      bz(nx+2,:,:) = 0._num
    END IF

    IF (proc_y_min == MPI_PROC_NULL .AND. ybc_min == BC_OTHER) THEN
      bx(:, 0,:) = 0._num
      bx(:,-1,:) = 0._num
      by(:,-1,:) = 0._num
      by(:,-2,:) = 0._num
      bz(:, 0,:) = 0._num
      bz(:,-1,:) = 0._num
    END IF

    IF (proc_y_max == MPI_PROC_NULL .AND. ybc_max == BC_OTHER) THEN
      bx(:,ny+1,:) = 0._num
      bx(:,ny+2,:) = 0._num
      by(:,ny+1,:) = 0._num
      by(:,ny+2,:) = 0._num
      bz(:,ny+1,:) = 0._num
      bz(:,ny+2,:) = 0._num
    END IF

    IF (proc_z_min == MPI_PROC_NULL .AND. zbc_min == BC_OTHER) THEN
      bx(:,:, 0) = 0._num
      bx(:,:,-1) = 0._num
      by(:,:, 0) = 0._num
      by(:,:,-1) = 0._num
      bz(:,:,-1) = 0._num
      bz(:,:,-2) = 0._num
    END IF

    IF (proc_z_max == MPI_PROC_NULL .AND. zbc_max == BC_OTHER) THEN
      bx(:,:,nz+1) = 0._num
      bx(:,:,nz+2) = 0._num
      by(:,:,nz+1) = 0._num
      by(:,:,nz+2) = 0._num
      bz(:,:,nz+1) = 0._num
      bz(:,:,nz+2) = 0._num
    END IF

  END SUBROUTINE bfield_bcs



  !****************************************************************************
  ! Boundary conditions for specific internal energy
  !****************************************************************************

  SUBROUTINE energy_bcs

    CALL energy_mpi

    IF (proc_x_min == MPI_PROC_NULL .AND. xbc_min == BC_OTHER) THEN
      energy( 0,:,:) = 0._num
      energy(-1,:,:) = 0._num
    END IF

    IF (proc_x_max == MPI_PROC_NULL .AND. xbc_max == BC_OTHER) THEN
      energy(nx+1,:,:) = 0._num
      energy(nx+2,:,:) = 0._num
    END IF

    IF (proc_y_min == MPI_PROC_NULL .AND. ybc_min == BC_OTHER) THEN
      energy(:, 0,:) = 0._num
      energy(:,-1,:) = 0._num
    END IF

    IF (proc_y_max == MPI_PROC_NULL .AND. ybc_max == BC_OTHER) THEN
      energy(:,ny+1,:) = 0._num
      energy(:,ny+2,:) = 0._num
    END IF

    IF (proc_z_min == MPI_PROC_NULL .AND. zbc_min == BC_OTHER) THEN
      energy(:,:, 0) = 0._num
      energy(:,:,-1) = 0._num
    END IF

    IF (proc_z_max == MPI_PROC_NULL .AND. zbc_max == BC_OTHER) THEN
      energy(:,:,nz+1) = 0._num
      energy(:,:,nz+2) = 0._num
    END IF

  END SUBROUTINE energy_bcs



  !****************************************************************************
  ! Boundary conditions for density
  !****************************************************************************

  SUBROUTINE density_bcs

    CALL density_mpi

    IF (proc_x_min == MPI_PROC_NULL .AND. xbc_min == BC_OTHER) THEN
      rho( 0,:,:) = 0._num
      rho(-1,:,:) = 0._num
    END IF

    IF (proc_x_max == MPI_PROC_NULL .AND. xbc_max == BC_OTHER) THEN
      rho(nx+1,:,:) = 0._num
      rho(nx+2,:,:) = 0._num
    END IF

    IF (proc_y_min == MPI_PROC_NULL .AND. ybc_min == BC_OTHER) THEN
      rho(:, 0,:) = 0._num
      rho(:,-1,:) = 0._num
    END IF

    IF (proc_y_max == MPI_PROC_NULL .AND. ybc_max == BC_OTHER) THEN
      rho(:,ny+1,:) = 0._num
      rho(:,ny+2,:) = 0._num
    END IF

    IF (proc_z_min == MPI_PROC_NULL .AND. zbc_min == BC_OTHER) THEN
      rho(:,:, 0) = 0._num
      rho(:,:,-1) = 0._num
    END IF

    IF (proc_z_max == MPI_PROC_NULL .AND. zbc_max == BC_OTHER) THEN
      rho(:,:,nz+1) = 0._num
      rho(:,:,nz+2) = 0._num
    END IF

  END SUBROUTINE density_bcs



  !****************************************************************************
  ! Full timestep velocity boundary conditions
  !****************************************************************************

  SUBROUTINE velocity_bcs

    CALL velocity_mpi

    IF (proc_x_min == MPI_PROC_NULL .AND. xbc_min == BC_OTHER) THEN
      vx(-2:0,:,:) = 0.0_num
      vy(-2:0,:,:) = 0.0_num
      vz(-2:0,:,:) = 0.0_num
    END IF

    IF (proc_x_max == MPI_PROC_NULL .AND. xbc_max == BC_OTHER) THEN
      vx(nx:nx+2,:,:) = 0.0_num
      vy(nx:nx+2,:,:) = 0.0_num
      vz(nx:nx+2,:,:) = 0.0_num
    END IF

    IF (proc_y_min == MPI_PROC_NULL .AND. ybc_min == BC_OTHER) THEN
      vx(:,-2:0,:) = 0.0_num
      vy(:,-2:0,:) = 0.0_num
      vz(:,-2:0,:) = 0.0_num
    END IF

    IF (proc_y_max == MPI_PROC_NULL .AND. ybc_max == BC_OTHER) THEN
      vx(:,ny:ny+2,:) = 0.0_num
      vy(:,ny:ny+2,:) = 0.0_num
      vz(:,ny:ny+2,:) = 0.0_num
    END IF

    IF (proc_z_min == MPI_PROC_NULL .AND. zbc_min == BC_OTHER) THEN
      vx(:,:,-2:0) = 0.0_num
      vy(:,:,-2:0) = 0.0_num
      vz(:,:,-2:0) = 0.0_num
    END IF

    IF (proc_z_max == MPI_PROC_NULL .AND. zbc_max == BC_OTHER) THEN
      vx(:,:,nz:nz+2) = 0.0_num
      vy(:,:,nz:nz+2) = 0.0_num
      vz(:,:,nz:nz+2) = 0.0_num
    END IF

  END SUBROUTINE velocity_bcs



  !****************************************************************************
  ! Half timestep velocity boundary conditions
  !****************************************************************************

  SUBROUTINE remap_v_bcs

    CALL remap_v_mpi

    IF (proc_x_min == MPI_PROC_NULL .AND. xbc_min == BC_OTHER) THEN
      vx1(-2:0,:,:) = 0.0_num
      vy1(-2:0,:,:) = 0.0_num
      vz1(-2:0,:,:) = 0.0_num
    END IF

    IF (proc_x_max == MPI_PROC_NULL .AND. xbc_max == BC_OTHER) THEN
      vx1(nx:nx+2,:,:) = 0.0_num
      vy1(nx:nx+2,:,:) = 0.0_num
      vz1(nx:nx+2,:,:) = 0.0_num
    END IF

    IF (proc_y_min == MPI_PROC_NULL .AND. ybc_min == BC_OTHER) THEN
      vx1(:,-2:0,:) = 0.0_num
      vy1(:,-2:0,:) = 0.0_num
      vz1(:,-2:0,:) = 0.0_num
    END IF

    IF (proc_y_max == MPI_PROC_NULL .AND. ybc_max == BC_OTHER) THEN
      vx1(:,ny:ny+2,:) = 0.0_num
      vy1(:,ny:ny+2,:) = 0.0_num
      vz1(:,ny:ny+2,:) = 0.0_num
    END IF

    IF (proc_z_min == MPI_PROC_NULL .AND. zbc_min == BC_OTHER) THEN
      vx1(:,:,-2:0) = 0.0_num
      vy1(:,:,-2:0) = 0.0_num
      vz1(:,:,-2:0) = 0.0_num
    END IF

    IF (proc_z_max == MPI_PROC_NULL .AND. zbc_max == BC_OTHER) THEN
      vx1(:,:,nz:nz+2) = 0.0_num
      vy1(:,:,nz:nz+2) = 0.0_num
      vz1(:,:,nz:nz+2) = 0.0_num
    END IF

  END SUBROUTINE remap_v_bcs



  !****************************************************************************
  ! Damped boundary conditions
  !****************************************************************************

  SUBROUTINE damp_boundaries

    ! These are not generic, always work damping options.
    ! Users should change the damping scheme for each problem

    REAL(num) :: a, d

    IF (.NOT.damping) RETURN

    IF (proc_x_min == MPI_PROC_NULL) THEN
      d = 0.7_num * x_min
      DO iz = -1, nz + 1
        DO iy = -1, ny + 1
          DO ix = -1, nx + 1
            IF (xb(ix) < d) THEN
              a = dt * (xb(ix) - d) / (x_min - d) + 1.0_num
              vx(ix,iy,iz) = vx(ix,iy,iz) / a
              vy(ix,iy,iz) = vy(ix,iy,iz) / a
              vz(ix,iy,iz) = vz(ix,iy,iz) / a
            END IF
          END DO
        END DO
      END DO
    END IF

    IF (proc_x_max == MPI_PROC_NULL) THEN
      d = 0.7_num * x_max
      DO iz = -1, nz + 1
        DO iy = -1, ny + 1
          DO ix = -1, nx + 1
            IF (xb(ix) > d) THEN
              a = dt * (xb(ix) - d) / (x_max - d) + 1.0_num
              vx(ix,iy,iz) = vx(ix,iy,iz) / a
              vy(ix,iy,iz) = vy(ix,iy,iz) / a
              vz(ix,iy,iz) = vz(ix,iy,iz) / a
            END IF
          END DO
        END DO
      END DO
    END IF

    IF (proc_y_min == MPI_PROC_NULL) THEN
      d = 0.7_num * y_min
      DO iz = -1, nz + 1
        DO iy = -1, ny + 1
          DO ix = -1, nx + 1
            IF (yb(iy) < d) THEN
              a = dt * (yb(iy) - d) / (y_min - d) + 1.0_num
              vx(ix,iy,iz) = vx(ix,iy,iz) / a
              vy(ix,iy,iz) = vy(ix,iy,iz) / a
              vz(ix,iy,iz) = vz(ix,iy,iz) / a
            END IF
          END DO
        END DO
      END DO
    END IF

    IF (proc_y_max == MPI_PROC_NULL) THEN
      d = 0.7_num * y_max
      DO iz = -1, nz + 1
        DO iy = -1, ny + 1
          DO ix = -1, nx + 1
            IF (yb(iy) > d) THEN
              a = dt * (yb(iy) - d) / (y_max - d) + 1.0_num
              vx(ix,iy,iz) = vx(ix,iy,iz) / a
              vy(ix,iy,iz) = vy(ix,iy,iz) / a
              vz(ix,iy,iz) = vz(ix,iy,iz) / a
            END IF
          END DO
        END DO
      END DO
    END IF

    IF (proc_z_min == MPI_PROC_NULL) THEN
      d = 0.7_num * z_min
      DO iz = -1, nz + 1
        DO iy = -1, ny + 1
          DO ix = -1, nx + 1
            IF (zb(iz) < d) THEN
              a = dt * (zb(iz) - d) / (z_min - d) + 1.0_num
              vx(ix,iy,iz) = vx(ix,iy,iz) / a
              vy(ix,iy,iz) = vy(ix,iy,iz) / a
              vz(ix,iy,iz) = vz(ix,iy,iz) / a
            END IF
          END DO
        END DO
      END DO
    END IF

    IF (proc_z_max == MPI_PROC_NULL) THEN
      d = 0.7_num * z_max
      DO iz = -1, nz + 1
        DO iy = -1, ny + 1
          DO ix = -1, nx + 1
            IF (zb(iz) > d) THEN
              a = dt * (zb(iz) - d) / (z_max - d) + 1.0_num
              vx(ix,iy,iz) = vx(ix,iy,iz) / a
              vy(ix,iy,iz) = vy(ix,iy,iz) / a
              vz(ix,iy,iz) = vz(ix,iy,iz) / a
            END IF
          END DO
        END DO
      END DO
    END IF

  END SUBROUTINE damp_boundaries

END MODULE boundary
