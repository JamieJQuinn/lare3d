!******************************************************************************
! Contains functions for dealing with anisotropic viscosity
!******************************************************************************

MODULE anisotropic_viscosity

  USE shared_data

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: add_braginskii_stress, calculate_heating

CONTAINS

  !****************************************************************************
  ! Add braginskii stress to qxx, qxy, etc
  !****************************************************************************

  SUBROUTINE add_braginskii_stress(&
    qxx, qxy, qxz, qyy, qyz, qzz,&
    sxx, sxy, sxz, syy, syz, szz,&
    bx, by, bz)

    REAL(num), INTENT(OUT) :: qxx, qxy, qxz, qyy, qyz, qzz
    REAL(num), INTENT(IN) :: sxx, sxy, sxz, syy, syz, szz
    REAL(num), INTENT(IN) :: bx, by, bz

    REAL(num) :: mB2, brag_visc1, brag_visc2
    REAL(num) :: wbdotb, a, b, c, d
    REAL(num) :: bsxx, bsxy, bsxz, bsyy, bsyz, bszz
    REAL(num) :: btxx, btxy, btxz, btyy, btyz, btzz

    ! Magnitude of B squared
    mB2 = bx**2 + by**2 + bz**2

    ! Viscosity parameters
    brag_visc1 = brag_visc_coeff(4.0_num*mB2)
    brag_visc2 = brag_visc_coeff(mB2)

    a = (3._num*visc3 + brag_visc1 - 4._num*brag_visc2)/(2._num*mB2**2)
    b = (brag_visc1 - visc3)/(2._num*mB2)
    c = (brag_visc2 - brag_visc1)/(mB2)
    d = brag_visc1

    ! Calculate B tensor product B
    btxx = bx**2
    btyy = by**2
    btzz = bz**2
    btxy = bx*by
    btxz = bx*bz
    btyz = by*bz

    ! Calculate WB.B
    wbdotb = calc_wbdotb(bx, by, bz, &
      sxx, sxy, sxz, syy, syz, szz)

    ! Calculate Braginskii stress
    bsxx = wbdotb*(a*btxx + b) + 2._num*d*sxx &
      + 4._num*c*(btxx*sxx + btxy*sxy + btxz*sxz)
    bsyy = wbdotb*(a*btyy + b) + 2._num*d*syy &
      + 4._num*c*(btyy*syy + btyz*syz + btxy*sxy)
    bszz = wbdotb*(a*btzz + b) + 2._num*d*szz &
      + 4._num*c*(btzz*szz + btxz*sxz + btyz*syz)

    bsxy = wbdotb*a*btxy + 2._num*d*sxy &
      + 2._num*c*(btxx* sxy + btxy* syy + btxz* syz &
                +  sxx*btxy +  sxy*btyy +  sxz*btyz)
    bsxz = wbdotb*a*btxz + 2._num*d*sxz &
      + 2._num*c*(btxx* sxz + btxy* syz + btxz* szz &
                +  sxx*btxz +  sxy*btyz +  sxz*btzz)
    bsyz = wbdotb*a*btyz + 2._num*d*syz &
      + 2._num*c*(btxy* sxz + btyy* syz + btyz* szz &
                +  sxy*btxz +  syy*btyz +  syz*btzz)

    qxx = qxx + bsxx
    qyy = qyy + bsyy
    qzz = qzz + bszz
    qxy = qxy + bsxy
    qxz = qxz + bsxz
    qyz = qyz + bsyz

    RETURN
  END SUBROUTINE add_braginskii_stress

  REAL(num) FUNCTION calc_wbdotb(bx, by, bz, sxx, sxy, sxz, syy, syz, szz)
    REAL(num), INTENT(IN) :: bx, by, bz, sxx, sxy, sxz, syy, syz, szz
    calc_wbdotb = &
      2._num*(bx*sxx + by*sxy + bz*sxz)*bx &
    + 2._num*(bx*sxy + by*syy + bz*syz)*by &
    + 2._num*(bx*sxz + by*syz + bz*szz)*bz
    RETURN
  END

  REAL(num) FUNCTION calc_wb2(bx, by, bz, sxx, sxy, sxz, syy, syz, szz)
    REAL(num), INTENT(IN) :: bx, by, bz, sxx, sxy, sxz, syy, syz, szz
    calc_wb2 = &
      4._num*(bx*sxx + by*sxy + bz*sxz)**2 &
    + 4._num*(bx*sxy + by*syy + bz*syz)**2 &
    + 4._num*(bx*sxz + by*syz + bz*szz)**2
    RETURN
  END

  REAL(num) FUNCTION brag_visc_coeff(mB2)
    REAL(num), INTENT(IN) :: mB2
    REAL(num) :: xi2
    xi2 = brag_alpha**2 * mB2
    brag_visc_coeff = visc3*(6._num/5._num*xi2 + 2.23_num)/(2.23_num + 4.03_num*xi2 + xi2**2)
    RETURN
  END

  SUBROUTINE calculate_heating(heating_array, isotropic)
    LOGICAL :: isotropic
    REAL(num), DIMENSION(:, :, :), INTENT(OUT) :: heating_array

    REAL(num) :: vxb, vxbm, vyb, vybm, vzb, vzbm
    REAL(num) :: dvxdx, dvydx, dvzdx
    REAL(num) :: dvxdy, dvydy, dvzdy
    REAL(num) :: dvxdz, dvydz, dvzdz
    REAL(num) :: dvxy, dvxz, dvyz
    REAL(num) :: sxx, syy, szz, sxy, sxz, syz
    REAL(num) :: traceW2
#ifdef BRAGINSKII_VISCOSITY
    REAL(num) :: mB2, brag_visc1, brag_visc2
    REAL(num) :: a, b, wbdotb, wb2
#endif

    heating_array = 0.0_num

    DO iz = 0, nz
      izm = iz - 1
      izp = iz + 1
      DO iy = 0, ny
        iym = iy - 1
        iyp = iy + 1
        DO ix = 0, nx
          ixm = ix - 1
          ixp = ix + 1

          ! vx at Bx(i,j,k)
          vxb  = (vx(ix ,iy ,iz ) + vx(ix ,iym,iz ) &
              +   vx(ix ,iy ,izm) + vx(ix ,iym,izm)) * 0.25_num
          ! vx at Bx(i-1,j,k)
          vxbm = (vx(ixm,iy ,iz ) + vx(ixm,iym,iz ) &
              +   vx(ixm,iy ,izm) + vx(ixm,iym,izm)) * 0.25_num
          ! vy at By(i,j,k)
          vyb  = (vy(ix ,iy ,iz ) + vy(ixm,iy ,iz ) &
              +   vy(ix ,iy ,izm) + vy(ixm,iy ,izm)) * 0.25_num
          ! vy at By(i,j-1,k)
          vybm = (vy(ix ,iym,iz ) + vy(ixm,iym,iz ) &
              +   vy(ix ,iym,izm) + vy(ixm,iym,izm)) * 0.25_num
          ! vz at Bz(i,j,k)
          vzb  = (vz(ix ,iy ,iz ) + vz(ixm,iy ,iz ) &
              +   vz(ix ,iym,iz ) + vz(ixm,iym,iz )) * 0.25_num
          ! vz at Bz(i,j,k-1)
          vzbm = (vz(ix ,iy ,izm) + vz(ixm,iy ,izm) &
              +   vz(ix ,iym,izm) + vz(ixm,iym,izm)) * 0.25_num

          dvxdx = (vxb - vxbm) / dxb(ix)
          dvydy = (vyb - vybm) / dyb(iy)
          dvzdz = (vzb - vzbm) / dzb(iz)

          ! vx at By(i,j,k)
          vxb  = (vx(ix ,iy ,iz ) + vx(ixm,iy ,iz ) &
              +   vx(ix ,iy ,izm) + vx(ixm,iy ,izm)) * 0.25_num
          ! vx at By(i,j-1,k)
          vxbm = (vx(ix ,iym,iz ) + vx(ixm,iym,iz ) &
              +   vx(ix ,iym,izm) + vx(ixm,iym,izm)) * 0.25_num
          ! vy at Bx(i,j,k)
          vyb  = (vy(ix ,iy ,iz ) + vy(ix ,iym,iz ) &
              +   vy(ix ,iy ,izm) + vy(ix ,iym,izm)) * 0.25_num
          ! vy at Bx(i-1,j,k)
          vybm = (vy(ixm,iy ,iz ) + vy(ixm,iym,iz ) &
              +   vy(ixm,iy ,izm) + vy(ixm,iym,izm)) * 0.25_num

          dvxdy = (vxb - vxbm) / dyb(iy)
          dvydx = (vyb - vybm) / dxb(ix)
          dvxy = dvxdy + dvydx

          sxy = dvxy * 0.5_num
          sxx = (2.0_num * dvxdx - dvydy - dvzdz) * third
          syy = (2.0_num * dvydy - dvxdx - dvzdz) * third
          szz = (2.0_num * dvzdz - dvxdx - dvydy) * third

          ! vx at Bz(i,j,k)
          vxb  = (vx(ix ,iy ,iz ) + vx(ixm,iy ,iz ) &
              +   vx(ix ,iym,iz ) + vx(ixm,iym,iz )) * 0.25_num
          ! vx at Bz(i,j,k-1)
          vxbm = (vx(ix ,iy ,izm) + vx(ixm,iy ,izm) &
              +   vx(ix ,iym,izm) + vx(ixm,iym,izm)) * 0.25_num
          ! vz at Bx(i,j,k)
          vzb  = (vz(ix ,iy ,iz ) + vz(ix ,iym,iz ) &
              +   vz(ix ,iy ,izm) + vz(ix ,iym,izm)) * 0.25_num
          ! vz at Bx(i-1,j,k)
          vzbm = (vz(ixm,iy ,iz ) + vz(ixm,iym,iz ) &
              +   vz(ixm,iy ,izm) + vz(ixm,iym,izm)) * 0.25_num

          dvxdz = (vxb - vxbm) / dzb(iz)
          dvzdx = (vzb - vzbm) / dxb(ix)
          dvxz = dvxdz + dvzdx

          sxz = dvxz * 0.5_num

          ! vy at Bz(i,j,k)
          vyb  = (vy(ix ,iy ,iz ) + vy(ixm,iy ,iz ) &
              +   vy(ix ,iym,iz ) + vy(ixm,iym,iz )) * 0.25_num
          ! vy at Bz(i,j,k-1)
          vybm = (vy(ix ,iy ,izm) + vy(ixm,iy ,izm) &
              +   vy(ix ,iym,izm) + vy(ixm,iym,izm)) * 0.25_num
          ! vz at By(i,j,k)
          vzb  = (vz(ix ,iy ,iz ) + vz(ixm,iy ,iz ) &
              +   vz(ix ,iy ,izm) + vz(ixm,iy ,izm)) * 0.25_num
          ! vz at By(i,j-1,k)
          vzbm = (vz(ix ,iym,iz ) + vz(ixm,iym,iz ) &
              +   vz(ix ,iym,izm) + vz(ixm,iym,izm)) * 0.25_num

          dvydz = (vyb - vybm) / dzb(iz)
          dvzdy = (vzb - vzbm) / dyb(iy)
          dvyz = dvydz + dvzdy

          syz = dvyz * 0.5_num

          IF (isotropic) THEN
            traceW2 = 4.0_num*(sxx**2 + syy**2 + szz**2 + 2*(sxy**2 + sxz**2 + syz**2))
#ifndef BRAGINSKII_VISCOSITY
#ifndef SWITCHING_VISCOSITY
            ! Heating array is offset, hence ix+1, etc
            heating_array(ix+1, iy+1, iz+1) = visc3/2.0_num*traceW2
#endif
#endif
#ifdef BRAGINSKII_VISCOSITY
            mB2 = bx(ix, iy, iz)**2 + by(ix, iy, iz)**2 + bz(ix, iy, iz)**2
            brag_visc1 = brag_visc_coeff(2*mB2)
            heating_array(ix+1, iy+1, iz+1) = brag_visc1/2.0_num*traceW2
#endif
          ELSE
#ifdef BRAGINSKII_VISCOSITY
            mB2 = bx(ix, iy, iz)**2 + by(ix, iy, iz)**2 + bz(ix, iy, iz)**2
            IF (mB2 /= 0.0_num) THEN
              brag_visc1 = brag_visc_coeff(2*mB2)
              brag_visc2 = brag_visc_coeff(mB2)
              a = (3._num*visc3 + brag_visc1 - 4._num*brag_visc2)/(4._num*mB2**2)
              b = (brag_visc2 - brag_visc1)/(mB2)
              wbdotb = calc_wbdotb(bx(ix, iy, iz), by(ix, iy, iz), bz(ix, iy, iz), &
                sxx, sxy, sxz, syy, syz, szz)
              wb2 = calc_wb2(bx(ix, iy, iz), by(ix, iy, iz), bz(ix, iy, iz), &
                sxx, sxy, sxz, syy, syz, szz)

              heating_array(ix+1, iy+1, iz+1) = a*wbdotb**2 + b*wb2
            END IF
#endif
          END IF
        END DO
      END DO
    END DO
  END SUBROUTINE calculate_heating

END MODULE anisotropic_viscosity
