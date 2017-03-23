!******************************************************************************
! Contains functions for dealing with anisotropic viscosity
!******************************************************************************

MODULE anisotropic_viscosity

  USE shared_data

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: add_braginskii_stress, calculate_viscous_heating, add_switching_stress

CONTAINS

  !****************************************************************************
  ! Add braginskii stress to qxx, qxy, etc
  !****************************************************************************

  SUBROUTINE add_braginskii_stress(&
    qxx, qxy, qxz, qyy, qyz, qzz,&
    sxx, sxy, sxz, syy, syz, szz,&
    bx, by, bz)

    REAL(num), INTENT(INOUT) :: qxx, qxy, qxz, qyy, qyz, qzz
    REAL(num), INTENT(IN) :: sxx, sxy, sxz, syy, syz, szz
    REAL(num), INTENT(IN) :: bx, by, bz

    REAL(num) :: mB2, brag_visc1, brag_visc2
    REAL(num) :: wbdotb, a, b, c, d
    REAL(num) :: bsxx, bsxy, bsxz, bsyy, bsyz, bszz
    REAL(num) :: btxx, btxy, btxz, btyy, btyz, btzz

    ! Magnitude of B squared
    mB2 = bx**2 + by**2 + bz**2
    mB2 = MAX(mB2, none_zero)

    ! Viscosity parameters
    brag_visc1 = brag_visc_coeff(4.0_num*mB2)
    brag_visc2 = brag_visc_coeff(mB2)

    ! Braginskii tensor coefficients
    a = (3._num*visc3 + brag_visc1 - 4._num*brag_visc2)/MAX(2._num*mB2**2, none_zero)
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
    wbdotb = calc_wbdotb(bx, by, bz, sxx, sxy, sxz, syy, syz, szz)

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

  SUBROUTINE add_switching_stress(&
    qxx, qxy, qxz, qyy, qyz, qzz,&
    sxx, sxy, sxz, syy, syz, szz,&
    bx, by, bz)

    REAL(num), INTENT(INOUT) :: qxx, qxy, qxz, qyy, qyz, qzz
    REAL(num), INTENT(IN) :: sxx, sxy, sxz, syy, syz, szz
    REAL(num), INTENT(IN) :: bx, by, bz

    REAL(num) :: mB2, s, wbdotb
    REAL(num) :: bsxx, bsxy, bsxz, bsyy, bsyz, bszz
    REAL(num) :: btxx, btxy, btxz, btyy, btyz, btzz

    mB2 = bx**2 + by**2 + bz**2

    btxx = bx**2
    btyy = by**2
    btzz = bz**2
    btxy = bx*by
    btxz = bx*bz
    btyz = by*bz

    wbdotb = calc_wbdotb(bx, by, bz, sxx, sxy, sxz, syy, syz, szz)

    s = calc_switching(mB2)

    bsxx = visc3*((1.0_num-s**2)*sxx + s**2/MAX(mB2**2, none_zero)*wbdotb*(3*btxx - mB2)/2.0_num)
    bsxy = visc3*((1.0_num-s**2)*sxy + s**2/MAX(mB2**2, none_zero)*wbdotb*(3*btxy - mB2)/2.0_num)
    bsxz = visc3*((1.0_num-s**2)*sxz + s**2/MAX(mB2**2, none_zero)*wbdotb*(3*btxz - mB2)/2.0_num)
    bsyy = visc3*((1.0_num-s**2)*syy + s**2/MAX(mB2**2, none_zero)*wbdotb*(3*btyy - mB2)/2.0_num)
    bsyz = visc3*((1.0_num-s**2)*syz + s**2/MAX(mB2**2, none_zero)*wbdotb*(3*btyz - mB2)/2.0_num)
    bszz = visc3*((1.0_num-s**2)*szz + s**2/MAX(mB2**2, none_zero)*wbdotb*(3*btzz - mB2)/2.0_num)

    qxx = qxx + bsxx
    qyy = qyy + bsyy
    qzz = qzz + bszz
    qxy = qxy + bsxy
    qxz = qxz + bsxz
    qyz = qyz + bsyz

    RETURN
  END SUBROUTINE

  REAL(num) FUNCTION calc_switching(mB2)
    ! evaluates interpolation function
    REAL(num), INTENT(IN) :: mB2
    REAL(num) :: a

    ! specify dependence of concentration param a on mag field
    a = brag_alpha**2 * mB2

    IF (a < 0.1_num) THEN
      calc_switching = 0.0_num
    ELSE IF (a > 0.1_num .AND. a < 2.0_num) THEN
      calc_switching = &
        -0.1053464561e-1_num*a**3 &
        +0.3160393684e-2_num*a**2 &
        +.3308402887_num*a &
        -0.3310509816e-1_num
    ELSE IF (a > 2.0_num .AND. a < 5.0_num) THEN
      calc_switching = &
         0.6347434291e-2_num*a**3 &
        -0.9813208574e-1_num*a**2 &
        +.5334252475_num*a &
        -.1681617374_num
    ELSE IF (a > 5.0_num .AND. a < 12.0_num) THEN
      calc_switching = &
         0.1251265283e-3_num*a**3 &
        -0.4797469303e-2_num*a**2 &
        +0.6675216532e-1_num*a &
        +.6096267329_num
    ELSE IF (a < 30.0_num) THEN
      calc_switching = &
         0.5424338620e-5_num*a**3 &
        -0.4881904758e-3_num*a**2 &
        +0.1504081939e-1_num*a &
        +.8164721166_num
    ELSE
      calc_switching = 1.0_num
    END IF

    ! Ensure switching stays between 0 & 1
    calc_switching = MIN(calc_switching, 1.0_num)
    calc_switching = MAX(calc_switching, 0.0_num)

    RETURN
  END

  REAL(num) FUNCTION calc_wbdotb(bx, by, bz, sxx, sxy, sxz, syy, syz, szz)
    ! Calculates (WB) dot B
    REAL(num), INTENT(IN) :: bx, by, bz, sxx, sxy, sxz, syy, syz, szz
    calc_wbdotb = 2._num*(&
      (bx*sxx + by*sxy + bz*sxz)*bx &
    + (bx*sxy + by*syy + bz*syz)*by &
    + (bx*sxz + by*syz + bz*szz)*bz)
    RETURN
  END

  REAL(num) FUNCTION calc_wb2(bx, by, bz, sxx, sxy, sxz, syy, syz, szz)
    ! Calculates |WB|^2
    REAL(num), INTENT(IN) :: bx, by, bz, sxx, sxy, sxz, syy, syz, szz
    calc_wb2 = 4._num*(&
      (bx*sxx + by*sxy + bz*sxz)**2 &
    + (bx*sxy + by*syy + bz*syz)**2 &
    + (bx*sxz + by*syz + bz*szz)**2)
    RETURN
  END

  REAL(num) FUNCTION brag_visc_coeff(mB2)
    ! Calculates viscosity parameter
    REAL(num), INTENT(IN) :: mB2
    REAL(num) :: xi2
    xi2 = (brag_alpha**2) * mB2
    brag_visc_coeff = visc3*(6._num/5._num*xi2 + 2.23_num)/(2.23_num + 4.03_num*xi2 + xi2**2)
    RETURN
  END

  SUBROUTINE calculate_viscous_heating(heating_array, isotropic)
    LOGICAL, INTENT(IN) :: isotropic
    REAL(num), DIMENSION(:, :, :), INTENT(OUT) :: heating_array

    REAL(num) :: sxx, syy, szz, sxy, sxz, syz
    REAL(num) :: traceW2, iso_heating_coeff
#ifdef BRAGINSKII_VISCOSITY
    REAL(num) :: mB2, brag_visc1, brag_visc2
    REAL(num) :: bx_cell, by_cell, bz_cell
    REAL(num) :: a, b, wbdotb, wb2
#endif
#ifdef SWITCHING_VISCOSITY
    REAL(num) :: mB2, wbdotb
    REAL(num) :: bx_cell, by_cell, bz_cell
#endif

    heating_array = 0.0_num

    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          CALL calculate_strain_rate(sxx, sxy, sxz, syy, syz, szz, ix, iy, iz)

          IF (isotropic) THEN
            iso_heating_coeff = visc3
#ifdef BRAGINSKII_VISCOSITY
            iso_heating_coeff = brag_visc_coeff(4.0_num*calculate_mB2(ix, iy, iz))
#endif
#ifdef SWITCHING_VISCOSITY
            iso_heating_coeff = (1.0_num - calc_switching(calculate_mB2(ix, iy, iz)))*visc3
#endif
            traceW2 = 4.0_num*(sxx**2 + syy**2 + szz**2 + 2._num*(sxy**2 + sxz**2 + syz**2))
            heating_array(ix, iy, iz) = iso_heating_coeff * traceW2 * 0.5_num
          ELSE
#ifdef BRAGINSKII_VISCOSITY
            bx_cell = (bx(ix,iy,iz) + bx(ix-1,iy  ,iz  )) * 0.5_num
            by_cell = (by(ix,iy,iz) + by(ix  ,iy-1,iz  )) * 0.5_num
            bz_cell = (bz(ix,iy,iz) + bz(ix  ,iy  ,iz-1)) * 0.5_num
            mB2 = calculate_mB2(ix, iy, iz)
            brag_visc1 = brag_visc_coeff(4.0_num*mB2)
            brag_visc2 = brag_visc_coeff(mB2)
            a = (3._num*visc3 + brag_visc1 - 4._num*brag_visc2) / MAX(4._num*mB2**2, none_zero)
            b = (brag_visc2 - brag_visc1) / MAX(mB2, none_zero)
            wbdotb = calc_wbdotb(bx_cell, by_cell, bz_cell, &
              sxx, sxy, sxz, syy, syz, szz)
            wb2 = calc_wb2(bx_cell, by_cell, bz_cell, &
              sxx, sxy, sxz, syy, syz, szz)

            heating_array(ix, iy, iz) = a*wbdotb**2 + b*wb2
#endif
#ifdef SWITCHING_VISCOSITY
            bx_cell = (bx(ix,iy,iz) + bx(ix-1,iy  ,iz  )) * 0.5_num
            by_cell = (by(ix,iy,iz) + by(ix  ,iy-1,iz  )) * 0.5_num
            bz_cell = (bz(ix,iy,iz) + bz(ix  ,iy  ,iz-1)) * 0.5_num
            mB2 = calculate_mB2(ix, iy, iz)
            wbdotb = calc_wbdotb(bx_cell, by_cell, bz_cell, &
              sxx, sxy, sxz, syy, syz, szz)

            heating_array(ix, iy, iz) = 0.75_num * visc3 * calc_switching(mB2)**4 &
              / MAX(mB2**2, none_zero) * wbdotb**2
#endif
          END IF
        END DO
      END DO
    END DO
  END SUBROUTINE calculate_viscous_heating

  REAL(num) FUNCTION calculate_mB2(ix, iy, iz)
    INTEGER, INTENT(IN) :: ix, iy, iz
    REAL(num) :: bx_cell, by_cell, bz_cell

    bx_cell = (bx(ix,iy,iz) + bx(ix-1,iy  ,iz  )) * 0.5_num
    by_cell = (by(ix,iy,iz) + by(ix  ,iy-1,iz  )) * 0.5_num
    bz_cell = (bz(ix,iy,iz) + bz(ix  ,iy  ,iz-1)) * 0.5_num
    calculate_mB2 = bx_cell**2 + by_cell**2 + bz_cell**2
    RETURN
  END

  SUBROUTINE calculate_strain_rate(sxx, sxy, sxz, syy, syz, szz, ix, iy, iz)
    REAL(num), INTENT(OUT) :: sxx, sxy, sxz, syy, syz, szz
    INTEGER, INTENT(IN) :: ix, iy, iz

    REAL(num) :: vxb, vxbm, vyb, vybm, vzb, vzbm
    REAL(num) :: dvxdx, dvydx, dvzdx
    REAL(num) :: dvxdy, dvydy, dvzdy
    REAL(num) :: dvxdz, dvydz, dvzdz
    REAL(num) :: dvxy, dvxz, dvyz
    INTEGER :: ixm, iym, izm

    izm = iz - 1
    iym = iy - 1
    ixm = ix - 1

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

    RETURN
  END SUBROUTINE calculate_strain_rate

END MODULE anisotropic_viscosity
