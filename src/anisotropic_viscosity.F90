!******************************************************************************
! Contains functions for dealing with anisotropic viscosity
!******************************************************************************

MODULE anisotropic_viscosity

  USE shared_data

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: add_braginskii_stress, add_switching_stress, &
    calc_iso_visc_heating_at, calc_aniso_visc_heating_at, &
    calc_max_visc_heating, calc_mB2, calc_switching2

CONTAINS

  !****************************************************************************
  ! Add braginskii stress to qxx, qxy, etc
  !****************************************************************************

  SUBROUTINE add_braginskii_stress(&
    qxx, qxy, qxz, qyy, qyz, qzz,&
    sxx, sxy, sxz, syy, syz, szz,&
    bx, by, bz,&
    rho)

    REAL(num), INTENT(INOUT) :: qxx, qxy, qxz, qyy, qyz, qzz
    REAL(num), INTENT(IN) :: sxx, sxy, sxz, syy, syz, szz
    REAL(num), INTENT(IN) :: bx, by, bz
    REAL(num), INTENT(IN) :: rho

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

    qxx = qxx + rho*bsxx
    qyy = qyy + rho*bsyy
    qzz = qzz + rho*bszz
    qxy = qxy + rho*bsxy
    qxz = qxz + rho*bsxz
    qyz = qyz + rho*bsyz

    RETURN
  END SUBROUTINE add_braginskii_stress

  !****************************************************************************
  ! Add switching viscous stress to qxx, qxy, etc
  !****************************************************************************

  SUBROUTINE add_switching_stress(&
    qxx, qxy, qxz, qyy, qyz, qzz,&
    sxx, sxy, sxz, syy, syz, szz,&
    bx, by, bz,&
    rho)

    REAL(num), INTENT(INOUT) :: qxx, qxy, qxz, qyy, qyz, qzz
    REAL(num), INTENT(IN) :: sxx, sxy, sxz, syy, syz, szz
    REAL(num), INTENT(IN) :: bx, by, bz
    REAL(num), INTENT(IN) :: rho

    REAL(num) :: mB2, s2, wbdotb
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

    s2 = calc_switching2(mB2)

    bsxx = visc3*((1.0_num-s2)*sxx*2.0_num + 1.5_num*s2/MAX(mB2**2, none_zero)*wbdotb*(btxx - mB2*third))
    bsyy = visc3*((1.0_num-s2)*syy*2.0_num + 1.5_num*s2/MAX(mB2**2, none_zero)*wbdotb*(btyy - mB2*third))
    bszz = visc3*((1.0_num-s2)*szz*2.0_num + 1.5_num*s2/MAX(mB2**2, none_zero)*wbdotb*(btzz - mB2*third))
    bsxy = visc3*((1.0_num-s2)*sxy*2.0_num + 1.5_num*s2/MAX(mB2**2, none_zero)*wbdotb*(btxy))
    bsxz = visc3*((1.0_num-s2)*sxz*2.0_num + 1.5_num*s2/MAX(mB2**2, none_zero)*wbdotb*(btxz))
    bsyz = visc3*((1.0_num-s2)*syz*2.0_num + 1.5_num*s2/MAX(mB2**2, none_zero)*wbdotb*(btyz))

    qxx = qxx + rho*bsxx
    qyy = qyy + rho*bsyy
    qzz = qzz + rho*bszz
    qxy = qxy + rho*bsxy
    qxz = qxz + rho*bsxz
    qyz = qyz + rho*bsyz

    RETURN
  END SUBROUTINE

  REAL(num) FUNCTION calc_switching2(mB2)
    ! evaluates interpolation function
    REAL(num), INTENT(IN) :: mB2
    REAL(num) :: xi2

    xi2 = (brag_alpha**2) * mB2

    calc_switching2 = xi2**2*(144.0_num*xi2**2 + 509.4_num*xi2**1 + 330.2964_num) &
      /(144.0_num*xi2**4 + 725.4_num*xi2**3 &
      + 925.8624_num*xi2**2 + 404.4105_num*xi2**1 + 44.7561_num)

    ! specify dependence of concentration param a on mag field
    !a = switching_param * mB2

    !IF (a < 0.5051_num) THEN
      !calc_switching2 = &
        !0.203207286957195e-1_num*a+0.843989591678155e-1_num*a**3
    !ELSE IF (a > 0.5051_num .AND. a < 1.54_num) THEN
      !calc_switching2 = &
        !0.159559056317913e-1_num-0.744480634040214e-1_num*a &
        !+.187623821223007_num*a**2-0.394206285197758e-1_num*a**3
    !ELSE IF (a > 1.54_num .AND. a < 3.011_num) THEN
      !calc_switching2 = &
        !-0.887980204109068e-1_num+.129618026285783_num*a &
        !+0.551133733724845e-1_num*a**2-0.107387134006150e-1_num*a**3
    !ELSE IF (a > 3.011_num .AND. a < 6.071_num) THEN
      !calc_switching2 = &
        !-0.498332905576744_num+.537656769397807_num*a &
        !-0.804026489164638e-1_num*a**2+0.426361387592081e-2_num*a**3
    !ELSE IF (a > 6.071_num .AND. a < 12.97_num) THEN
      !calc_switching2 = &
        !0.432999470914707_num+0.774365207593949e-1_num*a &
        !-0.459631575250207e-2_num*a**2+0.101403742282107e-3_num*a**3
    !ELSE IF (a > 12.97_num .AND. a < 22.59_num) THEN
      !calc_switching2 = &
        !0.609014168122852_num+0.367237920294194e-1_num*a &
        !-0.145732356052247e-2_num*a**2+0.207305941973063e-4_num*a**3
    !ELSE IF (a > 22.59_num .AND. a < 29.41_num) THEN
      !calc_switching2 = &
        !0.852464153339503_num+0.439311670189453e-2_num*a &
        !-0.261294335846660e-4_num*a**2-3.87808147946726e-7_num*a**3
    !ELSE IF (a < 29.41_num) THEN
      !calc_switching2 = &
        !0.816478868872687_num+0.806383596073046e-2_num*a &
        !-0.150941377101699e-3_num*a**2+0.102681208912720e-5_num*a**3
    !ELSE
      !calc_switching2 = 1.0_num
    !END IF

    ! Ensure switching stays between 0 & 1
    calc_switching2 = MIN(calc_switching2, 1.0_num)
    calc_switching2 = MAX(calc_switching2, 0.0_num)

    RETURN
  END FUNCTION

  REAL(num) FUNCTION calc_wbdotb(bx, by, bz, sxx, sxy, sxz, syy, syz, szz)
    ! Calculates (WB) dot B
    REAL(num), INTENT(IN) :: bx, by, bz, sxx, sxy, sxz, syy, syz, szz
    calc_wbdotb = 2._num*(&
      (bx*sxx + by*sxy + bz*sxz)*bx &
    + (bx*sxy + by*syy + bz*syz)*by &
    + (bx*sxz + by*syz + bz*szz)*bz)
    RETURN
  END FUNCTION

  REAL(num) FUNCTION calc_wb2(bx, by, bz, sxx, sxy, sxz, syy, syz, szz)
    ! Calculates |WB|^2
    REAL(num), INTENT(IN) :: bx, by, bz, sxx, sxy, sxz, syy, syz, szz
    calc_wb2 = 4._num*(&
      (bx*sxx + by*sxy + bz*sxz)**2 &
    + (bx*sxy + by*syy + bz*syz)**2 &
    + (bx*sxz + by*syz + bz*szz)**2)
    RETURN
  END FUNCTION

  REAL(num) FUNCTION brag_visc_coeff(mB2)
    ! Calculates viscosity parameter
    REAL(num), INTENT(IN) :: mB2
    REAL(num) :: xi2

    xi2 = (brag_alpha**2) * mB2
    brag_visc_coeff = visc3*(6._num/5._num*xi2 + 2.23_num)/(2.23_num + 4.03_num*xi2 + xi2**2)

    RETURN
  END FUNCTION

  REAL(num) FUNCTION calc_max_visc_heating(isotropic)
    LOGICAL, INTENT(IN) :: isotropic
    REAL(NUM) :: heating, max_heating
    LOGICAL :: anisotropic_viscosity_enabled
    anisotropic_viscosity_enabled = .FALSE.
#ifdef BRAGINSKII_VISCOSITY
    anisotropic_viscosity_enabled = .TRUE.
#endif
#ifdef SWITCHING_VISCOSITY
    anisotropic_viscosity_enabled = .TRUE.
#endif

    max_heating = 0.0_num
    heating = 0.0_num

    IF (isotropic) THEN
      DO iz = 1, nz
        DO iy = 1, ny
          DO ix = 1, nx
            heating = calc_iso_visc_heating_at(ix,iy,iz)
            IF (max_heating < heating) THEN
              max_heating = heating
            END IF
          END DO
        END DO
      END DO
    ELSE
      ! Only calculate if either anisotropic viscosity is enabled
      IF (anisotropic_viscosity_enabled) THEN
        DO iz = 1, nz
          DO iy = 1, ny
            DO ix = 1, nx
              heating = calc_aniso_visc_heating_at(ix,iy,iz)
              IF (max_heating < heating) THEN
                max_heating = heating
              END IF
            END DO
          END DO
        END DO
      END IF
    END IF

    calc_max_visc_heating = max_heating

    RETURN
  END FUNCTION

  REAL(num) FUNCTION calc_iso_visc_heating_at(ix, iy, iz)
    INTEGER, INTENT(IN) :: ix, iy, iz
    REAL(num) :: sxx, syy, szz, sxy, sxz, syz
    REAL(num) :: traceW2, iso_heating_coeff

    iso_heating_coeff = visc3
#ifdef BRAGINSKII_VISCOSITY
    iso_heating_coeff = brag_visc_coeff(4.0_num*calc_mB2(ix, iy, iz))
#endif
#ifdef SWITCHING_VISCOSITY
    iso_heating_coeff = (1.0_num - calc_switching2(calc_mB2(ix, iy, iz)))*visc3
#endif
    CALL calc_strain_rate(sxx, sxy, sxz, syy, syz, szz, ix, iy, iz)
    traceW2 = 4.0_num*(sxx**2 + syy**2 + szz**2 + 2._num*(sxy**2 + sxz**2 + syz**2))
    calc_iso_visc_heating_at = iso_heating_coeff * traceW2 * 0.5_num
    RETURN
  END FUNCTION

  REAL(num) FUNCTION calc_aniso_visc_heating_at(ix, iy, iz)
    INTEGER, INTENT(IN) :: ix, iy, iz
    REAL(num) :: sxx, syy, szz, sxy, sxz, syz
    REAL(num) :: mB2, wbdotb
    REAL(num) :: bx_cell, by_cell, bz_cell

#ifdef BRAGINSKII_VISCOSITY
    REAL(num) :: brag_visc1, brag_visc2
    REAL(num) :: a, b, wb2
#endif

    calc_aniso_visc_heating_at = 0.0_num

    CALL calc_strain_rate(sxx, sxy, sxz, syy, syz, szz, ix, iy, iz)

    bx_cell = (bx(ix,iy,iz) + bx(ix-1,iy  ,iz  )) * 0.5_num
    by_cell = (by(ix,iy,iz) + by(ix  ,iy-1,iz  )) * 0.5_num
    bz_cell = (bz(ix,iy,iz) + bz(ix  ,iy  ,iz-1)) * 0.5_num
    mB2 = calc_mB2(ix, iy, iz)
    wbdotb = calc_wbdotb(bx_cell, by_cell, bz_cell, &
      sxx, sxy, sxz, syy, syz, szz)
#ifdef BRAGINSKII_VISCOSITY
    brag_visc1 = brag_visc_coeff(4.0_num*mB2)
    brag_visc2 = brag_visc_coeff(mB2)
    a = (3._num*visc3 + brag_visc1 - 4._num*brag_visc2) / MAX(4._num*mB2**2, none_zero)
    b = (brag_visc2 - brag_visc1) / MAX(mB2, none_zero)
    wb2 = calc_wb2(bx_cell, by_cell, bz_cell, &
      sxx, sxy, sxz, syy, syz, szz)

    calc_aniso_visc_heating_at = a*wbdotb**2 + b*wb2
#endif
#ifdef SWITCHING_VISCOSITY
    calc_aniso_visc_heating_at = 0.75_num * visc3 * calc_switching2(mB2) &
      / MAX(mB2**2, none_zero) * wbdotb**2
#endif
    RETURN
  END FUNCTION

  REAL(num) FUNCTION calc_mB2(ix, iy, iz)
    INTEGER, INTENT(IN) :: ix, iy, iz
    REAL(num) :: bx_cell, by_cell, bz_cell

    bx_cell = (bx(ix,iy,iz) + bx(ix-1,iy  ,iz  )) * 0.5_num
    by_cell = (by(ix,iy,iz) + by(ix  ,iy-1,iz  )) * 0.5_num
    bz_cell = (bz(ix,iy,iz) + bz(ix  ,iy  ,iz-1)) * 0.5_num
    calc_mB2 = bx_cell**2 + by_cell**2 + bz_cell**2
    RETURN
  END FUNCTION

  SUBROUTINE calc_strain_rate(sxx, sxy, sxz, syy, syz, szz, ix, iy, iz)
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
  END SUBROUTINE calc_strain_rate

END MODULE anisotropic_viscosity
