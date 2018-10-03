!******************************************************************************
! Set up the control values required by the core code
!******************************************************************************

MODULE control

  USE shared_data
  USE normalise

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: user_normalisation, control_variables, set_output_dumps

CONTAINS

  !****************************************************************************
  ! Normalisation constants
  !****************************************************************************

  SUBROUTINE user_normalisation

    ! Set the normalising constants for LARE
    ! This is needed to allow the use of some physics modules which are coded
    ! in SI units

    ! Gamma is the ratio of specific heat capacities
    gamma = 5.0_num*third

    ! Average mass of an ion in proton masses
    ! The code assumes a single ion species with this mass
    mf = 1.2_num

    ! The equations describing the normalisation in LARE have three free
    ! parameters which must be specified by the end user. These must be the
    ! normalisation used for your initial conditions. Strictly only needed for
    ! non-ideal MHD terms.

    ! Normalisation constants from Bareford & Hood 2015
    ! DOI: 10.1098/rsta.2014.0266

    ! Magnetic field normalisation in Tesla
    B0 = 0.005_num

    ! Length normalisation in m
    L0 = 1.e6_num

    ! Density normalisation in kg / m^3
    RHO0 = 1.67e-12_num

  END SUBROUTINE user_normalisation



  !****************************************************************************
  ! General control variables. Commented in detail
  !****************************************************************************

  SUBROUTINE control_variables

    ! Set the number of gridpoints in x and y directions

    ! Low res
    !nx_global = 128
    !ny_global = 128
    !nz_global = 128

    ! Mid res
    !nx_global = 256
    !ny_global = 256
    !nz_global = 256

    ! High res
    nx_global = 512
    ny_global = 512
    nz_global = 512

    ! Set the maximum number of iterations of the core solver before the code
    ! terminates. If nsteps < 0 then the code will run until t = t_end
    nsteps = -1

    ! The maximum runtime of the code
    t_end = 10.0_num

    ! Shock viscosities as detailed in manual - they are dimensionless
    visc1 = 0.0_num
    visc2 = 0.0_num

    ! Real viscosity expressed as the inverse Reynolds number
    visc3 = 1.0e-4_num

    ! Set these constants to manually override the domain decomposition.
    ! If either constant is set to zero then the code will try to automatically
    ! decompose in this direction
    nprocx = 0
    nprocy = 0
    nprocz = 0

    ! The length of the domain in the x direction
    x_min = -3.5_num
    x_max =  3.5_num
    ! Should the x grid be stretched or uniform
    x_stretch = .FALSE.

    ! The length of the domain in the y direction
    y_min = -3.5_num
    y_max =  3.5_num
    ! Should the y grid be stretched or uniform
    y_stretch = .FALSE.

    ! The length of the domain in the z direction
    z_min = -0.25_num
    z_max =  0.25_num
    ! Should the z grid be stretched or uniform
    z_stretch = .FALSE.

    ! Turn on or off the resistive parts of the MHD equations
    resistive_mhd = .TRUE.

    ! The background resistivity expressed as the inverse Lundquist number
    eta_background = 1.0e-4_num

    ! The critical current for triggering anomalous resistivity
    ! and the resistivity when above the critical current.
    ! The resistivity is expressed as the inverse Lundquist number.
    j_max = 0.0_num
    eta0 = 0.0_num

    ! Turn on or off the Braginskii thermal conduction term in
    ! the MHD equations
    ! WARNING: this is not robust. It is known to have problems
    ! with steep temperature gradients and very hot regions with
    ! large thermal conductivity. For many problems it is however
    ! fine.
    conduction = .FALSE.
    ! Apply a flux limiter to stop heat flows exceeding free streaming limit
    heat_flux_limiter = .FALSE.
    ! Fraction of free streaming heat flux used in limiter
    flux_limiter = 0.05_num

    ! Use radiation as specified in SUBROUTINE rad_losses
    ! in src/core/conduct.f90
    radiation = .FALSE.
    ! Use coronal heating as specified in SUBROUTINE heating
    ! in src/core/conduct.f90
    coronal_heating = .FALSE.

    ! Remap kinetic energy correction. LARE does not perfectly conserve kinetic
    ! energy during the remap step. This missing energy can be added back into
    ! the simulation as a uniform heating. Setting rke to true turns on this
    ! addition.
    rke = .TRUE.

    ! The code to choose the initial conditions. The valid choices are
    ! IC_NEW     - Use set_initial_conditions in "initial_conditions.f90" to
    !              setup new initial conditions
    ! IC_RESTART - Load the output file with index restart_snapshot and use it
    !              as the initial conditions
    initial = IC_NEW
    restart_snapshot = 1

    ! If cowling_resistivity is true then the code calculates and
    ! applies the Cowling Resistivity to the MHD equations
    ! only possible if not EOS_IDEAL
    ! resistive_mhd must be TRUE for this to actaully be applied
    cowling_resistivity = .FALSE.

    ! Set the boundary conditions on the four edges of the simulation domain
    ! Valid constants are
    ! BC_PERIODIC - Periodic boundary conditions
    ! BC_OPEN     - Reimann far-field characteristic boundary conditions
    ! BC_OTHER    - Other boundary conditions specified in "boundary.f90"
    xbc_min = BC_OTHER
    xbc_max = BC_OTHER
    ybc_min = BC_OTHER
    ybc_max = BC_OTHER
    zbc_min = BC_OTHER
    zbc_max = BC_OTHER

    ! Set to true to turn on routine for damped boundaries.
    ! These routines are in boundary.f90 and you should check that they
    ! actually do what you want.
    damping = .FALSE.

    ! Set the equation of state. Valid choices are
    ! EOS_IDEAL - Simple ideal gas for perfectly ionised plasma
    ! EOS_PI    - Simple ideal gas for partially ionised plasma
    ! EOS_ION   - EOS_PI plus the ionisation potential
    ! N.B. read the manual for notes on these choices
    eos_number = EOS_IDEAL
    ! EOS_IDEAL also requires that you specific whether
    ! the gas is ionised or not. Some stratified atmospheres
    ! only work for neutral hydrogen even though using MHD
    ! For fully ionised gas set .FALSE.
    ! For neutral hydrogen set .TRUE.
    ! This flag is ignored for all other EOS choices.
    neutral_gas = .TRUE.

    ! Anisotropic viscosity parameters
    ! Typical coronal value calculated from cyc_freq * collision_time / B0
    ! collision time from Braginskii 1965
    ! Coronal temp taken to be 1e6 K
    brag_alpha = 1.201632375e-4_num/RHO0
    !switching_param = brag_alpha**2
    ! This value of the switching parameter will give null points a resolution of around 10gps
    switching_param = 150.0_num

    twisting_velocity = 0.05_num
    ramp_time = 2.0_num
  END SUBROUTINE control_variables



  !****************************************************************************
  ! Output controls.
  !****************************************************************************

  SUBROUTINE set_output_dumps

    ! The output directory for the code
    data_dir = 'Data'

    ! The interval between output snapshots.
    dt_snapshots = 1.0_num

    ! dump_mask is an array which specifies which quantities the code should
    ! output to disk in a data dump.
    ! The codes are
    ! 1  - rho
    ! 2  - energy
    ! 3  - vx
    ! 4  - vy
    ! 5  - vz
    ! 6  - bx
    ! 7  - by
    ! 8  - bz
    ! 9  - temperature
    ! 10 - pressure
    ! 11 - cs (sound speed)
    ! 12 - parallel_current
    ! 13 - perp_current
    ! 14 - neutral_faction
    ! 15 - eta_perp
    ! 16 - eta
    ! 17 - jx
    ! 18 - jy
    ! 19 - jz
    ! 20 - isotropic viscous heating
    ! 21 - anisotropic viscous heating
    ! 22 - switching function
    ! 23 - bx centered
    ! 24 - by centered
    ! 25 - bz centered
    ! If a given element of dump_mask is true then that field is dumped
    ! If the element is false then the field isn't dumped
    ! N.B. if dump_mask(1:8) not true then the restart will not work
    dump_mask = .FALSE.
    dump_mask(1:8) = .TRUE.
    IF (eos_number /= EOS_IDEAL) dump_mask(14) = .TRUE.
    IF (cowling_resistivity) dump_mask(15) = .TRUE.
    IF (resistive_mhd) dump_mask(16) = .TRUE.
    dump_mask(14:16) = .FALSE.
    dump_mask(22:25) = .TRUE.

  END SUBROUTINE set_output_dumps

END MODULE control
