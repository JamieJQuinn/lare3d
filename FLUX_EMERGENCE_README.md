# Flux Emergence

Initial temperature, energy and density profile, control parameters, boundary conditions, tube form, and arcade form are taken from Leake et al:

Leake, James E, Mark G Linton, and Tibor Török. ‘SIMULATIONS OF EMERGING MAGNETIC FLUX. I: THE FORMATION OF STABLE CORONAL FLUX ROPES’. Astrophysical Journal 778, no. 2 (2013). https://doi.org/10.1088/0004-637X/778/2/99.

This code can be run with an arcade-like or straight overlying field, see below.

## Running with straight overlying field

1. Ensure arcade isn't being generated in initial conditions
2. Set overlying field strength in initial conditions
3. Uncomment tube generation in initial conditions

## Running with arcade

The process for running with an arcade goes as follows:
1. Generate the arcade by uncommenting the arcade generation in the initial conditions and running for a while to get a stable initial condition (this needs to only be done once)
2. Comment out the arcade code
3. Set lare to restart from the stable arcade dump
4. The tube will be automatically added after the restart by calling `after_restart` in `diagnostics.f90`

## Density limiter

It might happen that the density eventually becomes so small that dt is practically too short. This may be overcome by enabling the `LIMIT_DENSITY` define in the Makefile and setting a minimum density in `control.f90`.
