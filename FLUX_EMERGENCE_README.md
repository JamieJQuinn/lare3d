# Flux Emergence with an Arcade

Initial temperature, energy and density profile, control parameters, boundary conditions, tube form, and arcade form are taken from Leake et al:

Leake, James E, Mark G Linton, and Tibor Török. ‘SIMULATIONS OF EMERGING MAGNETIC FLUX. I: THE FORMATION OF STABLE CORONAL FLUX ROPES’. Astrophysical Journal 778, no. 2 (2013). https://doi.org/10.1088/0004-637X/778/2/99.

## Running

The process for running with an arcade goes as follows:
- Generate the arcade by commenting out the tube insertion in the initial conditions and running for a while to get a stable initial condition
- Comment out the arcade code & uncomment the tube code
- Set lare to restart from the stable arcade dump

## Density limiter

It might happen that the density eventually becomes so small that dt is practically too short. This may be overcome by enabling the `LIMIT_DENSITY` define in the Makefile and setting a minimum density in `control.f90`.
