# Active Flux Method for 1D Advection Equation

This repository contains a Python implementation of the **Active Flux Method** for solving the 1D linear advection equation. The code is designed to demonstrate the high-resolution, characteristics of the active flux scheme, including convergence studies.

# Files

- `af.py`: Main script to run the active flux method simulation.
- `runcoverges.sh`: Shell script to perform convergence tests over a range of grid resolutions.
- `plotrate.py`: Script to compute and plot the Experimental Order of Convergence (EOC).
- `error_avg.txt`: Output file containing errors in cell averages.
- `error_int.txt`: Output file containing errors at cell interfaces.

---
# Running the Simulation

You can run the simulation with default parameters using:
$python af.py
If you wish to change the parameters, eg:
$python af.py -cfl 0.02

# EOC
To compute EOC run as follows

$sh runconverges.sh "40 80 160 320"

$python plotrate.py error_avg.txt

$python plotrate.py error_int.txt



