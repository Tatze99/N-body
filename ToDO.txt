/* To do:
 Programming part:
    -Write .csv/.txt -file with initialvalues for ~ 20-30 objects we want to consider with lines like: x, y, z, vx, vy, vz, m;
    -Test 1d integrators for different functions and determine order of convergence with log-log-plot
    -calculate E_tot, E_pot, E_kin in driver functions and output them for conservation plot
    -calculate L_tot in driver functions and output it for conservation plot
    -calculate Laplace integral in driver functions and output it for conservation plot
    -implement adaptive stepsize control in the driver functions

    -Write a script for plots:
        -Energy
        -L
        -Laplace integral
        -trajectories in position space
        -trajectories in phase space
        (-trajectories in momentum space)
    -mabye unify the driver functions, since they are pretty much the same

 Theoretical foundations:
    -datatypes and their representation
    -numerical errors: their occurances and reduction
    -idea of integrator schemes: from taylor-expansion to final formulas
    -representation of formulas in butcher tables and how to implement code from there
    -determining order of convergence (power laws) from log-log-plots

    -from newtons law / keplers law to final functions for acceleration calculations
    (-post-newtonian/relativistic corrections)
    -transformation of data/initial data in order to avoid poorly scaling prefactors like G (pow(10,-11))

 Other challenges:
    -how to call c++ programms / functions from a python script
    -how to unify the plotting python scrips to a complete script that does programm execution, numerical analysis and plotting
    -how to handle auxiliary functions of symplectic integrators if they appear in the output at all, since they are off by 0.5 time step
    -machine leaning approach
*/