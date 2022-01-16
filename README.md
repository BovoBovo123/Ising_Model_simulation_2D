# Ising_Model_simulation_2D

# Introduction

The subject of this project is the two dimensional Ising model. This is used to describe the interaction of spins, either spin up (+1) or spin down (-1), on a two-dimensional lattice (e.g. in a ferromagnet). At low temperatures the system is ordered, i.e. almost all spins are aligned; raising the temperature above a critical temperature Tc however causes thermal fluctuations to destroy the order and there is no preferred spin alignment any more. The task is to study this phase transition numerically using a Monte Carlo simulation in a Python code.

# Ising model

The Ising Hamiltonian is :

![equation1](https://latex.codecogs.com/svg.image?H&space;=&space;-J&space;\sum_{\langle&space;ij&space;\rangle}&space;\sigma_i&space;\sigma_j&space;&space;)

where:

- the spins <img src="https://latex.codecogs.com/svg.image?\sigma&space;" title="\sigma " /> can take values ±1;
- <ij> implies nearest-neighbor interactions only;
- J > 0 is the strength of exchange interaction; can be taken J = 1 to have dimensionless energy.  
  
This system undergoes a second order phase transition at the critical temperature Tc. For temperatures less than Tc, the system magnetizes, and the state is called the ferromagnetic or the ordered state. This amounts to a globally ordered state due to the presence of local interactions between the spin. For temperatures greater than Tc, the system is in the disordered or the paramagnetic state. In this case, there are no long-range correlations between the spins.

The order parameter for this system is the average magnetization: 
  
![equation2](https://latex.codecogs.com/svg.image?m&space;=&space;\frac{1}{N}&space;\sum_{i}&space;\sigma_i&space;)

and distinguishes the two phases realized by the systems: it is zero in the disordered paramagnetic state, while non-zero in the ordered ferromagnetic state.

# Monte Carlo simulation

The following code simulates the Ising model in 2D using the Metropolis algorithm. The main steps of Metropolis algorithm are:

- prepare an initial configuration of N spins;
- flip the spin of a randomly chosen lattice site;
- calculate the change in energy dE.
- if dE < 0, accept the move. Otherwise accept the move with probability exp^{-dE/T}. This satisfies the detailed balance condition, ensuring a final equilibrium state.

repeating the last three steps.
              
The code below can calculate and plot the energy and magnetization vs temperature and number of steps, as well as showing a representation of the lattice state at different "times" (i.e. number of steps).

# Structure of the project
## Python libraries used
The following libraries are used in the code:

- numpy 1.21.2
- matplotlib 3.4.3
- for an informative progress bar: tqdm 4.62.3
- logging, configparser and sys from python 3.9.7
- for testing: pytest 6.2.5
            
## Modules

In this project there are 2 modules for the definition of functions to study, plot and save relevant quantities of the system, and 1 module for simulation that calls them in an
example of execution. There is also a configuration file that is used in the example, and a testing module (functions that plot are not tested).
           
### functions_ising
            
Here a lattice of given dimensions can be created, with a spin configuration that can be random or polarized. The lattice can then be updated, simulating the Metropolis step at a certain inverse (dimensionless, putting the Boltzmann constant k = 1) temperature; energy and magnetization can be calculated. The lattice evolution configuration at certain time instants can be stored for later plotting.
Lattice parameters can be read from a configuration file, and energy and magnetization data can be saved in save files.
Logging is used to inform the user about some good practices for the functions.
            
### plots_ising
     
Here mean energy and magnetization can be plotted vs temperature to study the phase transition, as well as energy and magnetization at a specific temperature vs number of steps to study lattice thermalization and equilibrium. The lattice can be also visualized at specific time instants as a colored mesh.
     
### simulation            
     
Here all lattice parameters are read from the configuration file. A lattice if first created and then studied in a range of temperature, acquiring instantaneous data (i.e. step by step) as well as mean data vs temperature. Plots are then shown for the relevant quantities.
            
### configuration

Here lattice parameters can be changed to change the system conditions for the simulations. In particular in this example the user can specify: lattice dimensions and starting polarization, temperature values to study and visualize the system, number of steps to wait for thermalization and to average thermodinamical quantities, as well as visualize the system, and paths for save files or data load files.
            
### tests
            
Here properties that the functions must have are tested; all tests are passed.
            
## Example of execution and results
            
The current configuration and simulation files are used to study a 30 x 30 lattice, with spins randomly chosen, in 200 temperature points between 1 and 4; for each one, the lattice is updated 1000 times in order to go towards equilibrium, and then energy and magnetization are averaged over 1000 more steps. 
After a simulation time of about one hour, these are the resulting plots: 

![temperature_plot](https://user-images.githubusercontent.com/79457897/147859840-90a854ce-d55f-4612-9bdf-c5a397027e3b.png)

![steps_plot](https://user-images.githubusercontent.com/79457897/147859881-54921445-b503-475b-82bf-2c29fa8a61d5.png)

![evolution_plot](https://user-images.githubusercontent.com/79457897/147859884-433225d7-beec-4d92-bb14-d4a6cdcc8095.png)

It can be seen that:

- energy at low temperature shows that the system is in a more stable configuration than at high T;
- similarly, magnetization shows a phase transition between an ordered low-temperature state and a disordered high-temperature one with zero mean magnetization;
- the system goes to equilibrium for energy and magnetization vs number of steps; this plot can be used to optimize the number of steps in the simulation;
- from a random configuration, spin domain form; the times at which the lattice is shown may be varied to study the complete evolution.

The last two plots are for low temperature, more precisely for T = 1.
If the simulation were done changing only the parameter n_show from 0 to 199, the last two plots would be for T = 4:

![steps_plot](https://user-images.githubusercontent.com/79457897/147860011-92dd8731-33e8-487b-a9bf-633823616ad2.png)

![evolution_plot](https://user-images.githubusercontent.com/79457897/147860015-bb7b9062-6678-467b-a78f-dd469bed1fb8.png)

It can be seen now that the higher thermal energy makes energy and magnetization have big fluctuations around the mean values from their the plot vs temperature.
            
Note that intensive energy and magnetization are plotted vs temperature, while total ones are plotted vs number of steps; however, it can be seen that the values are compatible since at each step number the intensive one is the total one divided by the number of sites.


