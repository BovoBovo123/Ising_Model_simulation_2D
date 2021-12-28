# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 14:46:41 2021

@author: pietr
"""

import numpy as np
import matplotlib.pyplot as plt
import logging 
import configparser


def initialize_state(N, M, choice = False, spin_up_pol = 0.5, seed = 42, level = 20):
    """

    Parameters
    ----------
    N : TYPE
        DESCRIPTION.
    M : TYPE
        DESCRIPTION.
    choice : TYPE, optional
        DESCRIPTION. The default is False.
    spin_up_pol : TYPE, optional
        DESCRIPTION. The default is 0.5.
    seed : TYPE, optional
        DESCRIPTION. The default is 42.
    level : TYPE, optional
        DESCRIPTION. The default is 20.

    Returns
    -------
    None.

    """
    
    if level not in [0, 10, 20, 30, 40, 50]:
        raise ValueError('Logging level is expected to be 0, 10, 20, 30, 40 or 50 but is {0}\n'.format(level))
    
    logger = logging.getLogger(__name__)
    logging.basicConfig(level = level)
    
    np.random.seed(seed)
    
    if N < 1 or M < 1:
       raise ValueError('Both lattice dimensions must be >= 1, but are {0} and {1}\n'.format(N, M))
    
    size = (N, M)
    
    if isinstance(choice, bool) is False:
        raise TypeError('Expected a boolean to choose the spin up polarization percentage, but got {0}\n'.format(choice))
    
    elif choice is False:
        if spin_up_pol != 0.5:
         logger.warning('Spin up polarization percentage was set to {0}, but the choice switch is False, so the default value of 0.5 is being used\n'.format(spin_up_pol))

        #Generate the spin lattice randomly
        spin = [-1., 1.]
        initial_state = np.random.choice(spin, size)
        
    elif choice is True:
        if not 0 <= spin_up_pol <= 1:
            raise ValueError('Expected the percentage of spin up polarization (expressed between 0 and 1), but got {0}\n'.format(spin_up_pol))
        
        if spin_up_pol == 0.5:
         logger.info('The choice of setting the spin up polarization percentage was made, but its value has not been changed from the default 0.5\n')

        #Generate initial lattice with input spin up percentage polarization
        initial_random = np.random.random(size)
        initial_state = np.zeros(size)
        initial_state[initial_random >= spin_up_pol] = -1
        initial_state[initial_random < spin_up_pol] = 1
    
    return initial_state


def metropolis_move(lattice, beta):
    """

    Parameters
    ----------
    lattice : TYPE
        DESCRIPTION.
    beta : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    
    #Length and width for looping
    length = len(lattice)
    width = len(lattice[0])
    
    for i in range(length):
        for j in range(width):
            #Take a random lattice point 
            x = np.random.randint(0, length)
            y = np.random.randint(0, width)
            site_spin = lattice[x, y]
            
            #Nearest neighbours total spin, considering PBC
            neighbour_spin = lattice[(x+1)%length, y] + lattice[x, (y+1)%width] + lattice[(x-1)%length, y] + lattice[x, (y-1)%width]
            
            #Energy change due to spin flip
            energy_change = 2*site_spin*neighbour_spin
            
            #If the energy change is negative, accept the move and flip the spin, otherwise accept the move with probability exp(-cost*beta), and flip the spin. 
            if energy_change < 0:
                site_spin *= -1
            elif np.random.random() < np.exp(-energy_change*beta):
                site_spin *= -1
            
            #Update lattice with new spin state
            lattice[x, y] = site_spin
            
    return lattice


def calculate_energy(lattice):
    """

    Parameters
    ----------
    lattice : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    
    total_energy = 0.0
    
    #Length and width for looping
    length = len(lattice)
    width = len(lattice[0])
    
    #The energy is given by the product of the (x,y) spin and the 4 nearest neighbours spins.
    for x in range(length):
        for y in range(width):
            site_spin = lattice[x,y]
            neighbour_spin = lattice[(x+1)%length, y] + lattice[x, (y+1)%width] + lattice[(x-1)%length, y] + lattice[x, (y-1)%width]
            total_energy += -site_spin*neighbour_spin
    
    #Total energy with no double counting
    return total_energy/2


def calculate_magnetization(lattice):
    """

    Parameters
    ----------
    lattice : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    
    #Since spins are all +1 or -1
    total_magnetization = np.sum(lattice)
    
    return total_magnetization


def read_configuration(filename):
    """

    Returns
    -------
    None.

    """
    
    configuration = configparser.ConfigParser()
    configuration.read(filename)
    
    #Configuration file must be read correctly
    if len(configuration.sections()) == 0:
        raise FileNotFoundError('Cannot find a file named "{0}" with at least one section\n'.format(filename))
    
    return configuration


def save_steps_data(ene, mag, first_save, ene_path = 'ene_steps.txt', mag_path = 'mag_steps.txt', level = 20):
    #For energy
    if first_save[0] == False:
        try:
            with open(ene_path, 'a') as f:
                write_data = f.write('{0}\n'.format(ene))
        except PermissionError:
            logger = logging.getLogger(__name__)
            logging.basicConfig(level = level)
            logger.error('It appears you do not have the permission to create or open the file; if you want to save the data, try to create an empty file with the name of the save path, and ignore the logging\n')
            raise PermissionError('It appears you do not have the permission to create or open the file; if you want to save the data, try to create an empty file with the name of the save path, and ignore the logging\n')

    #First time write start of file and check if it already exists
    elif first_save[0] == True:
        try:
            with open(ene_path, 'x') as f:
                f.write('#Saving energy vs steps\n')
        except FileExistsError:
            logger = logging.getLogger(__name__)
            logging.basicConfig(level = level)
            logger.info('Save file for energy vs steps already exists; data will be appendend but do not make confusion with precedent one\n')
        except PermissionError:
            logger = logging.getLogger(__name__)
            logging.basicConfig(level = level)
            logger.error('It appears you do not have the permission to create or open the file; if you want to save the data, try to create an empty file with the name of the save path, and ignore the logging\n')
            raise PermissionError('It appears you do not have the permission to create or open the file; if you want to save the data, try to create an empty file with the name of the save path, and ignore the logging\n')
            
    else:
        raise TypeError('Expected True or False to check for the first call, but got {0}\n'.format(first_save[0]))

    #For magnetization
    if first_save[0] == False:
        try:
            with open(mag_path, 'a') as f:
                write_data = f.write('{0}\n'.format(mag))
        except PermissionError:
            logger = logging.getLogger(__name__)
            logging.basicConfig(level = level)
            logger.error('It appears you do not have the permission to create or open the file; if you want to save the data, try to create an empty file with the name of the save path, and ignore the logging\n')
            raise

    #First time write start of file and check if it already exists
    elif first_save[0] == True:
        try:
            with open(mag_path, 'x') as f:
                f.write('#Saving magnetization vs steps\n')
        except FileExistsError:
            logger = logging.getLogger(__name__)
            logging.basicConfig(level = level)
            logger.info('Save file for magnetization vs steps already exists; data will be appendend but do not make confusion with precedent one\n')
        except PermissionError:
            logger = logging.getLogger(__name__)
            logging.basicConfig(level = level)
            logger.error('It appears you do not have the permission to create or open the file; if you want to save the data, try to create an empty file with the name of the save path, and ignore the logging\n')
            raise PermissionError('It appears you do not have the permission to create or open the file; if you want to save the data, try to create an empty file with the name of the save path, and ignore the logging\n')
            
    else:
        raise TypeError('Expected True or False to check for the first call, but got {0}\n'.format(first_save[0]))
    
    first_save[0] = False


def save_temp_data(ene, mag, first_save, ene_path = 'ene_temp.txt', mag_path = 'mag_temp.txt', level = 20):
    #For energy
    if first_save[0] == False:
        try:
            with open(ene_path, 'a') as f:
                write_data = f.write('{0}\n'.format(ene))
        except PermissionError:
            logger = logging.getLogger(__name__)
            logging.basicConfig(level = level)
            logger.error('It appears you do not have the permission to create or open the file; if you want to save the data, try to create an empty file with the name of the save path, and ignore the logging\n')
            raise PermissionError('It appears you do not have the permission to create or open the file; if you want to save the data, try to create an empty file with the name of the save path, and ignore the logging\n')

    #First time write start of file and check if it already exists
    elif first_save[0] == True:
        try:
            with open(ene_path, 'x') as f:
                f.write('#Saving energy vs temperature\n')
        except FileExistsError:
            logger = logging.getLogger(__name__)
            logging.basicConfig(level = level)
            logger.info('Save file for energy vs temperature already exists; data will be appendend but do not make confusion with precedent one\n')
        except PermissionError:
            logger = logging.getLogger(__name__)
            logging.basicConfig(level = level)
            logger.error('It appears you do not have the permission to create or open the file; if you want to save the data, try to create an empty file with the name of the save path, and ignore the logging\n')
            raise PermissionError('It appears you do not have the permission to create or open the file; if you want to save the data, try to create an empty file with the name of the save path, and ignore the logging\n')
    
    else:
        raise TypeError('Expected True or False to check for the first call, but got {0}\n'.format(first_save[0]))

    #For magnetization
    if first_save[0] == False:
        try: 
            with open(mag_path, 'a') as f:
                write_data = f.write('{0}\n'.format(mag))
        except PermissionError:
            logger = logging.getLogger(__name__)
            logging.basicConfig(level = level)
            logger.error('It appears you do not have the permission to create or open the file; if you want to save the data, try to create an empty file with the name of the save path, and ignore the logging\n')
            raise PermissionError('It appears you do not have the permission to create or open the file; if you want to save the data, try to create an empty file with the name of the save path, and ignore the logging\n')

    #First time write start of file and check if it already exists
    elif first_save[0] == True:
        try:
            with open(mag_path, 'x') as f:
                f.write('#Saving magnetization vs temperature\n')
        except FileExistsError:
            logger = logging.getLogger(__name__)
            logging.basicConfig(level = level)
            logger.info('Save file for magnetization vs temperature already exists; data will be appendend but do not make confusion with precedent one\n')
        except PermissionError:
            logger = logging.getLogger(__name__)
            logging.basicConfig(level = level)
            logger.error('It appears you do not have the permission to create or open the file; if you want to save the data, try to create an empty file with the name of the save path, and ignore the logging\n')
            raise PermissionError('It appears you do not have the permission to create or open the file; if you want to save the data, try to create an empty file with the name of the save path, and ignore the logging\n')

    else:
        raise TypeError('Expected True or False to check for the first call, but got {0}\n'.format(first_save[0]))
    
    first_save[0] = False


def plots_T(T, energy, magnetization, saving = True, path = 'temperature_plot.png'):
    #Size set to fill well
    f = plt.figure(figsize=(24, 10));  

    #Energy plot vs T
    sub_f =  f.add_subplot(1, 2, 1);
    plt.scatter(T, energy, s = 50, marker = 'o', color = 'IndianRed')
    plt.xlabel("Temperature", fontsize  =22)
    plt.ylabel("Energy", fontsize = 22)         
    
    #Magnetization plot vs T
    sub_f =  f.add_subplot(1, 2, 2);
    plt.scatter(T, abs(magnetization), s = 50, marker = 'o', color = 'RoyalBlue')
    plt.xlabel("Temperature ", fontsize = 22)
    plt.ylabel("Magnetization ", fontsize = 22)   
    
    #Saving
    if saving == True:
        f.savefig(path)
    
    
def plots_steps(x_step, y_ene, y_mag, n_show, numb_T, saving = True, path = 'steps_plot.png'):
    if not 0 <= n_show <= numb_T - 1:
        raise ValueError('Must choose the index of the temperature to plot quantities vs steps between 0 and the number of temperature points - 1 (both included); got 0 <= {0} <= {1}\n'.format(n_show, numb_T))                     
    
    #Size set to fill well
    f = plt.figure(figsize=(24, 10));  

    #Energy plot vs steps
    sub_f =  f.add_subplot(1, 2, 1);
    plt.scatter(x_step, y_ene, s = 20, marker = 'o', color = 'IndianRed')
    plt.xlabel("Steps", fontsize=22)
    plt.ylabel("Energy ", fontsize=22)       

    #Magnetization plot vs steps
    sub_f =  f.add_subplot(1, 2, 2);
    plt.scatter(x_step, y_mag, s = 20, marker = 'o', color = 'RoyalBlue')
    plt.xlabel("Steps", fontsize=22)
    plt.ylabel("Magnetization ", fontsize=22) 
    
    #Saving
    if saving == True:
        f.savefig(path)
      

def simulate(lattice, beta, times = (5, 10, 50, 100, 1000), level = 20):
    """
    
    Parameters
    ----------
    lattice : TYPE
        DESCRIPTION.
    beta : TYPE
        DESCRIPTION.
    times : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    if len(times) != 5:
        raise ValueError('Must insert 5 value of times at which the lattice is shown; {0} were passed\n'.format(len(times)))
    
    for time in times:
        if time < 0:
            raise ValueError('Must insert non-negative values of times at which the lattice is shown\n')
        if time == 0:
            logger = logging.getLogger(__name__)
            logging.basicConfig(level = level)
            logger.info('THe initial lattice i.e. at time = 0 il always shown, so there will be a repetition in the evolution plot\n')

    initial_state = lattice.copy()
    evolution_steps = max(times)+ 1
    states_evolution = [initial_state]
    
    #Take data from selected points in evolution time
    for time in range(evolution_steps):
        evolved_state = metropolis_move(lattice, beta)
        if time in times:
            added_state = evolved_state.copy()
            states_evolution.append(added_state)
    
    return states_evolution


def plot_evolution(evolution_states, N, M, times = (5, 10, 50, 100, 1000), level = 20, saving = True, path = 'evolution_plot.png'):
    f = plt.figure(figsize=(30, 18), dpi=80)
    
    #Create a meshgrid; special case for N = M = 1
    if N*M == 1:
        x, y = np.meshgrid(range(2), range(2))
        logger = logging.getLogger(__name__)
        logging.basicConfig(level = level)
        logger.info('Time evolution cannot be displayed correctly for a single point lattice (i.e. M = N = 1)\n')
    else:
        x, y = np.meshgrid(range(N), range(M))
    
    #Add subplots at the different times
    for t in range(len(times)+1):
        sub_f = f.add_subplot(2, 3, t+1)
        plt.pcolormesh(x, y, evolution_states[t], cmap = plt.cm.RdBu, shading = 'nearest')
        plt.axis('off')
        
        if t != 0:
            plt.title('Time = 0', fontsize = 28)
            
        else:
            plt.title('Time = {0}'.format(times[t-1]), fontsize = 28)

    #Saving
    if saving == True:
        f.savefig(path)






