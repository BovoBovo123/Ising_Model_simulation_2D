# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 14:46:41 2021

@author: pietr
"""


import numpy as np
import logging 
import configparser


def initialize_state(N, M, spin_up_pol = None, seed = 42):
    """
    This function generate the spin lattice randomly, with a certain mean spin 
    polarization if given
    
    Parameters
    ----------
    N : int
        length of the lattice.
    M : int
        width of the lattice.
    spin_up_pol : float, optional
        mean spin up polarization. The default is None, that will generate a random lattice.
    seed : int, optional
        sets the seed using np.random.seed(). The default is 42.

    Returns
    -------
        the N*M spin lattice initial configuration.
    
    Raises
    ------
        ValueError is lattice dimensions are < 1.
        
    """
    
    np.random.seed(seed)
    
    if N < 1 or M < 1:
       raise ValueError('Both lattice dimensions must be >= 1, but are {0} and {1}\n'.format(N, M))
    
    size = (N, M)
    
    if spin_up_pol == None:
        #Generate the spin lattice randomly
        spin = [-1., 1.]
        initial_state = np.random.choice(spin, size)
        
    else:
        if not 0 <= spin_up_pol <= 1:
            logging.warning('Expected the percentage of spin up polarization (expressed between 0 and 1), but got {0}; the lattice will be completely polarized\n'.format(spin_up_pol))
        
        #Generate initial lattice with input spin up percentage polarization
        initial_random = np.random.random(size)
        initial_state = np.zeros(size)
        initial_state[initial_random >= spin_up_pol] = -1
        initial_state[initial_random < spin_up_pol] = 1
    
    return initial_state


def metropolis_move(lattice, beta):
    """
    This functions uses the Metropolis algorithm to update the lattice spins

    Parameters
    ----------
    lattice : 2D-like array
        lattice spin configuration.
    beta : float
        1/kT where T is the temperature and the Boltzmann constant k 
        is taken equal to 1.

    Returns
    -------
        the updated lattice spin configuration.

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
    This functions calculates the lattice energy (with PBC) using the Ising 
    Hamiltonian with the exchange constant J equal to 1

    Parameters
    ----------
    lattice : 2D-like array
        lattice spin configuration.

    Returns
    -------
        the lattice energy, considering periodic boundary conditions.

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
    This functions calculates the lattice magnetization

    Parameters
    ----------
    lattice : 2D-like array
        lattice spin configuration.

    Returns
    -------
        the lattice magnetization.

    """
    
    #Since spins are all +1 or -1
    total_magnetization = np.sum(lattice)
    
    return total_magnetization


def read_configuration(filename):
    """
    This function reads a configuration file

    Parameters
    ----------
    filename : string
        path of the configuration file to be read.

    Returns
    -------
        a configparser that reads the configuration file.
    
    Raises
    ------
        FileNotFoundError if the file that is read has no sections

    """
    
    configuration = configparser.ConfigParser()
    configuration.read(filename)
    
    #Configuration file must be read correctly
    if len(configuration.sections()) == 0:
        raise FileNotFoundError('Cannot find a file named "{0}" with at least one section\n'.format(filename))
    
    return configuration


def save_steps_data(ene, mag, ene_path = 'ene_steps.txt', mag_path = 'mag_steps.txt'):
    """
    This functions saves energy and magnetization points vs step number in files 
    with given path; to be called in a loop while calculating those quantities

    Parameters
    ----------
    ene : float
        energy value at a certain step number.
    mag : float
        magnetization value at a certain step number.
    ene_path : string, optional
        path for the energy save file. The default is 'ene_steps.txt'.
    mag_path : string, optional
        path for the magnetization save file. The default is 'mag_steps.txt'.

    Returns
    -------
        None.
    
    Raises
    ------
        PermissionError if the files cannot be created due to lack of permission

    """
    
    #For energy
    try:
        with open(ene_path, 'a') as f:
            write_data = f.write('{0}\n'.format(ene))
    except PermissionError:
        logging.error('It appears you do not have the permission to create or open the file; if you want to save the data, try to create an empty file with the name of the save path\n')
        raise PermissionError('It appears you do not have the permission to create or open the file; if you want to save the data, try to create an empty file with the name of the save path\n')

    #For magnetization
    try:
        with open(mag_path, 'a') as f:
            write_data = f.write('{0}\n'.format(mag))
    except PermissionError:
        logging.error('It appears you do not have the permission to create or open the file; if you want to save the data, try to create an empty file with the name of the save path\n')
        raise PermissionError('It appears you do not have the permission to create or open the file; if you want to save the data, try to create an empty file with the name of the save path\n')
 

def save_temp_data(ene, mag, ene_path = 'ene_temp.txt', mag_path = 'mag_temp.txt'):
    """
    This functions saves energy and magnetization points vs temperature in files 
    with given path; to be called in a loop while calculating those quantities

    Parameters
    ----------
    ene : float
        energy value at a certain temperature.
    mag : float
        magnetization value at a certain temperature.
    ene_path : string, optional
        path for the energy save file. The default is 'ene_temp.txt'.
    mag_path : string, optional
        path for the magnetization save file. The default is 'mag_temp.txt'.

    Returns
    -------
        None.
    
    Raises
    ------
        PermissionError if the files cannot be created due to lack of permission

    """
    
    #For energy
    try:
        with open(ene_path, 'a') as f:
            write_data = f.write('{0}\n'.format(ene))
    except PermissionError:
        logging.error('It appears you do not have the permission to create or open the file; if you want to save the data, try to create an empty file with the name of the save path\n')
        raise PermissionError('It appears you do not have the permission to create or open the file; if you want to save the data, try to create an empty file with the name of the save path\n')

    #For magnetization
    try: 
        with open(mag_path, 'a') as f:
            write_data = f.write('{0}\n'.format(mag))
    except PermissionError:
        logging.error('It appears you do not have the permission to create or open the file; if you want to save the data, try to create an empty file with the name of the save path\n')
        raise PermissionError('It appears you do not have the permission to create or open the file; if you want to save the data, try to create an empty file with the name of the save path\n')
      

def simulate(lattice, beta, times = (5, 10, 50, 100, 1000)):
    """
    This function simulates the lattice evolution for a given 
    number of steps (i.e. time)

    Parameters
    ----------
    lattice : 2D-like array
        lattice spin configuration to be studied.
    beta : float
        1/kT where T is the temperature and the Boltzmann constant k 
        is taken equal to 1.
    times : 1D-like array, optional
        five time instants when to store the evolved lattice spin configuration. 
        The default is (5, 10, 50, 100, 1000).

    Returns
    -------
        array containing the initial lattice spin configuration and the evolved ones.
    
    Raises
    ------
        ValueError if the number of time instants if different than five, or if
        any time instant is negative, or if time instants are repeated

    """
    
    if len(times) != 5:
        raise ValueError('Must insert 5 value of times at which the lattice is shown; {0} were passed\n'.format(len(times)))
    
    for time in times:
        if time < 0:
            raise ValueError('Must insert non-negative values of times at which the lattice is shown\n')
        if time == 0:
            logging.info('The initial lattice i.e. at time = 0 will always be shown, so there will be a repetition in the evolution plot\n')
    
    time_set = set(times)
    if len(times) != len(time_set):
        raise ValueError('Must insert time instants that are all different, but got only {0} that are unique'.format(time_set))
    
    sorted_times = sorted(times)
    if not np.array_equal(times, sorted_times):
        logging.info('Time instants are not sorted, so the future evolution plot might look strange and/or unclear')
    
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








