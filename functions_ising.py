# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 14:46:41 2021

@author: pietr
"""

import numpy as np
import matplotlib.pyplot as plt
import logging 

import configparser

logger = logging.getLogger(__name__)

def initialize_state(N, M, choice = False, spin_up_pol = 0.5, seed = 42):
    
    """DOCSTRING?
    a
    """
    
    np.random.seed(seed)
    
    if N < 1 or M < 1:
       raise ValueError('Both lattice dimensions must be > 1, but are {0} and {1}'.format(N,M))
    
    size = (N, M)
    
    if isinstance(choice, bool) == False:
        raise TypeError('Expected a boolean to choose the spin up polarization percentage, but got {0}'.format(choice))
    
    elif choice == False:
        if spin_up_pol != 0.5:
         logger.warning('''Spin up polarization percentage was set to {0}, but the choice switch is False, so the standard value of 0.5 is being used'''.format(spin_up_pol))

        #Generate the spin lattice randomly
        spin = [-1, 1]
        initial_state = np.random.choice(spin, size)
        
    elif choice == True:
        if not 0 <= spin_up_pol <= 1:
            raise ValueError('Expected the percentage of spin up polarization (expressed between 0 and 1), but got {0}'.format(spin_up_pol))
        
        if spin_up_pol == 0.5:
         logger.warning('''The choice of setting the spin up polarization percentage was made, but its value has not been changed from the standard 0.5''')

        #Generate initial lattice with input spin up percentage polarization
        initial_random = np.random.random(size)
        initial_state = np.zeros(size)
        initial_state[initial_random >= spin_up_pol] = -1
        initial_state[initial_random < spin_up_pol] = 1
    
    return initial_state


def metropolis_move(lattice, beta):
    """DOCSTRING?
    ciao
    usare beta o T?
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
            neighbour_spin = lattice[(x+1)%length, y] + lattice[x, (y+1)%width] 
            + lattice[(x-1)%length, y] + lattice[x, (y-1)%width]
            
            #Energy change due to flip
            energy_change = 2*site_spin*neighbour_spin
            
            #If the energy change is negative, accept the move and flip the spin, 
            #otherwise accept the move with probability exp(-cost*beta), and flip the spin. !?
            if energy_change < 0:
                site_spin *= -1
            elif np.random.random() < np.exp(-energy_change*beta):
                site_spin *= -1
            
            #Update lattice with new spin state
            lattice[x, y] = site_spin
            
    return lattice


def calculate_energy(lattice):
    """DOCSTRING?
    ciao2
    """
    
    total_energy = 0.0
    
    #Length and width for looping
    length = len(lattice)
    width = len(lattice[0])
    
    #The change in energy is given by the product of the (x,y) spin 
    #and the 4 nearest neighbours spins.
    for x in range(length):
        for y in range(width):
            site_spin = lattice[x,y]
            neighbour_spin = lattice[(x+1)%length, y] + lattice[x, (y+1)%width] 
            + lattice[(x-1)%length, y] + lattice[x, (y-1)%width]
            total_energy += - site_spin*neighbour_spin
    
    #Single site energy, considering the presence of four neighbours
    return total_energy/4


def calculate_magnetization(lattice):
    """DOCSTRING?
    lel
    """
    
    #Since spins are all +1 or -1
    total_magnetization = np.sum(lattice)
    
    return total_magnetization


def simulate(lattice, beta, times):
    """DOCSTRING
    serve beta o T? o nulla?
    """
    
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


def evolution_plot():
    """
    DOCSTRING?
    """
    
    figure = plt.figure(figsize=(30, 15), dpi=80)





