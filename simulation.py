# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 14:55:28 2021

@author: pietr
"""
import functions_ising as fi
import numpy as np
import matplotlib.pyplot as plt

import configparser

from tqdm import tqdm
from tqdm import trange

#Import configuration
filename = 'CONFIGURATION.txt'
configuration = fi.read_configuration(filename)

N = configuration.getint('SETTINGS', 'N')
M = configuration.getint('SETTINGS', 'M')

eq_steps = configuration.getint('SETTINGS', 'eq_steps')
mc_steps = configuration.getint('SETTINGS', 'mc_steps')

numb_T = configuration.getint('SETTINGS', 'numb_T')
T_init = configuration.getfloat('SETTINGS', 'T_init')
T_final = configuration.getfloat('SETTINGS', 'T_final')

choice = configuration.getboolean('SETTINGS', 'choice')
spin_up_pol = configuration.getfloat('SETTINGS', 'spin_up_pol')

level = configuration.getint('LOGGING', 'level')

seed = configuration.getint('SETTINGS', 'seed')

nT_show = configuration.getint('PLOTTING', 'nT_show')

t1 = configuration.getint('PLOTTING', 't1')
t2 = configuration.getint('PLOTTING', 't2')
t3 = configuration.getint('PLOTTING', 't3')
t4 = configuration.getint('PLOTTING', 't4')
t5 = configuration.getint('PLOTTING', 't5')
t6 = configuration.getint('PLOTTING', 't6')
times = (t1, t2, t3, t4, t5, t6)

#Get normalization factor by dividing by number of steps and system size to get intensive values
norm_intensive = 1.0/(mc_steps*N*M)

T = np.linspace(T_init, T_final, numb_T)
energy, magnetization = np.zeros(numb_T), np.zeros(numb_T)

T_show = T[nT_show]
beta_show = 1.0/T_show

x_step = []
y_ene = []
y_mag = []
    
for n_temp in trange(numb_T, desc = 'Loop over temperature values', position = 0):
    config = fi.initialize_state(N, M, choice, spin_up_pol, seed, level)        
    
    ene_count = mag_count = 0.0
    
    #Beta values, with Boltzmann constant k = 1
    beta = 1.0/T[n_temp]
    
    #Equilibrate the system
    for i in range(eq_steps):         
        fi.metropolis_move(config, beta)  
        
        #Data for plots vs steps
        if n_temp == nT_show:
            a = fi.calculate_energy(config)
            b = fi.calculate_magnetization(config)
            y_ene.append(a)
            y_mag.append(b)
            x_step.append(i)

    #Acquire energy and magnetization measurments
    for i in range(mc_steps):
        fi.metropolis_move(config, beta)          
        ene_step = fi.calculate_energy(config)     
        mag_step = fi.calculate_magnetization(config) 
        
        #Data for plots vs steps
        if n_temp == nT_show:
            c = fi.calculate_energy(config)
            d = fi.calculate_magnetization(config)
            y_ene.append(c)
            y_mag.append(d)
            x_step.append(i + eq_steps)

        ene_count += ene_step
        mag_count += mag_step

    #Divide by number of steps and system size to get intensive values    
    energy[n_temp] = norm_intensive*ene_count
    magnetization[n_temp] = norm_intensive*mag_count

#Plotting
fi.plots_T(T, energy, magnetization)
fi.plots_steps(x_step, y_ene, y_mag, nT_show, numb_T)

#Showing lattice evolution
initial_state = fi.initialize_state(N, M, choice, spin_up_pol, seed, level)  
evolution_states = fi.simulate(initial_state, beta_show, times)
fi.plot_evolution(evolution_states, times, N, M, level)


