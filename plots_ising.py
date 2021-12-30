#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 29 16:15:49 2021

@author: bovo123
"""


import numpy as np
import matplotlib.pyplot as plt
import logging 


def plots_T(T, energy, magnetization, saving = True, save_path = 'temperature_plot.png', load = False, load_path = ('ene_temp_path', 'mag_temp_path')):
    """
    This function plots energy and magnetization vs temperature, with data that is
    either given or loaded, and can save it 

    Parameters
    ----------
    T : 1D-like array
        temperature points.
    energy : 1D-like array
        energy points.
    magnetization : 1D-like array
        magnetization points.
    saving : bool, optional
        if True, the plot is saved. The default is True.
    save_path : string, optional
        path to the save file. The default is 'temperature_plot.png'.
    load : bool, optional
        if True, data is loaded. The default is False.
    load_path : 1D-like array, optional
        list of strings of two files from which to load data; the first should be
        for the energy. The default is ('ene_temp_path', 'mag_temp_path').

    Returns
    -------
        None.
    
    Raises
    ------
        TypeError if the load switch is not a boolean

    """
    
    #Data either given or loaded
    if load == True:
        if len(energy)*len(magnetization)!= 0:
            logging.warning('Non-empty arrays were given, but data is being loaded so they will be over-written\n')
        energy = np.read(load_path[0])
        magnetization = np.read(load_path[1])
    
    elif load == False:
        pass
    
    else:
        raise TypeError('Was expecting a boolean to decide if to load the data, but got {0}'.format(load))
    
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
        f.savefig(save_path)
    
    
def plots_steps(x_step, y_ene, y_mag, n_show, numb_T, saving = True, save_path = 'steps_plot.png', load = False, load_path = ('ene_steps_path', 'mag_steps_path')):
    """
    This function plots energy and magnetization vs temperature, with data that is
    either given or loaded, and can save it 

    Parameters
    ----------
    x_step : 1D-like array
        number of step points.
    y_ene : 1D-like array
        energy points.
    y_mag : 1D-like array
        magnetization points.
    n_show : int
        index of temperature for which plots are shown.
    numb_T : int
        number of temperature points simulated.
    saving : bool, optional
        if True, the plot is saved. The default is True.
    save_path : string, optional
        path to save file. The default is 'steps_plot.png'.
    load : bool, optional
        if True, data is loaded. The default is False.
    load_path : 1D-like array, optional
        list of strings of two files from which to load data; the first should be
        for the energy. The default is ('ene_steps_path', 'mag_steps_path').

    Returns
    -------
        None.
    
    Raises
    ------
        TypeError if the load switch is not a boolean
        ValueError if the temperature index is out of range

    """
    
    #Data either given or loaded
    if load == True:
        if len(y_ene)*len(y_mag)!= 0:
            logging.warning('Non-empty arrays were given, but data is being loaded so they will be over-written\n')
        y_ene = np.read(load_path[0])
        y_mag = np.read(load_path[1])
    
    elif load == False:
        pass
    
    else:
        raise TypeError('Was expecting a boolean to decide if to load the data, but got {0}'.format(load))

    #Correct choice of temperature to show plots    
    if not 0 <= n_show <= numb_T - 1:
        raise ValueError('Must choose the index of the temperature to plot quantities vs steps between 0 and the number of temperature points - 1 (both included); got 0 <= {0} <= {1}\n'.format(n_show, numb_T))                     
    
    #Size set to fill well
    f = plt.figure(figsize=(24, 10));  

    #Energy plot vs steps
    sub_f =  f.add_subplot(1, 2, 1);
    plt.scatter(x_step, y_ene, s = 50, marker = 'o', color = 'IndianRed')
    plt.xlabel("Steps", fontsize=22)
    plt.ylabel("Energy ", fontsize=22)       

    #Magnetization plot vs steps
    sub_f =  f.add_subplot(1, 2, 2);
    plt.scatter(x_step, y_mag, s = 50, marker = 'o', color = 'RoyalBlue')
    plt.xlabel("Steps", fontsize=22)
    plt.ylabel("Magnetization ", fontsize=22) 
    
    #Saving
    if saving == True:
        f.savefig(save_path)


def plot_evolution(evolution_states, N, M, times = (5, 10, 50, 100, 1000), saving = True, save_path = 'evolution_plot.png'):
    """
    This function plots the initial state of the lattice as well as its evolution
    at five time instants.

    Parameters
    ----------
    evolution_states : TYPE
        DESCRIPTION.
    N : int
        length of the lattice.
    M : int
        width of the lattice.
    times : 1D-like array, optional
        five time instants when to store the evolved lattice spin configuration. 
        The default is (5, 10, 50, 100, 1000).
    saving : bool, optional
        if True, the plot is saved. The default is True.
    save_path : string, optional
        path to save file. The default is 'evolution_plot.png'.

    Returns
    -------
        None.

    """
    
    f = plt.figure(figsize=(30, 18), dpi=80)
    
    #Create a meshgrid; special case for N = M = 1
    if N*M == 1:
        x, y = np.meshgrid(range(2), range(2))
        logging.info('Time evolution cannot be displayed correctly for a single point lattice (i.e. M = N = 1)\n')
    else:
        x, y = np.meshgrid(range(N), range(M))
    
    #Add subplots at the different times
    for t in range(len(times)+1):
        sub_f = f.add_subplot(2, 3, t+1)
        plt.pcolormesh(x, y, evolution_states[t], cmap = plt.cm.RdBu, shading = 'nearest')
        plt.axis('off')
        
        if t != 0:
            plt.title('Time = {0}'.format(times[t-1]), fontsize = 28)
            
        else:
            plt.title('Time = 0', fontsize = 28)

    #Saving
    if saving == True:
        f.savefig(save_path)





