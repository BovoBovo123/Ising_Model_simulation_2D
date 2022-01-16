# -*- coding: utf-8 -*-
"""
Created on Fri Dec 17 14:59:01 2021

@author: pietr
"""


import functions_ising as fi
import numpy as np
import pytest


#Test the lattice initialization function
def test_lattice_dimensions(N = 3, M = 2):
    lattice = fi.initialize_state(N, M)
    assert len(lattice) == 3
    assert len(lattice[0]) == 2
    
def test_spin_values(N = 3, M = 2):
    lattice = fi.initialize_state(N, M)
    abs_spin = np.abs(lattice)
    assert abs_spin.all() == 1
    
def test_spins_all_up(N = 3, M = 2, spin_up_pol = 1):
    lattice = fi.initialize_state(N, M, spin_up_pol)
    for i in range(N):
        for j in range(M):
            assert lattice[i][j] == 1
    
def test_spins_all_down(N = 3, M = 2, spin_up_pol = 0):
    lattice = fi.initialize_state(N, M, spin_up_pol)
    length = len(lattice)
    width = len(lattice[0])
    for i in range(length):
        for j in range(width):
            assert lattice[i][j] == -1    

def test_raises_error_lattice_dimensions(N = -1, M = -2):
    with pytest.raises(ValueError):
        lattice = fi.initialize_state(N, M)

#Test the function that updates the lattice
def test_evolution_shape(N = 2, M = 3, beta = 1.0):
    lattice = fi.initialize_state(N, M)
    lattice = fi.metropolis_move(lattice, beta)
    length = len(lattice)
    width = len(lattice[0])
    assert length == N
    assert width == M

def test_evolution_spins(N = 2, M = 3, beta = 1.0):
    lattice = fi.initialize_state(N, M)
    lattice = fi.metropolis_move(lattice, beta)
    abs_spin = np.abs(lattice)
    assert abs_spin.all() == 1 
    
def test_evolution_low_T(N = 2, M = 3, spin_up_pol = 1, beta = np.inf):
    lattice = fi.initialize_state(N, M, spin_up_pol)
    lattice = fi.metropolis_move(lattice, beta)
    length = len(lattice)
    width = len(lattice[0])
    for i in range(length):
        for j in range(width):
            assert lattice[i][j] == 1    
    
def test_evolution_high_T(N = 2, M = 3, spin_up_pol = 1, seed = 1, beta = 0.0):
    lattice = fi.initialize_state(N, M, spin_up_pol, seed)
    evolved_lattice = fi.metropolis_move(lattice, beta)
    simulated_lattice = ([1, -1, -1], [1, -1, -1])
    assert np.array_equal(evolved_lattice, simulated_lattice) == True

#Test the function that reads the configuration parameters
def test_read_configuration(filename = ''):
    with pytest.raises(FileNotFoundError):
        config = fi.read_configuration(filename)

#Test the function that simulated lattice evolution at certain time instants
def test_simulate_length(N = 2, M = 3, beta = 1.0, times = [1, 2, 3, 4, 5]):
    lattice = fi.initialize_state(N, M)
    evolved_states = fi.simulate(lattice, beta, times)
    assert len(evolved_states) == 6
    
def test_simulate_time_independent(N = 2, M = 3, seed = 1, beta = 1.0, times1 = [2, 4, 6, 8, 10], times2 = [2, 5, 7, 9, 10]):    
    lattice1 = fi.initialize_state(N, M, seed = seed)
    evolved_states1 = fi.simulate(lattice1, beta, times1)
    lattice2 = fi.initialize_state(N, M, seed = seed)
    evolved_states2 = fi.simulate(lattice2, beta, times2)
    assert np.array_equal(evolved_states1[1], evolved_states2[1]) == True
    assert np.array_equal(evolved_states1[5], evolved_states2[5]) == True













