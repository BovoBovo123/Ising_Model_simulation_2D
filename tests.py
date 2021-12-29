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
    
def test_lattice_dimensions_polarized(N = 3, M = 2, choice = True, spin_up_pol = 0.3):
    lattice = fi.initialize_state(N, M)
    assert len(lattice) == 3
    assert len(lattice[0]) == 2
    
def test_spin_values(N = 3, M = 2):
    lattice = fi.initialize_state(N, M)
    abs_spin = np.abs(lattice)
    assert abs_spin.all() == 1
    
def test_spin_values_polarized(N = 3, M = 2, choice = True, spin_up_pol = 0.3):
    lattice = fi.initialize_state(N, M)
    abs_spin = np.abs(lattice)
    assert abs_spin.all() == 1

def test_raises_error_lattice_dimensions(N = -1, M = -2):
    with pytest.raises(ValueError):
        lattice = fi.initialize_state(N, M)

def test_choice_not_boolean(N = 2, M = 3, choice = 'not a boolean'):
    with pytest.raises(TypeError):
        lattice = fi.initialize_state(N, M, choice) 

def test_spin_up_polarization_not_a_percentage(N = 2, M = 3, choice = True, spin_up_pol = 1.2):
    with pytest.raises(ValueError):
        lattice = fi.initialize_state(N, M, choice, spin_up_pol)   

#Test the function that updates the lattice
def test_metropolis_move(N = 2, M = 3, numb_T = 50, T_init = 1, T_final = 3):
    #Test that dimensions, spins, energy and magnetization are as expected 
    lattice = fi.initialize_state(N, M)
    T = np.linspace(T_init, T_final, numb_T)
    for temp in range(numb_T):
        beta = 1.0/T[temp]
        lattice = fi.metropolis_move(lattice, beta)
        abs_spin = np.abs(lattice)
        assert abs_spin.all() == 1 
        length = len(lattice)
        width = len(lattice[0])
        assert length == N
        assert width == M
        for x in range(length):
            for y in range(width):
                site_spin = lattice[x, y]
                assert abs(site_spin) == 1
                neighbour_spin = lattice[(x+1)%length, y] + lattice[x, (y+1)%width] + lattice[(x-1)%length, y] + lattice[x, (y-1)%width]
                assert -4 <= neighbour_spin <= 4
                energy_change = 2*site_spin*neighbour_spin
                assert -8 <= energy_change <= 8

#Test functions for energy and magnetization calculation
def test_correct_total_energy(N = 2, M = 3, choice = True, spin_up_pol = 1):
    lattice = fi.initialize_state(N, M, choice, spin_up_pol)
    assert fi.calculate_energy(lattice) == -N*M*4/2  
    
def test_correct_total_magnetization(N = 2, M = 3, choice = True, spin_up_pol = 1):
    lattice = fi.initialize_state(N, M, choice, spin_up_pol)
    assert fi.calculate_magnetization(lattice) == N*M

#Test the function that reads the configuration parameters
def test_read_configuration(filename = ''):
    with pytest.raises(FileNotFoundError):
        config = fi.read_configuration(filename)
        
#Test the function that simulates lattice evolution
def test_simulate_wrong_length():
    lattice = fi.initialize_state(N = 2, M = 3)
    beta = 1
    times = [1, 2, 3, 4]
    with pytest.raises(ValueError):
        states_evolution = fi.simulate(lattice, beta, times)

def test_simulate_correct_length():
    lattice = fi.initialize_state(N = 2, M = 3)
    beta = 1
    times = [1, 2, 3, 4, 5]
    six_states = fi.simulate(lattice, beta, times)
    assert len(six_states) == 6
    
def test_simulate_negative_times():
    lattice = fi.initialize_state(N = 2, M = 3)
    beta = 1
    times = [1, 2, 3, -4, -5]
    with pytest.raises(ValueError):
        states_evolution = fi.simulate(lattice, beta, times)
        
def test_simulate(N = 2, M = 3, beta = 1):
    #Test that dimensions, spins, energy and magnetization are as expected 
    lattice = fi.initialize_state(N, M)
    evolved_states = fi.simulate(lattice, beta)
    for state in evolved_states:
        abs_spin = np.abs(state)
        assert abs_spin.all() == 1 
        length = len(state)
        width = len(state[0])
        assert length == N
        assert width == M
        for x in range(length):
            for y in range(width):
               site_spin = lattice[x, y]
               assert abs(site_spin) == 1
               neighbour_spin = lattice[(x+1)%length, y] + lattice[x, (y+1)%width] + lattice[(x-1)%length, y] + lattice[x, (y-1)%width]
               assert -4 <= neighbour_spin <= 4
               energy_change = 2*site_spin*neighbour_spin
               assert -8 <= energy_change <= 8












