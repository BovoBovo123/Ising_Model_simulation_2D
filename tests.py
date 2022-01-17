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
    """
    Test that the shape of the lattice is the one expected.

    """
    
    lattice = fi.initialize_state(N, M)
    assert len(lattice) == 3
    assert len(lattice[0]) == 2
    
    
def test_spin_values(N = 3, M = 2):
    """
    Test that all spins are either +1 or -1; 
    this will also fix minimum and maximum energy and magnetization.

    """
    
    lattice = fi.initialize_state(N, M)
    abs_spin = np.abs(lattice)
    assert abs_spin.all() == 1
    
    
def test_spins_all_up(N = 3, M = 2, spin_up_pol = 1):
    """
    Test that the lattice can be generated with all spins up.

    """
    
    lattice = fi.initialize_state(N, M, spin_up_pol)
    length = len(lattice)
    width = len(lattice[0])
    for i in range(length):
        for j in range(width):
            assert lattice[i][j] == 1
    
    
def test_spins_all_down(N = 3, M = 2, spin_up_pol = 0):
    """
    Test that the lattice can be generated with all spins down.

    """
    
    lattice = fi.initialize_state(N, M, spin_up_pol)
    length = len(lattice)
    width = len(lattice[0])
    for i in range(length):
        for j in range(width):
            assert lattice[i][j] == -1    


def test_raises_error_lattice_dimensions(N = -1, M = -2):
    """
    Test that an error is raised if the lattice dimensions are negative.

    """
    
    with pytest.raises(ValueError):
        lattice = fi.initialize_state(N, M)


#Test the function that updates the lattice
def test_evolution_shape(N = 2, M = 3, beta = 1.0):
    """
    Test that the lattice maintains the correct shape when evolved.

    """
    
    lattice = fi.initialize_state(N, M)
    lattice = fi.metropolis_move(lattice, beta)
    length = len(lattice)
    width = len(lattice[0])
    assert length == N
    assert width == M


def test_evolution_spins(N = 2, M = 3, beta = 1.0):
    """
    Test that the lattice spins remain either +1 or -1 when evolved.

    """
    
    lattice = fi.initialize_state(N, M)
    lattice = fi.metropolis_move(lattice, beta)
    abs_spin = np.abs(lattice)
    assert abs_spin.all() == 1 
    
    
def test_evolution_low_T(N = 2, M = 3, spin_up_pol = 1, beta = np.inf):
    """
    Test that at zero temperature a fully polarized lattice does not change;
    this is because there are no thermal fluctuations to provide the energy to move
    from the energy minimum.

    """
    
    lattice = fi.initialize_state(N, M, spin_up_pol)
    lattice = fi.metropolis_move(lattice, beta)
    length = len(lattice)
    width = len(lattice[0])
    for i in range(length):
        for j in range(width):
            assert lattice[i][j] == 1    
    
    
def test_evolution_high_T(N = 2, M = 3, spin_up_pol = 1, seed = 1, beta = 0.0):
    """
    Test that a lattice evolves as simulated; in particular, that a fully
    polarized lattice at infinite temperature gets disordered easily.

    """
    
    lattice = fi.initialize_state(N, M, spin_up_pol, seed)
    evolved_lattice = fi.metropolis_move(lattice, beta)
    simulated_lattice = ([-1, -1, 1], [1, -1, -1])
    assert np.array_equal(evolved_lattice, simulated_lattice) == True


def test_low_T_ordering(N = 2, M = 3, spin_up_pol = 0.8, seed = 4, beta = np.inf):
    """
    Test that a lattice evolves as simulated; in particular, that an almost
    fully polarized lattice at zero temperature goes to fully polarized 
    equilibrium configuration.

    """
    
    lattice = fi.initialize_state(N, M, spin_up_pol, seed)
    evolved_lattice = fi.metropolis_move(lattice, beta)
    simulated_lattice = ([1, 1, 1], [1, 1, 1])
    assert np.array_equal(evolved_lattice, simulated_lattice) == True


def test_high_T_disorder(N = 2, M = 3, seed = 8, beta = 0.0):
    """
    Test that a lattice evolves as simulated; in particular, that a random
    lattice at infinite temperature does not order.

    """
    
    lattice = fi.initialize_state(N, M, seed = seed)
    evolved_lattice = fi.metropolis_move(lattice, beta)
    simulated_lattice = ([-1, 1, 1], [-1, 1, 1])
    assert np.array_equal(evolved_lattice, simulated_lattice) == True


#Test the functions that calculate energy and magnetization
def test_energy(N = 2, M = 3, seed = 2):
    """
    Test that the energy is calculated correctly.

    """
    
    lattice = fi.initialize_state(N, M, seed = seed)
    energy = fi.calculate_energy(lattice)
    calculated_energy = 0.0
    assert energy == calculated_energy


def test_mag(N = 2, M = 3, seed = 2):
    """
    Test that the magnetization is calculated correctly.

    """
    
    lattice = fi.initialize_state(N, M, seed = seed)
    mag = fi.calculate_magnetization(lattice)
    calculated_mag = 0.0
    assert mag == calculated_mag
    

#Test the function that reads the configuration parameters
def test_read_configuration(filename = ''):
    """
    Test that an error is raised if the configuration file is not found.

    """
    
    with pytest.raises(FileNotFoundError):
        config = fi.read_configuration(filename)


#Test the function that simulates lattice evolution at certain time instants
def test_simulate_length(N = 2, M = 3, beta = 1.0, times = [1, 2, 3, 4, 5]):
    """
    Test that the length of the array containing the lattice at various time instants 
    is one longer (due to the initial state added).

    """
    
    lattice = fi.initialize_state(N, M)
    evolved_states = fi.simulate(lattice, beta, times)
    assert len(evolved_states) == len(times) + 1
    
    
def test_simulate_times_independent(N = 2, M = 3, seed = 3, beta = 1.0, times1 = [2, 4, 6, 8, 10], times2 = [2, 5, 7, 9, 10]):    
    """
    Test that from the same initial state, the same evolved state is obtained at the same 
    time insant independently from the intermediate steps.

    """
    
    lattice = fi.initialize_state(N, M, seed = seed)
    evolved_states1 = fi.simulate(lattice, beta, times1)
    lattice = fi.initialize_state(N, M, seed = seed)
    evolved_states2 = fi.simulate(lattice, beta, times2)
    assert np.array_equal(evolved_states1[1], evolved_states2[1]) == True
    assert np.array_equal(evolved_states1[5], evolved_states2[5]) == True













