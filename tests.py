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

#Test functions for energy and magnetization calculation
def test_correct_total_energy(N = 2, M = 3, choice = True, spin_up_pol = 1):
    lattice = fi.initialize_state(N, M, choice, spin_up_pol)
    assert fi.calculate_energy(lattice) == -N*M*4/2  
    
def test_correct_total_magnetization(N = 2, M = 3, choice = True, spin_up_pol = 1):
    lattice = fi.initialize_state(N, M, choice, spin_up_pol)
    assert fi.calculate_magnetization(lattice) == N*M

















