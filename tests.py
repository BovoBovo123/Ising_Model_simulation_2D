# -*- coding: utf-8 -*-
"""
Created on Fri Dec 17 14:59:01 2021

@author: pietr
"""

import functions_ising
import numpy as np

import pytest

#Test for functions_ising.initialize_state function
def test_lattice_dimensions(N = 3, M = 2):
    lattice = functions_ising.initialize_state(N, M)
    assert len(lattice) == 3
    assert len(lattice[0]) == 2
    
    
def test_lattice_dimensions_polarized(N = 3, M = 2, choice = True, spin_up_pol = 0.3):
    lattice = functions_ising.initialize_state(N, M)
    assert len(lattice) == 3
    assert len(lattice[0]) == 2
    
def test_spin_values(N = 3, M = 2):
    lattice = functions_ising.initialize_state(N, M)
    abs_spin = np.abs(lattice)
    assert abs_spin.all() == 1
    
def test_spin_values_polarized(N = 3, M = 2, choice = True, spin_up_pol = 0.3):
    lattice = functions_ising.initialize_state(N, M)
    abs_spin = np.abs(lattice)
    assert abs_spin.all() == 1

def test_raises_error_lattice_dimensions(N = -1, M = -2):
    with pytest.raises(ValueError):
        lattice = functions_ising.initialize_state(N, M)

def test_choice_not_boolean(N = 2, M = 3, choice = 'not a boolean'):
    with pytest.raises(TypeError):
        lattice = functions_ising.initialize_state(N, M, choice) 

def test_spin_up_polarization_not_a_percentage(N = 2, M = 3, choice = True, spin_up_pol = 1.2):
    with pytest.raises(ValueError):
        lattice = functions_ising.initialize_state(N, M, choice, spin_up_pol)
        



























