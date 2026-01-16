#!/bin/bash
import numpy as np
import sys, os

Energy = 10 #eV

MASS = 30.9738 # Phosphorus
#Conversion
mass_unit   = 0.001 / (6.02214076 * 10**23) #g/mol -> kg 
eV_unit     = 1.60217646 * 10 **-19

velocity    = np.sqrt(2*Energy*eV_unit/(MASS*mass_unit)) / 100 # ==> m/s

print(velocity)
