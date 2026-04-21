#!/usr/bin/env python3

"""
This script calculates the calibration rates.
The equations and quantities are taken from DocDB-2904.

Eq (1) OR (2) are the same, but calculate different quantities 
(i.e., # of calibration gammas per crystal given calibration time 
OR time needed to calibrate each crystal given the number of calibration 
gammas per crystal, respectively).
"""

import matplotlib.pyplot as plt
import numpy as np
import math

# =============================================================================
# List of quantities:
# =============================================================================
sim_gam = 6.5e+7      # total number of photons simulated in GEANT4
n_xtal = 1348         # total number of crystals 
V_bath = 0.021        # Volume of bath at dt generator 
# n_16N = 1.4e+8      # Number of 16N in bath  
a = 0.40              # Attenuation from bath to furthest crystal
tao = 10.3            # Lifetime of 16N
e_geom = 0.3          # (n_xtal)/ (7e-5 * sim_gam ) #Geometric efficiency to enter crystal
f_16N = 0.3           # fractions of inelastic reactions producing 16N
w_bath = 0.16         # width of bath as seen by neutron
mfp = 3.0             # mean free path of neutron in fluid
e_bath = 0.88         # geometric efficiency of bath

# =============================================================================
# To account for changing neutron yield: 
# Change the neutron rate from the DT generator, R_DT, below and multiply 
# by the reduced yield.
# =============================================================================
R_DT = 3e+08    # neutron rate from DT generator

# =============================================================================
# Equations from paper:
# =============================================================================
T_prod_16N = f_16N * R_DT * (1 - math.exp(-w_bath/mfp)) * e_bath # rate to produce 16N in bath 
n_16N = T_prod_16N * tao                                         # Number of 16N in bath  
density_16N = n_16N / V_bath                                     # Number density of 16N in bath
V_xtal = 0.0019 / n_xtal                                         # Volume of fluid at a crystal
n_N_xtal = V_xtal * density_16N * a                              # Number of 16N at crystal


# =============================================================================
# Terminal Interface Logic
# =============================================================================
print("\n===========================================")
print(" Welcome to Calibration Rate Calculator")
print("===========================================\n")

try:
    option = int(input("What type of information will you input:\n (1) Running time (in minutes) \n (2) Number of calibration gammas (per crystal)\n\nNote: Enter integer value (1 or 2): "))
    print("") # Print a blank line for readability

    if option == 1:
        # =============================================================================
        # (1) Equation to calculate number of calib gammas per crystal, 
        # given calibration time (T_cal)
        # =============================================================================
        time_input = float(input("Please enter time in minutes: "))
        T_cal = time_input * 60       # time in seconds
        n_gam = (e_geom * T_cal * n_N_xtal) / tao  # number of calibration gammas per calibration per crystal

        print("\n--- RESULTS ---")
        print("Time you entered for calibration (in minutes): {:.3f}".format(T_cal/60))
        print("Time you entered for calibration (in seconds): {}".format(T_cal))
        print("Number of calibration gammas per crystal: {:,.0f}".format(n_gam))
        print("\n")

    elif option == 2:
        # =============================================================================
        # (2) Equation to calculate calibration time, 
        # given number of calibration gammas per crystal (calib_gammas)
        # =============================================================================
        gammas_input = float(input("Please enter number of calibration gammas: "))
        T_calib = (gammas_input * tao) / (e_geom * n_N_xtal) # Time duration of calibration (in seconds)

        print("\n--- RESULTS ---")
        print("Time needed for calibration provided {} calibration gammas: {:.3f} minutes".format(gammas_input, T_calib/60))
        print("Time needed for calibration provided {} calibration gammas: {:.2f} seconds".format(gammas_input, T_calib))
        print("\n")

    else:
        print("Invalid option. Please run the script again and type '1' or '2'.")

except ValueError:
    print("\nError: Invalid input! Please enter numbers only.")
