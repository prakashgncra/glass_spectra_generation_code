import numpy as np
import re


def generate_XI_forest_line_dict():
    # Change this routine if you want to add new specie 
    """
    --------------------- Dictionary for Line Transition  -------------------------

     Method to calculate i_alpha
     i_alpha   = (pi * e^2 * f_12) / (m_2 * c * nu_12)
     i_alpha_X = 8.83788e-21 * (f_12_X * lambda_rest_X) 

     HI   -->> Neutral Hydrogen
     HeI  -->> Neutral Helium
     HeII -->> Single ionized Helium
     lya  -->> Lyman alpha
     lyb  -->> Lyman beta
     lyg  -->> Lyman gamma
     lyd  -->> Lyman delta
    """

    line_dict = {\
                ("HI","lya")  :{"lambda_rest":1215.6701,"f_12":0.416,"damp_gamma":6.265e8, "i_alpha":4.469e-18,"n_proton":1,"n_neutron":0},\
                ("HI","lyb")  :{"lambda_rest":1025.7223,"f_12":0.079,"damp_gamma":1.897e8, "i_alpha":7.162e-19,"n_proton":1,"n_neutron":0},\
                ("HeII","lya"):{"lambda_rest":303.78220,"f_12":0.416,"damp_gamma":1.003e10,"i_alpha":1.117e-18,"n_proton":2,"n_neutron":2},\
                ("HeII","lyb"):{"lambda_rest":256.31700,"f_12":0.079,"damp_gamma":2.667e9, "i_alpha":1.790e-19,"n_proton":2,"n_neutron":2},\
                }
    return line_dict

def physical_constants():
    # No need to change this function 

    # Taken from http://www.astro.wisc.edu/~dolan/constants.html

    #--------------------- Physical Constants -------------------------

    Boltz_const_mag       = 1.380658          # Magnitude
    Boltz_const_pow_cgs   = -16               # Power in ergs / K
    Boltz_const_pow_si    = -23               # Power in J / K
    Boltz_const_cgs       = Boltz_const_mag * 10**Boltz_const_pow_cgs
    Boltz_const_si        = Boltz_const_mag * 10**Boltz_const_pow_si
    mass_proton_mag       = 1.6726231    # Magnitude In units of
    mass_proton_pow_cgs   = -24          # Power In units of 10^-24 gm
    mass_proton_pow_si    = -27          # Power In units of 10^-27 Kg
    mass_proton_cgs       = mass_proton_mag * 10**mass_proton_pow_cgs
    mass_proton_si        = mass_proton_mag * 10**mass_proton_pow_si
    mass_neutron_mag      = 1.6749286    # Magnitude In units of
    mass_neutron_pow_cgs  = -24          # Power In units of 10^-24 gm or 10^-27 Kg
    mass_neutron_pow_si   = -27          # Power In units of 10^-24 gm or 10^-27 Kg
    mass_neutron_cgs      = mass_neutron_mag * 10**mass_neutron_pow_cgs
    mass_neutron_si       = mass_neutron_mag * 10**mass_neutron_pow_si
    mass_electron_mag     = 9.1093897    # Magnitude In units of
    mass_electron_pow_cgs = -28          # Power In units of 10^-24 gm or 10^-27 Kg
    mass_electron_pow_si  = -31          # Power In units of 10^-24 gm or 10^-27 Kg
    mass_electron_cgs     = mass_electron_mag * 10**mass_electron_pow_cgs
    mass_electron_si      = mass_electron_mag * 10**mass_electron_pow_si

    h_planck_mag          = 6.6260755   # Magnitude of Planck's Constant
    h_planck_pow_cgs      = -27         # Power In units of 10^-27 erg s
    h_planck_pow_si       = -34         # Power In units of 10^-34 J s
    h_planck_cgs          = h_planck_mag * 10**h_planck_pow_cgs
    h_planck_si           = h_planck_mag * 10**h_planck_pow_si


    cp_by_cv        = 5.0 / 3.0         # Ratio of Specific Heats (Temperature - Internal Energy Conversion)
    mass_proton     = mass_proton_cgs   # In grams taken from http://www.astro.wisc.edu/~dolan/constants.html
    mass_neutron    = mass_neutron_cgs  # In grams taken from http://www.astro.wisc.edu/~dolan/constants.html
    mass_electron   = mass_electron_cgs # In grams taken from http://www.astro.wisc.edu/~dolan/constants.html
    Boltz_const     = Boltz_const_cgs   # In ergs / K taken from http://www.astro.wisc.edu/~dolan/constants.html
    h_planck        = h_planck_cgs      # In erg s taken from http://www.astro.wisc.edu/~dolan/constants.html

    #--------------------- Collect all parameters in phy_const dictonary -------------------------

    phy_const = vars()
    return phy_const
