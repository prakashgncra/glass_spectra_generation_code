import numpy as np
import sys
import matplotlib.pyplot as plt

from parameters import generate_XI_forest_line_dict,physical_constants
from sub_routines_tau import sub_routines_tau_primary


def generate_od_data_main():
    if len(sys.argv) == 2:
        flag_mpi = sys.argv[1]
    else:
        ECHO_ARGUMENT_REQUIRED(["flag_mpi 'y' for MPI and 'n' for Serial"])

    # Do not change variable names in this routine
    # You can change variable values
    n_XI_arr     = np.loadtxt("./test_data/n_HI_field.txt")
    T_XI_arr     = np.loadtxt("./test_data/T_HI_field.txt")
    v_XI_arr     = np.loadtxt("./test_data/v_HI_field.txt")
    z_arr_sim    = np.loadtxt("./test_data/z_arr.txt")
    hubble_param = 0.678
    L_box        = 40.0
    N_Grid       = 2048
    specie       = "HI"
    transition   = "lya"
    line_dict    = generate_XI_forest_line_dict()
    phy_const    = physical_constants()
    inp_dict     = vars()
    out_dict     = sub_routines_tau_primary(inp_dict)

    rank         = out_dict["rank"]
    if rank == 0:
        tau_XI_arr = out_dict["tau_XI_arr"]
        # Save or plot tau_XI_arr as needed
        print out_dict.keys()

def ECHO_ARGUMENT_REQUIRED(key_list):
    N_key_list = len(key_list)
    print
    print("Required Following " + str("%i"%N_key_list) + " Argument ... ")
    for key_indx in range(N_key_list):
        print("Argument " + str("%i"%(key_indx+1)) + " : " + key_list[key_indx])
    print("Exiting the code now !!!!")
    print
    exit()


generate_od_data_main()
