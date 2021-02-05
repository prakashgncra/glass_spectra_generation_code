import numpy as np
from scipy.special import wofz # Faddeeva function


def sub_routines_tau_primary(inp_dict):
    inp_dict  = STEP_0_CHECK_INPUT(inp_dict)
    flag_mpi  = inp_dict["flag_mpi"]
    if flag_mpi == "y":
        out_dict = MPI_LY_ALPHA_FOREST(inp_dict)
    else:
        out_dict = SERIAL_LY_ALPHA_FOREST(inp_dict)
    return out_dict

def STEP_0_CHECK_INPUT(inp_dict):
    key_list = ["n_XI_arr","v_XI_arr","T_XI_arr","L_box",\
                "N_Grid","hubble_param","line_dict","specie","transition"]
    CHECK_IF_KEY_PRESENT(key_list,inp_dict)
    inp_dict = CHECK_AND_INITIALIZE_KEY_VALUE(inp_dict,"flag_mpi","n")
    #inp_dict = CHECK_AND_INITIALIZE_KEY_VALUE(inp_dict,"specie","HI")
    #inp_dict = CHECK_AND_INITIALIZE_KEY_VALUE(inp_dict,"transition","lya")

    return inp_dict

def SERIAL_LY_ALPHA_FOREST(inp_dict):
    wavelength_arr,tau_XI_arr = GENERATE_SPECTRUM(inp_dict)
    out_dict              = {"tau_XI_arr":tau_XI_arr,"rank":0,"wavelength_arr":wavelength_arr}
    return out_dict

########################################################
######   Optical Depth Calculation (MPI Routines) ######
########################################################

def MPI_LY_ALPHA_FOREST(inp_dict):
    mpi_dict    = mpi_variable()
    comm        = mpi_dict["comm"]
    rank        = mpi_dict["rank"]
    N_proc      = mpi_dict["N_proc"]

    if rank == 0:
        n_XI_arr    = inp_dict["n_XI_arr"]
        T_XI_arr    = inp_dict["T_XI_arr"]
        v_XI_arr    = inp_dict["v_XI_arr"]
        z_arr_sim   = inp_dict["z_arr_sim"]
        N_los       = T_XI_arr.shape[1]
        od_los_dict = mpi_distribute_od_los_idx(N_los,1,N_proc)
        lb_arr      = od_los_dict["lb_arr"]
        ub_arr      = od_los_dict["ub_arr"]

        for proc_indx in xrange(1,N_proc):
            lb = lb_arr[proc_indx-1]
            ub = ub_arr[proc_indx-1]
            comm.send((n_XI_arr[:,lb:ub],v_XI_arr[:,lb:ub],T_XI_arr[:,lb:ub],z_arr_sim),dest=proc_indx,tag=proc_indx)

        print "Data sent to all processors ..."
        comm.Barrier()
        N_pixel   = z_arr_sim.shape[0]
        N_spectra = od_los_dict["N_spectra"]
        spec_lb   = od_los_dict["spec_lb"]
        spec_ub   = od_los_dict["spec_ub"]
        tau_all   = np.zeros((N_pixel,N_spectra))

        for proc_indx in xrange(1,N_proc):
            lb                 = spec_lb[proc_indx-1]
            ub                 = spec_ub[proc_indx-1]
            wavelength_arr,tau = comm.recv(source = proc_indx,tag = proc_indx)
            tau_all[:,lb:ub]   = tau[:,:]
            print "Data received from processor ",proc_indx
        out_dict    = {"tau_XI_arr":tau_all,"rank":rank,"wavelength_arr":wavelength_arr}
        comm.Barrier()
    else:
        n_XI_arr,v_XI_arr,T_XI_arr,z_arr_sim = comm.recv(source=0,tag=rank)
        inp_dict.update({"n_XI_arr":n_XI_arr,"v_XI_arr":v_XI_arr,"T_XI_arr":T_XI_arr})
        comm.Barrier()
        wavelength_arr,tau  = GENERATE_SPECTRUM(inp_dict,rank)
        comm.send((wavelength_arr,tau),dest=0,tag=rank)
        out_dict = {"rank":rank}
        comm.Barrier()
    return out_dict


##################################################################################
#######  CODES FOR SPECTRA GENERATION ALONG SINGLE AND MULTIPLE SIGHTLINE  #######
##################################################################################

def CALCULATE_PROFILE(z_val,velocity,doppler_b,n_XI_arr,z_arr_sim,lambda_rest,damp_gamma): # Used internally by SHOOT_ALONG_RANDOM_LOS GENERATE_SPECTRUM_LOG_NORMAL
    from scipy.constants import c as light_speed
    light_speed	   	   /= 1000.0 # Speed of Light in km/s
    voigt_factor_new	= np.zeros(n_XI_arr.shape)
    voigt_factor_new	= (velocity+(light_speed * (z_arr_sim - z_val) / (1 + z_val))) / doppler_b
    voigt_calc 			= voigt_faddeva(voigt_factor_new,doppler_b,lambda_rest,damp_gamma)           # Fastest and accurate method python
    tau_new 			= np.sum(n_XI_arr*voigt_calc/(doppler_b*(1.0+z_arr_sim)),axis=0) 
    return tau_new	

def CALCULATE_PROFILE_PERIODIC(z_val,velocity,doppler_b,n_XI_arr,z_arr_sim,lambda_rest,damp_gamma):
    from scipy.constants import c as light_speed
    light_speed	   	   /= 1000.0 # Speed of Light in km/s
    voigt_factor_new	= np.zeros(n_XI_arr.shape)

    # Following lines of code make sure that the tau is calculated assuming
    # periodic boundary condition. Such boundary conditions are used 
    # in Sherwood los and tau files LOS_extraction code.

    delta_z                = np.median(z_arr_sim[1:] - z_arr_sim[:-1])
    delta_z_by_2           = delta_z*0.5
    z_diff_arr             = z_arr_sim-z_val
    bool_arr_1             = z_diff_arr > delta_z_by_2
    bool_arr_2             = z_diff_arr < (delta_z_by_2*-1.0)
    z_diff_arr[bool_arr_2] = delta_z + z_diff_arr[bool_arr_2]

    # Rest of the code is same 

    voigt_factor_new	= (velocity+(light_speed * z_diff_arr / (1 + z_val))) / doppler_b
    voigt_calc 			= voigt_faddeva(voigt_factor_new,doppler_b,lambda_rest,damp_gamma)           # Fastest and accurate method python
    tau_all_arr         = n_XI_arr*voigt_calc/(doppler_b*(1.0+z_arr_sim))
    tau_new 			= np.sum(tau_all_arr,axis=0)

    return tau_new
        
def GENERATE_SPECTRUM(inp_dict,rank=0):
    """ line_dict has following structure
        line_dict = {("HI","lya"):{specie:"HI","transition":lya,"lambda_rest":1215.6701,"damp_gamma":6.265e8,"i_alpha":4.45e-18}}
        line_sub_dict -->> line_dict[("HI","lya")]
        atomic_line_info -->> {specie,transistion}
    """
    from scipy.constants import c as light_speed

    inp_dict        = COMPUTE_DOPPLER_B_PARAMETER(inp_dict)
    L_box_without_h = inp_dict["L_box"]
    N_Grid          = inp_dict["N_Grid"]
    n_XI_arr        = inp_dict["n_XI_arr"]
    z_arr_sim       = inp_dict["z_arr_sim"]
    doppler_b       = inp_dict["doppler_b"]
    v_XI_arr        = inp_dict["v_XI_arr"]
    specie          = inp_dict["specie"]
    transition      = inp_dict["transition"]
    hubble_param    = inp_dict["hubble_param"]
    line_dict       = inp_dict["line_dict"]
    light_speed	   /= 1000.0 # Speed of Light in km/s

    L_box           = L_box_without_h / hubble_param
    dL_box          = float(L_box) / (N_Grid)
    line_sub_dict   = line_dict[(specie,transition)]
    lambda_rest     = line_sub_dict["lambda_rest"]
    damp_gamma      = line_sub_dict["damp_gamma"]
    i_alpha         = line_sub_dict["i_alpha"]
    tau 			= np.zeros(n_XI_arr.shape)
    
    if n_XI_arr.ndim == 1: 
        N_res           = n_XI_arr.shape[0]
        z_arr_sim_nD    = z_arr_sim.copy() 
        for indx in xrange(N_res):
            z_val 		= z_arr_sim[indx]
            tau[indx]   = CALCULATE_PROFILE_PERIODIC(z_val,v_XI_arr,doppler_b,n_XI_arr,z_arr_sim_nD,lambda_rest,damp_gamma)
    else:
        N_res,N_spectra = n_XI_arr.shape[0],n_XI_arr.shape[1]
        z_arr_sim_nD    = np.array([z_arr_sim,]*N_spectra).T
        kounter         = 0
        perc_list       = [0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,110.0]
        for indx in xrange(N_res):
            z_val 		= z_arr_sim[indx]
            tau[indx,:] = CALCULATE_PROFILE_PERIODIC(z_val,v_XI_arr,doppler_b,n_XI_arr,z_arr_sim_nD,lambda_rest,damp_gamma)
            perc_val    = perc_list[kounter]

            if abs((indx * 100.0 / float(N_res)) - perc_val) < 1e-1:
                print str("%i"%perc_val) + "% Work completed on Rank: " + str("%i"%rank),"..."
                kounter += 1

    tau            = light_speed * i_alpha * (dL_box) * tau / np.sqrt(np.pi)
    tau            = tau * 3.0856e24
    wavelength_arr = (1.0 + z_arr_sim) * lambda_rest
    return wavelength_arr,tau

def COMPUTE_DOPPLER_B_PARAMETER(inp_dict):
    #from common_lib.parameter_files.physical_constants import phy_const
    specie        = inp_dict["specie"]
    transition    = inp_dict["transition"]
    line_dict     = inp_dict["line_dict"]
    T_XI_arr      = inp_dict["T_XI_arr"]
    Boltz_const   = phy_const["Boltz_const_mag"]   # Bolt_const ~ 1.38 (In atomic units)
    mass_proton   = phy_const["mass_proton_mag"]   # mass_proton ~ 1.67 (In atomic units)
    mass_neutron  = phy_const["mass_neutron_mag"]  # mass_proton ~ 1.67 (In atomic units)
    line_sub_dict = line_dict[(specie,transition)]
    n_proton      = line_sub_dict["n_proton"]
    n_neutron     = line_sub_dict["n_neutron"]
    mass_X        = (mass_proton * n_proton + mass_neutron * n_neutron)
    doppler_b 	  =	np.sqrt(2.0 * Boltz_const * T_XI_arr / mass_X)*1.e-1 # Note that the b_thermal depends on mass
    inp_dict.update({"doppler_b":doppler_b})
    return inp_dict

##################################################
######   Voigt Profile Calculation Routines ######
##################################################

def voigt_faddeva(voigt_factor,doppler_b,lambda_rest,damp_gamma):
    pi         = np.pi;
    alpha      = (damp_gamma * lambda_rest * 1e-8)/(4 * pi * doppler_b * 1.e5) #* (pi/4.e0)
    x          = voigt_factor
    y          = alpha
    z          = x + 1j*y
    I          = wofz(z).real
    return I

######################################
######   Miscellaneous Routines ######
######################################

def CHECK_IF_KEY_PRESENT(key_list,inp_dict):
    inp_dict_key_list = list(inp_dict.keys())
    for key_name_reqd in key_list:
        if key_name_reqd not in inp_dict_key_list:
            print("inp_dict Dictionary should contain key :",key_name_reqd)
            print("Exiting the code now !!!")
            exit()

def CHECK_AND_INITIALIZE_KEY_VALUE(inp_dict,key_name,default_key_value):
    if key_name not in list(inp_dict.keys()):
        inp_dict.update({key_name:default_key_value})
    return inp_dict

def mpi_variable():
    from mpi4py import MPI
    comm      = MPI.COMM_WORLD
    rank      = comm.Get_rank()
    N_proc    = comm.Get_size()
    status    = MPI.Status()
    proc_name = MPI.Get_processor_name()
    mpi_dict  = {"comm":comm,"rank":rank,"N_proc":N_proc,"status":status,"MPI":MPI,"proc_name":proc_name}
    return mpi_dict

def mpi_distribute_od_los_idx(N_los,N_box,N_proc):
    """
    Let N_los = 10002, N_box = 5 and N_proc = 84
    N_spectra = 2000
    N_spectra_per_proc = 25
    spec_lb_ub_arr ==>> 25, 25, 25, 25, ...., 25, 25 (84 times)
    N_spectra_extra = 2075-2000 = 75
    After the for loop output will be
    spec_lb_ub_arr ==>> 24, 24, 24, 24, ...., 25, 25 (24 -->> 75 times, 25 -->> 83-75 = 8 times)


    lb_arr ==>> [0  120  240  360  480 ... 9750 9875]
    ub_arr ==>> [120   240   360   480 ... 9875 10000]

    """
    N_spectra          = int(N_los / N_box)
    N_spectra_per_proc = int(N_spectra / (N_proc-1)) + 1
    spec_lb_ub_arr     = np.ones(N_proc-1,dtype=np.int64)*N_spectra_per_proc
    N_spectra_extra    = int(N_spectra_per_proc * (N_proc-1) - N_spectra)
    for indx in range(N_spectra_extra):
        spec_lb_ub_arr[indx] -= 1
    spec_ub     = np.cumsum(spec_lb_ub_arr)
    spec_lb     = np.zeros(spec_ub.shape[0],dtype = np.int64)
    spec_lb[1:] = spec_ub[:-1]

    lb_ub_arr   = spec_lb_ub_arr * N_box
    ub_arr      = np.cumsum(lb_ub_arr)
    lb_arr      = np.zeros(ub_arr.shape[0],dtype  = np.int64)
    lb_arr[1:]  = ub_arr[:-1]
    od_los_dict = {"lb_arr":lb_arr,"ub_arr":ub_arr,"N_spectra":N_spectra,"spec_lb":spec_lb,"spec_ub":spec_ub}
    return od_los_dict















