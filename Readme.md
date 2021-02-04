
This module allows to compute the optical depth for given specie (currently HI, HeII)
and transition (currently Ly-alpha and Ly-beta) along sightlines in simulation. 

============= Prerequisite ===============

numpy, scipy and mpi4py (if needs to run in parallel mode).
The code has been tested for python 2.7

============= Installation ===============

1) Copy the files in your machine / server.

2) You need to modify / supply 8 variables in file main_routines_tau.py

n_XI_arr     -->> Array containing number density of specie XI (in units of cm^-3)
T_XI_arr     -->> Array containing temperature of specie XI (in units of kelvin)
v_XI_arr     -->> Array containing line of sight peculiar velocity of specie XI (in units of km/s)
z_arr_sim    -->> Array containing redshifts along sightlines (dimensionless)
hubble_param -->> Hubble parameter (e.g., 0.678)
L_box        -->> Length of simulation box in h^-1 cMpc
N_Grid       -->> Number of grids on to which fields are gridded
specie       -->> Specie for which optical depth needs to be calculated e.g., "HI" or "HeII"
transition   -->> Transition of specie e.g., "lya" --> Ly-alpha, "lyb" -->> Ly-beta

3) Structure of array n_XI_arr, T_XI_arr, v_XI_arr:
If N_res number of pixels along sightlines and N_spec is number of different 
spectra/sightline/skewer then n_XI_arr, T_XI_arr, v_XI_arr should be of 
dimensions (N_res,N_spec)

z_arr_sim should be 1d array of dimension (N_res,)

4) If you need to add more specie and transition, you just need to add entries in 
variable line_dict (see parameters.py file and generate_XI_forest_line_dict function)
The code will automatically compute the optical depth for that specie and transition.
[Note: n_XI_arr, T_XI_arr, v_XI_arr now corresponds to number density of specie XI]

5) You need to modify main_routines_tau.py file as per your need to save tau array
in text, binary or hdf5 file

================ Usage ==================

6) Few sample fields are given in directory ./test_data/

7) You can directly copy directory and run the code on fields provided in ./test_data/ 

The syntax to run code in serial mode is:

python main_routines_tau.py n

"n" tells code to perform task serially.

The syntax to run code in serial mode is:
 
mpirun -np 5 python main_routines_tau.py y

"y" tells code to perform task in parallel.

================ Citation ==================

8) If you find this code useful in your research, please consider 
citing following work:

Link: https://ui.adsabs.harvard.edu/abs/2018MNRAS.474.2233G/abstract

Bibcode:

@ARTICLE{2018MNRAS.474.2233G,
       author = {{Gaikwad}, Prakash and {Choudhury}, Tirthankar Roy and {Srianand}, Raghunathan and {Khaire}, Vikram},
        title = "{Efficient adiabatic hydrodynamical simulations of the high-redshift intergalactic medium}",
      journal = {\mnras},
     keywords = {methods: numerical, intergalactic medium, quasars: absorption lines, large-scale structure of Universe, Astrophysics - Cosmology and Nongalactic Astrophysics, Astrophysics - Astrophysics of Galaxies},
         year = 2018,
        month = feb,
       volume = {474},
       number = {2},
        pages = {2233-2258},
          doi = {10.1093/mnras/stx2859},
archivePrefix = {arXiv},
       eprint = {1705.05374},
 primaryClass = {astro-ph.CO},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2018MNRAS.474.2233G},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
