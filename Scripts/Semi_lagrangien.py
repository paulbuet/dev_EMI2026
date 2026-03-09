import numpy as np
import xarray as xr

# On importe ici les classes extèrieures
from fonctions import InitialCond
from fonctions import Eq



class Model_sl():
    def __init__(self, number_stitches,time_step,number_particules,speed_max,esp,CFL, lenght_sim):   
        """
        Here we initialise the non-spatial fixed parameters and allow important variables 
        to travel between functions. We also call the initialisation.
        """

        self.nb_part = number_particules

        self.number_stitches = number_stitches

        self.length_sim = lenght_sim  # length of simulation in seconds

        self.delta_t = time_step # length of time step in seconds

        self.nb_step = self.length_sim // self.delta_t  # number of time step
   