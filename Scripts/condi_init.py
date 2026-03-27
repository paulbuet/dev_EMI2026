#!/usr/bin/env python3

# Copyright (c) Météo France (2025-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

# Le but de ce fichier est de mettre toutes les fonctions qu'on peut potentiellement appelés pour faire tourner nos futur codes
# Ensuite il suffira de faire un appelle du fichier et de la fonction qu nous interesse
# Exemple avec cette fonction test que j'appelle dans le code Sedimentation.py

### Imports ###

from math import *
import time
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy import integrate
from collections import deque
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import MaxNLocator
from pathlib import Path

from equations import Eq
from equations import Selection
from formatage import Formatage

def gaussienne (m, sigma,h):
        return (1/(sigma*np.sqrt(2*np.pi))) * exp (-1/2 * ((h-m)/sigma)**2)


class profil_rho_r:

    def __init__(self):

        # On définit les constantes que l'on va utiliser

        self.R = 287.058
        self.exp_baro = 9.80665/(self.R*0.0065)

        self.T0 = 288.150
        self.P0 = 101325

    def calcul(self,alt,r):
        alt = np.array(alt)
        r = np.array(r)
        
        T = self.T0 -0.0065 * alt

        P = self.P0 * (T/self.T0)**(self.exp_baro)

        rho = P/(self.R*T)

        rho_r = rho * r

        return rho_r,rho





class InitialCond :
    '''
    Class that poduces a dataset to initialize the sedimentation model, a type of initialization is used in a restricted area of the vertical gridding.
    The dataset contains the concentration profile of the specie for each bin of the sedimentation scheme.

    Parameters
    ----------
    nb_grid : number of grid cells in the model
    esp : specie considered "i", "s", "g", "r" or "c"
    mode : type of intial condition used, by default the last grid cell contains a relative concentration of 1
    Hmax : atltitude of the maximum of the gaussian fonction of the concentration if the type of initial condition is "gauss"
    sigma : value of the sigma parameter used for the gaussian function if the type of initial condition is "gauss"
    

    Attributes
    ----------
    self.levels_boundaries : list containing the altitudes of the boundaries of each grid mesh point
    self.grid : list containing the grid points
    self.data : DataArray having for variables
        concentration_bin_{i} : the concentrations of the specie in the i-th bin over the vertical grid
        diameter_bin_{i} : the diameter of the specie for the i-th bin
    and for coordinates 
        level : the levels of the grid points

    '''
    

    def __init__(self, nb_grid, esp, types, mode = "simple", Hmax = 5000, sigma = 2000, nb_classes = 300, r = 0.001) :
        eq=Eq(esp)
        form = Formatage(esp)
        selec = Selection(esp)



        if types == "bin":

            if nb_grid == "ARO" : 
                self.levels_boundaries = [5.00148256575414, 16.7609146275979, 31.9999856716034, 50.6506387418972, 72.6448134875948, 97.9144556307367, 126.391508411840, 158.007915671418, 192.695621996127, 230.386571302649, 271.012706897240, 314.505973414367, 360.798314810016, 409.821674828922, 461.507997847950, 515.890042408161, 573.093937517738, 633.231198491812, 696.398736766220, 762.678848044096, 832.139215500780, 904.832909766711, 980.798389785428, 1060.05950095623, 1142.62547636213, 1228.49093654830, 1317.63588971212, 1410.02573054417, 1505.73432230725, 1604.93984760156, 1707.79812299161, 1814.44258334521, 1924.98428740911, 2039.51193175523, 2280.76794378598, 2407.56182949566, 2538.47269089309, 2673.47735336959, 2812.53026865253, 2955.56351459630, 3102.48360541930, 3253.19498943504, 3407.62661862652, 3565.73195622143, 3727.48898443499, 3892.90019410005, 4061.99258757885, 4234.81767907587, 4411.45149612596, 4591.99457746971, 4776.57197432233, 4965.33325159684, 5158.45249292134, 5356.12828712931, 5558.65852457596, 5766.45816179863, 5980.00284181189, 6199.82890110527, 6426.53336977965, 6660.77397670184, 6903.26915353359, 7154.79804227343, 7416.20049682756, 7688.37709125546, 7972.28912861375, 8268.07660797884, 8575.22917035383, 8893.62412699775, 9223.52646472597, 9565.58886389478, 9920.85173976858, 10290.7432801852, 10677.0795176319, 11082.0546910758, 11510.3223903101, 11966.1795634141, 12454.7528314299, 12982.0128918015, 13554.6557230282, 14180.1025765563, 14866.4999945363, 15627.4731252487, 16485.0471290191, 17467.8152385289, 18610.9387773901, 19955.6070022098, 21550.0160473017, 23450.1477297912, 34461.1214067193]
                nb_grid = 89
            else :
                height_grid = 12e3
                self.levels_boundaries = np.linspace(1/nb_grid * height_grid, height_grid, nb_grid)
            nb_levels = len (self.levels_boundaries)

            boundaries = deque(self.levels_boundaries)
            boundaries.appendleft(0)
            
            self.grid = [(boundaries[i]+boundaries[i+1])/2 for i in range(len(boundaries)-1)]

            if mode == "simple" :
                relative_profile = [0 for i in range (nb_levels)]
                relative_profile [-1] = 1
            if mode == "gauss" :
                relative_profile = [gaussienne(Hmax, sigma, self.grid[i]) for i in range(nb_levels)]

            # Vérification : Affichage des profils relatifs (à décommenter)
            # x = [self.grid[i] for i in range(nb_levels)]
            # y = relative_profile
            # plt.plot(x,y)
            # plt.show()


            self.r_profile = np.array(relative_profile)*float(r)

            self.rho_r_profile,self.rho_profile = profil_rho_r().calcul(self.grid,self.r_profile)

            lam = eq.Liste_Lanbda_1_mom(self.rho_r_profile)
            dmin, dmax = 10e-6, 0.015
            
            N_profile = eq.contenu_to_conc(self.rho_r_profile)
            # N_profile = self.rho * r * (lam ** eq.b)/ (eq.a * eq.G(eq.b))

            levels_bin_concentrations_splittings= [selec.Classe_D_N (nb_classes, dmin, dmax, N_profile[ind_level], lam[ind_level]) if N_profile[ind_level] != 0 else " " for ind_level in range(nb_grid)] # List of lists containing the nb_grid classifications
            levels_bin_contents_splittings= [selec.Classe_D_rho_r (nb_classes, dmin, dmax, N_profile[ind_level], lam[ind_level]) if N_profile[ind_level] != 0 else " " for ind_level in range(nb_grid)] # List of lists containing the nb_grid classifications
            
            # Vérification : Affichage des profils des bins (à décommenter)
            # x = [self.grid[i] for i in range(nb_levels)]
            # y = self.bin_concentration[0][1]*np.array(relative_profile)
            # print(self.bin_concentration)
            # plt.plot(x,y)
            # plt.show()
            # print (bins)
            # print(levels_bin_concentrations_splittings)

            bins_concentrations_profiles = [[splitting[ind_bin][1] if splitting != " " else 0 for splitting in levels_bin_concentrations_splittings] for ind_bin in range(nb_classes)] # computinng of the n bin profiles
            data_vars1 = {f"concentration_bin_{ind_bin+1}" : ("level", bins_concentrations_profiles[ind_bin]) for ind_bin in range(nb_classes)}

            self.diameters = [[splitting for splitting in form.find_diameters_in_splittings(levels_bin_concentrations_splittings)][0][ind_bin][0] for ind_bin in range(nb_classes)]
            data_vars2 = {f"diameter_bin_{ind_bin+1}" : self.diameters[ind_bin] for ind_bin in range(nb_classes)} # addition of the diameters
            
            bins_rho_r_profiles = [[splitting[ind_bin][1] if splitting != " " else 0 for splitting in levels_bin_contents_splittings] for ind_bin in range(nb_classes)]
            data_vars3 = {f"rho_r_bin_{ind_bin+1}" : ("level",bins_rho_r_profiles[ind_bin][:]) for ind_bin in range(nb_classes)}
            data_vars1.update(data_vars2)
            data_vars1.update(data_vars3)

            # print (f" splitting {levels_bin_concentrations_splittings}")
            # print (f" relative profile {relative_profile}")
            # print (f" N profile {N_profile}")
            # print (f" rho_r profile {self.rho_r_profile}")
            # print (f" bins N profiles {bins_concentrations_profiles}")
            # print (f" bins rho_r profiles {bins_rho_r_profiles}")
            # print (f" rho {self.rho_profile}")
            # print (f" lambda {lam}")
            # print (f" diameters {self.diameters}")

            # x = [self.grid[i] for i in range(nb_levels)]
            # for i in range(nb_classes):
            #   y = bins_rho_r_profiles[i]
            #   plt.plot(x,y)
            # plt.show()
            # for i in range(nb_classes):
            #   y = bins_conentrations_profiles[i]
            #   plt.plot(x,y)
            # plt.show()

            self.data = xr.Dataset(data_vars= data_vars1, coords = {"level" : self.grid})

            # print (self.data["rho_r_bin_1"])
            # print ("OK")

        else:

            if nb_grid == "ARO" : 
                self.levels_boundaries = [5.00148256575414, 16.7609146275979, 31.9999856716034, 50.6506387418972, 72.6448134875948, 97.9144556307367, 126.391508411840, 158.007915671418, 192.695621996127, 230.386571302649, 271.012706897240, 314.505973414367, 360.798314810016, 409.821674828922, 461.507997847950, 515.890042408161, 573.093937517738, 633.231198491812, 696.398736766220, 762.678848044096, 832.139215500780, 904.832909766711, 980.798389785428, 1060.05950095623, 1142.62547636213, 1228.49093654830, 1317.63588971212, 1410.02573054417, 1505.73432230725, 1604.93984760156, 1707.79812299161, 1814.44258334521, 1924.98428740911, 2039.51193175523, 2280.76794378598, 2407.56182949566, 2538.47269089309, 2673.47735336959, 2812.53026865253, 2955.56351459630, 3102.48360541930, 3253.19498943504, 3407.62661862652, 3565.73195622143, 3727.48898443499, 3892.90019410005, 4061.99258757885, 4234.81767907587, 4411.45149612596, 4591.99457746971, 4776.57197432233, 4965.33325159684, 5158.45249292134, 5356.12828712931, 5558.65852457596, 5766.45816179863, 5980.00284181189, 6199.82890110527, 6426.53336977965, 6660.77397670184, 6903.26915353359, 7154.79804227343, 7416.20049682756, 7688.37709125546, 7972.28912861375, 8268.07660797884, 8575.22917035383, 8893.62412699775, 9223.52646472597, 9565.58886389478, 9920.85173976858, 10290.7432801852, 10677.0795176319, 11082.0546910758, 11510.3223903101, 11966.1795634141, 12454.7528314299, 12982.0128918015, 13554.6557230282, 14180.1025765563, 14866.4999945363, 15627.4731252487, 16485.0471290191, 17467.8152385289, 18610.9387773901, 19955.6070022098, 21550.0160473017, 23450.1477297912, 34461.1214067193]
                nb_grid = 89
            else :
                height_grid = 12e3
                self.levels_boundaries = np.linspace(1/nb_grid * height_grid, height_grid, nb_grid)
            nb_levels = len (self.levels_boundaries)

            boundaries = deque(self.levels_boundaries)
            boundaries.appendleft(0)
            
            self.grid = [(boundaries[i]+boundaries[i+1])/2 for i in range(len(boundaries)-1)]

            if mode in ["simple","ajout_continu"] :
                relative_profile = [ 0 for i in range(nb_levels)]
                relative_profile[-1] = 1
                #relative_profile[-2] = 1/2

            if mode == "gauss" :
                relative_profile = [gaussienne(Hmax, sigma, self.grid[i]) for i in range(nb_levels)]

            self.rho_r_profile,self.rho = profil_rho_r().calcul(self.grid,np.array(relative_profile)*float(r))

            eq=Eq(esp)

            N_profile = eq.contenu_to_conc(self.rho_r_profile)

            bulk_profile = np.array(N_profile)   # computinng of the n bin profiles
            data_vars1 = {"concentration" : ("level", bulk_profile), "rho_r": ("level",self.rho_r_profile)}

            if mode == "ajout_continu" :
                self.source_N = bulk_profile[-1]
                self.source_rho_r = self.rho_r_profile[-1]


            self.data = xr.Dataset(data_vars= data_vars1, coords = {"level" : self.grid})
    
    def ajout_continu_bulk (Data) : 
        Data["concentration"][-1] = self.source_N
        Data["rho_r"][-1] = self.source_rho_r

        


### Tests and verifications ###

#start = time.time()
#initial_conds = InitialCond(nb_grid = "ARO", esp = "r", mode = "gauss", nb_classes = 10, rho_r = 100, N = 10)

#print (initial_conds.levels_boundaries)
#print (initial_conds.grid)
#print (initial_conds.data)
#print (initial_conds.data["concentration_bin_1"])
#print (initial_conds.data["level"])

#z = initial_conds.grid
#for i in range(10): 
    #c = initial_conds.data[f"concentration_bin_{i+1}"]
    #plt.plot(c,z)
#plt.show()

#d = [float(initial_conds.data[f"diameter_bin_{i+1}"]) for i in range(10)]
#x = np.linspace(1,10,10)
#plt.plot(x,d)
#plt.show()
