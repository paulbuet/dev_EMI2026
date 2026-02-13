# Le but de ce fichier est de mettre toutes les fonctions qu'on peut potentiellement appelés pour faire tourner nos futur codes
# Ensuite il suffira de faire un appelle du fichier et de la fonction qu nous interesse
# Exemple avec cette fonction test que j'appelle dans le code Sedimentation.py

### Imports nécessaires ###

from math import *
import time
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

### Définitions des classes/fonctions ###

class eq :
    def __init__(self,esp):
        if esp=="i":
            self.a=0.82
            self.b=2.5
            self.c=800
            self.d=1
            self.alpha=3
            self.nu=3
        if esp=='s':
            self.a=0.02
            self.b=1.9
            self.c=5.1
            self.d=0.27
            self.alpha=1
            self.nu=1
        if esp=='g':
            self.a=19.6
            self.b=2.8
            self.c=124
            self.d=0.66
            self.alpha=1
            self.nu=1
        if esp=='r':
            self.a=524
            self.b=3
            self.c=842
            self.d=0.66
            self.alpha=1
            self.nu=1
        if esp=='c':
            self.a=524
            self.b=3
            self.c=3.2e7
            self.d=2
            self.alpha=1
            self.nu=3         

    def Gamma(self, diametre, lanbda) :
        return (self.alpha/gamma(self.nu))*(lanbda**(self.alpha*self.nu))*(diametre**(self.alpha*(self.nu-1)))*np.exp(-((lanbda*diametre)**self.alpha))
   
    def G(self, p) :
        return (gamma(self.nu+p/self.alpha)/gamma(self.nu))
    
    def Lanbda(self, rho_r, N):
        return (((rho_r)/(self.a*N*self.G(self.b)))**(-1/self.b))
    
    def Masse(self, diametre):
        return (self.a*(diametre**self.b))
    
    def Vitesse(self,diametre):
        return(self.c*(diametre**self.d))
        
    def Calcul_rho_r(self, m, n):
        return np.dot(m,n)

class InitialCond :
    '''
    Class that poduces a dataset to initialize the sedimentation model, a type of initialization is used in a restricted area of the vertical gridding.
    The dataset contains the concentration profile of the specie for each bin of the sedimentation scheme.

    Attributes
    ----------
    nb_grid : number of grid cells in the model
    mode : type of intial condition used, by default the last grid cell contains a relative concentration of 1
    area_init : area in which the relative specie concentration is initialized, by default the initialization area is [30,50]
    bin_concentration : concentrations of the specie in each bin, by default only one bin is considered and its concentration is 1

    '''
    def __init__(self, nb_grid, mode = "simple", area_init = [30, 50], bin_concentration = [1]) :

        if nb_grid == "ARO" : 
            self.vertical_levels_grid = np.array([5.00148256575414, 16.7609146275979, 31.9999856716034, 50.6506387418972, 72.6448134875948, 97.9144556307367, 126.391508411840, 158.007915671418, 192.695621996127, 230.386571302649, 271.012706897240, 314.505973414367, 360.798314810016, 409.821674828922, 461.507997847950, 515.890042408161, 573.093937517738, 633.231198491812, 696.398736766220, 762.678848044096, 832.139215500780, 904.832909766711, 980.798389785428, 1060.05950095623, 1142.62547636213, 1228.49093654830, 1317.63588971212, 1410.02573054417, 1505.73432230725, 1604.93984760156, 1707.79812299161, 1814.44258334521, 1924.98428740911, 2039.51193175523, 2280.76794378598, 2407.56182949566, 2538.47269089309, 2673.47735336959, 2812.53026865253, 2955.56351459630, 3102.48360541930, 3253.19498943504, 3407.62661862652, 3565.73195622143, 3727.48898443499, 3892.90019410005, 4061.99258757885, 4234.81767907587, 4411.45149612596, 4591.99457746971, 4776.57197432233, 4965.33325159684, 5158.45249292134, 5356.12828712931, 5558.65852457596, 5766.45816179863, 5980.00284181189, 6199.82890110527, 6426.53336977965, 6660.77397670184, 6903.26915353359, 7154.79804227343, 7416.20049682756, 7688.37709125546, 7972.28912861375, 8268.07660797884, 8575.22917035383, 8893.62412699775, 9223.52646472597, 9565.58886389478, 9920.85173976858, 10290.7432801852, 10677.0795176319, 11082.0546910758, 11510.3223903101, 11966.1795634141, 12454.7528314299, 12982.0128918015, 13554.6557230282, 14180.1025765563, 14866.4999945363, 15627.4731252487, 16485.0471290191, 17467.8152385289, 18610.9387773901, 19955.6070022098, 21550.0160473017, 23450.1477297912, 34461.1214067193])
        else :
            height_grid = 12e3
            self.vertical_levels_grid = np.linspace(1/nb_grid * height_grid, height_grid, nb_grid)
        nb_levels = len (self.vertical_levels_grid)

        concentration_profile = [0 for i in range (nb_levels)]
        if mode == "simple" :
            concentration_profile [-1] = 1
        if mode == "gauss" :
            concentration_profile [area_init[0],area_init[1]] = [1 for i in range(len[area_init])]

        bin_profile = [np.array(concentration_profile) * bin_concentration[i] for i in range(len(bin_concentration))]
        data_vars = {f"concentration_bin{ind_bin+1}" : ("level", bin_profile[ind_bin]) for ind_bin in range(len(bin_concentration))}
        self.data = xr.Dataset(data_vars=data_vars, coords = {"level" : self.vertical_levels_grid})

initial_conds = InitialCond(nb_grid = 50, bin_concentration = [1,2,3,4,5])
print (initial_conds.vertical_levels_grid)
print (initial_conds.data)
print (initial_conds.data["concentration_bin1"])
