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

### Classes and functions definition ###

class Eq :
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
            self.C=5
            self.x=1
        if esp=='g':
            self.a=19.6
            self.b=2.8
            self.c=124
            self.d=0.66
            self.alpha=1
            self.nu=1
            self.C=5e5
            self.x=-0.5
        if esp=='r':
            self.a=524
            self.b=3
            self.c=842
            self.d=0.66
            self.alpha=1
            self.C=8e6
            self.x=-1            
            self.nu=1
        if esp=='c':
            self.a=524
            self.b=3
            self.c=3.2e7
            self.d=2
            self.alpha=1
            self.nu=3         

    def Gamma(self, diametre, lam) :
        return (self.alpha/gamma(self.nu))*(lam**(self.alpha*self.nu))*(diametre**(self.alpha*self.nu-1))*np.exp(-((lam*diametre)**self.alpha))
   
    def Gamma_fois_masse(self, diametre, lam) :
        return self.Gamma(diametre, lam)*(self.a*(diametre**self.b))

    def G(self, p) :
        return (gamma(self.nu+p/self.alpha)/gamma(self.nu))
    
    def Lanbda_2M(self, rho_r, N):
        return (((rho_r)/(self.a*N*self.G(self.b)))**(-1/self.b))
    
    def Lanbda_1M(self, rho_r):
        return (rho_r/(self.a*self.C*self.G(self.b)))**(1/(self.x-self.b))
    
    def Liste_Lanbda(self, rho_r, Liste_N):
        Liste_N=np.array(Liste_N)
        Ga=self.G(self.b)
        Liste_N = np.array([np.nan if x == 0 else x for x in Liste_N])
        Liste_lanbda=(rho_r/(Liste_N*self.a*Ga))**(-1/self.b)
        Liste_lanbda = np.array([0 if x == np.nan else x for x in Liste_lanbda])
        return Liste_lanbda
    
    def Liste_Dmin_Dmax(self, Liste_Lanbda):
        Liste_Lanbda=np.array(Liste_Lanbda)
        indices_nan = np.array([i for i, x in enumerate(Liste_Lanbda) if np.isnan(x)])
        Liste_Lanbda_sans_nan = np.array([x for x in Liste_Lanbda if not np.isnan(x)])
        Liste_Dm = [self.Dmin_Dmax(elem)[::-1] for elem in Liste_Lanbda_sans_nan]
        Liste_Dm_avec_nan = Liste_Dm.copy()
        for i in indices_nan:
            Liste_Dm_avec_nan.insert(i, (np.nan, np.nan))
        Liste_Dm_avec_nan=np.array(Liste_Dm_avec_nan)
        return Liste_Dm_avec_nan

    def Liste_Vitesse_Concentration(self, Liste_Lanbda):

        Liste_Dm=self.Liste_Dmin_Dmax(Liste_Lanbda)
        Vitesse= self.c*(Liste_Dm**self.d)
        return Vitesse

    def Gamma_Masse(self, M, lam):
        return (((M/self.a)**((1/self.b)-1))*self.Gamma(((M/self.a)**(1/self.b)), lam)/(self.a*self.b))

    def Vitesse_Masse(self, M):
        return (self.c*((M/self.a)**(self.d/self.b)))

    def Massemin_Massemax(self, lam):
        def Fct(M):
            return quad(self.Gamma_Masse, 0, M, args=(lam))[0]
 
        M_high = 1/lam
        while Fct(M_high) < 0.999:
            M_high *= 2
        Massemin = brentq(lambda M: Fct(M) - 0.01, 0, M_high)
        Massemax = brentq(lambda M: Fct(M) - 0.999, 0, M_high)
 
        return Massemin, Massemax
        
    def Liste_Massemin_Massemax(self, Liste_Lanbda):
        Liste_Lanbda=np.array(Liste_Lanbda)
        indices_nan = np.array([i for i, x in enumerate(Liste_Lanbda) if np.isnan(x)])
        Liste_Lanbda_sans_nan = np.array([x for x in Liste_Lanbda if not np.isnan(x)])
        Liste_M=[self.Massemin_Massemax(elem)[::-1] for elem in Liste_Lanbda_sans_nan]
        Liste_M_avec_nan = Liste_M.copy()
        for i in indices_nan:
            Liste_M_avec_nan.insert(i, (np.nan, np.nan))
        Liste_M_avec_nan=np.array(Liste_M_avec_nan)
        return Liste_M_avec_nan
    
    def Liste_Vitesse_Masse(self, Liste_lanbda):
        Liste_M=np.array(self.Liste_Massemin_Massemax(Liste_lanbda))
        Vitesse=self.Vitesse_Masse(Liste_M)
    

        return Vitesse

    def Masse(self, diametre):
        return (self.a*(diametre**self.b))
    
    def Vitesse(self,diametre):

        return (self.c*(diametre**self.d))
    
    def Calcul_rho_r(self, m, n):
        return np.dot(m,n)
    
    def Dmin_Dmax(self, lam):
        def F(D):
            return quad(self.Gamma, 0, D, args=(lam))[0]
 
        D_high = 1.0 / lam # échelle naturelle
        while F(D_high) < 0.999:
            D_high *= 2

        Dmin = brentq(lambda D: F(D) - 0.01, 0, D_high)
        Dmax = brentq(lambda D: F(D) - 0.999, 0, D_high)

        print("les infos que je veux : ", lam, Dmin, Dmax)

        return Dmin, Dmax
    
    def Classe_D(self, nb_classes, Dmin, Dmax, N, lam):
        Result=[]
        Intervalle=(Dmax-Dmin)/nb_classes
        for i in range(nb_classes):
            Di=(1+2*i)*Intervalle/2 + Dmin
            
            P_i=quad(self.Gamma, Dmin+i*Intervalle, Dmin+(i+1)*Intervalle, args=(lam))[0]/0.98
            Ni=N*P_i
            Result.append([Di, Ni]) #Liste de deux paramètres : diamètre moyen, quantité associé par rapport au nombre total de particule.
        return Result
    
    def Liste_rho_r(self, Liste_nb,nb_part,masse,dz):
        return np.array(Liste_nb)*masse/nb_part/dz
        

    
    def Calcul_Masse_Tot(self,liste_rho_r, hauteur_col):
        liste_rho_r_sans_nan=np.array([x for x in liste_rho_r if not np.isnan(x)])
        rho_tot=np.sum(liste_rho_r_sans_nan)
        Masse_tot=rho_tot*hauteur_col
        return Masse_tot
    
    def Epaiss_to_diam(self,h_interface):
        
        nb_interf = len(h_interface)

        tab_diam = [np.pad([h_interface[stitch_dep] - h_interface[stitch_arr] for stitch_arr in range(stitch_dep+1)],(0,nb_interf-stitch_dep-1))for stitch_dep in range(nb_interf)]
        
        return tab_diam
    
    def Diameters_conc(self, p_list, lambda_list, Dmin=0, Dmax=1, nD=20000):

        D = np.linspace(Dmin, Dmax, nD)

        # densité pour tous les lambda
        P = self.Gamma(D[:, None], lambda_list[None, :])

        # position du maximum pour chaque lambda
        imax = np.argmax(P, axis=0)

        result = []

        for i in range(len(lambda_list)):
            p = p_list[i]

            # branche gauche
            Dl = D[:imax[i]+1]
            Pl = P[:imax[i]+1, i]

            # branche droite
            Dr = D[imax[i]:]
            Pr = P[imax[i]:, i]

            D_left = np.interp(p, Pl, Dl)
            D_right = np.interp(p, Pr[::-1], Dr[::-1])

            result.append([D_left, D_right])

        return result

    def calcul_liste_vitesse_conc(self, liste_lam):
        liste_lam=np.array(liste_lam)
        Liste_Dmax=(((self.alpha * self.nu - 1) / self.alpha)**(1/self.alpha)) / liste_lam
        Liste_prob_max=self.Gamma(Liste_Dmax, liste_lam)
        un_pourcent_prob_max=0.01*Liste_prob_max
        return self.Diameters_conc(un_pourcent_prob_max, liste_lam)
    
    def Liste_Lanbda_1_mom(self, liste_rho_r):
        liste_rho_r=np.array(liste_rho_r)
        Gam=self.G(self.b)
        l_lam_1_mom=(liste_rho_r/(self.a*self.C*Gam))**(1/(self.x-self.b))
        return l_lam_1_mom
    


    #Fonction pour le nouveau schéma de Sébastien Riette

    def calcul_diametre(self, liste_h, dt):
        liste_h=np.array(liste_h)
        liste_d=(liste_h/(self.c*dt))**(1/self.d)
        return liste_d
    
    def Calcul_integrale_conc(self, liste_d, lam):
        integrale=[]
        for i in range(len(liste_d)-1):
            integrale_entre_d_i_et_d_i_plus_1, err = integrate.quad(self.Gamma, liste_d[i], liste_d[i+1], args=(lam,))
            integrale.append(integrale_entre_d_i_et_d_i_plus_1)
        return integrale
    
    def Calcul_integrale_mass(self, liste_d, lam):
        integrale=[]
        for i in range(len(liste_d)-1):
            integrale_entre_d_i_et_d_i_plus_1, err = integrate.quad(self.Gamma_Masse, liste_d[i], liste_d[i+1], args=(lam,))
            integrale.append(integrale_entre_d_i_et_d_i_plus_1)
        return integrale
    
    def content_to_conc_bulk_1M(self,rho_r):
        return self.C * self.Lanbda_1M(rho_r)**(self.x)
    


"""

eq_rain=Eq("r")

print(eq_rain.Liste_Lanbda_1_mom([263, 241, 72]))


"""




class Affichage :

    def Affichage_Concentration(Concentration, typ, model, path_fig): #type="concentration" ou "masse"
        Temps_simu=len(Concentration)
        nb_boites=len(Concentration[0])
        Concentration=np.array(Concentration)
        Transpose=Concentration.T
        plt.figure(figsize=(Temps_simu, nb_boites))
        orig_map=plt.cm.get_cmap('gist_ncar')
        reversed_map = orig_map.reversed()
        plt.pcolormesh(Transpose,cmap=reversed_map)
        plt.title(f"Evolution de la {typ} de particules en fonction du temps", fontsize=22)
        plt.xlabel("Temps", fontsize=18)
        plt.ylabel("Mailles du modèle", fontsize=18)
        plt.colorbar()
        file_location = Path(path_fig) / Path(model) / Path(typ) 
        plt.savefig(str(file_location))


    def Affichage_Precipitation(Precip, model, path_fig):
        Precip=np.array(Precip)
        liste=np.zeros(len(Precip))
        Cumul=[]

        for i in range(len(Precip)):
            liste[i]=1
            Cumul.append(np.dot(Precip, liste))
        fig, ax1 = plt.subplots()
        ax = plt.gca()
        ax.xaxis.set_major_locator(MultipleLocator(1))
        ax.yaxis.set_major_locator(MultipleLocator(1))
        ax.set_xlim(left=0)
        x_max = len(Precip)-1
        ticks = np.linspace(0, x_max, 11)
        ax.set_xticks(ticks)
        ax.set_xticklabels([f"{t:.2f}" for t in ticks])
        ax.yaxis.set_major_locator(MaxNLocator(10))

        time=np.linspace(1, len(Precip), len(Precip))
        time2 = time -(time[1]-time[0])/2
        ax1.bar(time2, Precip, color="blue", label="Cumul")
        ax1.set_xlabel("temps", fontsize=18)
        ax1.set_ylabel("Précip par pas de temps", fontsize=18)
        ax2 = ax1.twinx()
        ax2.plot(time, Cumul, '--', color="red", label="Précip par pas de temps")
        ax2.set_ylabel('Cumul', fontsize=18)
        plt.title("Evolution des précipitations par pas de temps et cumulée", fontsize=22)
        plt.grid(axis='x', which='major', markevery=[1,2,3],lw=2, ls=':')
        fig.legend(loc=2)
        file_location = Path(path_fig) / Path(model) / "Précipitations.png"
        plt.savefig(file_location)

    def Afficher () :
        plt.show()





def gaussiennem (m, sigma,h):
        return  1 * exp (-1/2 * ((h-m)/sigma)**2)


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
    

    def __init__(self, nb_grid, esp, types, mode = "simple", Hmax = 5000, sigma = 2000, nb_classes = 10, r = 0.001) :

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
                relative_profile = [gaussiennem(Hmax, sigma, self.grid[i]) for i in range(nb_levels)]

            # Vérification : Affichage des profils relatifs (à décommenter)
            # x = [self.grid[i] for i in range(nb_levels)]
            # y = relative_profile
            # plt.plot(x,y)
            # plt.show()


            self.r_profile = np.array(relative_profile)*r

            self.rho_r_profile,self.rho_profile = profil_rho_r().calcul(self.grid,self.r_profile)

            eq=Eq(esp)
            lam = eq.Lanbda_1M(self.rho_profile)
            print (lam)
            # print(f"rho {self.rho}")

            N_profile = eq.content_to_conc_bulk_1M(self.rho_profile)
            # N_profile = self.rho * r * (lam ** eq.b)/ (eq.a * eq.G(eq.b))
            print(N_profile)

            dminmax = [eq.Dmin_Dmax(lam[ind_level]) for ind_level in range(nb_levels)]
            print (dminmax)
            levels_bin_splittings= [eq.Classe_D (nb_classes, dminmax[ind_level][0], dminmax[ind_level][1], N_profile[ind_level], lam[ind_level]) for ind_level in range(nb_grid)] # List of lists containing the nb_grid classifications
            # print (f" levels_bins {levels_bins}")


            # Vérification : Affichage des profils des bins (à décommenter)
            # x = [self.grid[i] for i in range(nb_levels)]
            # y = self.bin_concentration[0][1]*np.array(relative_profile)
            # print(self.bin_concentration)
            # plt.plot(x,y)
            # plt.show()
            # print (bins)


            bins_concentrations_profiles = [[splitting[ind_bin][1] for splitting in levels_bin_splittings] for ind_bin in range(nb_classes)] # computinng of the n bin profiles
            print(bins_concentrations_profiles)
            data_vars1 = {f"concentration_bin_{ind_bin+1}" : ("level", bins_concentrations_profiles[ind_bin]) for ind_bin in range(nb_classes)}

            self.diameters = [bins[ind_bin][0] for ind_bin in range(nb_classes)]
            data_vars2 = {f"diameter_bin_{ind_bin+1}" : self.diameters[ind_bin] for ind_bin in range(nb_classes)} # addition of the diameters
            
            bins_rho_r_profiles = [bins_concentrations_profiles[ind_bin] * eq.Masse(np.array(bins[ind_bin][0])) for ind_bin in range(nb_classes)]
            
            data_vars3 = {f"rho_r_bin_{ind_bin+1}" : bins_rho_r_profiles[ind_bin][:] for ind_bin in range(nb_classes)}

            data_vars1.update(data_vars2)
            data_vars1.update(data_vars3)

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

            print (self.data["rho_r_bin_1"])

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

            if mode == "simple" :
                relative_profile = [ 0 for i in range(nb_levels)]
                relative_profile[-1] = 1
                #relative_profile[-2] = 1/2


            if mode == "gauss" :
                relative_profile = [gaussiennem(Hmax, sigma, self.grid[i]) for i in range(nb_levels)]

            self.rho_r_profile,self.rho = profil_rho_r().calcul(self.grid,relative_profile*r)

            eq=Eq(esp)

            relative_profile = eq.content_to_conc_bulk_1M(self.rho_r_profile)

            bulk_profile = np.array(relative_profile)   # computinng of the n bin profiles
            data_vars1 = {"concentration" : ("level", bulk_profile), "rho_r": ("level",self.rho_r_profile)}


            self.data = xr.Dataset(data_vars= data_vars1, coords = {"level" : self.grid})


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
