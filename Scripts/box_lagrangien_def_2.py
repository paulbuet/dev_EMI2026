#!/usr/bin/env python3

# Copyright (c) Météo France (2025-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

# On utilisera dans la description du Box-Lagrangien les notations de Kato(1995)

# On importe les librairies nécessaires
import numpy as np
import xarray as xr
import xarray_regrid
from fipy.meshes.uniformGrid1D import UniformGrid1D
from fonctions import Affichage
# On importe les classes nécessaires

from fonctions  import InitialCond as IC
from fonctions import Eq

class Model_bl_def_2:

    def __init__(self,nb_stitches,esp,r,N,CFL,time_step, duree_simu):
        """
        nb_stitches (int): nombre de mailles
        esp (str): nom de la particule que l'on sédimente
        r (float): rapport de mélange max 
        N (int): nombre de particules
        """

        # On configure le modèle avec les données d'entrée
        condi_config = IC(nb_stitches, esp, "Bulk", mode = "simple", r = r, N = N)
        self.Eq_config = Eq(esp)

        # On récupère les concentrations et le profil de rho_r et de rho
        self.grid_0 = condi_config.data
        rho = condi_config.rho

        # On calcule le facteur mutiplicatif pour la vitesse
        rho_form = sorted(np.concatenate(([1.225],np.concatenate((rho,rho)))))[:-1]
        self.coeff_mult = (rho_form/rho_form[0])**0.4


        # On récupère les valeurs de hauteurs des interfaces et on ajoute le sol

        h_interfaces = np.concatenate(([0],condi_config.levels_boundaries))

        # On calcule l'épaisseur de chaque maille
        self.epaiss_maille = [h_interfaces[stitche+1]-h_interfaces[stitche] for stitche in range(nb_stitches)]

        # On s'occupe ici de définir les variables temporelles en secondes

        duree_simu = 500
        self.duree_simu = duree_simu
        self.time_step = time_step

        self.nb_time_step = self.duree_simu // time_step

        # On configure la vitesse max de chute d'une particule

        #Il y a pas de problème de cond CFL ici ???

        if CFL == "Yes":
            self.vitesse_max = epaiss_maille[0] / self.time_step
        else:
            self.vitesse_max = float('inf')

        # On formate les hauteurs d'interfaces afin de créer les mailles déformées

        #Pareil je comprend pas pourquoi tu fais ça @Paul_buet ???

        self.h_interfaces_form = sorted(np.concatenate((h_interfaces,h_interfaces)))[1:-1]
        self.z_top = self.h_interfaces_form[-1]


    def regridage(self,h_interface,stitche,variable):

        # On note la quantité de variable de départ
        qt_dep = self.grid_t[variable].values[stitche]

        # On prends les hauteurs de interfaces de la maille déformée et on calcule son épaisseur
        h_interface_sup = h_interface[stitche*2+1]
        h_interface_inf = h_interface[stitche*2]

        dz_def = h_interface_sup - h_interface_inf

        # On calcule le centre de la maille déformée et on crée une grille de maille déformées régulières
        centre_stitche = (h_interface_sup + h_interface_inf)/2

        # On calcule le nombre de maille à rajouter pour couvrir toute la maille initiale
        nb_stitche_sup = int((self.z_top-h_interface_sup) // dz_def + 2)
        nb_stitche_inf = int(h_interface_inf // dz_def +2)

        # On crée alors la nouvelle grille, l'origine est une interface, pas le centre d'une maille
        grille_sup = UniformGrid1D(dx = dz_def,nx = nb_stitche_sup,  origin = (h_interface_sup,))
        grille_inf = UniformGrid1D(dx = -dz_def,nx = nb_stitche_inf,  origin = (h_interface_sup,))

        grille_tot = np.concatenate((np.flip(grille_inf.cellCenters[0]) ,grille_sup.cellCenters[0]))

        # On remets la valeur de la variable dans la bonne maille déformée de la nouvelle grille et
        data_var = np.zeros(len(grille_tot)-1)
        data_var = np.insert(data_var,nb_stitche_inf-1,qt_dep)
        
        # rgo_r étant une variable intensive, on la change avec la déformation
        if variable == "rho_r":
            data_var *= self.epaiss_maille[stitche]/dz_def

        # On crée le dataset avec les mailles déformée et leur valeur de variables associées
        dataset_def = xr.Dataset(data_vars= dict(variable = (("level"),data_var)), coords = {"level" : grille_tot})

        # On regrid sur la maille non déformée
        grille_regrid = dataset_def.regrid.conservative(self.grid_0)

        # On densifie les variables regridder
        grille_regrid[variable] = grille_regrid["variable"].copy(data=grille_regrid["variable"].data.todense())

        # De la même manière que tout à l'heure mais à l'inverse, on a rho une valeur intensive
        if variable == "rho_r":
            grille_regrid[variable] *= dz_def/self.epaiss_maille[stitche]

        # On applique ensuite un correctif pour la conservation, il n'est pas le même selon si on a touché le sol
        if h_interface_inf<0:

            proportion_tombe = -h_interface_inf/dz_def

        else:

            proportion_tombe = 0

        qt_voulu_fin = qt_dep*(1-proportion_tombe)
        qt_obtenu = sum(grille_regrid[variable].values)

        # Les mailles étant régulières, on a rho_r qui se comporte comme une valeur extensive donc conservative

        grille_regrid[variable].values *= (qt_voulu_fin)/qt_obtenu
        grille_regrid = grille_regrid.drop_vars("variable")

        return grille_regrid,proportion_tombe*qt_dep

    def init_grid_dt(self):

        rho_r = np.zeros(len(self.grid_0["level"].values))

        concentration = np.zeros(len(self.grid_0["level"].values))
        data_var = {"concentration" : ("level", concentration), "rho_r": ("level",rho_r)}

        grid_dt = xr.Dataset(data_vars= data_var, coords = {"level" : self.grid_0["level"].values})

        return grid_dt



        






    def run(self):
        # On rentre nos profil de rho_r et de la concentration et on les enregistre pour les afficher

        self.grid_t = self.grid_0

        rho_r_profil_conc = self.grid_0["rho_r"].values
        concentration_profil = self.grid_0["concentration"].values

        rho_r_profil_mass = self.grid_0["rho_r"].values
        masse_profil = self.grid_0["concentration"].values

        Liste_rho_r = [rho_r_profil_mass]
        Liste_concentration = [concentration_profil]
        Liste_precip = [0]

        for t in range(self.nb_time_step):

            # On initialise à chaque pas de temps un nouvel état
            grid_dt = self.init_grid_dt()
            somme_tomb = 0
            
            # On calcule les lambda
            Lambda_concentration = self.Eq_config.Liste_Lanbda(rho_r_profil_conc,concentration_profil)
            Lambda_masse = self.Eq_config.Liste_Lanbda(rho_r_profil_mass,masse_profil)

            # On calcule les vitesses
            vitesse_concentration = self.Eq_config.Liste_Vitesse_Concentration(Lambda_concentration)
            vitesse_masse = self.Eq_config.Liste_Vitesse_Masse(Lambda_masse)

            # On formate pour appliquer aux interfaces
            vitesse_concentration = np.nan_to_num(vitesse_concentration.flatten())
            vitesse_masse = np.nan_to_num(vitesse_masse.flatten())

            # On applique le coefficient multiplicateur
            vitesse_concentration *= self.coeff_mult
            vitesse_masse *= self.coeff_mult

            # On seuil par la vitesse maximale
            vitesse_concentration = np.clip(vitesse_concentration,0,self.vitesse_max)
            vitesse_masse = np.clip(vitesse_masse,0,self.vitesse_max)

            # On calcule la hauteur des interfaces des mailles déformées

            h_interface_def_conc = self.h_interfaces_form - vitesse_concentration * self.time_step
            h_interface_def_mass = self.h_interfaces_form - vitesse_masse * self.time_step

            # On parcourt les mailles où il y a des particules
            for stitche in np.where(np.logical_or(Lambda_masse > 0, Lambda_concentration > 0))[0]:
                
                nouv_conc_prof,qt_tomb_conc = self.regridage(h_interface_def_conc,stitche,"concentration")
                nouv_mass_prof,qt_tomb_mass = self.regridage(h_interface_def_mass,stitche,"rho_r")

                print(nouv_conc_prof)

                grid_dt["concentration"] = grid_dt["concentration"] + nouv_conc_prof["concentration"]
                grid_dt["rho_r"] = grid_dt["rho_r"] + nouv_mass_prof["rho_r"]

                somme_tomb += qt_tomb_mass

            self.grid_t = grid_dt

            rho_r_profil_conc = self.grid_t["rho_r"].values
            concentration_profil = self.grid_t["concentration"].values

            rho_r_profil_mass = self.grid_t["rho_r"].values
            masse_profil = self.grid_t["concentration"].values
        
            Liste_rho_r.append(self.grid_t["rho_r"].values)
            Liste_concentration.append(self.grid_t["concentration"].values)
            Liste_precip.append(somme_tomb)

            print(grid_dt)
        return Liste_concentration,Liste_precip,Liste_rho_r 



            










dodo = Model_bl_def_2(20,'r',0.001,60,'No',100, 2000)

profil = dodo.run()
concentration_formate = np.array(profil[0])

mass_form = np.array(profil[2])
Affichage.Affichage_Concentration(concentration_formate, "concentration", "Box_Lagrangien", )
Affichage.Affichage_Concentration(mass_form, "masse", "Box_Lagrangien")
Affichage.Affichage_Precipitation(profil[1], "Box_Lagrangien")