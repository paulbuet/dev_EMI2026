#!/usr/bin/env python3

# Copyright (c) Météo France (2025-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

# On utilisera dans la description du Box-Lagrangien les notations de Kato(1995)

# On importe les librairies nécessaires
import numpy as np

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

        # On récupère les concentrations et le profil de rho_r
        self.grid_0 = condi_config.data

        # On récupère les valeurs de hauteurs des interfaces

        h_interfaces = condi_config.levels_boundaries

        # On calcule l'épaisseur de chaque maille
        epaiss_maille = [h_interfaces[stitche+1]-h_interfaces[stitche] for stitche in range(nb_stitches-1)]

        # On s'occupe ici de définir les variables temporelles en secondes

        self.duree_simu = duree_simu

        self.nb_time_step = self.duree_simu // time_step

        # On configure la vitesse max de chute d'une particule

        #Il y a pas de problème de cond CFL ici ???

        if CFL == "Yes":
            vitesse_max = epaiss_maille / time_step
        else:
            vitesse_max = float('inf')

        # On formate les hauteurs d'interfaces afin de créer les mailles déformées

        #Pareil je comprend pas pourquoi tu fais ça @Paul_buet ???

        h_interfaces_form = np.concatenate(([0],sorted(np.concatenate((h_interfaces,h_interfaces)))[:-1]))


    def run(self):
        # On rentre nos profil de rho_r et de la concentration et on les enregistre pour les afficher
        rho_r_profil = self.grid_0["rho_r"].values
        concentration_profil = self.grid_0["concentration"].values

        Liste_rho_r = [rho_r_profil]
        Liste_concentration = [concentration_profil]

        Lambda = self.Eq_config.Liste_Lanbda(rho_r_profil,concentration_profil)   # On trouve lambda

        # On calcule les vitesses
        vitesse_concentration = self.Eq_config.Liste_Vitesse_Concentration(Lambda) 
        vitesse_masse = self.Eq_config.Liste_Vitesse_Masse(Lambda) 

        print(vitesse_concentration,vitesse_masse)


        print(Lambda)






dodo = model_bl_def_2(20,'i',0.001,60,'No',100)

dodo.run()