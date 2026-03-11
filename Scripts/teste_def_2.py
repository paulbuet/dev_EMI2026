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
from fonctions import Affichage

class model_bl_def_3:

    def __init__(self,nb_stitches,esp,r,N,CFL,time_step):
        
        # On configure le modèle avec les données d'entrée
        condi_config = IC(nb_stitches, esp, "Bulk", mode = "simple", r = r, N = N)
        self.Eq_config = Eq(esp)

        # On récupère les concentrations et le profil de rho_r et de rho
        self.grid_0 = condi_config.data
        rho = condi_config.rho

        # On récupère les valeurs de hauteurs des interfaces et on ajoute le sol

        h_interfaces = np.concatenate(([0],condi_config.levels_boundaries))

        # On s'occupe du temps

        duree_sim = 10000

        print(time_step)

        self.nb_time_step = duree_sim // time_step

        self.time_step = time_step

        # On calcule la liste des diamètre necessaire pour la sédimentation des mailles
        self.epaiss_maille = [h_interfaces[stitche+1]-h_interfaces[stitche] for stitche in range(nb_stitches)]

        dist_maille = np.array(self.Eq_config.Epaiss_to_diam(h_interfaces))

        self.diam_dist = [self.Eq_config.calcul_diametre(dist_maille[stitch],time_step) for stitch in range(len(dist_maille))]


        # On configure la vitesse max de chute d'une particule

        if CFL == "Yes":
            self.vitesse_max = self.epaiss_maille[0] / time_step
        else:
            self.vitesse_max = float('inf')

    def run(self):

        # On rentre nos profil de rho_r et de la concentration et on les enregistre pour les afficher
        
        self.grid_t = self.grid_0

        rho_r_profil = self.grid_0["rho_r"].values
        concentration_profil = self.grid_0["concentration"].values

        Liste_rho_r = [rho_r_profil]
        Liste_concentration = [concentration_profil]
        Liste_precip = [0]
        chute_conc_dt = 0
        for t in range(self.nb_time_step):

            # On calcule lambda

            Lambda = self.Eq_config.Liste_Lanbda(rho_r_profil,concentration_profil)


            concentration_profil_dt = np.zeros(len(self.epaiss_maille))
            rho_r_profil_dt = np.zeros(len(self.epaiss_maille))

            print("somme avant:", sum(rho_r_profil), sum(concentration_profil)+chute_conc_dt)
            print("rho_r=",rho_r_profil,"C=", concentration_profil)            
            

            for stitch in range(len(self.epaiss_maille)):
                
                concentration_profil_intermed = - concentration_profil[stitch] * np.array(self.Eq_config.Calcul_integrale_conc(self.diam_dist[stitch+1][:stitch+2],Lambda[stitch]))
                chute_conc = concentration_profil[stitch] * np.array(self.Eq_config.Calcul_integrale_conc([self.diam_dist[stitch+1][0],np.inf],Lambda[stitch]))

                

                if sum(np.array(self.Eq_config.Calcul_integrale_conc([0,np.inf],Lambda[stitch])))<0.99:
                    print(sum(np.array(self.Eq_config.Calcul_integrale_conc([0,np.inf],Lambda[stitch]))), "ça beugue parce que lambda = ",Lambda[stitch])
                    print("C'est possiblement parce que rho_r et C sont:",rho_r_profil[stitch], concentration_profil[stitch])
                rho_r_profil_intermed = - rho_r_profil[stitch] * np.array(self.Eq_config.Calcul_integrale_mass(self.diam_dist[stitch+1][:stitch+2],Lambda[stitch]))

                concentration_profil_intermed = np.pad(concentration_profil_intermed,(0,len(self.epaiss_maille)-stitch-1))
                rho_r_profil_intermed = np.pad(rho_r_profil_intermed,(0,len(self.epaiss_maille)-stitch-1))

                

                concentration_profil_dt += np.nan_to_num(concentration_profil_intermed)
                rho_r_profil_dt += np.nan_to_num(rho_r_profil_intermed)
                chute_conc_dt += np.nan_to_num(chute_conc)

                print("ce qui tombe:", chute_conc_dt)

            
            print("somme après", sum(rho_r_profil_dt), sum(concentration_profil_dt)+chute_conc_dt)

            Liste_rho_r.append(rho_r_profil_dt)
            Liste_concentration.append(concentration_profil_dt)

            rho_r_profil = rho_r_profil_dt
            concentration_profil = concentration_profil_dt



        Affichage.Affichage_Concentration(Liste_concentration,"concentration","Box_Lagrangien","./fig")
        Affichage.Affichage_Concentration(Liste_rho_r,"masse","Box_Lagrangien","./fig")

        
        


dodo = model_bl_def_3(40,"r",0.001,10000,"Yes",100)

print(dodo.run())