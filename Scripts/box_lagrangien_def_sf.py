#!/usr/bin/env python3

# Copyright (c) Météo France (2025-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

# On utilisera dans la description du Box-Lagrangien les notations de Kato(1995)

# On importe les librairies nécessaires
import numpy as np
from tqdm import tqdm

# On importe les classes nécessaires
from condi_init  import InitialCond as IC
from equations import Eq
from formatage import Formatage 


class model_bl_def_sf:

    def __init__(self,nb_stitches,time_step,esp,r,mode):
        
        # On configure le modèle avec les données d'entrée
        condi_config = IC(nb_stitches, esp, "Bulk", mode = mode, r = r)
        self.Eq_config = Eq(esp)

        # On récupère les concentrations et le profil de rho_r et de rho
        self.grid_0 = condi_config.data
        rho = condi_config.rho

        # On récupère les valeurs de hauteurs des interfaces et on ajoute le sol

        h_interfaces = np.concatenate(([0],condi_config.levels_boundaries))
        self.h_interfaces = h_interfaces

        # On s'occupe du temps

        duree_sim = 2000


        self.nb_time_step = duree_sim // time_step

        self.time_step = time_step

        # On calcule la liste des diamètre necessaire pour la sédimentation des mailles
        self.epaiss_maille = np.array([h_interfaces[stitche+1]-h_interfaces[stitche] for stitche in range(nb_stitches)])


        dist_maille = np.array(Formatage(esp).Epaiss_to_diam(h_interfaces))

        for i in range(1,len(dist_maille)):

            dist_maille[i,:i] -= dist_maille[i,i-1]/2
                       
        self.diam_dist = [self.Eq_config.calcul_diametre(dist_maille[stitch],time_step) for stitch in range(len(dist_maille))]

    def run(self):

        # On rentre nos profil de rho_r et de la concentration et on les enregistre pour les afficher
        
        self.grid_t = self.grid_0

        rho_r_profil = self.grid_0["rho_r"].values
        concentration_profil = self.grid_0["concentration"].values


        Liste_rho_r = [rho_r_profil]
        Liste_concentration = [concentration_profil]
        Liste_precip = [0]
        chute_conc_dt = 0
        chute_mass_dt = 0
        chute_tot =[]

        print (" ")
        print ("---------------------------------------")

        format_affichage_time = "{l_bar}|{bar}| {n_fmt}/{total_fmt} pas de temps | temps écoulé : {elapsed} < temps restant {remaining} |" # definition of the format of the tqdm bar
        format_affichage_level = "-> {desc} |{bar}| {n_fmt}/{total_fmt} niveaux  | {elapsed}<{remaining} |                                                                                           "

        # On calcule lambda
        Lambda = self.Eq_config.Liste_Lanbda(rho_r_profil,concentration_profil)

        for t_time in tqdm(range(self.nb_time_step), bar_format = format_affichage_time, desc = f"Avancement total Box Lagrangien déformable Step_Forward : ", position = 0, colour = 'blue'):
            
            chute_tot_int=0
        

            concentration_profil_dt = np.zeros(len(self.epaiss_maille))
            rho_r_profil_dt = np.zeros(len(self.epaiss_maille))
            

            for stitch_dep in tqdm(np.where(concentration_profil!=0)[0], bar_format = format_affichage_level, desc = f"Calculs t = {t_time * self.time_step} / {self.nb_time_step * self.time_step} s : ", leave = False, colour = "green"):
                
                
                conc_profil_int = []
                rho_r_profil_intermed = []
                for stitch_arr in range(stitch_dep+1):
                    conc_sed_i_to_j = self.Eq_config.calcul_maille_arrivee(self.h_interfaces[stitch_arr], self.h_interfaces[stitch_arr+1], self.h_interfaces[stitch_dep], self.h_interfaces[stitch_dep+1], self.grid_0["concentration"].values[stitch_dep], "concentration", self.time_step*(t_time+1), Lambda[stitch_dep])
                    rho_r_sed = self.Eq_config.calcul_maille_arrivee(self.h_interfaces[stitch_arr], self.h_interfaces[stitch_arr+1], self.h_interfaces[stitch_dep], self.h_interfaces[stitch_dep+1], self.grid_0["concentration"].values[stitch_dep], "masse", self.time_step*(t_time+1), Lambda[stitch_dep])
                    chute_int = self.Eq_config.calcul_maille_arrivee(-10000, self.h_interfaces[0], self.h_interfaces[stitch_dep], self.h_interfaces[stitch_dep+1], self.grid_0["concentration"].values[stitch_dep], "masse", self.time_step*(t_time+1), Lambda[stitch_dep])*10000/self.epaiss_maille[stitch_dep]
                    conc_profil_int.append(conc_sed_i_to_j)
                    rho_r_profil_intermed.append(rho_r_sed)
                    chute_tot_int+= chute_int
                    # print(stitch_arr)
                    # print("test, conc_sed_i_to_j : ", conc_sed_i_to_j)

                rho_r_profil_intermed = np.pad(rho_r_profil_intermed,(0,len(self.epaiss_maille)-stitch_dep-1))
                conc_profil_int = np.pad(conc_profil_int, (0, len(self.epaiss_maille)-stitch_dep-1))

                concentration_profil_dt += np.nan_to_num(conc_profil_int)
                rho_r_profil_dt += np.nan_to_num(rho_r_profil_intermed)

                #print("test, conc_profil_int : ", conc_profil_int, len(conc_profil_int))
                
                chuté = concentration_profil[stitch_dep]


                #concentration_profil_intermed =- concentration_profil[stitch] * np.array(self.Eq_config.Calcul_integrale_conc(self.diam_dist[stitch+1][:stitch+2],Lambda[stitch])) #self.Eq_config.calcul_maille_arrivee(h1, h2, h3, h4, concentration_profil[stitch], "concentration", self.time_step, Lambda[stitch])
                #chute_conc = concentration_profil[stitch] * np.array(self.Eq_config.Calcul_integrale_conc([self.diam_dist[stitch+1][0],np.inf],Lambda[stitch]))

                #rho_r_profil_intermed = - concentration_profil[stitch] * np.array(self.Eq_config.Calcul_integrale_mass(self.diam_dist[stitch+1][:stitch+2],Lambda[stitch]))
                
                #chute_mass = concentration_profil[stitch] * np.array(self.Eq_config.Calcul_integrale_mass([self.diam_dist[stitch+1][0],np.inf],Lambda[stitch]))

                #concentration_profil_intermed = np.pad(concentration_profil_intermed,(0,len(self.epaiss_maille)-stitch-1))
                #print("test 2 : ", t, stitch, concentration_profil_intermed)
                #rho_r_profil_intermed = np.pad(rho_r_profil_intermed,(0,len(self.epaiss_maille)-stitch-1))


                #concentration_profil_dt += np.nan_to_num(concentration_profil_intermed)
                #rho_r_profil_dt += np.nan_to_num(rho_r_profil_intermed)
                #chute_conc_dt += np.nan_to_num(chute_conc[0])
                #chute_mass_dt += np.nan_to_num(chute_mass[0])

            chute_tot.append(chute_tot_int)
            Liste_rho_r.append(rho_r_profil_dt)
            Liste_concentration.append(concentration_profil_dt)
            Liste_precip.append(chute_mass_dt)

            rho_r_profil = rho_r_profil_dt
            concentration_profil = concentration_profil_dt


        return Liste_concentration, chute_tot, Liste_rho_r

        
        
