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

    def __init__(self,nb_mailles,delta_t,esp,r,mode,duree_sim):

        self.nb_mailles = nb_mailles
        self.delta_t = delta_t
        self.esp = esp

        self.nb_step = duree_sim // self.delta_t

        # On configure les classes extèrieures
        self.Eq_config = Eq(esp)
        
        # On configure le modèle avec les données d'entrée
        condi_config = IC(nb_mailles, esp, "Bulk", mode = mode, r = r)
        

        # On récupère les concentrations et le profil de rho_r
        self.grid_0 = condi_config.data
        self.hauteur_interf = np.concatenate(([0],condi_config.levels_boundaries))

        # On calcule la profil des diamètre necessaire pour la sédimentation des mailles
        self.epaiss_maille = np.array([self.hauteur_interf[stitche+1]-self.hauteur_interf[stitche] for stitche in range(nb_mailles)])

        """
        dist_maille = np.array(Formatage(esp).Epaiss_to_diam(h_interfaces))

        for i in range(1,len(dist_maille)):

            dist_maille[i,:i] -= dist_maille[i,i-1]/2
                       
        self.diam_dist = [self.Eq_config.calcul_diametre(dist_maille[stitch],delta_t) for stitch in range(len(dist_maille))]
        """
    def run(self):

        # On rentre nos profil de rho_r et de la concentration et on les enregistre pour les afficher

        profil_rho_r_0 = self.grid_0["rho_r"].values
        profil_concentration_0 = self.grid_0["concentration"].values


        profil_rho_r = [profil_rho_r_0]
        profil_concentration = [profil_concentration_0]
        liste_precip = [0]

        print (" ")
        print ("---------------------------------------")

        format_affichage_time = "{l_bar}|{bar}| {n_fmt}/{total_fmt} pas de temps | temps écoulé : {elapsed} < temps restant {remaining} |" # definition of the format of the tqdm bar
        format_affichage_level = "-> {desc} |{bar}| {n_fmt}/{total_fmt} niveaux  | {elapsed}<{remaining} |                                                                                           "

        # On calcule lambda
        Lambda = self.Eq_config.Liste_Lanbda(profil_rho_r_0,profil_concentration_0)

        self.lam_init = Lambda[-1]
        self.conc_tot_init = profil_concentration_0[-1]

        for t_time in tqdm(range(self.nb_step), bar_format = format_affichage_time, desc = f"Avancement total Box Lagrangien déformable Step_Forward : ", position = 0, colour = 'blue'):
            
            # On initialise les champs au temps dt
            precip_dt =0
            profil_concentration_dt = np.zeros(len(self.epaiss_maille))
            profil_rho_r_dt = np.zeros(len(self.epaiss_maille))
            

            for maille_dep in tqdm(np.where(profil_concentration_0!=0)[0], bar_format = format_affichage_level, desc = f"Calculs t = {t_time * self.delta_t} / {self.nb_step * self.nb_step} s : ", leave = False, colour = "green"):
                
                
                profil_conc_maill =  []
                profil_rho_r_maill = []

                # On regarde ce qui tombe au sol et on l'enregistre
                precip_maill = self.Eq_config.calcul_precip(self.hauteur_interf[maille_dep], self.hauteur_interf[maille_dep+1], profil_concentration_0[maille_dep],self.delta_t*(t_time+1), Lambda[maille_dep])
                precip_dt+= np.nan_to_num(precip_maill)

                # On fait ensuite sédimenter sur les mailles inférieures
                for maille_arr in range(maille_dep+1):
                    conc_sed_i_to_j = self.Eq_config.calcul_maille_arrivee(self.hauteur_interf[maille_arr], self.hauteur_interf[maille_arr+1], self.hauteur_interf[maille_dep], self.hauteur_interf[maille_dep+1],profil_concentration_0[maille_dep], "concentration", self.delta_t*(t_time+1), Lambda[maille_dep])
                    rho_r_sed = self.Eq_config.calcul_maille_arrivee(self.hauteur_interf[maille_arr], self.hauteur_interf[maille_arr+1], self.hauteur_interf[maille_dep], self.hauteur_interf[maille_dep+1], profil_concentration_0[maille_dep], "masse", self.delta_t*(t_time+1), Lambda[maille_dep])
                    
                    # Tu créer le profil sédimenté  en faisant apparaître 
                    profil_conc_maill.append(conc_sed_i_to_j)
                    profil_rho_r_maill.append(rho_r_sed)
                    
                    #print("test, conc_sed_i_to_j : ", conc_sed_i_to_j)
                

                profil_rho_r_maill = np.pad(profil_rho_r_maill,(0,len(self.epaiss_maille)-maille_dep-1))
                profil_conc_maill = np.pad(profil_conc_maill, (0, len(self.epaiss_maille)-maille_dep-1))

                profil_concentration_dt += np.nan_to_num(profil_conc_maill)
                profil_rho_r_dt += np.nan_to_num(profil_rho_r_maill)

                #print("test, conc_profil_int : ", conc_profil_int, len(conc_profil_int))


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

            profil_concentration.append(profil_concentration_dt)
            profil_rho_r.append(profil_rho_r_dt)
            liste_precip.append(precip_dt-sum(liste_precip))


        return profil_concentration, liste_precip, profil_rho_r

        
        
