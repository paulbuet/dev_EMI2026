#!/usr/bin/env python3

# Copyright (c) Météo France (2025-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

# On utilisera dans la description du Box-Lagrangien les notations de Kato(1995)

# On importe les librairies nécessaires
import numpy as np
import xarray as xr
import xarray_regrid
from tqdm import tqdm

# On importe ici les classes extèrieures
from condi_init import InitialCond
from equations import Eq

class Model_bl_sf():
   
   def __init__(self,nb_mailles,nb_bins,rapp_mel,delta_t,vit_max,esp,CFL,type_init,duree_sim):
        """
        On initialise le modèle Box-Lagrangien indéformable en prenant en ompte les paramètres d'entrée du modèle
        """

        self.nb_mailles = nb_mailles
        self.delta_t = delta_t
        self.nb_bins = nb_bins
        self.esp = esp
        
        # On calcule le nombre de pas de temps
        self.nb_step = duree_sim // delta_t

        # ---Configuration des classes extérieures---
        self.Eq_config = Eq(self.esp)

        # --- # On configure le modèle avec les données d'entrée ---
   
        condi_init = InitialCond(self.nb_mailles,self.esp,"bin",nb_classes = self.nb_bins,mode=type_init,r=rapp_mel)

        #---Récupération des variables---
   
        self.grid0 = condi_init.data
        self.hauteur_interf = np.concatenate(([0],condi_init.levels_boundaries))              # On ajoute le sol
        self.taille_diam = np.array(condi_init.diameters+[condi_init.diameters[-1]*2])


        # Initialisation des variables

        self.eau_au_sol = [0]   

        self.epaiss_maille = np.array([self.hauteur_interf[stitche+1]-self.hauteur_interf[stitche] for stitche in range(self.nb_mailles)])

        if CFL == "Yes":
            # On prend la condition CFL la plus contraignante soit celle où la vitesse max est la plus faible
            self.vit_max = min(self.epaiss_maille) / self.delta_t

        else:
            self.vit_max = vit_max


        self.h_max = self.hauteur_interf[-1]

        self.profil_conc_bin = [[self.grid0[f"concentration_bin_{diam}"].values for diam in range(1,self.nb_bins+1)]] #profil des valeurs par bin et par pas de temps
        self.profil_cont_bin = [[self.grid0[f"rho_r_bin_{diam}"].values for diam in range(1,self.nb_bins+1)]]

        # speed is calculated

        self.vit_bin = [self.Eq_config.Vitesse(self.taille_diam[n_diam]) * int(self.Eq_config.Vitesse(self.taille_diam[n_diam]) <= self.vit_max) + self.vit_max * int(self.Eq_config.Vitesse(self.taille_diam[n_diam]> self.vit_max)) for n_diam in range(self.nb_bins)]

        # On enregistre ces données pour la courbe de référence
        self.lam_init = self.Eq_config.Liste_Lanbda(sum(self.profil_cont_bin[0]),sum(self.profil_conc_bin[0]))[-1]
        self.conc_tot_init = sum(self.profil_conc_bin[0])[-1]
        

   def conc_to_mass(self,grid,var, diam):
       """
       Prend en entrée les concentrations d'un bin par maille sur la colonne et renvoie sa masse totale 
       à partir du diamètre du bin
       """
       nb_part_col = sum(grid[var].values*self.epaiss_maille) 

       masse_part = self.Eq_config.Masse((self.taille_diam[diam-1]+self.taille_diam[diam])/2)
       return nb_part_col * masse_part
   

   def sedim(self,grid, V, dt):
        """
        Cette fonction décale les box au temps t d'une vitesse propre claculée en fonction du diamètre du bin
        """

        dep = -min(V,self.vit_max) * dt   # On calcule le futur mouvement verticale

        # On applique au centre des mailles le déplacement
        grid_dep = grid.assign_coords(level=grid["level"] + dep)
        return grid_dep
   
   def ajout_maille(self,grid_dep,diam):
    
        """ 
        On veut que la grille déplacée recouvre totalement la grille initiale pour utiliser la fonction regrid
        Ainsi, on ajoute les mailles nécessaires
        """

        h_interf_haute = grid_dep["level"].values[-1]    # hauteur du du centre de la maille supérieure de la colonne déplacée

        # On ajoute une maille à la colonne déplacée et on ajoute cette nouvelle maille aux données tel que la concentration dans le bin soit nulle
        grid_augmente = np.sort(np.unique(np.append(grid_dep["level"].values,h_interf_haute+self.epaiss_maille[-1])))
        grid_dep = grid_dep.reindex(level = grid_augmente)
        grid_dep[f"concentration_bin_{diam}"] = grid_dep[f"concentration_bin_{diam}"].copy(data=grid_dep[f"concentration_bin_{diam}"].data).fillna(0)

        h_interf_haute = grid_dep["level"].values[-1] # On recalcule la hauteur de cette colonne déplacée avec la maille ajoutée

        if h_interf_haute < self.h_max:  # Tant que l'on ne dépasse pas la hauteur de la colonne initiale

            """
            ###SAME PROCESS###
            """

            grid_augmente = np.sort(np.unique(np.append(grid_dep["level"].values,self.h_max)))
            grid_dep = grid_dep.reindex(level = grid_augmente)
            grid_dep[f"concentration_bin_{diam}"] = grid_dep[f"concentration_bin_{diam}"].copy(data=grid_dep[f"concentration_bin_{diam}"].data).fillna(0)
        return grid_dep
   
   def run(self):
        
        """
        On fait tourner le modèle pour chaque pas de temps, pour chaque mailles, pour chaque diamètres.
        On enregistre le profils des concentration en fonction de la hauteur à chaque pas de temps
        """ 

        print (" ")
        print ("---------------------------------------")

        format_affichage_time = "{l_bar}|{bar}| {n_fmt}/{total_fmt} pas de temps | temps écoulé : {elapsed} < temps restant {remaining} |" # definition of the format of the tqdm bar
        format_affichage_diam = "-> {desc} |{bar}| {n_fmt}/{total_fmt} bins  | {elapsed}<{remaining} |                                                                                           "

        for t_time in tqdm(range(self.nb_step), bar_format = format_affichage_time, desc = f"Avancement total Box Lagrangien Step_Forward à {self.nb_bins} bins : ", position = 0, colour = 'blue'):
            
            profil_conc_bin_dt = []
            profil_cont_bin_dt = []
            eau_au_sol_dt=[]


            for diam in tqdm(range(1,self.nb_bins+1), bar_format = format_affichage_diam, desc = f"Calculs t = {t_time * self.delta_t} / {self.nb_step * self.delta_t} s : ", leave = False, colour = "green"):  

                # Sedimentation is processed
                grid_dep= self.sedim(self.grid0,self.vit_bin[diam-1],(t_time+1) * self.delta_t)

                # we add the sticthes above
                grid_dep= self.ajout_maille(grid_dep, diam)

                # We put the stagerd grid on the good grid
                grid_regrid = grid_dep.regrid.conservative(self.grid0, time_dim=None)

                
                
                # Here, we transform a sparse vector to a dense one
                grid_regrid[f"concentration_bin_{diam}"] = grid_regrid[f"concentration_bin_{diam}"].copy(data=grid_regrid[f"concentration_bin_{diam}"].data.todense())

               
                # We want to see how much rain touch the ground during the time step
                M0 = round(self.conc_to_mass(self.grid0,f"concentration_bin_{diam}",diam),3)
                M1 = round(self.conc_to_mass(grid_regrid,f"concentration_bin_{diam}",diam),3)
                
                eau_au_sol_bin = M0-M1


                profil_conc_bin_dt.append(grid_regrid[f"concentration_bin_{diam}"].values)
                profil_cont_bin_dt.append(grid_regrid[f"concentration_bin_{diam}"].values*Eq(self.esp).Masse(self.taille_diam[diam-1]))
                eau_au_sol_dt.append(eau_au_sol_bin) 
                

            self.profil_conc_bin.append(profil_conc_bin_dt)
            self.profil_cont_bin.append(profil_cont_bin_dt)
            self.eau_au_sol.append(sum(eau_au_sol_dt)-sum(self.eau_au_sol))

        print ("---------------------------------------")
        print(" ")
            
        return self.profil_conc_bin,self.eau_au_sol,self.profil_cont_bin 


