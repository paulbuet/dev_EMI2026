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

# On importe ici les classes extèrieures
from fonctions import InitialCond
from fonctions import Eq

class Model_bl_def():
   
   def __init__(self, number_stitches,time_step,speed_max,esp,CFL):
        """
        Here we initialise the non-spatial fixed parameters and allow important variables 
        to travel between functions. We also call the initialisation.
        """

        self.number_stitches = number_stitches

        self.length_sim = 200  # length of simulation in seconds

        self.delta_t = time_step # length of time step in seconds

        self.nb_step = self.length_sim // self.delta_t  # number of time step
   
        # Ajouter le calcul des différents diamètres dans une liste self.diameter

        """
        grid is in the form of a numpy array of size nb_mesh*(variables*2):

        concentration_bin_1 : number of particles with diameter 1 in each level: N0(1)  N1(1)  N2 (1)  N3(1) ...
        ...
        concentration_bin_i : number of particles with diameter i in each level: N0(i)  N1(i)  N2 (i)  N3(i) ...
        ...
        concentration_bin_n : number of particles with diameter n in each level: N0(n)  N1(n)  N2 (n)  N3(n) ...

        level : height of upper mesh interfaces: Z0  Z1  Z2  Z3  ...
        
        """
   
        condi_init = InitialCond(self.number_stitches,esp,nb_classes = 1)
   
        self.grid0 = condi_init.data

        self.vertical_boundaries = condi_init.levels_boundaries

        self.rho_r = condi_init.rho_r

        self.vec_bound = sorted(np.concatenate((self.vertical_boundaries,self.vertical_boundaries)))

        # Initialisation of variables

        dz = self.grid0["level"].values[2]-self.grid0["level"].values[1]   # length of a stitch

        if CFL == "Yes":
            self.speed_max = dz / self.delta_t
        else:
            self.speed_max = speed_max

        # La classe Eq étant appelé plusieurs fois pour la même espèce, on la configure ici
        self.eq = Eq(esp)
        
        # On cazlcule la hauteur de l'interface la plus haute
        self.z_top_ref = self.grid0["level"].values[-1] + self.dz/2


   def regridage(self,grid,variable,vec_bound,h_bot,stitch):
        """
        Cette fonction prends en entrée des variables qui lui permettent de regridder une maille déformée sur
        le maillage régulier de départ de manière conservatif.

        grid (type = dataset): dataset contenant les level en coordonées et la variable (masse/concentration) en data variables
        variable (type = str): chaine de caractère pouvant être masse ou concentration suivant le critère de sédimentation
        vec_bound (type  = list): liste d'interface supérieurs des mailles déformée 
        h_bot (type = float): hauteur minimale atteinte une particule au sol avec la vitesse de chute maximale
        stitch (type = int): numéro de la maille déformée observée

        retourne la grille régulière avec des 0 partout sauf sur les mailles affectées par la maille déformée:

        grid_i (type = dataset)
        """
        # On cherche l'épaisseur de la maille déformée
        dz = (vec_bound[stitch*2+1]  - vec_bound[stitch*2]) /2

        # Pour éviter les erreurs de regridage, on crée des mailles de même épaisseurs en dessous et au dessus
        # On calcule ici leurs nombres
        nb_stit_inf = (abs(h_bot)+ vec_bound[stitch*2]) // dz + 1
        nb_stit_sup = (self.z_top_ref - vec_bound[stitch*2+1]) // dz + 1

        # On crée ensuite la maille déformée totale
        mesh_inf = UniformGrid1D(dz = dz,nx = nb_stit_inf,  origin = (vec_bound[i*2+1],))
        mesh_sup = UniformGrid1D(dz = dz,nx = nb_stit_sup,  origin = (vec_bound[i*2+1],))

        mesh_tot = mesh_inf + mesh_sup

        # On insère la valeur portée par la maille déformée dans la grille déformée
        data = np.zeros(nb_stit_inf+nb_stit_sup)
        data = np.insert(data,grid[variable].values[stitch],np.where(mesh_tot.cellCenters == vec_bound[stitch*2+1]-dz/2))

        # On crée le dataset avec les mailles déformée et leur valeur de variables associées
        data_i = xr.Dataset(data_vars= dict(variable = ((level),data)), coords = {"level" : mesh_tot.cellCenters})

        # On regrid sur la maille non déformée
        grid_i = data_i.regrid.conservative(self.grid0,time_dim=None)

        # On densifie le vecteur 
        grid_i[variable] = grid_i[variable].copy(data=grid_i[variable].data.todense())

        return grid_i

   
   def run(self):
        
        """
        On fait tourner le modèle pour chaque pas de temps, pour chaque mailles, pour chaque diamètres.
        On enregistre le profils des concentration en fonction de la hauteur à chaque pas de temps
        """ 
        
        # On initialise la concentration, la masse ainsi que les variables de stockage

        grid_t_conc = self.grid0

        grid_t_mass = self.grid0

        rho_r_conc = self.rho_r
        rho_r_mass = self.rho_r

        # Les masses
        Masse= round(self.eq.Calcul_Masse_Tot(rho_r_mass,self.z_top_ref),4)  
        list_mass_tot = [Masse]
        
        #Les précip
        water_on_floor = 0               
        wat_flo_tot = [water_on_floor]
        wat_flo_on_time = [sum(wat_flo_tot)]

        #Les data
        list_data = [grid_t_conc["concentration"].values]
        list_mass= [grid_t_mass["masse"].values]
        
        

        for t in range(self.nb_step):
            # On initialise aussi concentration et masse au pas de temps suivant

            grid_dt_conc = xr.Dataset(data_vars={}, coords = {"level" : grid_t["level"]})
            grid_dt_mass = xr.Dataset(data_vars={}, coords = {"level" : grid_t["level"]})

            # On calcule les lambda de chaque mailles 

            list_lamb_conc = self.eq.Liste_Lanbda(rho_r_conc, grid_t_conc["concentration"].values)
            list_lamb_mass = self.eq.Liste_Lanbda(rho_r_mass, grid_t_conc["masse"].values)

            # On initialise les lignes de valeurs à 0

            grid_dt_conc = grid_dt_conc.assign(**{"concentration":(("level",),np.zeros(self.number_stitches))})
            grid_dt_mass = grid_dt_mass.assign(**{"masse":(("level",),np.zeros(self.number_stitches))})

            # On calcule les vitesses minimales et maximales de chaques mailles
            self.speed_conc = self.eq.Liste_Vitesse_Concentration(list_lamb_conc)
            self.speed_mass = self.eq.Liste_Vitesse_Masse(list_lamb_mass)

            # On applique la vitesse respective des particules en fonction du critère pris (concentration/ masse)
            new_vec_bound_conc = new_vec_bound_conc - self.speed_conc * self.delta_t
            new_vec_bound_mass = new_vec_bound_mass - self.speed_mass * self.delta_t

            # On enregistre la hauteur minimale que peux atteindre un particule en concentration et en masse
            self.speed_max_conc = np.max(self.speed_conc)
            height_bot_conc = -self.speed_max_conc * self.delta_t
            self.speed_max_mass = np.max(self.speed_mass)
            height_bot_mass = -self.speed_max_mass * self.delta_t


            # On itère sur le nombre de mailles
            for stitch in range(self.number_stitches):
                
                # On regride la maille déformée sur la maille de départ puis on ajoute ces valeurs à celles déjà calculées
                grid_stit_conc = self.regridage(grid_t_conc,"concentration",new_vec_bound_conc,height_bot_conc,stitch)
                grid_stit_mass = self.regridage(grid_t_mass,"masse",new_vec_bound_mass,height_bot_mass,stitch)


                grid_dt_conc = grid_dt_conc + grid_stit_conc
                grid_dt_mass = grid_dt_mass + grid_stit_mass

            # On recalcule notre profil de rho_r sur nos deux sédimentations
            rho_r_conc = self.eq.Liste_rho_r(grid_dt_conc["concentration"].values,list_lamb_conc)
            rho_r_mass = self.eq.Liste_rho_r(grid_dt_conc["masse"].values,list_lamb_mass)

            # On calcule la masse de nos nouvelles grilles et on l'enregistre
            Masse_dt= round(self.eq.Calcul_Masse_Tot(rho_r_mass,self.z_top_ref),4)
            list_mass_tot.append(Masse_dt)

            #On s'occupe des précip
            water_on_floor = list_mass_tot[-2] - Masse_dt
            wat_flo_tot.append(water_on_floor)
            wat_flo_on_time.append(sum(wat_flo_tot))

            #On enregistre les dataset
            list_data.append(grid_dt_conc["concentration"].values)
            list_mass.append(grid_dt_mass["masse"].values)
            
            # On prends le nouveau dataset comme l'ancien
            grid_t_conc = grid_dt_conc.copy()
            grid_t_mass = grid_dt_mass.copy()
        
        # On return les profils stockés pour l'affichage

        return list_data,wat_flo_on_time,list_mass 








        