#!/usr/bin/env python3

# Copyright (c) Météo France (2025-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

# On utilisera dans la description du Box-Lagrangien les notations de Kato(1995)

# On importe les librairies nécessaires
import numpy as np
import xarray as xr
import xarray_regrid

# On importe ici les classes extèrieures
from fonctions import InitialCond
from fonctions import Eq

class Model_bl():
   
   def __init__(self, number_stitches,number_bin,number_particules,delta_t,speed_max,esp,CFL,type_init):
        """
        Here we initialise the non-spatial fixed parameters and allow important variables 
        to travel between functions. We also call the initialisation.
        """

        self.number_stitches = number_stitches

        N = number_particules

        self.length_sim = 200  # length of simulation in seconds

        self.delta_t = delta_t # length of time step in seconds

        self.nb_step = self.length_sim // self.delta_t  # number of time step

        self.nb_diam = number_bin # number of type of diameter
        self.esp = esp
   
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
   
        condi_init = InitialCond(self.number_stitches,'i',nb_classes = self.nb_diam,N=N,mode=type_init)
   
        self.grid0 = condi_init.data

        self.vertical_boundaries = condi_init.levels_boundaries

        self.size_diam = np.array(condi_init.bin_concentration)[:,0]

        # Initialisation of variables

        self.water_on_floor=0   # Mass of water wich has touched the ground

        self.wat_flo_on_time=[0]  # List of this mass in time, here, time = 0 soit no precipitation

        self.dz = self.grid0["level"].values[2]-self.grid0["level"].values[1]   # length of a stitch

        if CFL == "Yes":
            self.speed_max = self.dz / self.delta_t
        else:
            self.speed_max = speed_max

        self.z_top_ref = self.grid0["level"].values[-1] + self.dz/2

        self.list_data = [[self.grid0[f"concentration_bin_{diam}"].values for diam in range(1,self.nb_diam+1)]] #liste des valeurs par bin et par pas de temps
        self.list_mass = [[self.grid0[f"concentration_bin_{diam}"].values*Eq(self.esp).Masse(self.size_diam[diam-1]) for diam in range(1,self.nb_diam+1)]]
        
   def mass(self,grid,var, diam):
       return sum(grid[var].values)*Eq(self.esp).Masse(self.size_diam[diam-1])

   def advect_down(self,ds, V, dt):
        """
        Cette fonction décale les box au temps t d'une vitesse propre claculée en fonction du diamètre du bin
        """

        shift = -V * dt   # On calcule le futur mouvement verticale

        # On applique au centre des mailles le déplacement
        ds = ds.assign_coords(level=ds["level"] + shift)
        return ds
   
   def add_stitche(self,new_grid,diam):
    
        # We create new stitches above the new_grid to have a initially grid tottaly cover

        nb_stitche_create = 1

        z_mid_last_i = new_grid["level"].values[-1]    # height of the highest center of stitch

        grid_add = np.sort(np.unique(np.append(new_grid["level"].values,z_mid_last_i+self.dz)))

        # We add this new stitch with 0 concentration to the dataset

        new_grid = new_grid.reindex(level = grid_add)

        new_grid[f"concentration_bin_{diam}"] = new_grid[f"concentration_bin_{diam}"].copy(data=new_grid[f"concentration_bin_{diam}"].data).fillna(0)

        nb_stitche_create  += 1   # We add 1 to the number of stitch created to see if we need an other one


        while new_grid["level"].values[-1] <=self.z_top_ref:  # If we are above the maximum level, we stop

            z_mid_last_i = new_grid["level"].values[-1]  # New evaluation of the height of the highest center of stitch

            """
            ###SAME PROCESS###
            """

            grid_add = np.sort(np.unique(np.append(new_grid["level"].values,z_mid_last_i+self.dz)))

            new_grid = new_grid.reindex(level = grid_add)

            new_grid[f"concentration_bin_{diam}"] = new_grid[f"concentration_bin_{diam}"].copy(data=new_grid[f"concentration_bin_{diam}"].data).fillna(0)

            nb_stitche_create  += 1

        return new_grid
   
   def run(self):
        
        """
        On fait tourner le modèle pour chaque pas de temps, pour chaque mailles, pour chaque diamètres.
        On enregistre le profils des concentration en fonction de la hauteur à chaque pas de temps
        """ 
        
        grid_t = self.grid0

        for t_time in range(self.nb_step):
            
            grid_dt = xr.Dataset(data_vars={}, coords = {"level" : grid_t["level"]})

            list_data_bin = []
            list_mass_bin = []
            wat_flo_tot=[]


            for diam in range(1,self.nb_diam+1):  

                # speed is calculated

                speed = Eq(self.esp).Vitesse(self.size_diam[diam-1])



                if speed> self.speed_max:

                    speed = self.speed_max


                # Sedimentation is processed

                new_grid = self.advect_down(grid_t,speed,self.delta_t)

                # we add the sticthes above
                new_grid = self.add_stitche(new_grid, diam)

                # We put the stagerd grid on the good grid
                grid_on_old = new_grid.regrid.conservative(grid_t, time_dim=None)

                
                
                # Here, we transform a sparse vector to a dense one
                grid_on_old[f"concentration_bin_{diam}"] = grid_on_old[f"concentration_bin_{diam}"].copy(data=grid_on_old[f"concentration_bin_{diam}"].data.todense())

               
                # We want to see how much rain touch the ground during the time step
                M0 = round(self.mass(grid_t,f"concentration_bin_{diam}",diam),3)
                M1 = round(self.mass(grid_on_old,f"concentration_bin_{diam}",diam),3)
                
                self.water_on_floor = M0-M1


                list_data_bin.append(grid_on_old[f"concentration_bin_{diam}"].values)
                list_mass_bin.append(grid_on_old[f"concentration_bin_{diam}"].values*Eq(self.esp).Masse(self.size_diam[diam-1]))
                wat_flo_tot.append(self.water_on_floor)
                
                grid_dt = grid_dt.assign(**{f"concentration_bin_{diam}":(("level",),grid_on_old[f"concentration_bin_{diam}"].values)})

            self.list_mass.append(list_mass_bin)
            self.wat_flo_on_time.append(sum(wat_flo_tot))
            self.list_data.append(list_data_bin)
            grid_t = grid_dt.copy()


            
        return self.list_data,self.wat_flo_on_time,self.list_mass 










        