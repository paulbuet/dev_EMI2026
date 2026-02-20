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

class Model_bl():
   
   def __init__(self, number_stitches,time_step,speed_max,esp,CFL):
        """
        Here we initialise the non-spatial fixed parameters and allow important variables 
        to travel between functions. We also call the initialisation.
        """

        self.number_stitches = number_stitches

        self.length_sim = 200  # length of simulation in seconds

        self.delta_t = time_step # length of time step in seconds

        self.nb_step = self.length_sim // self.delta_t  # number of time step

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
   
        condi_init = InitialCond(self.number_stitches,self.esp,nb_classes = 1,N=N)
   
        self.grid0 = condi_init.data

        self.vertical_boundaries = condi_init.levels_boundaries

        self.vec_bound = sorted(np.concatenate((self.vertical_boundaries,self.vertical_boundaries)))

        # Initialisation of variables

        self.water_on_floor=0   # Mass of water wich has touched the ground

        self.wat_flo_on_time=[0]  # List of this mass in time, here, time = 0 soit no precipitation

        self.dz = self.grid0["level"].values[2]-self.grid0["level"].values[1]   # length of a stitch

        if CFL == "Yes":
            self.speed_max = self.dz / self.delta_t
        else:
            self.speed_max = speed_max

        self.z_top_ref = self.grid0["level"].values[-1] + self.dz/2

        self.list_mass = [[self.grid0[f"concentration_bin_{diam}"].values*Eq(self.esp).Masse(self.size_diam[diam-1]) for diam in range(1,self.nb_diam+1)]]

        # speed is calculated

        self.speeds = [Eq(self.esp).Vitesse(self.size_diam[n_diam]) * int(self.size_diam[n_diam] <= self.speed_max) + self.speed_max * int(self.size_diam[n_diam] > self.speed_max) for n_diam in range(self.nb_diam)]
        

   def mass(self,grid,var, diam):
       return sum(grid[var].values)*Eq(self.esp).Masse(self.size_diam[diam-1])

  
   
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

        for t in range(self.nb_step):

            grid_dt = xr.Dataset(data_vars={}, coords = {"level" : grid_t["level"]})

            grid_dt = grid_dt.assign(**{"concentration":(("level",),np.zeros())})

            self.new_vec_bound = self.vec_bound - self.speed * self.delta_t

            self.speed_max = np.max(self.speed)

            height_bot = -self.speed_max * self.delta_t

            for i in range(self.number_stitches):
                dz = (self.new_vec_bound[i*2+1]  - self.new_vec_bound[i*2]) /2
                nb_stit_inf = (abs(height_bot)+self.new_vec_bound[i*2]) // dz + 1
                nb_stit_sup = (self.z_top_ref - self.new_vec_bound[i*2+1]) // dz + 1

                mesh_inf = UniformGrid1D(dz = dz,nx = nb_stit_inf,  origin = (self.new_vec_bound[i*2+1],))
                mesh_sup = UniformGrid1D(dz = dz,nx = nb_stit_sup,  origin = (self.new_vec_bound[i*2+1],))

                mesh_tot = mesh_inf + mesh_sup

                data = np.zeros(nb_stit_inf+nb_stit_sup)

                data = np.insert(data,concentration,np.where(mesh_tot.cellCenters==self.new_vec_bound[i*2+1]-dz))

                self.data_i = xr.Dataset(data_vars= data, coords = {"level" : mesh_tot.cellCenters})

                grid_i = self.data_i.regrid.conservative(grid0,time_dim=None)

                grid_i["concentration"] = grid_i["concentration"].copy(data=grid_i["concentration"].data.todense())


                grid_dt = grid_dt + grid_i


            M0 = round(self.mass(grid_dt),3)
            M1 = round(self.mass(grid_t),3)
            
            self.water_on_floor = M0-M1

            list_data.append(grid_dt[f"concentration"].values)
            wat_flo_tot.append(self.water_on_floor)
            self.list_mass.append(list_mass_bin)
            self.wat_flo_on_time.append(sum(wat_flo_tot))
            self.list_data.append(list_data_bin)
            grid_t = grid_dt.copy()








        