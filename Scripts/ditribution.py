import numpy as np
from box_lagrangien_vectorised import Model_bl
from box_lagrangien_sf_vectorised import Model_bl_sf
from box_lagrangien import Model_bl as Model_bl_old
from box_lagrangien_sf import Model_bl_sf as Modelbl_sf_old
from box_lagrangien_def_1 import Model_bl_def 
from box_lagrangien_def_2 import Model_bl_def_2
"""
from phyex import Eule, Eule2, Stat
"""
import time
from pathlib import Path
import os, sys
import matplotlib.pyplot as plt

from fonctions import Affichage

class distribution:
    def __init__(self,model,type_advance,number_stitches,deformable,number_bin,number_particules,time_step,speed_max,esp,CFL, efficiency_test, type_init, path_phyex, path_fig):

        path_to_phyex = Path(path_phyex) / "PHYEX"
        sys.path.append(str(path_to_phyex))
        from phyex import Eule, Eule2, Stat

        if efficiency_test == "Yes" : 
            
            if model == "Box_Lagrangien":
                if type_advance == "Step_By_Step":

                    if deformable == "No":
                        a = time.time()
                        model_config = Model_bl_old(number_stitches,number_bin,number_particules,time_step,speed_max,esp,CFL,type_init)
                        
                        results = model_config.run()
                        concentration_formate = np.array(results[0]).sum(axis=1)
                        b = time.time()
                        print (f"temps : {b-a} s")

                        a = time.time()
                        model_config = Model_bl(number_stitches,number_bin,number_particules,time_step,speed_max,esp,CFL,type_init)

                        results = model_config.run()
                        concentration_formate = np.array(results[0]).sum(axis=1)
                        
                        b = time.time()
                        b = time.time()
                        print (f"temps v : {b-a} s")

                        mass_form = np.array(results[2]).sum(axis=1)
                        Affichage.Affichage_Concentration(concentration_formate, "Concentration", model, path_fig)
                        Affichage.Affichage_Concentration(mass_form, "Masse", model, path_fig)
                        Affichage.Affichage_Precipitation(results[1], model, path_fig)
                        Affichage.Afficher()
                else:
                    if deformable == "No":
                        model_config = Model_bl_sf(number_stitches,number_bin,number_particules,time_step,speed_max,esp,CFL,type_init)

                        results = model_config.run()
                        concentration_formate = np.array(results[0]).sum(axis=1)

                        mass_form = np.array(results[2]).sum(axis=1)
                        Affichage.Affichage_Concentration(concentration_formate, "Concentration", model, path_fig)
                        Affichage.Affichage_Concentration(mass_form, "Masse", model, path_fig)
                        Affichage.Affichage_Precipitation(results[1], model, path_fig)
                        Affichage.Afficher()
        else : 
            
            if model == "Box_Lagrangien":
                if type_advance == "Step_By_Step":                    
                    if deformable == "No":     #  Par défaut on arrive ici.

                        model_config = Model_bl(number_stitches,number_bin,number_particules,time_step,speed_max,esp,CFL,type_init)

                        results = model_config.run()
                        concentration_formate = np.array(results[0]).sum(axis=1)
                        print(results[0])

                        mass_form = np.array(results[2]).sum(axis=1)
                        #Affichage.Affichage_Concentration(concentration_formate, "Concentration", model, path_fig)
                        #Affichage.Affichage_Concentration(mass_form, "Masse", model, path_fig)
                        #Affichage.Affichage_Precipitation(results[1], model, path_fig)
                        #Affichage.Afficher()
                    else:
                        model_config = Model_bl_def(number_stitches,time_step,number_particules,speed_max,esp,CFL,type_init)

                        profil = model_config.run()
                        concentration_formate = np.array(profil[0])

                        mass_form = np.array(profil[2])
                        Affichage.Affichage_Concentration(concentration_formate, "Concentration", model)
                        Affichage.Affichage_Concentration(mass_form, "Masse", model)
                        Affichage.Affichage_Precipitation(profil[1], model)
                        Affichage.Afficher()

                else:
                    if deformable == "No":
                        model_config = Model_bl_sf(number_stitches,number_bin,number_particules,time_step,speed_max,esp,CFL,type_init)

                        results = model_config.run()
                        concentration_formate = np.array(results[0]).sum(axis=1)

                        mass_form = np.array(results[2]).sum(axis=1)
                        Affichage.Affichage_Concentration(concentration_formate, "concentration", model, path_fig)
                        Affichage.Affichage_Concentration(mass_form, "masse", model, path_fig)
                        Affichage.Affichage_Precipitation(results[1], model, path_fig)

            elif model in ('EULE', 'EULE2', 'STAT'):
                cls = {'EULE': Eule, 'EULE2': Eule2, 'STAT': Stat}[model]
                #print (cls)
                model_config = cls(number_stitches,number_bin,number_particules,time_step,speed_max,esp,CFL,type_init)
                #model_config_bl = Model_bl_sf(number_stitches,number_bin,number_particules,time_step,speed_max,esp,CFL,type_init)

                results = model_config.run()

                concentration_formate = np.array(results[3])
                Affichage.Affichage_Concentration(concentration_formate, "Concentration", model, path_fig)
                concentration_formate = np.array(results[2])
                Affichage.Affichage_Concentration(concentration_formate, "Masse", model, path_fig)
                #Affichage.Affichage_Precipitation(results[0], model, path_fig)
                Affichage.Afficher()

                print (f"{model_config.__dict__}")
                #print (f"{model_config_bl.__dict__}")

                #plt.plot(model_config.levels, model_config.rho_r_profile)
                #plt.show()
