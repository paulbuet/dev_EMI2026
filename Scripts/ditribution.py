import numpy as np
from box_lagrangien_vectorised import Model_bl
from box_lagrangien_sf_vectorised import Model_bl_sf
from box_lagrangien import Model_bl as Model_bl_old
from box_lagrangien_sf import Model_bl_sf as Modelbl_sf_old
from box_lagrangien_def import model_bl_def_3 
"""
from phyex import Eule, Eule2, Stat
"""
import time
from pathlib import Path
import os, sys
import matplotlib.pyplot as plt

from affichage import Affichage

class distribution:
    def __init__(self,model,type_advance,number_stitches,deformable,number_bin,mixing_ratio,time_step,speed_max,esp,CFL, efficiency_test, type_init, path_phyex, path_fig):
        
        chemin= os.path.join("./results", f"Model_{model}", f"Type_{type_advance}", f"Stitches_number_{number_stitches}", f"Time_step_{time_step}")
        os.makedirs(chemin, exist_ok=True)

        path_to_phyex = Path(path_phyex) / "PHYEX"
        sys.path.append(str(path_to_phyex))
        from phyex import Eule, Eule2, Stat

        if efficiency_test == "Yes" : 
            
            if model == "Box_Lagrangien":
                if type_advance == "Step_By_Step":

                    if deformable == "No":
                        a = time.time()
                        model_config = Model_bl_old(number_stitches,number_bin,mixing_ratio,time_step,speed_max,esp,CFL,type_init)
                        
                        results = model_config.run()
                        concentration_formate = np.array(results[0]).sum(axis=1)
                        b = time.time()
                        print (f"temps : {b-a} s")

                        a = time.time()
                        model_config = Model_bl(number_stitches,number_bin,mixing_ratio,time_step,speed_max,esp,CFL,type_init)

                        results = model_config.run()
                        concentration_formate = np.array(results[0]).sum(axis=1)
                        
                        b = time.time()
                        print (f"temps v : {b-a} s")

                        mass_form = np.array(results[2]).sum(axis=1)
                        Affichage.Affichage_Concentration(concentration_formate, "Concentration", chemin, b-a)
                        Affichage.Affichage_Concentration(mass_form, "Masse", chemin, b-a)
                        Affichage.Affichage_Precipitation(results[1], chemin, b-a)
                        Affichage.Afficher()

                else:
                    if deformable == "No":
                        a = time.time()
                        model_config = Model_bl_sf(number_stitches,number_bin,mixing_ratio,time_step,speed_max,esp,CFL,type_init)

                        results = model_config.run()
                        concentration_formate = np.array(results[0]).sum(axis=1)
                        b = time.time()
                        mass_form = np.array(results[2]).sum(axis=1)
                        Affichage.Affichage_Concentration(concentration_formate, "Concentration", chemin, b-a)
                        Affichage.Affichage_Concentration(mass_form, "Masse", chemin, b-a)
                        Affichage.Affichage_Precipitation(results[1], chemin, b-a)
                        Affichage.Afficher()
        else : 
            
            if model == "Box_Lagrangien":
                if type_advance == "Step_By_Step":                    
                    if deformable == "No":     #  Par défaut on arrive ici.
                        a = time.time()
                        model_config = Model_bl(number_stitches,number_bin,mixing_ratio,time_step,speed_max,esp,CFL,type_init)

                        results = model_config.run()
                        concentration_formate = np.array(results[0]).sum(axis=1)
                        # print(results[0])
                        b = time.time()
                        mass_form = np.array(results[2]).sum(axis=1)
                        Affichage.Affichage_Concentration(concentration_formate, "Concentration", chemin, b-a)
                        Affichage.Affichage_Concentration(mass_form, "Masse", chemin,b-a)
                        Affichage.Affichage_Precipitation(results[1], chemin, b-a)
                        Affichage.Afficher()
                    else:
                        a = time.time()
                        model_config = model_bl_def_3(number_stitches,time_step,esp,mixing_ratio,type_init)

                        profil = model_config.run()
                        concentration_formate = np.array(profil[0])
                        b = time.time()
                        mass_form = np.array(profil[2])
                        Affichage.Affichage_Concentration(concentration_formate, "Concentration", chemin, b-a)
                        Affichage.Affichage_Concentration(mass_form, "Masse", chemin, b-a)
                        Affichage.Affichage_Precipitation(profil[1], model, chemin, b-a)
                        Affichage.Afficher()

                else:
                    if deformable == "No":
                        a = time.time()
                        model_config = Model_bl_sf(number_stitches,number_bin,mixing_ratio,time_step,speed_max,esp,CFL,type_init)

                        results = model_config.run()
                        concentration_formate = np.array(results[0]).sum(axis=1)
                        b = time.time()
                        mass_form = np.array(results[2]).sum(axis=1)
                        Affichage.Affichage_Concentration(concentration_formate, "concentration", chemin, b-a)
                        Affichage.Affichage_Concentration(mass_form, "masse", chemin, b-a)
                        Affichage.Affichage_Precipitation(results[1], chemin, b-a) 
                        Affichage.Afficher()

            elif model in ('EULE', 'EULE2', 'STAT'):
                a = time.time()
                cls = {'EULE': Eule, 'EULE2': Eule2, 'STAT': Stat}[model]
                #print (cls)
                model_config = cls(number_stitches,number_bin,time_step,speed_max,esp,CFL,type_init)
                #model_config_bl = Model_bl_sf(number_stitches,number_bin,time_step,speed_max,esp,CFL,type_init)

                results = model_config.run()
                b = time.time()
                concentration_formate = np.array(results[1])
                Affichage.Affichage_Concentration(concentration_formate, "Concentration", chemin, b-a)
                concentration_formate = np.array(results[2])
                Affichage.Affichage_Concentration(concentration_formate, "Masse", chemin, b-a)
                #Affichage.Affichage_Precipitation(results[0], model, path_fig)
                Affichage.Afficher()

                print (f"{model_config.__dict__}")
                #print (f"{model_config_bl.__dict__}")

                #plt.plot(model_config.levels, model_config.rho_r_profile)
                #plt.show()
