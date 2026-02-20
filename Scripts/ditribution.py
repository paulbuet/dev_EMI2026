import numpy as np
from box_lagrangien_vectorised import Model_bl
from box_lagrangien_sf_vectorised import Model_bl_sf
from box_lagrangien import Model_bl as Model_bl_old
from box_lagrangien_sf import Model_bl_sf as Modelbl_sf_old
import time

from fonctions import Affichage

class distribution:
    def __init__(self,model,type_advance,number_stitches,deformable,number_bin,number_particules,time_step,speed_max,esp,CFL, efficiency_test, type_init):
        if efficiency_test == "Yes" : 
            
            if model == "Box_Lagrangien":
                if type_advance == "Step_By_Step":

                    if deformable == "No":
                        a = time.time()
                        model_config = Model_bl_old(number_stitches,number_bin,number_particules,time_step,speed_max,esp,CFL,type_init)

                        profil = model_config.run()
                        concentration_formate = np.array(profil[0]).sum(axis=1)
                        b = time.time()
                        print (f"temps : {b-a} s")

                        a = time.time()
                        model_config = Model_bl(number_stitches,number_bin,number_particules,time_step,speed_max,esp,CFL,type_init)

                        profil = model_config.run()
                        concentration_formate = np.array(profil[0]).sum(axis=1)
                        
                        b = time.time()
                        b = time.time()
                        print (f"temps v : {b-a} s")

                        mass_form = np.array(profil[2]).sum(axis=1)
                        Affichage.Affichage_Concentration(concentration_formate, "concentration", model)
                        Affichage.Affichage_Concentration(mass_form, "masse", model)
                        Affichage.Affichage_Precipitation(profil[1], model)
                else:
                    if deformable == "No":
                        model_config = Model_bl_sf(number_stitches,number_bin,number_particules,time_step,speed_max,esp,CFL,type_init)

                        profil = model_config.run()
                        concentration_formate = np.array(profil[0]).sum(axis=1)

                        mass_form = np.array(profil[2]).sum(axis=1)
                        Affichage.Affichage_Concentration(concentration_formate, "concentration", model)
                        Affichage.Affichage_Concentration(mass_form, "masse", model)
                        Affichage.Affichage_Precipitation(profil[1], model)
        else : 
            
            if model == "Box_Lagrangien":
                if type_advance == "Step_By_Step":
                    
                    if deformable == "No":

                        model_config = Model_bl(number_stitches,number_bin,number_particules,time_step,speed_max,esp,CFL,type_init)

                        profil = model_config.run()
                        concentration_formate = np.array(profil[0]).sum(axis=1)

                        mass_form = np.array(profil[2]).sum(axis=1)
                        Affichage.Affichage_Concentration(concentration_formate, "concentration", model)
                        Affichage.Affichage_Concentration(mass_form, "masse", model)
                        Affichage.Affichage_Precipitation(profil[1], model)
                else:
                    if deformable == "No":
                        model_config = Model_bl_sf(number_stitches,number_bin,number_particules,time_step,speed_max,esp,CFL,type_init)

                        profil = model_config.run()
                        print (profil)
                        concentration_formate = np.array(profil[0]).sum(axis=1)

                        mass_form = np.array(profil[2]).sum(axis=1)
                        Affichage.Affichage_Concentration(concentration_formate, "concentration", model)
                        Affichage.Affichage_Concentration(mass_form, "masse", model)
                        Affichage.Affichage_Precipitation(profil[1], model) 