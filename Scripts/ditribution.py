import numpy as np
from box_lagrangien_vectorised import Model_bl
from box_lagrangien_sf_vectorised import Model_bl_sf
from box_lagrangien import Model_bl as Model_bl_old
from box_lagrangien_sf import Model_bl_sf as Modelbl_sf_old
from phyex import Eule, Eule2, Stat
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
                        concentration_formate = np.array(profil[0]).sum(axis=1)

                        mass_form = np.array(profil[2]).sum(axis=1)
                        Affichage.Affichage_Concentration(concentration_formate, "concentration", model)
                        Affichage.Affichage_Concentration(mass_form, "masse", model)
                        Affichage.Affichage_Precipitation(profil[1], model)

            elif model in ('EULE', 'EULE2', 'STAT'):
                cls = {'EULE': Eule, 'EULE2': Eule2, 'STAT': Stat}[model]
                print (cls)
                model_config = cls(number_stitches,number_bin,number_particules,time_step,speed_max,esp,CFL)
                model_config_bl = Model_bl_sf(number_stitches,number_bin,number_particules,time_step,speed_max,esp,CFL,type_init)

                model_config.run()

                print (f"{model_config.__dict__}")
                print (f"{model_config_bl.__dict__}")
                print(model_config.rho_r_profile)

