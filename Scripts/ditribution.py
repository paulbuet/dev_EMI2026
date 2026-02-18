import numpy as np
from box_lagrangien_vectorised import Model_bl
from box_lagrangien_sf_vectorised import Model_bl_sf

from fonctions import Affichage

class distribution:
    def __init__(self,model,type_advance,number_stitches,deformable,number_bin,number_particules,time_step,speed_max,esp,CFL):
        
        if model == "Box_Lagrangien":
            if type_advance == "Step_By_Step":
                
                if deformable == None:
                    model_config = Model_bl(number_stitches,number_bin,number_particules,time_step,speed_max,esp,CFL)

                    profil = model_config.run()
                    concentration_formate = np.array(profil[0]).sum(axis=1)

                    mass_form = np.array(profil[2]).sum(axis=1)
                    Affichage.Affichage_Concentration(concentration_formate, "concentration")
                    Affichage.Affichage_Concentration(mass_form, "masse")
                    Affichage.Affichage_Precipitation(profil[1])
            else:
                if deformable == None:
                    model_config = Model_bl_sf(number_stitches,number_bin,number_particules,time_step,speed_max,esp,CFL)

                    profil = model_config.run()
                    concentration_formate = np.array(profil[0]).sum(axis=1)

                    mass_form = np.array(profil[2]).sum(axis=1)
                    Affichage.Affichage_Concentration(concentration_formate, "concentration")
                    Affichage.Affichage_Concentration(mass_form, "masse")
                    Affichage.Affichage_Precipitation(profil[1])
