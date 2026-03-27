
from box_lagrangien_vectorised import Model_bl
from box_lagrangien_def import model_bl_def
from box_lagrangien_sf_vectorised import Model_bl_sf
from box_lagrangien_def_sf import model_bl_def_sf


from affichage import Affichage
from equations import Eq


"""
from phyex import Eule, Eule2, Stat
"""

import numpy as np
import time
from pathlib import Path
import sys
import matplotlib.pyplot as plt


class distribution:
    def __init__(self,model,type_advance,number_stitches,deformable,number_bin,mixing_ratio,time_step,speed_max,esp,CFL, type_init, path_phyex, path_fig, diag):
        h_tot=12000
        duree_sim = 2000
        chemin= os.path.join(f"./fig",f"{model}", f"{type_advance}", f"{deformable}" )
        os.makedirs(chemin, exist_ok=True)
        
        path_to_phyex = Path(path_phyex) / "PHYEX"
        sys.path.append(str(path_to_phyex))
        #from phyex import Eule, Eule2, Stat
        if model == "Box_Lagrangien":
            if type_advance == "Step_By_Step":                    
                if deformable  == "No":     #  Par défaut on arrive ici.
                    a = time.time()
                    model_config = Model_bl(number_stitches,number_bin,mixing_ratio,time_step,speed_max,esp,CFL,type_init,duree_sim)
                    
                    results = model_config.run()

                    lam=model_config.lam_init
                    N=model_config.conc_tot_init
                    Quantiles=Eq(esp).sedimentation_times(N, lam, h_tot,number_stitches)

                    concentration_formate = np.array(results[0]).sum(axis=1)
                    # print(results[0])
                    b = time.time()
                    mass_form = np.array(results[2]).sum(axis=1)
                    Affichage.Affichage_Concentration(concentration_formate, "Concentration", model, path_fig,type_advance,deformable)
                    Affichage.Affichage_Concentration(mass_form, "Masse", model, path_fig,type_advance,deformable)
                    Affichage.Affichage_Precipitation(results[1], model,path_fig,type_advance,deformable, Quantiles, duree_sim)
                    #Affichage.Afficher()
                else:
                    a = time.time()
                    model_config = model_bl_def(number_stitches,time_step,esp,mixing_ratio,type_init,duree_sim)

                    results = model_config.run()

                    lam=model_config.lam_init
                    N=model_config.conc_tot_init
                    Quantiles=Eq(esp).sedimentation_times(N, lam, h_tot,number_stitches)

                    concentration_formate = np.array(results[0])
                    b = time.time()
                    mass_form = np.array(results[2])
                    Affichage.Affichage_Concentration(concentration_formate, "Concentration", model, path_fig,type_advance,deformable)
                    Affichage.Affichage_Concentration(mass_form, "Masse", model, path_fig,type_advance,deformable)
                    Affichage.Affichage_Precipitation(results[1], model,path_fig,type_advance,deformable, Quantiles, duree_sim)
                    #Affichage.Afficher()
                        results = model_config.run()
                        concentration_formate = np.array(results[0])
                        mass_form = np.array(results[2])

                        fig_config = Affichage(model, path_fig,type_advance,deformable)
                        fig_config.afficher(concentration_formate,mass_form,results[1])

                else:
                    if deformable == "No":
                        a = time.time()
                        model_config = Model_bl_sf(number_stitches,number_bin,mixing_ratio,time_step,speed_max,esp,CFL,type_init)
                        b = time.time()

                        results = model_config.run()
                        concentration_formate = np.array(results[0]).sum(axis=1)
                        mass_form = np.array(results[2]).sum(axis=1)

                        fig_config = Affichage(model, path_fig,type_advance,deformable)
                        fig_config.afficher(concentration_formate,mass_form,results[1])

                    else:
                        a = time.time()
                        model_config = model_bl_def_sf(number_stitches,time_step,esp,mixing_ratio,type_init)
                        b = time.time()

                        results = model_config.run()
                        concentration_formate = np.array(results[0])
                        mass_form = np.array(results[2])
                        
                        fig_config = Affichage(model, path_fig,type_advance,deformable)
                        fig_config.afficher(concentration_formate,mass_form,results[1])

            elif model in ('EULE', 'EULE2', 'STAT'):
                a = time.time()
                cls = {'EULE': Eule, 'EULE2': Eule2, 'STAT': Stat}[model]
                model_config = cls(number_stitches,number_bin,time_step,speed_max,esp,CFL,type_init)
            else:
                if deformable == "No":
                    a = time.time()
                    model_config = Model_bl_sf(number_stitches,number_bin,mixing_ratio,time_step,speed_max,esp,CFL,type_init,duree_sim)
                    
                    results = model_config.run()

                    lam=model_config.lam_init
                    N=model_config.conc_tot_init
                    Quantiles=Eq(esp).sedimentation_times(N, lam, h_tot,number_stitches)

                    concentration_formate = np.array(results[0]).sum(axis=1)
                    b = time.time()
                    mass_form = np.array(results[2]).sum(axis=1)
                    Affichage.Affichage_Concentration(concentration_formate, "Concentration", model, path_fig,type_advance,deformable)
                    Affichage.Affichage_Concentration(mass_form, "Masse", model, path_fig,type_advance,deformable)
                    Affichage.Affichage_Precipitation(results[1], model,path_fig,type_advance,deformable, Quantiles, duree_sim)
                    #Affichage.Afficher()
                else:
                    a = time.time()
                    model_config = model_bl_def_sf(number_stitches,time_step,esp,mixing_ratio,type_init,duree_sim)
                    results = model_config.run()

                    lam=model_config.lam_init
                    N=model_config.conc_tot_init
                    Quantiles=Eq(esp).sedimentation_times(N, lam, h_tot,number_stitches)

                    concentration_formate = np.array(results[0])
                    b = time.time()
                    mass_form = np.array(results[2])
                    Affichage.Affichage_Concentration(concentration_formate, "Concentration", model, path_fig,type_advance,deformable)
                    Affichage.Affichage_Concentration(mass_form, "Masse", model, path_fig,type_advance,deformable)
                    Affichage.Affichage_Precipitation(results[1], model,path_fig,type_advance,deformable, Quantiles, duree_sim)
                    #Affichage.Afficher()
        elif model in ('EULE', 'EULE2', 'STAT'):
            a = time.time()
            cls = {'EULE': Eule, 'EULE2': Eule2, 'STAT': Stat}[model]
            #print (cls)
            model_config = cls(number_stitches,number_bin,time_step,speed_max,esp,CFL,type_init)
            #model_config_bl = Model_bl_sf(number_stitches,number_bin,time_step,speed_max,esp,CFL,type_init)

            results = model_config.run()

            lam=model_config.lam_init
            N=model_config.conc_tot_init
            Quantiles=Eq(esp).sedimentation_times(N, lam, h_tot,number_stitches)

            concentration_formate = np.array(results[2])
            Affichage.Affichage_Concentration(concentration_formate, "Concentration", model, path_fig, type_advance, deformable)
            mass_form = np.array(results[1])
            Affichage.Affichage_Concentration(mass_form, "Masse", model, path_fig, type_advance, deformable)
            Affichage.Affichage_Precipitation(results[0], model, path_fig, type_advance, deformable, Quantiles, model_config.duree_sim)
            Affichage.Afficher()
                results = model_config.run()
                concentration_formate = np.array(results[2])
                mass_form = np.array(results[1])

                fig_config = Affichage(model, path_fig,type_advance,deformable)
                fig_config.afficher(concentration_formate,mass_form,results[0])

            #print (f"{model_config.__dict__}")
            #print (f"{model_config_bl.__dict__}")

            #plt.plot(model_config.levels, model_config.rho_r_profile)
            #plt.show()
