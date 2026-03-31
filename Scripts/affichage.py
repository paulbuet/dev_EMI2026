#!/usr/bin/env python3

# Copyright (c) Météo France (2025-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

# Le but de ce fichier est de mettre toutes les fonctions qu'on peut potentiellement appelés pour faire tourner nos futur codes
# Ensuite il suffira de faire un appelle du fichier et de la fonction qu nous interesse
# Exemple avec cette fonction test que j'appelle dans le code Sedimentation.py

### Imports ###

from math import *
import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import MaxNLocator
from pathlib import Path
from matplotlib.colors import Normalize, LogNorm



class Figure :
    def sedimentation_time2(Concentration, typ, chemin, params):

        Concentration = np.array(Concentration)
        Transpose = Concentration.T

        # 🔥 Création figure + grille
        fig = plt.figure(figsize=(10, 6))
        gs = fig.add_gridspec(2, 1, height_ratios=[4, 1])

        # --- Graphe principal ---
        ax = fig.add_subplot(gs[0])

        orig_map = plt.cm.get_cmap('gist_ncar')
        reversed_map = orig_map.reversed()
        norm = Normalize(vmin=0, vmax=max(Concentration[0]))

        im = ax.pcolormesh(Transpose, cmap=reversed_map, norm=norm)

        ax.set_title(f"{typ} : évolution dans le temps", fontsize=22)
        ax.set_xlabel("Temps", fontsize=18)
        ax.set_ylabel("Mailles du modèle", fontsize=18)

        fig.colorbar(im, ax=ax)

        # --- Texte (zone dédiée) ---
        ax_caption = fig.add_subplot(gs[1])
        ax_caption.axis("off")

        if params["model"] in ('EULE', 'EULE2', 'STAT'):
            #Dans ce cas : params = model,path_fig, number_stitches, time_step, esp
            ax_caption.text(0.5, 0, f"Précipitation au sol, model {params["model"]},\n {params["number_stitches"]} mailles, pas de temps de {params["time_step"]} s, durée de la simulation {params["duree_sim"]},\n espece {params["esp"]}, temps de calculs : {params["time_run"]} s ", ha='center', va='bottom', wrap=True, fontsize=10)
        else :
            if params["deformability"] == "Yes" :
                params["deformability"]= "déformable"
                ax_caption.text(0.5, 0, f"Précipitation au sol, model {params["model"]}, {params["type_advance"]}, {params["deformability"]},\n {params["number_stitches"]} mailles, pas de temps de {params["time_step"]} s, durée de la simulation {params["duree_sim"]},\n espece {params["esp"]}, temps de calculs : {params["time_run"]} s ", ha='center', va='bottom', wrap=True, fontsize=10)
            else :
                params["deformability"] = "indéformable"
                ax_caption.text(0.5, 0, f"Sédimentation de la {typ}, model {params["model"]}, {params["type_advance"]}, {params["deformability"]},\n {params["number_stitches"]} mailles, nombre de bin : {params["number_bin"]}, pas de temps de {params["time_step"]} s, durée de la simulation {params["duree_sim"]},\n espece {params["esp"]}, temps de calculs : {params["time_run"]} s ", ha='center', wrap=True, fontsize=12)


        # --- Sauvegarde ---

        file_location = chemin / Path(typ)
        fig.savefig(str(file_location))
        plt.close()




    def precipitation(Precip, chemin, params, Quantiles):


        # --- Données ---
        Precip = np.array(Precip, dtype=float)

        reference = np.array(Quantiles[0], dtype=float) # Courbe du Cumul théorique (courbe verte)
        time3 = np.array(Quantiles[1], dtype=float)

        # Ajout du point initial (0,0) (pour que la courbe de reférence parte de (0,0))
        reference = np.insert(reference, 0, 0)
        time3 = np.insert(time3, 0, 0)

        # --- Cumul correct ---
        Cumul = np.cumsum(Precip)

        # --- Temps ---
        n = len(Precip)
        time = np.linspace(0, params["duree_sim"], n) #params -2 c'est durée simulation
        dt = time[1] - time[0]
        time2 = time - dt / 2 # pour centrer les barres

        # --- Figure ---
        fig, ax1 = plt.subplots(figsize=(10, 6))

        # --- Création d'un espace texte ---

        gs = fig.add_gridspec(2, 1, height_ratios=[4, 1])

        # --- Axe principal (barres) ---
        ax1.bar(time2, Precip, width=dt, color="blue", label="Précip par pas de temps")
        ax1.set_xlabel("Temps", fontsize=14)
        ax1.set_ylabel("Précipitation", fontsize=14)

        # --- Axe secondaire (cumul) ---
        ax2 = ax1.twinx()
        ax2.plot(time, Cumul, '--', color="red", linewidth=2, label="Cumul")
        ax2.plot(time3, reference, '--', color="green", linewidth=2, label="Cumul théorique")
        ax2.set_ylabel("Cumul", fontsize=14)

        # --- Axes propres ---
        ax1.set_xlim(0, max(time3.max(), time.max()) * 1.05)

        ax1.set_ylim(0, Precip.max() * 1.2 if Precip.max() > 0 else 1)
        ax2.set_ylim(0, max(Cumul.max(), reference.max()) * 1.2 if max(Cumul.max(), reference.max()) > 0 else 1)

        ax1.xaxis.set_major_locator(MaxNLocator(10))
        ax1.yaxis.set_major_locator(MaxNLocator(10))
        ax2.yaxis.set_major_locator(MaxNLocator(10))

        # --- Grille ---
        ax1.grid(True, axis='x', linestyle=':', linewidth=1)

        # --- Légende combinée ---
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        fig.legend(lines1 + lines2, labels1 + labels2, loc="upper left")

        # --- Titre ---
        plt.title("Évolution des précipitations et cumul", fontsize=16)

        # --- Texte ---

        ax_caption = fig.add_subplot(gs[1])
        ax_caption.axis("off")

        if params["model"] in ('EULE', 'EULE2', 'STAT'):
            #Dans ce cas : params = model,path_fig, number_stitches, time_step, esp
            ax_caption.text(0.5, -1.5, f"Précipitation au sol, model {params["model"]},\n {params["number_stitches"]} mailles, pas de temps de {params["time_step"]} s, durée de la simulation {params["duree_sim"]},\n espece {params["esp"]}, temps de calculs : {params["time_run"]} s ", ha='center', va='bottom', wrap=True, fontsize=10)
        else :
            if params["deformability"] == "Yes" :
                params["deformability"]= "déformable"
                ax_caption.text(0.5, -1.5, f"Précipitation au sol, model {params["model"]}, {params["type_advance"]}, {params["deformability"]},\n {params["number_stitches"]} mailles, pas de temps de {params["time_step"]} s, durée de la simulation {params["duree_sim"]},\n espece {params["esp"]}, temps de calculs : {params["time_run"]} s ", ha='center', va='bottom', wrap=True, fontsize=10)
            else :
                params["deformability"] = "indéformable"
                ax_caption.text(0.5, -1.5, f"Précipitation au sol, model {params["model"]}, {params["type_advance"]}, {params["deformability"]},\n {params["number_stitches"]} mailles, nombre de bin : {params["number_bin"]}, pas de temps de {params["time_step"]} s, durée de la simulation {params["duree_sim"]},\n espece {params["esp"]}, temps de calculs : {params["time_run"]} s ", ha='center', va='bottom',wrap=True, fontsize=10)

        plt.subplots_adjust(bottom=0.2)

        # --- Sauvegarde ---
        file_location = chemin / Path("Précipitation")
        fig.savefig(str(file_location))
        plt.close()

    
class Affichage :

    def __init__(self,param_en_plus):

        # Param_en plus contient : model, path_fig, type_advance, deformable, number_stitches, time_step, esp, number_bin


        if param_en_plus["model"] in ('EULE', 'EULE2', 'STAT'): #Si le modèle c'est EULE, EULE2 ou STAT
            self.chemin = "." / Path(param_en_plus["path_fig"]) / Path(param_en_plus["model"]) / Path(f"Number_stitches_{param_en_plus["number_stitches"]}") / Path(f"duree_simu_{param_en_plus["duree_sim"]}") / Path(f"time_step_{param_en_plus["time_step"]}") / Path(f"espece_{param_en_plus["esp"]}") / Path (f"type_init_{param_en_plus["type_init"]}")
        else :
            if param_en_plus["deformability"] == "No" : #deformable = NO
                self.chemin = "."/Path(param_en_plus["path_fig"]) / Path(param_en_plus["model"]) / Path(param_en_plus["type_advance"])/Path(f"déformable_{param_en_plus["deformability"]}") / Path(f"Number_stitches_{param_en_plus["number_stitches"]}") / Path(f"duree_simu_{param_en_plus["duree_sim"]}") / Path(f"Number_bin_{param_en_plus["number_bin"]}") / Path(f"time_step_{param_en_plus["time_step"]}") / Path(f"espece_{param_en_plus["esp"]}") / Path (f"type_init_{param_en_plus["type_init"]}")
            else : #deformable = Yes
                self.chemin = "."/Path(param_en_plus["path_fig"]) / Path(param_en_plus["model"]) / Path(param_en_plus["type_advance"])/Path(f"déformable_{param_en_plus["deformability"]}") / Path(f"Number_stitches_{param_en_plus["number_stitches"]}") / Path(f"duree_simu_{param_en_plus["duree_sim"]}") / Path(f"time_step_{param_en_plus["time_step"]}") / Path(f"espece_{param_en_plus["esp"]}") / Path (f"type_init_{param_en_plus["type_init"]}")
        
        os.makedirs(self.chemin, exist_ok=True)
        self.params=param_en_plus
    
    def afficher (self,Concentration,Contenu,Precip, Quantiles) :

        
        Figure.sedimentation_time2(Concentration,"Concentration", self.chemin, self.params)
        Figure.sedimentation_time2(Contenu, "Contenu", self.chemin, self.params)
        Figure.precipitation(Precip, self.chemin, self.params, Quantiles)

        plt.show()
