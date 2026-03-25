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
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import MaxNLocator
from pathlib import Path
from matplotlib.colors import Normalize, LogNorm



class Affichage :

    def Affichage_Concentration(Concentration, typ, model, path_fig, type_advance,deformable): #type="concentration" ou "masse"
        Temps_simu=len(Concentration)
        nb_boites=len(Concentration[0])
        Concentration=np.array(Concentration)
        Transpose=Concentration.T
        plt.figure() # (figsize=(Temps_simu, nb_boites))
        orig_map=plt.cm.get_cmap('gist_ncar')
        reversed_map = orig_map.reversed()
        norm = Normalize(vmin=0,vmax = max(Concentration[0]))
        plt.pcolormesh(Transpose,cmap=reversed_map, norm=norm)
        plt.title(f"{typ} : évolution dans le temps", fontsize=22)
        plt.xlabel("Temps", fontsize=18)
        plt.ylabel("Mailles du modèle", fontsize=18)
        plt.colorbar()
        file_location = "."/Path(path_fig) / Path(model) / Path(type_advance)/Path(deformable) / Path(typ)
        plt.savefig(str(file_location))


    def Affichage_Precipitation(Precip, model, path_fig, type_advance,deformable, Quantiles, duree_sim):
        import numpy as np
        import matplotlib.pyplot as plt
        from matplotlib.ticker import MultipleLocator, MaxNLocator
        from pathlib import Path

        # --- Données ---
        Precip = np.array(Precip, dtype=float)

        reference = np.array(Quantiles[0], dtype=float)
        time3 = np.array(Quantiles[1], dtype=float)

        # Ajout du point initial (0,0)
        reference = np.insert(reference, 0, 0)
        time3 = np.insert(time3, 0, 0)

        # --- Cumul correct ---
        Cumul = np.cumsum(Precip)

        # --- Temps ---
        n = len(Precip)
        time = np.linspace(0, duree_sim, n)
        dt = time[1] - time[0]
        time2 = time - dt / 2 # pour centrer les barres

        # --- Figure ---
        fig, ax1 = plt.subplots(figsize=(10, 6))

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

        # --- Sauvegarde ---
        file_location = Path(path_fig) / model / type_advance / deformable / "Precipitation"
        file_location.mkdir(parents=True, exist_ok=True)

        fig.savefig(file_location / "precipitation.png", dpi=300, bbox_inches='tight')

    def Afficher () :
        plt.show()
