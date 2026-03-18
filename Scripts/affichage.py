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

    def Affichage_Concentration(Concentration, typ, model, path_fig): #type="concentration" ou "masse"
        Temps_simu=len(Concentration)
        nb_boites=len(Concentration[0])
        Concentration=np.array(Concentration)
        Transpose=Concentration.T
        plt.figure(figsize=(Temps_simu, nb_boites))
        orig_map=plt.cm.get_cmap('gist_ncar')
        reversed_map = orig_map.reversed()
        norm = Normalize(vmin=0,vmax = max(Concentration[0]))
        plt.pcolormesh(Transpose,cmap=reversed_map, norm=norm)
        plt.title(f"{typ} : évolution dans le temps", fontsize=22)
        plt.xlabel("Temps", fontsize=18)
        plt.ylabel("Mailles du modèle", fontsize=18)
        plt.colorbar()
        file_location = Path(path_fig) / Path(model) / Path(typ) 
        plt.savefig(str(file_location))


    def Affichage_Precipitation(Precip, model, path_fig):
        Precip=np.array(Precip)
        liste=np.zeros(len(Precip))
        Cumul=[]

        for i in range(len(Precip)):
            liste[i]=1
            Cumul.append(np.dot(Precip, liste))
        fig, ax1 = plt.subplots()
        ax = plt.gca()
        ax.xaxis.set_major_locator(MultipleLocator(1))
        ax.yaxis.set_major_locator(MultipleLocator(1))
        ax.set_xlim(left=0)
        x_max = len(Precip)-1
        ticks = np.linspace(0, x_max, 11)
        ax.set_xticks(ticks)
        ax.set_xticklabels([f"{t:.2f}" for t in ticks])
        ax.yaxis.set_major_locator(MaxNLocator(10))

        time=np.linspace(1, len(Precip), len(Precip))
        time2 = time -(time[1]-time[0])/2
        ax1.bar(time2, Precip, color="blue", label="Précip par pas de temps")
        ax1.set_xlabel("temps", fontsize=18)
        ax1.set_ylabel("Précip par pas de temps", fontsize=18)
        ax2 = ax1.twinx()
        ax2.plot(time, Cumul, '--', color="red", label="Cumul")
        ax2.set_ylabel('Cumul', fontsize=18)
        plt.title("Evolution des précipitations par pas de temps et cumulée", fontsize=22)
        plt.grid(axis='x', which='major', markevery=[1,2,3],lw=2, ls=':')
        fig.legend(loc=2)
        file_location = Path(path_fig) / Path(model) / "Précipitations.png"
        plt.savefig(file_location)

    def Afficher () :
        plt.show()