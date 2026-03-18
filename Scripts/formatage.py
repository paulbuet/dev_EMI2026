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
from scipy import integrate

from equations import Eq

class Formatage:
    
    def __init__(self,esp):
        

        self.Eq_config = Eq(esp)

        # On initialise les constantes requises suivant l'espèce

        if esp != 'i' and esp!= 'c':

            self.a=self.Eq_config.a
            self.b=self.Eq_config.b
            self.c=self.Eq_config.c
            self.d=self.Eq_config.d
            self.alpha=self.Eq_config.alpha
            self.C=self.Eq_config.C
            self.x=self.Eq_config.x            
            self.nu=self.Eq_config.nu

        else:
            self.a=self.Eq_config.a
            self.b=self.Eq_config.b
            self.c=self.Eq_config.c
            self.d=self.Eq_config.d
            self.alpha=self.Eq_config.alpha         
            self.nu=self.Eq_config.nu


    def Epaiss_to_diam(self,h_interface):
        
        nb_interf = len(h_interface)

        tab_diam = [np.pad([h_interface[stitch_dep] - h_interface[stitch_arr] for stitch_arr in range(stitch_dep+1)],(0,nb_interf-stitch_dep-1))for stitch_dep in range(nb_interf)]
        
        return tab_diam
    
    def find_diameters_in_splittings (self,splitting_list):
        '''
        Returns an iterable with the first-found splitting in the splitting list given. 
        It aims at providing a single splitting list of [Di,Ni] sub-lists to acces the diameter values 
        without accessing other splitting lists containing basically the same values of diameter.
        '''
        i = 0
        found = False
        while i <= len(splitting_list) and found == False :
            if splitting_list[i] != " " :
                found = True
                yield splitting_list[i]
            i+= 1