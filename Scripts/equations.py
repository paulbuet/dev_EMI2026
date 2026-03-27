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
from scipy.integrate import quad, dblquad
from scipy.optimize import brentq
from scipy import integrate

### Classes and functions definition ###

class Eq :
    def __init__(self,esp):

        #On initialise les constantes en fonction de l'espèce 
    
        if esp=="i":
            self.a=0.82
            self.b=2.5
            self.c=800
            self.d=1
            self.alpha=3
            self.nu=3

        if esp=='s':
            self.a=0.02
            self.b=1.9
            self.c=5.1
            self.d=0.27
            self.alpha=1
            self.nu=1
            self.C=5
            self.x=1

        if esp=='g':
            self.a=19.6
            self.b=2.8
            self.c=124
            self.d=0.66
            self.alpha=1
            self.nu=1
            self.C=5e5
            self.x=-0.5

        if esp=='r':
            self.a=524
            self.b=3
            self.c=842
            self.d=0.66
            self.alpha=1
            self.C=8e6
            self.x=-1            
            self.nu=1

        if esp=='c':
            self.a=524
            self.b=3
            self.c=3.2e7
            self.d=2
            self.alpha=1
            self.nu=3         

    def Gamma(self, diametre, lam) :
        """
        fonction qui retourne la proportion de particule (float, sans unité) d'un certain diamètre (float, en m) en fonction du lambda (float, sans unité)
        """

        return (self.alpha/gamma(self.nu))*(lam**(self.alpha*self.nu))*(diametre**(self.alpha*self.nu-1))*np.exp(-((lam*diametre)**self.alpha))
   
    def Gamma_fois_masse(self, diametre, lam) :
        """
        fonction qui retourne la masse d'une proportion de particule (float, en kg) d'un certain diamètre (float,en m) en fonction du lambda (float) 
        """
        propor = self.Gamma(diametre, lam)   # proportion de particules du diamètre

        masse_part = self.Masse(diametre)   # masse d'une de ces particules

        return propor * masse_part     # Soit en kg

    def G(self, p) :
        """
        fonction qui retourne en fonction de p la fonction G(p), soit un rapport de loi gamma
        """

        gamma1 = gamma(self.nu+p/self.alpha)

        gamma2 = gamma(self.nu)

        return gamma1 / gamma2      # Sans unité
    
    def Lanbda(self, rho_r, N):
        """
        fonction qui retourne Lambda (float, sans unité) en fonction du contenu (float, en kg.m-3) et de la concentration (float, en m-3)
        """

        return (((rho_r)/(self.a*N*self.G(self.b)))**(-1/self.b))  # On est ici en 2 moments
    
    def Liste_Lanbda(self, rho_r, Liste_N):
        """
        fonction qui retourne Lambda (np.array, sans unité) sur chaque maille en fonction du contenu (np.array, en kg.m-3)
            et de la concentration (np.array, en m-3) sur toute la colonne
        """

        Liste_N=np.array(Liste_N)
        Ga=self.G(self.b)

        Liste_N = np.array([np.nan if x == 0 else x for x in Liste_N])
        Liste_lanbda=(rho_r/(Liste_N*self.a*Ga))**(-1/self.b)
        Liste_lanbda = np.array([0 if x == np.nan else x for x in Liste_lanbda])   # On mets np.nan si il n'y a pas de concentration sur la maille

        return Liste_lanbda
    

    def Masse(self, diametre):
        """
        fonction qui en fonction du diamètre (float, en m) renvoie la masse de la particule (float, en kg)
        """
        return (self.a*(diametre**self.b))
    
    def Vitesse(self,diametre):
        """
        fonction qui en fonction du diamètre (float, en m) renvoie la vitesse de chute (float, en m.s-1)
        """

        return (self.c*(diametre**self.d))
    
    def calcul_diametre(self, liste_h, dt):
        """
        fonction qui retourne les diamètres d'une particule (np.array, en m) si elle veut parcourir les distances de liste_h (np.array, en m)
            dans un temps dt (float, en s)
        """

        liste_h=np.array(liste_h)

        # On calcule la vitesse pour ensuite calculer les diamètres
        vitesse = liste_h/dt
        liste_d = ( vitesse/self.c)**(1/self.d)

        return liste_d
    
    def Calcul_integrale_conc(self, liste_d, lam):
        """
        fonction qui retourne les proportions de concentrations de particules en fonction des intervalles de diamètres
        et de la PDF qui dépends de Lambda
        """
        propor_tot=[]

        for i in range(len(liste_d)-1):

            propor, err = integrate.quad(self.Gamma, liste_d[i], liste_d[i+1], args=(lam,))
            propor_tot.append(propor)

        return propor_tot
    
    def Calcul_integrale_mass(self, liste_d, lam):
        """
        semblable à la fonction précédente, ressort les masses respectives de ces proportions
        """
        masse_tot=[]

        for i in range(len(liste_d)-1):

            masse, err = integrate.quad(self.Gamma_fois_masse, liste_d[i], liste_d[i+1], args=(lam,))
            masse_tot.append(masse)

        # Attention, ce n'est pas la masse totale par maille, c'est en kg.m-3 mais pour une proportion normée, 
        # pas pour le nombre véritable par tranche de diamètre
        # Il faut alors multiplier par le nombre de particule initiale pour obtenir la masse totale par tranche de diamètre
    
        return masse_tot  
    
    def contenu_to_conc(self,rho_r):
        """
        fonction qui renvoie à l'aide du 1er moment la concentration de particule par maille en fonction de son contenu
        """

        return self.C * (rho_r/(self.a*self.C*self.G(self.b)))**(self.x/(self.x-self.b))
    
    def Liste_Lanbda_1_mom(self, liste_rho_r):
        liste_rho_r=np.array(liste_rho_r)
        Gam=self.G(self.b)
        l_lam_1_mom=(liste_rho_r/(self.a*self.C*Gam))**(1/(self.x-self.b))
        return l_lam_1_mom
    
    def calcul_maille_arrivee(self, h1, h2, h3, h4, N, type, dt, lam) :
        #h1 : hauteur de l'interface du haut de la maille de départ
        #h2 : hauteur de l'interface du bas de la maille de départ
        #h3 : hauteur de l'interface du haut de la maille de d'arrivée
        #h4 : hauteur de l'interface du bas de la maille de d'arrivée
        #N : concentration actuelle dans la maille
        #type : choisir "concentration" pour faire tomber la concentration ou "masse" pour la masse
        #dt : pas de temps
        f = lambda D,h,sigma,beta,N : N*sigma*D**beta*self.alpha/gamma(self.nu)*lam**(self.alpha*self.nu)*D**(self.alpha*self.nu-1)  * np.exp(-(lam*D)**self.alpha)
            
        if h2==h4:
            if type == "concentration" : 
                integrale = -dblquad(f, h3, h4, lambda h : (((h-h1)/(dt*self.c))**(1/self.d)), lambda h : 0,args=(1,0,N))[0]
                #print('val int : ', integrale)
                new_val = integrale/(h2-h1)

                return new_val

            if type == "masse" :
                integrale = -dblquad(f, h3, h4, lambda h : (((h-h1)/(dt*self.c))**(1/self.d)), lambda h : 0,args=(self.a,self.b,N))[0]
                #print('val int : ', integrale)
                new_val = integrale/(h2-h1)

                return new_val

        if type == "concentration" : 
            integrale = -dblquad(f, h3, h4, lambda h : (((h-h1)/(dt*self.c))**(1/self.d)), lambda h : (((h-h2)/(dt*self.c))**(1/self.d)),args=(1,0,N))[0]
            #print('val int : ', integrale)
            new_val = integrale/(h2-h1)

            return new_val

        if type == "masse" :
            integrale = -dblquad(f, h3, h4, lambda h : (((h-h1)/(dt*self.c))**(1/self.d)), lambda h : (((h-h2)/(dt*self.c))**(1/self.d)),args=(self.a,self.b,N))[0]
            #print('val int : ', integrale)
            new_val = integrale/(h2-h1)

            return new_val
        
    def calcul_precip(self,h3, h4, N, dt, lam) :
        #h3 : hauteur de l'interface du haut de la maille de d'arrivée
        #h4 : hauteur de l'interface du bas de la maille de d'arrivée
        #N : concentration actuelle dans la maille
        #dt : pas de temps
        f = lambda D,h,sigma,beta,N : N*sigma*D**beta*self.alpha/gamma(self.nu)*lam**(self.alpha*self.nu)*D**(self.alpha*self.nu-1)  * np.exp(-(lam*D)**self.alpha)
        
        integrale = -dblquad(f, h3, h4, lambda h : (((np.inf)/(dt*self.c))**(1/self.d)), lambda h : (((h)/(dt*self.c))**(1/self.d)),args=(self.a,self.b,N))[0]
        #print('val int : ', integrale)

        return integrale
      
    def calcul_percentil_chute(self, lam, h_tot): #lam : valeur de lambda dans la maille, h_tot : hauteur totale de la colonne d'eau
        
        


        D=np.linspace(1e-5, 1, 10000)
        f=self.Gamma_fois_masse(D, lam)
        cdf = np.cumsum((f[:-1]+f[1:])/2 * np.diff(D))
        cdf = np.insert(cdf, 0, 0)
        cdf/=cdf[-1]
        q=np.arange(0.02, 1, 0.02)
        q = np.append(q, 1)
        D_q = np.interp(q, cdf, D)

        q = np.insert(q, 0, 0)
        D_q = np.insert(D_q, 0, 0)
        
        print("ça c'est q : ", q, " et ça c'est D_q : ", D_q)
        Liste_masse=[]
        for i in range(len(q)-1):
            result, error =  integrate.quad(self.Gamma_fois_masse, D_q[i], D_q[i+1], args=lam)
            fact = - result/(D_q[i]-D_q[i+1])    
            Liste_masse.append(fact)
        Liste_masse = np.array(Liste_masse)
        print("ça c'est la somme de ma liste de masse : ",np.sum(Liste_masse))

        Liste_vitesse=[]
        for i in range(len(q)-1):
            result, error = integrate.quad(self.Vitesse, D_q[i], D_q[i+1])
            fact = - result/(D_q[i]-D_q[i+1])
            Liste_vitesse.append(fact)
        Liste_vitesse = np.array(Liste_vitesse)

        Liste_temps_chute = h_tot/Liste_vitesse
        Liste_temps_chute = Liste_temps_chute[::-1]

        return q, Liste_temps_chute
    

    def sedimentation_times(self, N, lam, h_tot,number_stitches):

        """
        Calcule le temps théorique pour lequel une fraction p de la masse
        a sédimenté.

        Paramètres
        ----------
        p_values : array-like
            Fractions de masse (ex: np.linspace(0.02, 0.98, 49))
        N : float
            Concentration initiale
        lam : float
            Paramètre lambda de la loi gamma
        alpha, nu : float
            Paramètres de la loi gamma
        a, b : float
            Loi de masse : m(D) = a D^b
        c, d : float
            Loi de vitesse : v(D) = c D^d
        H : float
            Hauteur de la colonne
        D_min, D_max : float
            Bornes d'intégration

        Retour
        ------
        t_values : array
            Temps associés à chaque p
        D_values : array
            Diamètres seuils associés
        """

        p_values = np.linspace(0.02, 0.99, 50)
        D_min=1e-6
        D_max=0.015
        print("lambda : ", lam)
        # --- Loi gamma généralisée ---
        

        # --- Masse totale ---
        rho_r_tot = N * quad(self.Gamma_fois_masse, D_min, D_max, args=(lam))[0]
        print("rho_r_tot:",rho_r_tot)

        t_values = []
        D_values = []
        mass_in_time=[]

        for p in p_values:


            # Fonction à annuler : masse au-dessus de Dx - p*Mtot
            def func(Dx):
                integral = N * quad(self.Gamma_fois_masse, Dx, D_max,args=(lam))[0]
                return integral - p * rho_r_tot

            # Résolution de Dx
            Dx = brentq(func, D_min, D_max)
            D_values.append(Dx)

            mass_in_time.append(p*rho_r_tot*h_tot/number_stitches)

        D_values = np.array(D_values)
        print("D_values:",D_values)

        # Vitesse et temps
        v = self.c * D_values**self.d
        t_values = h_tot / v
        

        

        print(np.array(mass_in_time), np.array(t_values))

        return np.array(mass_in_time), np.array(t_values)


class Selection:

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

    
    def Dmin_Dmax(self, lam):

        """ fonction qui en fonction de lambda (float, sans unité) retourne le diamètre minimale au dessous duquel il reste seulement 1% de partucles
            et le diamètre maximale au dessus duquel il reste seulement 1% des particules.
        """

        def F(D):
            return quad(self.Eq_config.Gamma, 0, D, args=(lam))[0]
 
        D_high = 1.0 / lam # échelle naturelle
        while F(D_high) < 0.999:
            D_high *= 2

        Dmin = brentq(lambda D: F(D) - 0.01, 0, D_high)       # On cherche D tel que l'intégrale de 0 à D fasse 1% de 1 soit 0.01 car normé
        Dmax = brentq(lambda D: F(D) - 0.999, 0, D_high)        # Même chose pour Dmax mais tel que 99% des particules soient plus petites

        return Dmin, Dmax
    
    def Classe_D_N(self, nb_classes, Dmin, Dmax, N, lam):
        list_interv = np.linspace(0,nb_classes,nb_classes+2)
        list_interv = [sum(list_interv[:i+1]) for i in range(nb_classes+1)]
        Result=[]
        Intervalle=(Dmax-Dmin)/((nb_classes+1)*nb_classes/2)
        for i in range(nb_classes):
            Di=Dmin+list_interv[i]*Intervalle
            
            P_i=quad(self.Eq_config.Gamma, Dmin+list_interv[i]*Intervalle,  Dmin+list_interv[i+1]*Intervalle, args=(lam))[0]
            Ni=N*P_i
            Result.append([Di, Ni]) #Liste de deux paramètres : diamètre moyen, quantité associé par rapport au nombre total de particule.
        return Result
    
    def Classe_D_rho_r(self, nb_classes, Dmin, Dmax, N, lam):
        list_interv = np.linspace(0,nb_classes,nb_classes+2)
        list_interv = [sum(list_interv[:i+1]) for i in range(nb_classes+1)]
        Result=[]
        Intervalle=(Dmax-Dmin)/((nb_classes+1)*nb_classes/2)
        for i in range(nb_classes):
            Di=Dmin+list_interv[i]*Intervalle
            
            P_i=quad(self.Eq_config.Gamma_fois_masse, Dmin+list_interv[i]*Intervalle, Dmin+list_interv[i+1]*Intervalle, args=(lam))[0]
            rho_r_i=N*P_i
            Result.append([Di, rho_r_i]) #Liste de deux paramètres : diamètre moyen, quantité associé par rapport au nombre total de particule.
        return Result