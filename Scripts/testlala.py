from scipy.integrate import dblquad
from math import *
import numpy as np

def calcul_maille_arrivee(h1, h2, h3, h4, N, type, dt, lam) :
        a=19.6
        b=2.8
        c=124
        d=0.66
        alpha=1
        nu=1
        C=5e5
        x=-0.5
        #h1 : hauteur de l'interface du haut de la maille de départ
        #h2 : hauteur de l'interface du bas de la maille de départ
        #h3 : hauteur de l'interface du haut de la maille de d'arrivée
        #h4 : hauteur de l'interface du bas de la maille de d'arrivée
        #N : concentration actuelle dans la maille
        #type : choisir "concentration" pour faire tomber la concentration ou "masse" pour la masse
        #dt : pas de temps
        f = lambda D,h,sigma,beta,N : N*sigma*D**beta*alpha/gamma(nu)*lam**(alpha*nu)*D**(alpha*nu-1)  * np.exp(-(lam*D)**alpha)


        if type == "concentration" : 
            integrale = -dblquad(f, h3, h4, lambda h : (((h-h1)/(dt*c))**(1/d)), lambda h : (((h-h2)/(dt*c))**(1/d)),args=(1,0,N))[0]
            print('val int : ', integrale)
            new_val = integrale/(h2-h1)

            return new_val

        if type == "masse" :
            integrale = -dblquad(f, h3, h4, lambda h : (((h-h1)/(dt*c))**(1/d)), lambda h : (((h-h2)/(dt*c))**(1/d)),args=(a,b,N))[0]
            print('val int : ', integrale)
            new_val = integrale/(h2-h1)

            return new_val

print(calcul_maille_arrivee(0,10,15,20,1000,"concentration",10,300))