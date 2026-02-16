from math import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import brentq

class eq :
    def __init__(self,esp):
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
        if esp=='g':
            self.a=19.6
            self.b=2.8
            self.c=124
            self.d=0.66
            self.alpha=1
            self.nu=1
        if esp=='r':
            self.a=524
            self.b=3
            self.c=842
            self.d=0.66
            self.alpha=1
            self.nu=1
        if esp=='c':
            self.a=524
            self.b=3
            self.c=3.2e7
            self.d=2
            self.alpha=1
            self.nu=3         

    def Gamma(self, diametre, lam) :
        return (self.alpha/gamma(self.nu))*(lam**(self.alpha*self.nu))*(diametre**(self.alpha*self.nu-1))*np.exp(-((lam*diametre)**self.alpha))
   
    def G(self, p) :
        return (gamma(self.nu+p/self.alpha)/gamma(self.nu))
    
    def Lanbda(self, rho_r, N):
        return (((rho_r)/(self.a*N*self.G(self.b)))**(-1/self.b))
    
    def Masse(self, diametre):
        return (self.a*(diametre**self.b))
    
    def Vitesse(self,diametre):
        return(self.c*(diametre**self.d))
    
    def Calcul_rho_r(self, m, n):
        return np.dot(m,n)
    
    def Dmin_Dmax(self, lam):

 
        # fonction de répartition
        def F(D):
            return quad(self.Gamma, 0, D, args=(lam))[0]
 
        # équations à résoudre
        def f_min(D):
            return F(D) - 0.01
 
        def f_max(D):
            return F(D) - 0.99
 
        Dmin = brentq(f_min, 0, 500)
        Dmax = brentq(f_max, 0, 500)
 
        return Dmin, Dmax
    
    def Classe_D(self, nb_classes, Dmin, Dmax, N, lam):
        Result=[]
        Intervalle=(Dmax-Dmin)/nb_classes
        for i in range(nb_classes):
            Di=(1+2*i)*Intervalle/2 + Dmin
            
            P_i=quad(self.Gamma, Dmin+i*Intervalle, Dmin+(i+1)*Intervalle, args=(lam))[0]/0.98
            print(i, P_i)
            Ni=N*P_i
            Result.append([Di, Ni]) #Liste de deux paramètres : diamètre moyen, quantité associé par rapport au nombre total de particule.
        return Result
    





eq_rain = eq("i")

lam=5
N=1000

Pi=quad(eq_rain.Gamma, 0, 500, args=(lam))[0]
print("Intégrale entre 0 et 500 : ", Pi)
dmin, dmax = eq_rain.Dmin_Dmax(lam)
print("Integrale entre 0 et Dmax : ", quad(eq_rain.Gamma, 0, dmax, args=(lam))[0])
print("Integrale entre Dmin et 500 : ", quad(eq_rain.Gamma, dmin, 500, args=(lam))[0])
print("Integrale entre Dmin et Dmax : ", quad(eq_rain.Gamma, dmin, dmax, args=(lam))[0])


Resultat=eq_rain.Classe_D(6, dmin, dmax, N, lam)
print("Pour la grèle, avec lambda=0.5, en fixant 6 différentes classes et un nombres totales de particules à 1000 on obtient la répartition : ", Resultat)
somme=0
for i in range(6):
    somme+=Resultat[i][1]
print("Le résultat de la somme totale des particules est : ", somme)






