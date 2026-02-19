from math import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import brentq
from matplotlib.ticker import MultipleLocator


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
    
    def Liste_Lanbda(self, rho_r, Liste_N):
        Liste_N=np.array(Liste_N)
        Ga=self.G(self.b)
        Liste_N = np.array([np.nan if x == 0 else x for x in Liste_N])
        Liste_lanbda=(rho_r/(Liste_N*self.a*Ga))**(-1/self.b)
        Liste_lanbda = np.array([0 if x == np.nan else x for x in Liste_lanbda])
        print("Liste lambda : ", Liste_lanbda)
        return Liste_lanbda
    
    def Liste_Dmin_Dmax(self, rho_r, Liste_N):
        Liste_Lanbda=self.Liste_Lanbda(rho_r, Liste_N)
        Liste_Lanbda=np.array(Liste_Lanbda)
        indices_nan = np.array([i for i, x in enumerate(Liste_Lanbda) if np.isnan(x)])
        Liste_Lanbda_sans_nan = np.array([x for x in Liste_Lanbda if not np.isnan(x)])
        Liste_Dm = [self.Dmin_Dmax(elem) for elem in Liste_Lanbda_sans_nan]
        Liste_Dm_avec_nan = Liste_Dm.copy()
        for i in sorted(indices_nan, reverse=True):
            Liste_Dm_avec_nan.insert(i, (np.nan, np.nan))
        print(Liste_Dm_avec_nan)
        Liste_Dm_avec_nan=np.array(Liste_Dm_avec_nan)
        return Liste_Dm_avec_nan

    def Liste_Vitesse(self, rho_r, Liste_N):
        Liste_Dm=self.Liste_Dmin_Dmax(rho_r, Liste_N)
        Vitesse= self.a*(Liste_Dm**self.b)
        return Vitesse

    def Dmin_Dmax(self, lam):
 
        def F(D):
            return quad(self.Gamma, 0, D, args=(lam))[0]
 
        D_high = 1.0 / lam # échelle naturelle
        while F(D_high) < 0.999:
            D_high *= 2

        Dmin = brentq(lambda D: F(D) - 0.1, 0, D_high)
        Dmax = brentq(lambda D: F(D) - 0.9, 0, D_high)
 
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
    
class Affichage :

    def Affichage_Concentration(Concentration, typ): # type = "concentration" ou "masse"
        Temps_simu=len(Concentration)
        nb_boites=len(Concentration[0])
        #time=np.linspace(1, Temps_simu, Temps_simu)
        Concentration=np.array(Concentration)
        Transpose=Concentration.T
        print(Temps_simu, nb_boites)
        plt.figure(figsize=(Temps_simu, nb_boites))
        plt.pcolormesh(Transpose,cmap='binary')
        plt.title(f"Evolution de la {typ} de particules dans le temps")
        plt.ylabel("Mailles du modèle")
        plt.xlabel("Temps")
        plt.savefig(f"results/{typ}.png")
        plt.show()

    def Affichage_Precipitation(Precip):
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

        time=np.linspace(1, len(Precip), len(Precip))
        time2 = time -(time[1]-time[0])/2
        ax1.bar(time2, Precip, color="blue", label="Précip_horaire")
        ax1.set_xlabel("temps")
        ax1.set_ylabel("Cumul")
        ax2 = ax1.twinx()
        ax2.plot(time, Cumul, '--', color="red", label="Cumul")
        ax2.set_ylabel('Précipitations horaires')
        plt.title("Evolution de la précip en horaire et cumulée")
        plt.grid(axis='x', which='major', markevery=[1,2,3],lw=2, ls=':')
        fig.legend()
        plt.savefig("results/Précipitations.png")
        plt.show()



    


Concentration=[[3, 2, 3, 5], [6, 7, 2, 28], [1, 8, 1, 1], [1, 5, 2, 6], [1, 4, 7, 6], [1, 4, 3, 5], [1, 5, 2, 7]]
Precip=[0, 0, 0, 0, 1, 1, 2, 3, 9, 2, 1, 1, 1, 0, 0]

#Affichage.Affichage_Concentration(Concentration, "concentration")
#Affichage.Affichage_Precipitation(Precip)


a=[0, 1, 1, 3, 3, 6, 3, 7, 0, 7 , 1, 2]
a=[np.nan if x==0 else x for x in a]
print(a)
a=np.array(a)
b=0.2*a
print(b)
c=4/a
print(c)


eq_rain = eq("r")

lam=980.6
N=1000
rho_r=1

Liste_N=[0, 7, 6, 3, 6, 12, 24, 23, 47, 48, 59, 3, 89, 2, 6, 75, 1]
print(eq_rain.Liste_Vitesse(rho_r, Liste_N))











