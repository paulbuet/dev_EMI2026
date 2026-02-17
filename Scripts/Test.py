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
    
    def Dmin_Dmax(self, lam):
 
        def F(D):
            return quad(self.Gamma, 0, D, args=(lam))[0]
 
        D_high = 1.0 / lam # échelle naturelle
        while F(D_high) < 0.999:
            D_high *= 2

        Dmin = brentq(lambda D: F(D) - 0.01, 0, D_high)
        Dmax = brentq(lambda D: F(D) - 0.99, 0, D_high)
 
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

Affichage.Affichage_Concentration(Concentration, "concentration")
Affichage.Affichage_Precipitation(Precip)


eq_rain = eq("c")

lam=0.00001
N=151515

Pi=quad(eq_rain.Gamma, 0, 500, args=(lam))[0]
print("Intégrale entre 0 et 500 : ", Pi)
dmin, dmax = eq_rain.Dmin_Dmax(lam)
print("Integrale entre 0 et Dmax : ", quad(eq_rain.Gamma, 0, dmax, args=(lam))[0])
print("Integrale entre Dmin et 500 : ", quad(eq_rain.Gamma, dmin, 500, args=(lam))[0])
print("Integrale entre Dmin et Dmax : ", quad(eq_rain.Gamma, dmin, dmax, args=(lam))[0])


Resultat=eq_rain.Classe_D(10, dmin, dmax, N, lam)
print("Pour la grèle, avec lambda=0.5, en fixant 6 différentes classes et un nombres totales de particules à 1000 on obtient la répartition : ", Resultat)
somme=0
for i in range(10):
    somme+=Resultat[i][1]
print("Le résultat de la somme totale des particules est : ", somme)






