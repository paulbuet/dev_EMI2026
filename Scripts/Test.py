from math import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

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

    def Gamma(self, diametre, lanbda) :
        return (self.alpha/gamma(self.nu))*(lanbda**(self.alpha*self.nu))*(diametre**(self.alpha*(self.nu-1)))*np.exp(-((lanbda*diametre)**self.alpha))
   
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


eq_rain = eq("r")

n_rain = np.linspace(10, 10, 100)
D_rain = np.linspace(0.001, 0.1, 100)
m_rain = eq_rain.Masse(D_rain)
N=sum(n_rain)
print(N)

rho_r = 1
lanbda = eq_rain.Lanbda(rho_r,N)

def integrand(D):
    return eq_rain.Gamma(D, lanbda)

result, error = quad (integrand, 0, 100)
print('résultat : ', result)

#a=quad(eq_rain.Gamma, 0, 50, args=(lanbda))
#print(a)



