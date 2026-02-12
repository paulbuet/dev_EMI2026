# Le but de ce fichier est de mettre toutes les fonctions qu'on peut potentiellement appelés pour faire tourner nos futur codes
# Ensuite il suffira de faire un appelle du fichier et de la fonction qu nous interesse
# Exemple avec cette fonction test que j'appelle dans le code Sedimentation.py
from math import *
import time
import numpy as np
import matplotlib.pyplot as plt
class eq :
    def test() :
        return 'toto'
    def Gamma(diametre, alpha, lanbda, nu) :
        return (alpha/gamma(diametre))*(lanbda**(alpha*nu))*(diametre**(alpha*(nu-1)))*exp(-((lanbda*diametre)**alpha))

param1=1
param2=1
param3=1



a=np.linspace(1, 20, 100)
b=[]
for i in range(100):
    b.append(eq.Gamma(a[i], param1, param2, param3))


plt.plot(a, b, '--', color='r')


param1=1
param2=1
param3=1.2

a=np.linspace(1, 20, 100)
b=[]
for i in range(100):
    b.append(eq.Gamma(a[i], param1, param2, param3))


plt.plot(a, b, '--', color='b')


param1=5
param2=0.5
param3=5

a=np.linspace(1, 20, 100)
b=[]
for i in range(100):
    b.append(eq.Gamma(a[i], param1, param2, param3))


plt.plot(a, b, '--', color='y')
plt.show()
