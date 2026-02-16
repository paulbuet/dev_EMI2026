
# On utilisera dans la description du Box-Lagrangien les notations de Kato(1995)

from fonctions import Initial_Cond
from fonctions import eq

class Model_bl():
   
   def __init__(self, type_advance,number_stitches,deformable):
        """
        Here we initialise the non-spatial fixed parameters and allow important variables 
        to travel between functions. We also call the initialisation.
        """

        self.type_advance = type_advance
        self.number_stitches = number_stitches
        self.deformable = deformable

        self.length_sim = 30  # length of simulation in seconds

        self.speed_max = 6  # max of speed in m/s

        self.delta_t = 10 # length of time step in seconds

        self.nb_step = self.length_sim // self.delta_t  # number of time step

        self.nb_diam = 1 # number of type of diameter
   
        # Ajouter le calcul des différents diamètres dans une liste self.diameter

        """
        grid is in the form of a numpy array of size nb_mesh*(variables*2):

        concentration_bin_1 : number of particles with diameter 1 in each level: N0(1)  N1(1)  N2 (1)  N3(1) ...
        ...
        concentration_bin_i : number of particles with diameter i in each level: N0(i)  N1(i)  N2 (i)  N3(i) ...
        ...
        concentration_bin_n : number of particles with diameter n in each level: N0(n)  N1(n)  N2 (n)  N3(n) ...

        level : height of upper mesh interfaces: Z0  Z1  Z2  Z3  ...
        
        """
   
        condi_init = Initial_Cond(self.number_stitches)
   
        self.grid = condi_init.data()

    def calc_height(self,diameter:int,stitche:int):

        """
        Cette fonction prend en entrée le diamètre des particules, ainsi que leur maille,
        on calcule alors leur vitesse et on fait descendre la boite à cette vitesse la.
        diameter (type = int): indice du diamètre des particules que l'on descend
        stitches (type=int) : indice de la maille que l'on descend
        On return enfin les deux hauteurs.
        """

        speed = eq.Vitesse(self.diameter[diameter])   # On calcule la vitesse en m/s

        low_height = self.grid[stitche-1][1] - speed * self.delta_t  

        high_height = self.grid[stitche][1] - speed * self.delta_t

        return [low_height,high_height]
   
    def boundaries(self,list_height,n_start):

        """
        Cette fonction renvoie les indices des interfaces des mailles où ZTk et ZBk sont sous la forme L2,L2_1,L1,L1_1 
        Pour cela, elle itère les mailles en partant de celle qu'on a descendu afin de trouver le plus rapidement les interfaces.
        list_height (type = lis): liste comprenant ZBk et ZTk, ici k est en fait t_time.
        n_start (type = int): indice de la maille de laquelle on est descendu
        """
        
        L2_1 = n_start
        while self.grid[L2_1][1] >= list_height[1]:
            L2_1 -=1
        L2 = L2_1+1

        L1_1 = n_start
        while self.grid[L1_1][1] >= list_height[0]:
            L1_1 -=1
        L1 = L1_1+1

        return L2,L2_1,L1,L1_1


        



    def run():
        
        """
        On fait tourner le modèle pour chaque pas de temps, pour chaque mailles, pour chaque diamètres.
        On enregistre le profils des concentration en fonction de la hauteur à chaque pas de temps
        """ 
        

        for t_time in range(self.nb_step):

            for n_stitches in range(self.number_stitches):

                for diam in range(self.nb_diam):

                    height_fall = self.calc_height(diam,n_stitches)

                    L2,L2_1,L1,L1_1 = self.boundaries(height_fall,n_stitches)

                    self.grid[L2][2*(diam+1)] += (height_fall[1] - self.grid[L2_1][1])/(grid[L2][1]-grid[L2_1][1])* self.grid[n_stitches][2*(diam+1)]

                    self.grid[L1][2*(diam+1)] += (self.grid[L1][1] - height_fall[0] )/(grid[L1][1]-grid[L1_1][1])* self.grid[n_stitches][2*(diam+1)]












        