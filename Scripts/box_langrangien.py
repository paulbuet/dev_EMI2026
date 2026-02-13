
from fonctions import Initial_Cond

class Model_bl():

    def __init__(self, type_advance,number_stitches,deformable):
        """
        Here we initialise the non-spatial fixed parameters and allow important variables 
        to travel between functions.
        """

        self.type_advance = type_advance
        self.number_stitches = number_stitches
        self.deformable = deformable

        self.length_sim = 1000  # in seconds

        self.speed_max = 6  # in m/s



    def run():
        
        condi_init = Initial_Cond(self.number_stitches)

        grid = condi_init.grid()



        