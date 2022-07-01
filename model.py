import numpy as np
import cell
import bacterium as bact


class Model:
    """Class for the modelisation
    It is composed of :
    - a list of all the bacterias
    - time step : dt
    
    Process : 
    - Launch the calculus of the movement of the bacteria
    - integration of the velocities of the bacterias"""

    def __init__(self):
        """Constructor for the model class"""

        #Creation of the list of bacteria
        self.bacteria = []
        self.N_bacteria = 0
        self.rest_spring_l = 2
        self.theta = np.pi
        self.ray = 0.5
    
    def generate_bacterium(self):
        Cell1 = cell.Cell(X=np.array([7,5]),mass=0.1,ray=self.ray)
        self.bacteria = [bact.Bacterium(N=1,Cells=[Cell1],Torques=[],l=self.rest_spring_l,theta=self.theta)]
        self.N_bacteria +=1
        Cell2 = cell.Cell(X=np.array([Cell1.X[0]+self.bacteria[0].spring_rest_l,Cell1.X[1]]),mass=0.1,ray=self.ray)
        self.bacteria[0].Cells.append(Cell2)
        self.bacteria[0].N+=1
