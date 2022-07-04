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
        self.ray = 0.5
        self.rest_spring_l = 2*self.ray
        self.theta = np.pi
        self.dt = 0.05
    
    def generate_bacterium(self):
        Cell1 = cell.Cell(X=np.array([1,5]),mass=0.1,ray=self.ray)
        self.bacteria = [bact.Bacterium(N=1,Cells=[Cell1],Torques=[],l=self.rest_spring_l,theta=self.theta,color=(139,0,0))]
        self.N_bacteria +=1
        Cell2 = cell.Cell(X=np.array([Cell1.X[0]+self.bacteria[0].spring_rest_l+10,Cell1.X[1]]),mass=0.1,ray=self.ray)
        self.bacteria[0].Cells.append(Cell2)
        self.bacteria[0].N+=1
        Cell3 = cell.Cell(X=np.array([Cell2.X[0]+self.bacteria[0].spring_rest_l,Cell2.X[1]+5]),mass=0.1,ray=self.ray)
        self.bacteria[0].Cells.append(Cell3)
        self.bacteria[0].N+=1

        Cell1 = cell.Cell(X=np.array([2,4]),mass=0.1,ray=self.ray)
        Cell2 = cell.Cell(X=np.array([Cell1.X[0],Cell1.X[1]+self.rest_spring_l+2]),mass=0.1,ray=self.ray)
        Cell3 = cell.Cell(X=np.array([Cell2.X[0],Cell2.X[1]+self.rest_spring_l]),mass=0.1,ray=self.ray)
        Cell4 = cell.Cell(X=np.array([Cell3.X[0]+self.rest_spring_l+5,Cell3.X[1]]),mass=0.1,ray=self.ray)
        Cell5 = cell.Cell(X=np.array([Cell4.X[0],Cell4.X[1]+self.rest_spring_l+4]),mass=0.1,ray=self.ray)
        self.bacteria.append (bact.Bacterium(N=5,Cells=[Cell1,Cell2,Cell3,Cell4,Cell5],Torques=[],l=self.rest_spring_l,theta=self.theta,color=(0,128,0)))

    def Move_bacteria(self):
        """Calculates the velocities of each cells of all the bacteria
        and Integrates them to determine the positions"""

        for bact in self.bacteria:
            bact.spring_velocity()
            bact.Euler_explicit(self.dt)