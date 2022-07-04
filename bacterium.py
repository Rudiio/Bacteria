import numpy as np
import cell 
from numpy.linalg import norm


class Bacterium:
    """ Class representing a bacterium 
    It composed of :
    - the number of cells in the Bacterium: N
    - a list of cells: Cells
    - a list of torques: Torques 
    - the rest torque theta
    - the rest length of the springs l 
    - the stiffness constants ks,kt1, kt2"""

    def __init__(self,N=0,Cells : list[cell.Cell] = [],Torques = [],l = 1.0,theta = np.pi,color = (0,0,0)):
        """Constructor for the class Bacterium"""
        self.N = N
        self.Cells = Cells
        self.Torques = Torques
        self.theta = theta
        self.spring_rest_l = l
        self.ks = 1
        self.kt1 = 1
        self.kt2 = 1
        self.color = color

    def __str__(self):
        """Display the bacterium informations whit the print function"""
        print(f"N = {self.N}")
        for i in self.Cells :
            print(i)
        
        for i in self.Torques:
            print(i)
        return("")
    
    def spring_velocity(self):
        """Calculate the velocity that comes from the spring forces/torques
        of each cell of the bacterium"""

        # Loop on all the cells
        for k in range(self.N):
            self.Cells[k].V = self.linear_spring(k)

    def linear_spring(self,k):
        """Calculates the velocity created by the linear springs"""
        if k==0:
            Xj = self.Cells[k].X
            Xj1 = self.Cells[k+1].X

            return (self.ks/(self.spring_rest_l**2)*(norm(Xj1 - Xj) - self.spring_rest_l)*
                    (Xj1 - Xj)/(norm(Xj1 - Xj))) 
        elif k==self.N-1:
            Xj = self.Cells[k].X
            Xj_1 = self.Cells[k-1].X

            return (-self.ks/(self.spring_rest_l**2)*(norm(Xj - Xj_1) - self.spring_rest_l)*
                    (Xj - Xj_1)/norm(Xj-Xj_1))
        else:
            Xj = self.Cells[k].X
            Xj1 = self.Cells[k+1].X
            Xj_1 = self.Cells[k-1].X

            return (self.ks/(self.spring_rest_l**2)*(norm(Xj1 - Xj) - self.spring_rest_l)*
                    (Xj1 - Xj)/(norm(Xj1 - Xj)) -self.ks/(self.spring_rest_l**2)*(norm(Xj - Xj_1) - self.spring_rest_l)*
                    (Xj - Xj_1)/norm(Xj-Xj_1))

    def torsion_spring(self,k):
        """Calculates the velocity created by the linear springs"""

        return 0

    def Euler_explicit(self,dt):
        """Calculates the position of all the cells by using 
        an Euler explicit integrator"""

        for i in range(self.N):
            self.Cells[i].X = self.Cells[i].X + dt*self.Cells[i].V
        

if __name__=="__main__":
    print("test")