import numpy as np
import cell 


class Bacterium:
    """ Class representing a bacterium 
    It composed of :
    - the number of cells in the Bacterium: N
    - a list of cells: Cells
    - a list of torques: Torques """

    def __init__(self,N=0,Cells = [],Torques = []):
        """Constructor for the class Bacterium"""
        self.N = N
        self.Cells = Cells
        self.Torques = Torques

    def __str__(self):
        """Display the bacterium informations whit the print function"""
        print(f"N = {self.N}")
        for i in self.Cells :
            print(i)
        
        for i in self.Torques:
            print(i)
        return("")

if __name__=="__main__":
    print("test")