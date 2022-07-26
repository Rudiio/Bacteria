import numpy as np

class Disk:
    """ Class reprsenting a Disk a bacterium 
    It is composed of :
    - An array for its center : X in micrometers
    - An array for its velocity : V  micrometers/s
    - its mass in kg
    - its ray in micrometers
    """

    def __init__(self,X = np.zeros(2),V = np.zeros(2),mass=1,ray = 1):
        """Constructor for the class Disk"""
        self.X = X
        self.V = V
        self.mass = mass
        self.radius = ray
        self.next_bact = -1
        self.next_disk = -1

    def __str__(self):
        """Display the Disks information with the print fuction"""
        return (f"{self.X}")
    
    def vplus(self,v):
        """Modify the velocity of the bacterium"""
        self.V += v
        

if __name__=="__main__":
    print("salut")