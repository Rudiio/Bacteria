import numpy as np
import disk
import bacterium as bact
import random
from numpy.linalg import norm

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

        # Bacteria parameters and variables
        self.bacteria = []
        self.N= 0
        
        # disks parameters
        self.radius = 0.5

        # Springs parameters
        self.rest_spring_l = self.radius
        self.theta = np.pi

        # Simulation parameters
        self.disk_add_method = 2 # decide the position of the new disk
        self.dt = 0.5
        self.time = 0
    
    ### ---------------- Bacteria generators ---------------------------

    def generate_bacterium(self):
        # disk1 = disk.Disk(X=np.array([0,0]),mass=0.1,ray=self.radius)
        # self.bacteria = [bact.Bacterium(N=1,Disks=[disk1],l=self.rest_spring_l,theta=self.theta,color=(139,0,0))]
        # self.N+=1
        # disk2 = disk.Disk(X=np.array([disk1.X[0]+self.bacteria[0].spring_rest_l+10,disk1.X[1]]),mass=0.1,ray=self.radius)
        # self.bacteria[0].Disks.append(disk2)
        # self.bacteria[0].p_i+=1
        # disk3 = disk.Disk(X=np.array([disk2.X[0]+self.bacteria[0].spring_rest_l,disk2.X[1]+5]),mass=0.1,ray=self.radius)
        # self.bacteria[0].Disks.append(disk3)
        # self.bacteria[0].p_i+=1

        disk1 = disk.Disk(X=np.array([1,5]),mass=0.1,ray=self.radius)
        self.bacteria = [bact.Bacterium(N=1,Disks=[disk1],l=self.rest_spring_l,theta=self.theta,color=(139,0,0))]
        self.N+=1
        disk2 = disk.Disk(X=np.array([disk1.X[0]+self.bacteria[0].spring_rest_l+10,disk1.X[1]]),mass=0.1,ray=self.radius)
        self.bacteria[0].Disks.append(disk2)
        self.bacteria[0].p_i+=1
        disk3 = disk.Disk(X=np.array([disk2.X[0]+self.bacteria[0].spring_rest_l,disk2.X[1]+5]),mass=0.1,ray=self.radius)
        self.bacteria[0].Disks.append(disk3)
        self.bacteria[0].p_i+=1

        disk1 = disk.Disk(X=np.array([4,7]),mass=0.1,ray=self.radius)
        disk2 = disk.Disk(X=np.array([disk1.X[0],disk1.X[1]-self.rest_spring_l]),mass=0.1,ray=self.radius)
        disk3 = disk.Disk(X=np.array([disk2.X[0],disk2.X[1]-self.rest_spring_l]),mass=0.1,ray=self.radius)
        disk4 = disk.Disk(X=np.array([disk3.X[0]-self.rest_spring_l,disk3.X[1]]),mass=0.1,ray=self.radius)
        disk5 = disk.Disk(X=np.array([disk4.X[0]-self.rest_spring_l,disk4.X[1]]),mass=0.1,ray=self.radius)
        self.bacteria.append (bact.Bacterium(N=5,Disks=[disk1,disk2,disk3,disk4,disk5],l=self.rest_spring_l,theta=self.theta,color=(0,128,0)))
        
        # disk1 = disk.Disk(X=np.array([4,7]),mass=0.1,ray=self.radius)
        # self.bacteria.append (bact.Bacterium(N=1,Disks=[disk1],l=self.rest_spring_l,theta=self.theta,color=(0,128,0)))
        
        self.N_bacteria()

    def generate_random_bacteria(self,N=1):
        """Generate a random bacteria of N disks"""

        disks = []
        d = 4*self.radius #disks distance

        #Generation of the head
        x1 = random.uniform(5,15)
        x2 = random. uniform(5,15)
        X = np.array([x1,x2])
        disks.append(disk.Disk(X,ray=self.radius))

        ok = True
        i=1

        #Generation of the other disks
        while i <N:
            # Position generation
            ok= True
            alpha = random.uniform(0,2*np.pi)
            # if alpha >np.pi/2:
            #     alpha += np.pi
            x1 = disks[i-1].X[0] + d*np.cos(alpha)
            x2 = disks[i-1].X[1] + d*np.sin(alpha)

            # Overlapping checking
            for j in range(0,i):
                if(norm([x1-disks[j].X[0],x2-disks[j].X[1]])<=2*self.radius):
                    ok = False

            # Ok position
            if(ok==True):
                disks.append(disk.Disk(np.array([x1,x2]),ray=self.radius))
                i+=1
        
        # Generation of a color
        color= (random.randint(0,255),random.randint(0,255),random.randint(0,255))

        # Creation of the bacterium 
        new_bacteria = bact.Bacterium(len(disks),disks,l=self.radius,t_i=self.time,color=color)

        #Adding it into the list of badcteria
        self.bacteria.append(new_bacteria)
        self.N_bacteria()

    ### ----------------- Mechanical processes ------------------

    def bacteria_processes(self):
        """ Apply all the processes of the bacteria """

        for bact in self.bacteria:
            # Calculation of the velocity
            bact.spring_velocity()

            # Numerical intergration
            bact.Euler_explicit(self.dt)

            # Bacterium growth
            bact.growth(self.time,self.disk_add_method)
    
    ### ----------------- Simulation informations updaters ---------------
    def N_bacteria(self):
        """Update the number of bacteria"""
        self.N = len(self.bacteria)

    def increment_time(self):
        """Increment the total time of simulation by dt"""
        self.time +=self.dt