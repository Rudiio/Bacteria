from turtle import distance
import numpy as np
import disk
import bacterium as bact
import random

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
        self.N= 0
        self.ray = 0.5
        self.rest_spring_l = self.ray
        self.theta = np.pi
        self.dt = 0.5
        self.time = 0
    
    def generate_bacterium(self):
        disk1 = disk.Disk(X=np.array([1,5]),mass=0.1,ray=self.ray)
        self.bacteria = [bact.Bacterium(N=1,Disks=[disk1],l=self.rest_spring_l,theta=self.theta,color=(139,0,0))]
        self.N+=1
        disk2 = disk.Disk(X=np.array([disk1.X[0]+self.bacteria[0].spring_rest_l+10,disk1.X[1]]),mass=0.1,ray=self.ray)
        self.bacteria[0].Disks.append(disk2)
        self.bacteria[0].p_i+=1
        disk3 = disk.Disk(X=np.array([disk2.X[0]+self.bacteria[0].spring_rest_l,disk2.X[1]+5]),mass=0.1,ray=self.ray)
        self.bacteria[0].Disks.append(disk3)
        self.bacteria[0].p_i+=1

        # disk1 = disk.Disk(X=np.array([1,5]),mass=0.1,ray=self.ray)
        # self.bacteria = [bact.Bacterium(N=1,Disks=[disk1],l=self.rest_spring_l,theta=self.theta,color=(139,0,0))]
        # self.N+=1
        # disk2 = disk.Disk(X=np.array([disk1.X[0]+self.bacteria[0].spring_rest_l+10,disk1.X[1]]),mass=0.1,ray=self.ray)
        # self.bacteria[0].Disks.append(disk2)
        # self.bacteria[0].p_i+=1
        # disk3 = disk.Disk(X=np.array([disk2.X[0]+self.bacteria[0].spring_rest_l,disk2.X[1]+5]),mass=0.1,ray=self.ray)
        # self.bacteria[0].Disks.append(disk3)
        # self.bacteria[0].p_i+=1

        # disk1 = disk.Disk(X=np.array([4,7]),mass=0.1,ray=self.ray)
        # disk2 = disk.Disk(X=np.array([disk1.X[0],disk1.X[1]-self.rest_spring_l]),mass=0.1,ray=self.ray)
        # disk3 = disk.Disk(X=np.array([disk2.X[0],disk2.X[1]-self.rest_spring_l]),mass=0.1,ray=self.ray)
        # disk4 = disk.Disk(X=np.array([disk3.X[0]-self.rest_spring_l,disk3.X[1]]),mass=0.1,ray=self.ray)
        # disk5 = disk.Disk(X=np.array([disk4.X[0]-self.rest_spring_l,disk4.X[1]]),mass=0.1,ray=self.ray)
        # self.bacteria.append (bact.Bacterium(N=5,Disks=[disk1,disk2,disk3,disk4,disk5],l=self.rest_spring_l,theta=self.theta,color=(0,128,0)))
        
        self.N_bacteria()

    def generate_random_bacteria(self,N=1):
        """Generate a random bacteria of N disks"""
        disks = []
        d = 4*self.ray #disks distance
        #Generation of the first disk
        x1 = random.uniform(5,10)
        x2 = random. uniform(5,10)
        X = np.array([x1,x2])
        disks.append(disk.Disk(X,ray=self.ray))

        #Generation of the other disks
        for i in range(1,N):
            alpha = random.uniform(0,2*np.pi)
            x1 = disks[i-1].X[0] + d*np.cos(alpha)
            x2 = disks[i-1].X[1] + d*np.sin(alpha)
            disks.append(disk.Disk(np.array([x1,x2]),ray=self.ray))
        
        color= (random.randint(0,255),random.randint(0,255),random.randint(0,255))
        new_bacteria = bact.Bacterium(len(disks),disks,color=color)
        self.bacteria.append(new_bacteria)
        self.N_bacteria()


    def N_bacteria(self):
        """Update the number of bacteria"""
        self.N = len(self.bacteria)

    def Move_bacteria(self):
        """Calculates the velocities of each disks of all the bacteria
        and Integrates them to determine the positions"""

        for bact in self.bacteria:
            bact.spring_velocity()
            bact.Euler_explicit(self.dt)
    
    def increment_time(self):
        """Increment the total time of simulation by dt"""
        self.time +=self.dt