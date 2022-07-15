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
        
        # Simulation parameters
        self.disk_add_method =  2 # decide the position of the new disk
        self.dt = 0.5
        self.time = 0

        # Bacteria parameters and variables
        self.bacteria = []
        self.N= 0
        
        # disks parameters
        self.radius = 0.5

        # Growth parameters
        self.k = 0.01 #0.025 # 0.001    # Growth constant
        self.max_length = 10
        self.max_disks = 20 #1/(self.k*self.dt)

        # print(self.max_disks)

        # Springs parameters
        self.rest_spring_l = self.radius
        self.theta = np.pi

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
        self.bacteria = [bact.Bacterium(N=1,Disks=[disk1],l=self.rest_spring_l,theta=self.theta,growth_k=self.k,color=(139,0,0))]
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
        self.bacteria.append (bact.Bacterium(N=5,Disks=[disk1,disk2,disk3,disk4,disk5],l=self.rest_spring_l,theta=self.theta,growth_k=self.k,color=(0,128,0)))
        
        # disk1 = disk.Disk(X=np.array([4,7]),mass=0.1,ray=self.radius)
        # self.bacteria.append (bact.Bacterium(N=1,Disks=[disk1],l=self.rest_spring_l,theta=self.theta,color=(0,128,0)))
        
        self.N_bacteria()

    def gen_cell_pos(self,x,y,method):
        """Generates a cell bacterium"""
        cell = disk.Disk(np.array([x,y]),ray=self.radius)
        color= (random.randint(0,255),random.randint(0,255),random.randint(0,255))
        self.bacteria.append (bact.Bacterium(N=1,Disks=[cell],l=self.rest_spring_l,theta=self.theta,gm=method,growth_k=self.k,color=color))
        self.N +=1

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
        new_bacteria = bact.Bacterium(len(disks),disks,l=self.radius,t_i=self.time,gm=self.disk_add_method,growth_k=self.k,color=color)

        #Adding it into the list of badcteria
        self.bacteria.append(new_bacteria)
        self.N_bacteria()

    def generate_random_bacterium_no_collision(self,N=1):
        """Generate a random bacteria of N disks without collision"""

        disks = []
        d = 4*self.radius #disks distance

        #Generation of the head
        x1 = random.uniform(5,10)
        x2 = random. uniform(5,10)
        X = np.array([x1,x2])

        while(self.detect_collision(X,self.bacteria)):
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
            if(ok==True and self.detect_collision(np.array([x1,x2]),self.bacteria)==False):
                disks.append(disk.Disk(np.array([x1,x2]),ray=self.radius))
                i+=1
        
        # Generation of a color
        color= (random.randint(0,255),random.randint(0,255),random.randint(0,255))

        # Creation of the bacterium 
        new_bacteria = bact.Bacterium(len(disks),disks,l=self.radius,t_i=self.time,gm=self.disk_add_method,growth_k=self.k,color=color)

        #Adding it into the list of badcteria
        self.bacteria.append(new_bacteria)
        self.N_bacteria()
    
    def generate_random_bacterium_no_collision_method(self,N=1,method=1):
        """Generate a random bacteria of N disks without collision"""

        disks = []
        d = 4*self.radius #disks distance

        #Generation of the head
        x1 = random.uniform(5,25)
        x2 = random. uniform(5,25)
        X = np.array([x1,x2])

        while(self.detect_collision(X,self.bacteria)):
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
            if(ok==True and self.detect_collision(np.array([x1,x2]),self.bacteria)==False):
                disks.append(disk.Disk(np.array([x1,x2]),ray=self.radius))
                i+=1
        
        # Generation of a color
        color= (random.randint(0,255),random.randint(0,255),random.randint(0,255))

        # Creation of the bacterium 
        new_bacteria = bact.Bacterium(len(disks),disks,l=self.radius,t_i=self.time,gm=method,growth_k=self.k,color=color)

        #Adding it into the list of badcteria
        self.bacteria.append(new_bacteria)
        self.N_bacteria()

    ### ----------------- Collision detection -------------------

    def detect_collision(self,X : np.array,bact: list[bact.Bacterium] =[]):
        """ Detect if the position X is in collision with another bacterium
        return True if it detects a collision"""

        for b in bact:
            for d in b.Disks:
                if( norm(X - d.X ) <= 2*d.radius):
                    return True
        return False
                
    ### ---------------- Division -------------------------------
    def division(self,i,bacte: bact.Bacterium):
        """ Divide the bacterium if the maximum number of disks is reached.
        Create two daughters bacteria"""

        # Checking if the bacterium should divide
        if(bacte.p_i >= self.max_disks):
            # Extracting the lists
            L1,L2 = bacte.division()
            
            # Noises
            self.division_noise(L1,L2)

            # Oppose the direction of the daughters
            # L2.reverse()

            # Creation of the daughters
            color= (random.randint(0,255),random.randint(0,255),random.randint(0,255))
            D1 = bact.Bacterium(N=len(L1),Disks=L1,l=self.rest_spring_l,t_i=self.time,gm=self.disk_add_method,
                growth_k=self.k,color=color)
            color= (random.randint(0,255),random.randint(0,255),random.randint(0,255))
            D2 = bact.Bacterium(N=len(L2),Disks=L2,l=self.rest_spring_l,t_i=self.time,gm=self.disk_add_method,
                growth_k=self.k,color=color)
            
            # Deleting the mother
            self.bacteria.remove(bacte)
            self.N -=1

            #Adding the daughter
            self.bacteria.append(D1)
            self.N +=1
            self.bacteria.append(D2)
            self.N +=1

    def division_noise(self,L1:list[disk.Disk],L2 : list[disk.Disk]):
        """ Add angular noise to the lists of disks after division"""

        # Noises
        theta = 1.e-1 #np.pi/2 #1.e-1
        dtheta1 = random.uniform(-theta/2,theta/2)
        dtheta2 = -dtheta1

        # Length of the lists
        l1= len(L1)
        l2= len(L2)

        # Rotation matrix
        M1 = np.array(([np.cos(dtheta1),-np.sin(dtheta1)],[np.sin(dtheta1),np.cos(dtheta1)]))
        M2 = np.array(([np.cos(dtheta2),-np.sin(dtheta2)],[np.sin(dtheta2),np.cos(dtheta2)]))

        # Rotation centers 
        if l1%2==0:
            rc1 = np.array([(L1[l1//2].X[0]+L1[l1//2+1].X[0])/2,(L1[l1//2].X[1]+L1[l1//2+1].X[1])/2])
        else:
            rc1 = L1[l1//2].X
        
        if l2%2==0:
            rc2 = np.array([(L2[l2//2].X[0]+L2[l2//2+1].X[0])/2,(L2[l2//2].X[1]+L2[l2//2+1].X[1])/2])
        else:
            rc2 = L2[l2//2].X
        
        # Rotations
        for j in range(0,max(l1,l2)):
            if j<l1 or (l1%2==0 and j!=l1//2 and j<l1):
                L1[j].X = np.dot(M1,L1[j].X  - rc1 )+ rc1
            if j<l2 or (l2%2==0 and j!=l2//2 and j<l2):
                L2[j].X = np.dot(M2,L2[j].X  - rc2 )+ rc2
                
    ### ----------------- All Processes ------------------

    def bacteria_processes(self):
        """ Apply all the processes of the bacteria """

        for i in range(0,self.N):
            bact = self.bacteria[i]

            # Calculation of the velocity
            bact.spring_velocity(i,self.bacteria)

        for bact in self.bacteria:
            # Numerical integration
            bact.Euler_explicit(self.dt)

            # Bacterium growth
            # bact.growth(self.time,self.disk_add_method)
            bact.growth(self.time,bact.growth_method,self.max_disks)

            #Bacterium division
            self.division(i,bact)
            
    ### ----------------- Simulation informations updaters ---------------
    def N_bacteria(self):
        """Update the number of bacteria"""
        self.N = len(self.bacteria)

    def increment_time(self):
        """Increment the total time of simulation by dt"""
        self.time +=self.dt