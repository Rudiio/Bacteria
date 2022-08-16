import time

# Classes
import disk
import bacterium as bact

# Calculations
import numpy as np
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
       
        # Stiffnesses with micrometer
        self.ks = 1       # Springs linear stiffness
        self.kt_par = 10   # Springs torsion parallel stiffness
        self.kt_bot = 10   # Springs torsion bot stiffness
        self.kc = 10       # Collision stiffness

        self.stiffness = (self.ks,self.kt_par,self.kt_bot,self.kc)

        # Bacteria parameters and variables
        self.bacteria = []
        self.N= 0
        
        # disks parameters in micrometer
        self.radius = 0.7 # 0.45  #2.15 #0.5  

        # Initial length in micrometer
        self.l_ini = 4.4  # 2.41

        # Growth parameters
        self.gk = 0.02843 #0.0168 #0.025 # 0.001    # Growth rate
        self.max_length = 15     # division length (useless)
        self.max_disks = 20 #1/(self.k*self.dt)

        # Springs parameters
        self.rest_spring_l = self.radius
        self.theta = np.pi

        # friction parameter
        self.mu = 0.2

        # Simulation parameters
        self.disk_add_method = 4 # decide the position of the new disk
        dtt = np.array([self.radius**2/(self.ks),self.radius/(self.kt_par),self.radius/(self.kt_bot),4*self.radius**2/(self.kc)])
        self.dt = self.mu*self.l_ini*dtt.min()*0.1 # 0.01  # in min
        if(self.ks<10):
            self.dt/=2
        self.time = 0.0
        self.tmax = 200

        # Mesh parameters
        # Size of the mesh
        self.Nx = 60
        self.Ny = 60
        self.dx = 2*self.radius

        # Position of the down left corner
        self.xmin = -self.dx*self.Nx/2 
        self.ymin = -self.dx*self.Ny/2 

        # Generating the first bacterium
        self.generate_bacterium()
        
    ### ---------------- Bacteria generators ---------------------------

    def generate_bacterium(self):
        l = 0
        disk1 = disk.Disk(X=np.array([-self.l_ini/2,0]),mass=1,ray=self.radius)
        D = [disk1]
        l+=2*self.radius
        i=0
        while l<self.l_ini:
            new_disk = disk.Disk(X=np.array([D[i].X[0]+self.rest_spring_l,D[i].X[1]]),mass=1,ray=self.radius)
            D.append(new_disk)
            l+=self.rest_spring_l
            i+=1
        self.bacteria.append(bact.Bacterium(N=len(D),Disks=D,l=self.radius,t_i = self.time,gm=self.disk_add_method,theta=self.theta,stiffness=self.stiffness,growth_k=self.gk,color=(0,128,0)))
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
        x1 = 5
        x2 = 5
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

    def daughter_cell_evolution(self):
        """" """
        if(self.N==4):
            self.bacteria[0].color=(139,0,0)
            self.bacteria[1].color= (0,0,128)
            self.bacteria[2].color=(255,140,0)
            self.bacteria[3].color=(34,139,34)
            
    ### ----------------- Collision detection -------------------

    def reset_bacteria_pointers(self):
        """Reset all the pointers from all disks of the simulation"""

        for bact in self.bacteria:
            bact.reset_pointers()

    def detect_collision(self,X : np.array,bact=[]):
        """ Detect if the position X is in collision with another bacterium
        return True if it detects a collision"""

        for b in bact:
            for d in b.Disks:
                if( norm(X - d.X ) <= 2*d.radius):
                    return True
        return False

    def space_mesh(self):
        """ Constructs an abstractive mesh of the space.
        Returns of an array containing the first cell of each case.
        Links all the corresponding disks"""

        # Resentting all the pointers
        self.reset_bacteria_pointers()

        # Constructing the array 
        first_cell = -1*np.ones((2,self.Nx*self.Ny))

        for i in range(self.N):
            bacte : bact.Bacterium = self.bacteria[i] 
            for j in range(bacte.p_i):
                d : disk.Disk = bacte.Disks[j]
                
                # Calculating the space index
                k1  = np.floor((d.X[0]-self.xmin)/self.dx)
                k2  = np.floor((d.X[1]-self.ymin)/self.dx)
                l = int(k1 +self.Nx*k2)
                
                # Adding the cell into the first_cell array if case
                if first_cell[0,l]==-1 and first_cell[1,l]==-1:
                    first_cell[0,l]=i
                    first_cell[1,l]=j
                
                # Linking the corresponding disks
                else:
                    ni = int(first_cell[0,l])
                    nj = int(first_cell[1,l])

                    while(self.bacteria[ni].Disks[nj].next_bact!=-1 and self.bacteria[ni].Disks[nj].next_disk!=-1 ):
                        temp_i =ni
                        temp_j =nj
                        ni = self.bacteria[temp_i].Disks[temp_j].next_bact
                        nj = self.bacteria[temp_i].Disks[temp_j].next_disk
                    
                    self.bacteria[ni].Disks[nj].next_bact = i
                    self.bacteria[ni].Disks[nj].next_disk = j
                
        return first_cell

    ### ---------------- Division -------------------------------
    def division(self,bacte: bact.Bacterium):
        """ Divide the bacterium if the maximum number of disks is reached.
        Create two daughters bacteria"""

        # Checking if the bacterium should divide
        if(bacte.L >= bacte.max_length):
            # Extracting the lists
            L1,L2 = bacte.division()
            
            # Noises
            self.division_noise(L1,L2)

            # Make the bacteria collide
            if(self.disk_add_method==7):
                L1.reverse()

            # Oppose the direction of the daughters
            # L2.reverse()

            # Creation of the daughters
            color= (random.randint(0,255),random.randint(0,255),random.randint(0,255))
            D1 = bact.Bacterium(N=len(L1),Disks=L1,l=self.rest_spring_l,t_i=self.time,gm=self.disk_add_method,
                growth_k=self.gk,stiffness=self.stiffness,gen=bacte.gen+1,color=color)
            color= (random.randint(0,255),random.randint(0,255),random.randint(0,255))
            D2 = bact.Bacterium(N=len(L2),Disks=L2,l=self.rest_spring_l,t_i=self.time,gm=self.disk_add_method,
                growth_k=self.gk,stiffness=self.stiffness,gen=bacte.gen+1,color=color)
            
            # Deleting the mother
            self.bacteria.remove(bacte)
            self.N -=1

            #Adding the daughter
            self.bacteria.append(D1)
            self.N +=1
            self.bacteria.append(D2)
            self.N +=1

    def division_noise(self,L1,L2):
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
            rc1 = np.array([(L1[l1//2].X[0]+L1[l1//2-1].X[0])/2,(L1[l1//2].X[1]+L1[l1//2-1].X[1])/2])
        else:
            rc1 = L1[l1//2].X
        
        if l2%2==0:
            rc2 = np.array([(L2[l2//2].X[0]+L2[l2//2-1].X[0])/2,(L2[l2//2].X[1]+L2[l2//2-1].X[1])/2])
        else:
            rc2 = L2[l2//2].X
        
        # Rotations
        for j in range(0,max(l1,l2)):
            if j<l1 or (l1%2==0 and j!=l1//2 and j<l1):
                L1[j].X = np.dot(M1,L1[j].X  - rc1 )+ rc1
            if j<l2 or (l2%2==0 and j!=l2//2 and j<l2):
                L2[j].X = np.dot(M2,L2[j].X  - rc2 )+ rc2
                
    ### ----------------- All Processes ------------------
    
    def interaction_calculation(self):
        """ Calculate the intern and extern interactions of the bacteria wich modify the velocity
        of  the disks"""

        # Construction of the space mesh
        first_cell = self.space_mesh()

        for i in range(0,self.N):
            bact = self.bacteria[i]

            # Calculation of the velocity
            mesh_param = (self.xmin,self.ymin,self.dx,self.Nx)
            bact.spring_velocity(i,self.bacteria,first_cell,mesh_param,self.mu)

    def bacteria_processes(self):
        """ Apply all the processes of the bacteria """
        self.interaction_calculation()

        for bact in self.bacteria:
            # Numerical integration
            bact.Euler_explicit(self.dt)

            # Bacterium growth
            bact.growth(self.time,bact.growth_method,self.max_length)
            bact.update_length()
            
            #Bacterium division
            self.division(bact)
            
    ### ----------------- Simulation informations updaters ---------------
    def N_bacteria(self):
        """Update the number of bacteria"""
        self.N = len(self.bacteria)

    def increment_time(self):
        """Increment the total time of simulation by dt"""
        self.time +=self.dt
    
    def write_columns(self,s):
        file = open(s,"a")
        file.write("N disks\t")
        file.write("Disks radius\t")
        file.write("time\t")
        file.write("Gen\t")
        file.write("ks\t")
        file.write("kt_bot\t")
        file.write("kt_par\t")
        file.write("kc\t")
        file.write("color\t")
        file.write("current length\t")
        file.write("div_length\t")
        file.write("X")
        file.write("\n")
        file.close()

    def write_txt(self,s):
        """" Writes the data into a txt file"""
        file = open(s,"a")
        for bact in self.bacteria:
            file.write(f"{bact.p_i}\t")
            file.write(f"{self.radius}\t")
            file.write(f"{self.time:.3f}\t")
            file.write(f"{bact.gen}\t")
            file.write(f"{self.ks}\t")
            file.write(f"{self.kt_bot}\t")
            file.write(f"{self.kt_par}\t")
            file.write(f"{self.kc}\t")
            file.write(f"{bact.color}\t")
            file.write(f"{bact.L}\t")
            file.write(f"{bact.max_length}\t")
            for j in range(bact.p_i):
                file.write(f"{bact.Disks[j].X[0]} {bact.Disks[j].X[1]} ")
            file.write("\n")
        file.close()
    
    def mainloop(self):
        """Main loop of the application/simulation"""

        s = r"./simulations/simuc11.txt"
        self.write_columns(s)
        start = time.time()
        last_int = 0
        while self.time <= self.tmax :
            if(int(self.time)%3==0 and int(self.time)!=last_int):
                last_int = int(self.time)
                self.write_txt(s)

            #Move the bacteria
            self.bacteria_processes()

            # Increasing the time
            self.increment_time()
            print("time=",self.time)


        end = time.time()
        print(end-start)
        
        #Writing the final data
        self.write_txt(s)


if __name__=="__main__":
    simu = Model()
    simu.mainloop()
