# class
import scipy
import disk 

# Calculations
import numpy as np
from numpy.linalg import norm
import random
from scipy import stats
import joblib

def g(x,mu,sigma):
    """ Density of a Normal law"""
    return 1/(sigma*np.sqrt(2*np.pi))*np.exp(-0.5*((x-mu)/sigma)**2)

def increment_reject():
    """ Simulates the law of the increment distribution
    Only the increment of E.coli can be simulated here.
    A new hull is needed for the pseudomonas increment. """

    # Density to fit
    data = np.load("./increment.npy")
    f = stats.gaussian_kde(data)

    # Parameters of our gaussian hull
    mu=3.23
    sigma=1.4
    c=2
    t = False
    
    while not t:
        x = random.gauss(mu,sigma)
        u = random.random()
        
        if(c*g(x,mu,sigma)*u < f(x)):
            t = True
            return x 

class Bacterium:
    """ Class representing a bacterium 
    It composed of :
    - the number of Disks in the Bacterium: N
    - a list of Disks: Disks
    - the rest torque theta
    - the rest length of the springs l 
    - the stiffness constants ks,kt1, kt2"""

    def __init__(self,N=0,Disks = [],l = 1.0,theta = np.pi,t_i = 0,gm=1,growth_k = 0.01,gen=1,stiffness=(1,1,1,1), color = (0,0,0)):
        """Constructor for the class Bacterium"""

        # Disks variables
        self.p_i = N        # number of disks
        self.Disks = Disks  # List of disks
        self.color = color  # disks colors

        # Generation 
        self.gen = gen

        # bacteria length
        self.L = 1
        self.update_length()
        self.max_length = self.L + increment_reject()   # length of  divison in micrometer
        # print(self.max_length)

        # Springs parameters
        self.theta = theta  # Rest torque
        self.spring_rest_l = l  # Rest length

        self.ks = stiffness[0] # Springs linear stiffness
        self.kt_par = stiffness[1] # Springs torsion parallel stiffness
        self.kt_bot = stiffness[2] # Springs torsion bot stiffness

        # Growth variable
        self.gk = growth_k    # Growth constant
        self.t_i = t_i    # time of the last addition
        self.max_disks = 20 # Number maximal of disk contained by a bacterium

        # softening parameter (to avoid division by zero)
        self.eps = 0.0001 

        # Growth method
        self.growth_method = gm
        self.side = -1
        
        # Collision constant
        self.kc = stiffness[3]  # Collision stiffness

    def __str__(self):
        """Display the bacterium informations whit the print function"""
        print(f"N = {self.p_i}")
        for i in self.Disks :
            print(i)
        return("")
    
    def get_segment_length(self):
        """ return an array with the length of all the segments """
        if(self.p_i>1):
            L = np.zeros(self.p_i -1 )
            for i in range(0,self.p_i-1):
                L[i] = norm(self.Disks[i].X - self.Disks[i+1].X)
            return L
        else:
            return np.zeros(0)

    def update_length(self):
        """Returns the length of a bacterium"""
        self.L = self.get_segment_length().sum() + 2*self.Disks[0].radius

    def points(self):
        """Returns an array of the points of the bacterium"""
        P =[]
        for p in self.Disks:
            P.append([p.X[0],p.X[1]])
        return P

    def reset_pointers(self):
        """Reset the pointers of the disks of the bacteria"""

        for d in self.Disks:
            d.next_bact = -1
            d.next_disk = -1

    ###-----------------  Velocity calculation -----------------------

    def spring_velocity(self,ci,bacteria,first_cell,mesh_param,mu):
        """Calculate the velocity that comes from the spring forces/torques
        of each cell of the bacterium"""

        # Loop on all the Disks
        for k in range(self.p_i):
            if(self.p_i >1):
                self.Disks[k].V = self.linear_spring(k) + self.torsion_spring_par(k) +self.torsion_spring_bot(k)
            # print("no opti",self.repulsion(ci,bacteria,k))
            self.repulsion_opti(ci,bacteria,k,first_cell,mesh_param)
            self.Disks[k].V *=1/(mu*self.L)

    def linear_spring(self,k):
        """Calculates the velocity created by the linear springs"""

        # case : the disk is the head of the bacterium
        if k==0:
            Xj = self.Disks[k].X
            Xj1 = self.Disks[k+1].X

            return (self.ks/(self.spring_rest_l**2)*(norm(Xj1 - Xj) - self.spring_rest_l)*
                    (Xj1 - Xj)/(norm(Xj1 - Xj))) 
        
        # case the disk is the tail of the bacterium
        elif k==self.p_i-1:
            Xj = self.Disks[k].X
            Xj_1 = self.Disks[k-1].X

            return (-self.ks/(self.spring_rest_l**2)*(norm(Xj - Xj_1) - self.spring_rest_l)*
                    (Xj - Xj_1)/norm(Xj-Xj_1))
        
        # all the others disks 
        else:
            Xj = self.Disks[k].X
            Xj1 = self.Disks[k+1].X
            Xj_1 = self.Disks[k-1].X

            return (self.ks/(self.spring_rest_l**2)*(norm(Xj1 - Xj) - self.spring_rest_l)*
                    (Xj1 - Xj)/(norm(Xj1 - Xj)) -self.ks/(self.spring_rest_l**2)*(norm(Xj - Xj_1) - self.spring_rest_l)*
                    (Xj - Xj_1)/norm(Xj-Xj_1))

    def torsion_spring_par(self,k):
        """Calculates the velocity created by the torsion springs for the parallel component"""
        V = 0

        # case : the disk is the head of the bacterium
        if k==0 and self.p_i>=3 :  
            Xj = self.Disks[k].X
            Xj1 = self.Disks[k+1].X
            Xj2 = self.Disks[k+2].X

            # term in j,j+1,j+2 with cos
            V -= self.kt_par/(norm(Xj2 - Xj1)*norm(Xj - Xj1) +self.eps)*((np.dot((Xj2 - Xj1),(Xj - Xj1)))/(norm(Xj2 - Xj1)*norm(Xj - Xj1)+self.eps) - np.cos(self.theta))*\
                ((Xj2 - Xj1)-np.dot((Xj2 - Xj1),(Xj - Xj1))*(Xj - Xj1)/(norm(Xj-Xj1)**2+self.eps))
        
        # case : the disk is the second of the bacterium
        elif k==1 :
            Xj_1 = self.Disks[k-1].X
            Xj = self.Disks[k].X

            # term in j-1,j,j+1 with cos
            if(self.p_i >=3):
                Xj1 = self.Disks[k+1].X
                V += self.kt_par/(norm(Xj1 - Xj)*norm(Xj_1 - Xj)+self.eps)*((np.dot((Xj1 - Xj),(Xj_1 - Xj)))/(norm(Xj1 - Xj)*norm(Xj_1 - Xj)+self.eps) - np.cos(self.theta))*\
                    ((Xj_1 - Xj)-np.dot((Xj1 - Xj),(Xj_1 - Xj))*(Xj1 - Xj)/(norm(Xj1 -Xj)**2+self.eps)+\
                        (Xj1 - Xj)-np.dot((Xj1 - Xj),(Xj_1 - Xj))*(Xj_1-Xj)/(norm(Xj_1 - Xj)**2+self.eps))

            # term in j,j+1,j+2 with cos
            if(self.p_i >3):
                Xj2 = self.Disks[k+2].X
                V -= self.kt_par/(norm(Xj2 - Xj1)*norm(Xj - Xj1)+self.eps)*((np.dot((Xj2 - Xj1),(Xj - Xj1)))/(norm(Xj2 - Xj1)*norm(Xj - Xj1)+self.eps) - np.cos(self.theta))*\
                    ((Xj2 - Xj1)-np.dot((Xj2 - Xj1),(Xj - Xj1))*(Xj - Xj1)/(norm(Xj-Xj1)**2+self.eps))

        # case : the disk is the tail of the bacterium
        elif k==self.p_i-1 and self.p_i>=3:
            Xj_2 = self.Disks[k-2].X
            Xj_1 = self.Disks[k-1].X
            Xj = self.Disks[k].X

            # term in j-2,j-1,j with cos
            V -= self.kt_par/(norm(Xj -Xj_1)*norm(Xj_2 - Xj_1)+self.eps)*((np.dot((Xj-Xj_1),(Xj_2 - Xj_1)))/(norm(Xj -Xj_1)*norm(Xj_2 - Xj_1)+self.eps) - np.cos(self.theta))*\
                ((Xj_2 - Xj_1)-np.dot((Xj-Xj_1),(Xj_2 - Xj_1))*(Xj - Xj_1)/(norm(Xj - Xj_1)**2+self.eps))

        # case : the disk is just before the tail of the bacterium
        elif k==self.p_i-2:
            Xj = self.Disks[k].X

            # term in j-2,j-1,j with cos
            if(self.p_i>3):
                Xj_2 = self.Disks[k-2].X
                Xj_1 = self.Disks[k-1].X
                V -= self.kt_par/(norm(Xj -Xj_1)*norm(Xj_2 - Xj_1)+self.eps)*((np.dot((Xj-Xj_1),(Xj_2 - Xj_1)))/(norm(Xj -Xj_1)*norm(Xj_2 - Xj_1)+self.eps) - np.cos(self.theta))*\
                    ((Xj_2 - Xj_1)-np.dot((Xj-Xj_1),(Xj_2 - Xj_1))*(Xj - Xj_1)/(norm(Xj - Xj_1)**2+self.eps))
            
            # term in j-1,j,j+1 with cos
            if(self.p_i>=3):
                Xj_1 = self.Disks[k-1].X
                Xj1 = self.Disks[k+1].X
                V += self.kt_par/(norm(Xj1 - Xj)*norm(Xj_1 - Xj)+self.eps)*((np.dot((Xj1 - Xj),(Xj_1 - Xj)))/(norm(Xj1 - Xj)*norm(Xj_1 - Xj)+self.eps) - np.cos(self.theta))*\
                    ((Xj_1 - Xj)-np.dot((Xj1 - Xj),(Xj_1 - Xj))*(Xj1 - Xj)/(norm(Xj1 -Xj)**2+self.eps)+\
                        (Xj1 - Xj)-np.dot((Xj1 - Xj),(Xj_1 - Xj))*(Xj_1-Xj)/(norm(Xj_1 - Xj)**2+self.eps))
        
        # case : all the other disks
        else :
            Xj_2 = self.Disks[k-2].X
            Xj_1 = self.Disks[k-1].X
            Xj = self.Disks[k].X
            Xj1 = self.Disks[k+1].X
            Xj2 = self.Disks[k+2].X

            # term in j-2,j-1,j with cos
            V -= self.kt_par/(norm(Xj -Xj_1)*norm(Xj_2 - Xj_1)+self.eps)*((np.dot((Xj-Xj_1),(Xj_2 - Xj_1)))/(norm(Xj -Xj_1)*norm(Xj_2 - Xj_1)+self.eps) - np.cos(self.theta))*\
                ((Xj_2 - Xj_1)-np.dot((Xj-Xj_1),(Xj_2 - Xj_1))*(Xj - Xj_1)/(norm(Xj - Xj_1)**2+self.eps))
            
            # term in j-1,j,j+1 with cos
            V += self.kt_par/(norm(Xj1 - Xj)*norm(Xj_1 - Xj)+self.eps)*((np.dot((Xj1 - Xj),(Xj_1 - Xj)))/(norm(Xj1 - Xj)*norm(Xj_1 - Xj)+self.eps) - np.cos(self.theta))*\
                ((Xj_1 - Xj)-np.dot((Xj1 - Xj),(Xj_1 - Xj))*(Xj1 - Xj)/(norm(Xj1 -Xj)**2+self.eps)+\
                    (Xj1 - Xj)-np.dot((Xj1 - Xj),(Xj_1 - Xj))*(Xj_1-Xj)/(norm(Xj_1 - Xj)**2+self.eps))
            
            # term in j,j+1,j+2 with cos
            V -= self.kt_par/(norm(Xj2 - Xj1)*norm(Xj - Xj1)+self.eps)*((np.dot((Xj2 - Xj1),(Xj - Xj1)))/(norm(Xj2 - Xj1)*norm(Xj - Xj1)+self.eps) - np.cos(self.theta))*\
                ((Xj2 - Xj1)-np.dot((Xj2 - Xj1),(Xj - Xj1))*(Xj - Xj1)/(norm(Xj-Xj1)**2+self.eps))

        return V

    def torsion_spring_bot(self,k):
        """Calculates the velocity created by the torsion springs for the parallel component"""
        V = 0

        # case : the disk is the head of the bacterium
        if k==0 and self.p_i>=3:
            Xj = self.Disks[k].X
            Xj1 = self.Disks[k+1].X
            Xj2 = self.Disks[k+2].X

            # term in j,j+1,j+2 with sin
            V -= self.kt_bot/(norm(Xj2 - Xj1)*norm(Xj - Xj1) + self.eps)*((np.cross((Xj2 - Xj1),(Xj -Xj1)))/(norm(Xj2 - Xj1)*norm(Xj - Xj1)+ self.eps) - np.sin(self.theta))*\
                np.array([-(Xj2-Xj1)[1]+(((Xj-Xj1)[0])/(norm(Xj-Xj1)**2 + self.eps))*np.cross((Xj-Xj1),(Xj2-Xj1)),(Xj2 - Xj1)[0] + ((Xj - Xj1)[1])/(norm(Xj -Xj1)**2 + self.eps)*np.cross((Xj-Xj1),(Xj2-Xj1))])
        
        # case : the disk is the second of the bacterium
        elif k==1:
            Xj = self.Disks[k].X
            Xj_1 = self.Disks[k-1].X

            # term in j-1,j,j+1 with sin
            if(self.p_i>=3):
                Xj1 = self.Disks[k+1].X

                # term in j-1,j,j+1
                V += self.kt_bot/(norm(Xj1 - Xj)*norm(Xj_1 - Xj) + self.eps)*((np.cross((Xj1 - Xj),(Xj_1 - Xj)))/(norm(Xj1 - Xj)*norm(Xj_1 - Xj) + self.eps) - np.sin(self.theta))*\
                    (np.array([(Xj_1 - Xj)[1] - (((Xj1-Xj)[0])/(norm(Xj1 - Xj)**2 + self.eps))*np.cross((Xj1 - Xj),(Xj_1 -Xj)),-(Xj_1 - Xj)[0] - (((Xj1-Xj)[1])/(norm(Xj1 - Xj)**2 + self.eps))*np.cross((Xj1 - Xj),(Xj_1 -Xj))]) +\
                    np.array([-(Xj1 -Xj)[1] + (((Xj_1 -Xj)[0])/(norm(Xj_1 -Xj)**2 + self.eps))*np.cross((Xj_1 - Xj),(Xj1 - Xj)),(Xj1 -Xj)[0] + (((Xj_1 -Xj)[1])/(norm(Xj_1 -Xj)**2 + self.eps))*np.cross((Xj_1 - Xj),(Xj1 - Xj))]))
            
            if(self.p_i>3):
                Xj2 = self.Disks[k+2].X

                # term in j,j+1,j+2 with sin
                V -= self.kt_bot/(norm(Xj2 - Xj1)*norm(Xj - Xj1) + self.eps)*((np.cross((Xj2 - Xj1),(Xj -Xj1)))/(norm(Xj2 - Xj1)*norm(Xj - Xj1)+ self.eps) - np.sin(self.theta))*\
                np.array([-(Xj2-Xj1)[1]+(((Xj-Xj1)[0])/(norm(Xj-Xj1)**2 + self.eps))*np.cross((Xj-Xj1),(Xj2-Xj1)),(Xj2 - Xj1)[0] + ((Xj - Xj1)[1])/(norm(Xj -Xj1)**2+self.eps)*np.cross((Xj-Xj1),(Xj2-Xj1))])

        # case : the disk is the tail of the bacterium
        elif k==self.p_i - 1 and self.p_i>=3:
            Xj_2 = self.Disks[k-2].X
            Xj_1 = self.Disks[k-1].X
            Xj = self.Disks[k].X

            # term in j-2,j-1,j with sin
            V -= self.kt_bot/(norm(Xj - Xj_1)*norm(Xj_2 - Xj_1) + self.eps)*((np.cross((Xj - Xj_1),(Xj_2 - Xj_1)))/(norm(Xj - Xj_1)*norm(Xj_2 - Xj_1)+ self.eps) -np.sin(self.theta))*\
                np.array([(Xj_2 - Xj_1)[1] - ((Xj - Xj_1)[0])/(norm(Xj-Xj_1)**2 + self.eps)*np.cross((Xj - Xj_1),(Xj_2 -Xj_1)),-(Xj_2 - Xj_1)[0] - ((Xj - Xj_1)[1])/(norm(Xj-Xj_1)**2 + self.eps)*np.cross((Xj - Xj_1),(Xj_2 -Xj_1))])
        
        # case : the disk is just before the tail of the bacterium
        elif k==self.p_i -2 :
            Xj = self.Disks[k].X

            if self.p_i >3:
                Xj_2 = self.Disks[k-2].X
                Xj_1 = self.Disks[k-1].X

                # term in j-2,j-1,j with sin
                V -= self.kt_bot/(norm(Xj - Xj_1)*norm(Xj_2 - Xj_1) + self.eps)*((np.cross((Xj - Xj_1),(Xj_2 - Xj_1)))/(norm(Xj - Xj_1)*norm(Xj_2 - Xj_1)+ self.eps) -np.sin(self.theta))*\
                    np.array([(Xj_2 - Xj_1)[1] - ((Xj - Xj_1)[0])/(norm(Xj-Xj_1)**2 + self.eps)*np.cross((Xj - Xj_1),(Xj_2 -Xj_1)),-(Xj_2 - Xj_1)[0] - ((Xj - Xj_1)[1])/(norm(Xj-Xj_1)**2 + self.eps)*np.cross((Xj - Xj_1),(Xj_2 -Xj_1))])

            if self.p_i >=3:
                Xj1 = self.Disks[k+1].X
                Xj_1 = self.Disks[k-1].X

                # term in j-1,j,j+1
                V += self.kt_bot/(norm(Xj1 - Xj)*norm(Xj_1 - Xj) + self.eps)*((np.cross((Xj1 - Xj),(Xj_1 - Xj)))/(norm(Xj1 - Xj)*norm(Xj_1 - Xj) + self.eps) - np.sin(self.theta))*\
                    (np.array([(Xj_1 - Xj)[1] - (((Xj1-Xj)[0])/(norm(Xj1 - Xj)**2 + self.eps))*np.cross((Xj1 - Xj),(Xj_1 -Xj)),-(Xj_1 - Xj)[0] - (((Xj1-Xj)[1])/(norm(Xj1 - Xj)**2 + self.eps))*np.cross((Xj1 - Xj),(Xj_1 -Xj))]) +\
                    np.array([-(Xj1 -Xj)[1] + (((Xj_1 -Xj)[0])/(norm(Xj_1 -Xj)**2 + self.eps))*np.cross((Xj_1 - Xj),(Xj1 - Xj)),(Xj1 -Xj)[0] + (((Xj_1 -Xj)[1])/(norm(Xj_1 -Xj)**2 + self.eps))*np.cross((Xj_1 - Xj),(Xj1 - Xj))]))

        # case : all the other disks
        else :
            Xj_2 = self.Disks[k-2].X
            Xj_1 = self.Disks[k-1].X
            Xj = self.Disks[k].X
            Xj1 = self.Disks[k+1].X
            Xj2 = self.Disks[k+2].X

            # term in j-2,j-1,j with sin
            V -= self.kt_bot/(norm(Xj - Xj_1)*norm(Xj_2 - Xj_1) + self.eps)*((np.cross((Xj - Xj_1),(Xj_2 - Xj_1)))/(norm(Xj - Xj_1)*norm(Xj_2 - Xj_1)+ self.eps) -np.sin(self.theta))*\
                np.array([(Xj_2 - Xj_1)[1] - ((Xj - Xj_1)[0])/(norm(Xj-Xj_1)**2 + self.eps)*np.cross((Xj - Xj_1),(Xj_2 -Xj_1)),-(Xj_2 - Xj_1)[0] - ((Xj - Xj_1)[1])/(norm(Xj-Xj_1)**2 + self.eps)*np.cross((Xj - Xj_1),(Xj_2 -Xj_1))])
                    
            # term in j-1,j,j+1
            V += self.kt_bot/(norm(Xj1 - Xj)*norm(Xj_1 - Xj) + self.eps)*((np.cross((Xj1 - Xj),(Xj_1 - Xj)))/(norm(Xj1 - Xj)*norm(Xj_1 - Xj) + self.eps) - np.sin(self.theta))*\
                (np.array([(Xj_1 - Xj)[1] - (((Xj1-Xj)[0])/(norm(Xj1 - Xj)**2 + self.eps))*np.cross((Xj1 - Xj),(Xj_1 -Xj)),-(Xj_1 - Xj)[0] - (((Xj1-Xj)[1])/(norm(Xj1 - Xj)**2 + self.eps))*np.cross((Xj1 - Xj),(Xj_1 -Xj))]) +\
                np.array([-(Xj1 -Xj)[1] + (((Xj_1 -Xj)[0])/(norm(Xj_1 -Xj)**2 + self.eps))*np.cross((Xj_1 - Xj),(Xj1 - Xj)),(Xj1 -Xj)[0] + (((Xj_1 -Xj)[1])/(norm(Xj_1 -Xj)**2 + self.eps))*np.cross((Xj_1 - Xj),(Xj1 - Xj))]))

            # term in j,j+1,j+2 with sin
            V -= self.kt_bot/(norm(Xj2 - Xj1)*norm(Xj - Xj1) + self.eps)*((np.cross((Xj2 - Xj1),(Xj -Xj1)))/(norm(Xj2 - Xj1)*norm(Xj - Xj1)+ self.eps) - np.sin(self.theta))*\
                np.array([-(Xj2-Xj1)[1]+(((Xj-Xj1)[0])/(norm(Xj-Xj1)**2 + self.eps))*np.cross((Xj-Xj1),(Xj2-Xj1)),(Xj2 - Xj1)[0] + ((Xj - Xj1)[1])/(norm(Xj -Xj1)**2 + self.eps)*np.cross((Xj-Xj1),(Xj2-Xj1))])
        
        return V

    def repulsion(self,ci :int,bacteria, j):
        """Calculate the non-overlapping forces for all the disk of a bacterium
        toward the other disks of the bacteria of the simulation"""

        N = len(bacteria)
        v=0

        # Current disk
        Xj = self.Disks[j].X

        # Itering on the bacteria
        for i in range(0,N):
            cbact = bacteria[i]
            p_i = cbact.p_i

            if(i!=ci):
                for l in range(0,p_i):
                    Xl = cbact.Disks[l].X
                    
                    # checking the overlapping condition
                    if(norm(Xj - Xl) <= 2*self.Disks[0].radius):
                        v += self.kc/((2*self.Disks[0].radius)**2)*(1-(2*self.Disks[0].radius)/(norm(Xj-Xl)+self.eps))*(Xj-Xl)
                        self.Disks[j].V -= v

    def repulsion_2disks(self,cdisk :disk.Disk,odisk:disk.Disk):
        """" Calculate the repulsion force applied by a odisk on cdisk"""
        Xj = cdisk.X
        Xl = odisk.X

        if(norm(Xj - Xl) <= 2*self.Disks[0].radius):
            return -self.kc/((2*self.Disks[0].radius)**2)*(1-(2*self.Disks[0].radius)/(norm(Xj-Xl)+self.eps))*(Xj-Xl)

        return 0          

    def repulsion_opti(self,ci: int,bacteria,j,first_cell,mesh_param):
        """ Calculate the disks  inter bacteria interactions with the use of a spatial mesh"""

        N = len(bacteria)
        v=0

        # mesh param extraction
        xmin =mesh_param[0]
        ymin = mesh_param[1]
        dx = mesh_param[2]
        Nx = mesh_param[3]

        # Current disk
        Xj = self.Disks[j].X

        # Calculating the space index
        k1  = np.floor((Xj[0]-xmin)/dx)
        k2  = np.floor((Xj[1]-ymin)/dx)
        l = int(k1 +Nx*k2)
        
        # Itering on the cases around
        for k1 in [-1,0,1]:
            for k2 in [-1,0,1]:
                # cases number
                k = l + k1 + k2*Nx

                # head
                ni = int(first_cell[0,k])
                nj = int(first_cell[1,k])

                while(ni!=-1 and ni!=-1):
                    # force calculation
                    if(ci!=ni):
                        self.Disks[j].V += self.repulsion_2disks(self.Disks[j],bacteria[ni].Disks[nj]) 

                    temp_i = ni
                    temp_j = nj

                    # next disk
                    ni = bacteria[temp_i].Disks[temp_j].next_bact
                    nj = bacteria[temp_i].Disks[temp_j].next_disk

    ###------------------ Model bacterium processes -------------------------
    
    ## GROWTH
    
    def growth(self,t,method,max_length):
        """ Handle the growth process of thre bacterium 
        add a new disk if t_i > t. The position of the new bacteria depends on method"""

        rate_i = 1/(self.gk*self.p_i) 

        #updating the length of the bacteria
        self.update_length()

        if(method==2 or method==5):
            rate_i = 2/(self.gk*self.p_i) 

        if( t - self.t_i >= rate_i and self.L< self.max_length and self.p_i <self.max_disks):
            # Saving the moment we added a new disk
            self.t_i = rate_i + self.t_i
            
            # One side not equilibrum
            if method==1 :
                self.add_disk1()
            
            # Two sides equilibrum
            elif method==2:
                self.add_disk2()
            
            # At the center
            elif method==3:
                self.add_disk32()
            
            # One random side equilibrum
            elif method==4:
                self.add_disk4()
            
            # Two sides equilibrum
            elif method==5:
                self.add_disk5()
            
            # One fixed side equilibrum
            elif method==6:
                self.add_disk6()
                # One fixed side equilibrum
            elif method==7:
                self.add_disk7()

    def add_disk1(self):
        """Add a new disk is added not at equilibrum into the bacterium on one side of the bacterium
        chosen randomly"""
        
        side = random.random()

        # the disk is added in the head
        if(side > 0.5):
            l = random.uniform(0.01,self.spring_rest_l-0.01)

            # calculation of the angle
            if(self.p_i==1):
                alpha = random.uniform(0,2*np.pi)
            else :
                d_head = self.Disks[0]
                d_tail = self.Disks[1]
                vec1 = np.array([d_tail.X[0] - d_head.X[0],d_tail.X[1] - d_head.X[1]])
                vec2 = np.array([2,0])

                if(d_tail.X[1] - d_head.X[1]<=0):
                    alpha = 2*np.pi - np.arccos(np.dot(vec1,vec2)/(norm(vec1)*norm(vec2)+self.eps))
                else:
                    alpha = np.arccos(np.dot(vec1,vec2)/(norm(vec1)*norm(vec2)+self.eps))
                
            X = np.array([self.Disks[0].X[0] - l*np.cos(alpha),self.Disks[0].X[1] - l*np.sin(alpha)])
            self.Disks.insert(0,disk.Disk(X,ray=self.Disks[0].radius))
        
        # the disk is added in the tail
        else:
            l = random.uniform(0.01,self.spring_rest_l-0.01)

            # Calculation of the angle
            if(self.p_i==1):
                alpha = random.uniform(0,2*np.pi)
            else :
                d_head = self.Disks[self.p_i -2]
                d_tail = self.Disks[self.p_i -1]
                vec1 = np.array([d_tail.X[0] - d_head.X[0],d_tail.X[1] - d_head.X[1]])
                vec2 = np.array([2,0])

                if(d_tail.X[1] - d_head.X[1]<=0):
                    alpha = 2*np.pi - np.arccos(np.dot(vec1,vec2)/(norm(vec1)*norm(vec2)+self.eps))
                else:
                    alpha = np.arccos(np.dot(vec1,vec2)/(norm(vec1)*norm(vec2)+self.eps))
            
            X = np.array([self.Disks[self.p_i-1].X[0] + l*np.cos(alpha),self.Disks[self.p_i-1].X[1] + l*np.sin(alpha)])
            self.Disks.append(disk.Disk(X,ray=self.Disks[0].radius))
        
        # Updating the number of disks
        self.p_i +=1
    
    def add_disk2(self):
        """Add a new disk not at equilibrum into the bacterium on both side of the bacterium"""

        # the disk is added in the head
        l = random.uniform(0.01,self.spring_rest_l-0.01)
        alpha = 0

        #Calculation of the first angle
        if(self.p_i==1):
            alpha = random.uniform(0,2*np.pi)
        else :
            d_head = self.Disks[0]
            d_tail = self.Disks[1]
            vec1 = np.array([d_tail.X[0] - d_head.X[0],d_tail.X[1] - d_head.X[1]])
            vec2 = np.array([2,0])

            if(d_tail.X[1] - d_head.X[1]<=0):
                alpha = 2*np.pi - np.arccos(np.dot(vec1,vec2)/(norm(vec1)*norm(vec2)+self.eps))
            else:
                alpha = np.arccos(np.dot(vec1,vec2)/(norm(vec1)*norm(vec2)+self.eps))
            
        X = np.array([self.Disks[0].X[0] - l*np.cos(alpha),self.Disks[0].X[1] - l*np.sin(alpha)])
        self.Disks.insert(0,disk.Disk(X,ray=self.Disks[0].radius))
        
        # Updating the number of disks
        self.p_i +=1

        # Calculation of the second angle
        if(self.p_i-1!=1):
            d_head = self.Disks[self.p_i-2]
            d_tail = self.Disks[self.p_i-1]
            
            vec1 = np.array([d_tail.X[0] - d_head.X[0],d_tail.X[1] - d_head.X[1]])
            vec2 = np.array([2,0])

            if(d_tail.X[1] - d_head.X[1]<=0):
                alpha = 2*np.pi - np.arccos(np.dot(vec1,vec2)/(norm(vec1)*norm(vec2)+self.eps))
            else:
                alpha = np.arccos(np.dot(vec1,vec2)/(norm(vec1)*norm(vec2)+self.eps))

        # the disk is added in the tail
        Y = np.array([self.Disks[self.p_i-1].X[0] + l*np.cos(alpha),self.Disks[self.p_i-1].X[1] + l*np.sin(alpha)])
        self.Disks.append(disk.Disk(Y,ray=self.Disks[0].radius))
        
        # Updating the number of disks
        self.p_i +=1

    def add_disk3(self):
        """Add a new disk at the center of the bacterium"""

        if self.p_i == 1:
            l = random.uniform(0.01,self.spring_rest_l-0.01)
            alpha = random.uniform(0,2*np.pi)
            X = np.array([self.Disks[0].X[0] + l*np.cos(alpha),self.Disks[0].X[1] + l*np.sin(alpha)])
            self.Disks.insert(0,disk.Disk(X,ray=self.Disks[0].radius))
            self.p_i +=1
        else:
            if(self.p_i%2==0):
                mid = self.p_i//2 - 1
                X = np.array([(self.Disks[mid].X[0] + self.Disks[mid+1].X[0])/2,(self.Disks[mid].X[1] + self.Disks[mid+1].X[1])/2])
                self.Disks.insert(mid+1,disk.Disk(X,ray=self.Disks[0].radius))
            else:
                mid = int(self.p_i/2)
                X = np.array([(self.Disks[mid-1].X[0] + self.Disks[mid].X[0])/2,(self.Disks[mid-1].X[1] + self.Disks[mid].X[1])/2])
                self.Disks.insert(mid,disk.Disk(X,ray=self.Disks[0].radius))
            self.p_i +=1

    def add_disk32(self):
        """Add a new disk at the center of the bacterium"""

        if self.p_i == 1:
            l = random.uniform(0.01,self.spring_rest_l-0.01)
            alpha = random.uniform(0,2*np.pi)
            X = np.array([self.Disks[0].X[0] + l*np.cos(alpha),self.Disks[0].X[1] + l*np.sin(alpha)])
            self.Disks.insert(0,disk.Disk(X,ray=self.Disks[0].radius))
            self.p_i +=1
        else:
            if(self.p_i%2==0):
                mid = self.p_i//2 - 1
                X = np.array([(self.Disks[mid].X[0] + self.Disks[mid+1].X[0])/2,(self.Disks[mid].X[1] + self.Disks[mid+1].X[1])/2])
                self.Disks.insert(mid+1,disk.Disk(X,ray=self.Disks[0].radius))
            else:
                mid = self.p_i//2
                L = self.get_segment_length()
                S = L.sum()
                CS = np.cumsum(L)
                i=0
                # print(CS)
                while(CS[i]<S/2 and i<self.p_i):
                    i+=1
                
                X = np.array([(self.Disks[i].X[0] + self.Disks[i+1].X[0])/2,(self.Disks[i].X[1] + self.Disks[i+1].X[1])/2])
                self.Disks.insert(i+1,disk.Disk(X,ray=self.Disks[0].radius))
                    
            self.p_i +=1

    def add_disk42(self):
        """Add a new disk at equilibrum into the bacterium on one side of the bacterium
        chosen randomly"""
        side = random.random()
        alpha = 0

        # the disk is added in the head
        if(side < 0.5 ):
            l = self.spring_rest_l
            if(self.p_i==1):
                alpha = random.uniform(0,2*np.pi)
            else :
                d_head = self.Disks[0]
                d_tail = self.Disks[self.p_i -1]
                vec1 = np.array([d_tail.X[0] - d_head.X[0],d_tail.X[1] - d_head.X[1]])
                vec2 = np.array([2,0])

                if(d_tail.X[1] - d_head.X[1]<=0):
                    alpha = 2*np.pi - np.arccos(np.dot(vec1,vec2)/(norm(vec1)*norm(vec2)+self.eps))
                else:
                    alpha = np.arccos(np.dot(vec1,vec2)/(norm(vec1)*norm(vec2)+self.eps))

                
            X = np.array([self.Disks[0].X[0] - l*np.cos(alpha),self.Disks[0].X[1] - l*np.sin(alpha)])
            self.Disks.insert(0,disk.Disk(X,ray=self.Disks[0].radius))
        
        # the disk is added in the tail
        else:
            l = self.spring_rest_l
            if(self.p_i==1):
                alpha = random.uniform(0,2*np.pi)
            else :
                d_head = self.Disks[0]
                d_tail = self.Disks[self.p_i -1]
                vec1 = np.array([d_tail.X[0] - d_head.X[0],d_tail.X[1] - d_head.X[1]])
                vec2 = np.array([2,0])

                if(d_tail.X[1] - d_head.X[1]<=0):
                    alpha = 2*np.pi - np.arccos(np.dot(vec1,vec2)/(norm(vec1)*norm(vec2)+self.eps))
                else:
                    alpha = np.arccos(np.dot(vec1,vec2)/(norm(vec1)*norm(vec2)+self.eps))
            
            X = np.array([self.Disks[self.p_i-1].X[0] + l*np.cos(alpha),self.Disks[self.p_i-1].X[1] + l*np.sin(alpha)])
            self.Disks.append(disk.Disk(X,ray=self.Disks[0].radius))
        
        # Updating the number of disks
        self.p_i +=1
    
    def add_disk4(self):
        """Add a new disk is added not at equilibrum into the bacterium on one side of the bacterium
        chosen randomly"""
        
        side = random.random()

        # the disk is added in the head
        if(side > 0.5):
            l = self.spring_rest_l

            # calculation of the angle
            if(self.p_i==1):
                alpha = random.uniform(0,2*np.pi)
            else :
                d_head = self.Disks[0]
                d_tail = self.Disks[1]
                vec1 = np.array([d_tail.X[0] - d_head.X[0],d_tail.X[1] - d_head.X[1]])
                vec2 = np.array([2,0])

                if(d_tail.X[1] - d_head.X[1]<=0):
                    alpha = 2*np.pi - np.arccos(np.dot(vec1,vec2)/(norm(vec1)*norm(vec2)+self.eps))
                else:
                    alpha = np.arccos(np.dot(vec1,vec2)/(norm(vec1)*norm(vec2)+self.eps))
                
            X = np.array([self.Disks[0].X[0] - l*np.cos(alpha),self.Disks[0].X[1] - l*np.sin(alpha)])
            self.Disks.insert(0,disk.Disk(X,ray=self.Disks[0].radius))
        
        # the disk is added in the tail
        else:
            l = self.spring_rest_l

            # Calculation of the angle
            if(self.p_i==1):
                alpha = random.uniform(0,2*np.pi)
            else :
                d_head = self.Disks[self.p_i -2]
                d_tail = self.Disks[self.p_i -1]
                vec1 = np.array([d_tail.X[0] - d_head.X[0],d_tail.X[1] - d_head.X[1]])
                vec2 = np.array([2,0])

                if(d_tail.X[1] - d_head.X[1]<=0):
                    alpha = 2*np.pi - np.arccos(np.dot(vec1,vec2)/(norm(vec1)*norm(vec2)+self.eps))
                else:
                    alpha = np.arccos(np.dot(vec1,vec2)/(norm(vec1)*norm(vec2)+self.eps))
            
            X = np.array([self.Disks[self.p_i-1].X[0] + l*np.cos(alpha),self.Disks[self.p_i-1].X[1] + l*np.sin(alpha)])
            self.Disks.append(disk.Disk(X,ray=self.Disks[0].radius))
        
        # Updating the number of disks
        self.p_i +=1

    def add_disk52(self):
        """Add a new disk at equilibrum into the bacterium on both side of the bacterium"""
        alpha = 0

        # the disk is added in the head
        l = self.spring_rest_l
        if(self.p_i==1):
            alpha = random.uniform(0,2*np.pi)
        else :
            d_head = self.Disks[0]
            d_tail = self.Disks[self.p_i -1]
            vec1 = np.array([d_tail.X[0] - d_head.X[0],d_tail.X[1] - d_head.X[1]])
            vec2 = np.array([2,0])

            if(d_tail.X[1] - d_head.X[1]<=0):
                alpha = 2*np.pi - np.arccos(np.dot(vec1,vec2)/(norm(vec1)*norm(vec2)+self.eps))
            else:
                alpha = np.arccos(np.dot(vec1,vec2)/(norm(vec1)*norm(vec2)+self.eps))
            
        X = np.array([self.Disks[0].X[0] - l*np.cos(alpha),self.Disks[0].X[1] - l*np.sin(alpha)])
        self.Disks.insert(0,disk.Disk(X,ray=self.Disks[0].radius))
        
        # Updating the number of disks
        self.p_i +=1

        # the disk is added in the tail
        Y = np.array([self.Disks[self.p_i-1].X[0] + l*np.cos(alpha),self.Disks[self.p_i-1].X[1] + l*np.sin(alpha)])
        self.Disks.append(disk.Disk(Y,ray=self.Disks[0].radius))
        
        # Updating the number of disks
        self.p_i +=1

    def add_disk5(self):
        """Add a new disk not at equilibrum into the bacterium on both side of the bacterium"""

        # the disk is added in the head
        l = self.spring_rest_l
        alpha = 0

        #Calculation of the first angle
        if(self.p_i==1):
            alpha = random.uniform(0,2*np.pi)
        else :
            d_head = self.Disks[0]
            d_tail = self.Disks[1]
            vec1 = np.array([d_tail.X[0] - d_head.X[0],d_tail.X[1] - d_head.X[1]])
            vec2 = np.array([2,0])

            if(d_tail.X[1] - d_head.X[1]<=0):
                alpha = 2*np.pi - np.arccos(np.dot(vec1,vec2)/(norm(vec1)*norm(vec2)+self.eps))
            else:
                alpha = np.arccos(np.dot(vec1,vec2)/(norm(vec1)*norm(vec2)+self.eps))
            
        X = np.array([self.Disks[0].X[0] - l*np.cos(alpha),self.Disks[0].X[1] - l*np.sin(alpha)])
        self.Disks.insert(0,disk.Disk(X,ray=self.Disks[0].radius))
        
        # Updating the number of disks
        self.p_i +=1

        # Calculation of the second angle
        if(self.p_i-1!=1):
            d_head = self.Disks[self.p_i-2]
            d_tail = self.Disks[self.p_i-1]
            
            vec1 = np.array([d_tail.X[0] - d_head.X[0],d_tail.X[1] - d_head.X[1]])
            vec2 = np.array([2,0])

            if(d_tail.X[1] - d_head.X[1]<=0):
                alpha = 2*np.pi - np.arccos(np.dot(vec1,vec2)/(norm(vec1)*norm(vec2)+self.eps))
            else:
                alpha = np.arccos(np.dot(vec1,vec2)/(norm(vec1)*norm(vec2)+self.eps))

        # the disk is added in the tail
        Y = np.array([self.Disks[self.p_i-1].X[0] + l*np.cos(alpha),self.Disks[self.p_i-1].X[1] + l*np.sin(alpha)])
        self.Disks.append(disk.Disk(Y,ray=self.Disks[0].radius))
        
        # Updating the number of disks
        self.p_i +=1

    def add_disk6(self):
        """Add a new disk is added not at equilibrum into the bacterium on one side of the bacterium
        chosen randomly"""
        
        if self.side==-1:
            self.side = random.random()

        # the disk is added in the head
        if(self.side > 0.5):
            l = self.spring_rest_l

            # calculation of the angle
            if(self.p_i==1):
                alpha = random.uniform(0,2*np.pi)
            else :
                d_head = self.Disks[0]
                d_tail = self.Disks[1]
                vec1 = np.array([d_tail.X[0] - d_head.X[0],d_tail.X[1] - d_head.X[1]])
                vec2 = np.array([2,0])

                if(d_tail.X[1] - d_head.X[1]<=0):
                    alpha = 2*np.pi - np.arccos(np.dot(vec1,vec2)/(norm(vec1)*norm(vec2)+self.eps))
                else:
                    alpha = np.arccos(np.dot(vec1,vec2)/(norm(vec1)*norm(vec2)+self.eps))
                
            X = np.array([self.Disks[0].X[0] - l*np.cos(alpha),self.Disks[0].X[1] - l*np.sin(alpha)])
            self.Disks.insert(0,disk.Disk(X,ray=self.Disks[0].radius))
        
        # the disk is added in the tail
        else:
            l = self.spring_rest_l

            # Calculation of the angle
            if(self.p_i==1):
                alpha = random.uniform(0,2*np.pi)
            else :
                d_head = self.Disks[self.p_i -2]
                d_tail = self.Disks[self.p_i -1]
                vec1 = np.array([d_tail.X[0] - d_head.X[0],d_tail.X[1] - d_head.X[1]])
                vec2 = np.array([2,0])

                if(d_tail.X[1] - d_head.X[1]<=0):
                    alpha = 2*np.pi - np.arccos(np.dot(vec1,vec2)/(norm(vec1)*norm(vec2)+self.eps))
                else:
                    alpha = np.arccos(np.dot(vec1,vec2)/(norm(vec1)*norm(vec2)+self.eps))
            
            X = np.array([self.Disks[self.p_i-1].X[0] + l*np.cos(alpha),self.Disks[self.p_i-1].X[1] + l*np.sin(alpha)])
            self.Disks.append(disk.Disk(X,ray=self.Disks[0].radius))
        
        # Updating the number of disks
        self.p_i +=1

    def add_disk7(self):
        """Add a new disk is added not at equilibrum into the bacterium on one side of the bacterium
        chosen randomly"""
        
        # the disk is added in the head
        l = self.spring_rest_l

        # calculation of the angle
        if(self.p_i==1):
            alpha = random.uniform(0,2*np.pi)
        else :
            d_head = self.Disks[0]
            d_tail = self.Disks[1]
            vec1 = np.array([d_tail.X[0] - d_head.X[0],d_tail.X[1] - d_head.X[1]])
            vec2 = np.array([2,0])

            if(d_tail.X[1] - d_head.X[1]<=0):
                alpha = 2*np.pi - np.arccos(np.dot(vec1,vec2)/(norm(vec1)*norm(vec2)+self.eps))
            else:
                alpha = np.arccos(np.dot(vec1,vec2)/(norm(vec1)*norm(vec2)+self.eps))
            
        X = np.array([self.Disks[0].X[0] - l*np.cos(alpha),self.Disks[0].X[1] - l*np.sin(alpha)])
        self.Disks.insert(0,disk.Disk(X,ray=self.Disks[0].radius))
    
        
        # Updating the number of disks
        self.p_i +=1

    ## DIVISION

    def division(self):
        """Divide the bactierium in two : return 2 lists containing the disks of the bacterium
        If p_i is even then the lists have the same length, else one will be longer"""

        mid = self.p_i//2

        return self.Disks[0:mid],self.Disks[mid:self.p_i]
        
    ###-----------------  Position calculation -----------------------

    def Euler_explicit(self,dt):
        """Calculates the position of all the Disks by using 
        an Euler explicit integrator"""

        for i in range(self.p_i):
            self.Disks[i].X = self.Disks[i].X + dt*self.Disks[i].V      # Calculation in meter and s (Si)
        
if __name__=="__main__":
    print("test")