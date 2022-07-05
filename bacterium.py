import numpy as np
from pyparsing import NoMatch
import disk 
from numpy.linalg import norm


class Bacterium:
    """ Class representing a bacterium 
    It composed of :
    - the number of Disks in the Bacterium: N
    - a list of Disks: Disks
    - the rest torque theta
    - the rest length of the springs l 
    - the stiffness constants ks,kt1, kt2"""

    def __init__(self,N=0,Disks : list[disk.Disk] = [],l = 1.0,theta = np.pi,color = (0,0,0)):
        """Constructor for the class Bacterium"""
        self.p_i = N
        self.Disks = Disks
        self.theta = theta
        self.spring_rest_l = l
        self.ks = 1.e-2
        self.kt_par = 1.e-2
        self.kt_bot = 1.e-2
        self.color = color
        self.eps = 0.01

    def __str__(self):
        """Display the bacterium informations whit the print function"""
        print(f"N = {self.p_i}")
        for i in self.Disks :
            print(i)
        return("")
    
    def spring_velocity(self):
        """Calculate the velocity that comes from the spring forces/torques
        of each cell of the bacterium"""

        # Loop on all the Disks
        for k in range(self.p_i):
            # print(self.torsion_spring_par(k))
            self.Disks[k].V = self.linear_spring(k) + self.torsion_spring_par(k) + self.torsion_spring_bot(k)

    def linear_spring(self,k):
        """Calculates the velocity created by the linear springs"""
        if k==0:
            Xj = self.Disks[k].X
            Xj1 = self.Disks[k+1].X

            return (self.ks/(self.spring_rest_l**2)*(norm(Xj1 - Xj) - self.spring_rest_l)*
                    (Xj1 - Xj)/(norm(Xj1 - Xj))) 
        elif k==self.p_i-1:
            Xj = self.Disks[k].X
            Xj_1 = self.Disks[k-1].X

            return (-self.ks/(self.spring_rest_l**2)*(norm(Xj - Xj_1) - self.spring_rest_l)*
                    (Xj - Xj_1)/norm(Xj-Xj_1))
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
        if k==0 and self.p_i>=3 :  
            Xj = self.Disks[k].X
            Xj1 = self.Disks[k+1].X
            Xj2 = self.Disks[k+2].X

            # term in j,j+1,j+2 with cos
            V -= self.kt_par/(norm(Xj2 - Xj1)*norm(Xj - Xj1) +self.eps)*((np.dot((Xj2 - Xj1),(Xj - Xj1)))/(norm(Xj2 - Xj1)*norm(Xj - Xj1)+self.eps) - np.cos(self.theta))*\
                ((Xj2 - Xj1)-np.dot((Xj2 - Xj1),(Xj - Xj1))*(Xj - Xj1)/(norm(Xj-Xj1)**2+self.eps))
        
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

        elif k==self.p_i-1 and self.p_i>=3:
            Xj_2 = self.Disks[k-2].X
            Xj_1 = self.Disks[k-1].X
            Xj = self.Disks[k].X

            # term in j-2,j-1,j with cos
            V -= self.kt_par/(norm(Xj -Xj_1)*norm(Xj_2 - Xj_1)+self.eps)*((np.dot((Xj-Xj_1),(Xj_2 - Xj_1)))/(norm(Xj -Xj_1)*norm(Xj_2 - Xj_1)+self.eps) - np.cos(self.theta))*\
                ((Xj_2 - Xj_1)-np.dot((Xj-Xj_1),(Xj_2 - Xj_1))*(Xj - Xj_1)/(norm(Xj - Xj_1)**2+self.eps))

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

        if k==0 and self.p_i>=3:
            Xj = self.Disks[k].X
            Xj1 = self.Disks[k+1].X
            Xj2 = self.Disks[k+2].X

            # term in j,j+1,j+2 with sin
            V -= self.kt_bot/(norm(Xj2 - Xj1)*norm(Xj - Xj1) + self.eps)*((np.cross((Xj2 - Xj1),(Xj -Xj1)))/(norm(Xj2 - Xj1)*norm(Xj - Xj1)+ self.eps) - np.sin(self.theta))*\
                np.array([-(Xj2-Xj1)[1]+(((Xj-Xj1)[0])/(norm(Xj-Xj1)**2 + self.eps))*np.cross((Xj-Xj1),(Xj2-Xj1)),(Xj2 - Xj1)[0] + ((Xj - Xj1)[1])/(norm(Xj -Xj1)**2 + self.eps)*np.cross((Xj-Xj1),(Xj2-Xj1))])
        
        elif k==1:
            Xj = self.Disks[k].X
            Xj_1 = self.Disks[k-1].X

            # term in j-1,j,j+1 with sin
            if(self.p_i>=3):
                Xj1 = self.Disks[k+1].X

                V += self.kt_bot/(norm(Xj1 - Xj)*norm(Xj_1 - Xj) + self.eps)*((np.cross((Xj1 - Xj),(Xj_1 - Xj)))/(norm(Xj1 - Xj)*norm(Xj_1 - Xj) + self.eps) - np.sin(self.theta))*\
                    (np.array([(Xj_1 - Xj)[1] - ((Xj1-Xj)[0])/(norm(Xj1 - Xj)**2 + self.eps)*np.cross((Xj1 - Xj),(Xj_1 -Xj)),-(Xj_1 - Xj)[0] - ((Xj1-Xj)[1])/(norm(Xj1 - Xj)**2)*np.cross((Xj1 - Xj),(Xj_1 -Xj))]) +\
                        np.array([-(Xj1 -Xj)[1] + ((Xj_1 -Xj)[0])/(norm(Xj_1 -Xj)**2 + self.eps)*np.cross((Xj1 - Xj),(Xj_1 - Xj)),(Xj1 -Xj)[0] + ((Xj_1 -Xj)[1])/(norm(Xj_1 -Xj)**2)*np.cross((Xj1 - Xj),(Xj_1 - Xj))]))
            
            if(self.p_i>3):
                Xj2 = self.Disks[k+2].X

                # term in j,j+1,j+2 with sin
                V -= self.kt_bot/(norm(Xj2 - Xj1)*norm(Xj - Xj1) + self.eps)*((np.cross((Xj2 - Xj1),(Xj -Xj1)))/(norm(Xj2 - Xj1)*norm(Xj - Xj1)+ self.eps) - np.sin(self.theta))*\
                np.array([-(Xj2-Xj1)[1]+(((Xj-Xj1)[0])/(norm(Xj-Xj1)**2))*np.cross((Xj-Xj1),(Xj2-Xj1)),(Xj2 - Xj1)[0] + ((Xj - Xj1)[1])/(norm(Xj -Xj1)**2)*np.cross((Xj-Xj1),(Xj2-Xj1))])

        elif k==self.p_i - 1 and self.p_i>=3:
            Xj_2 = self.Disks[k-2].X
            Xj_1 = self.Disks[k-1].X
            Xj = self.Disks[k].X

            # term in j-2,j-1,j with sin
            V -= self.kt_bot/(norm(Xj - Xj_1)*norm(Xj_2 - Xj_1) + self.eps)*((np.cross((Xj - Xj_1),(Xj_2 - Xj_1)))/(norm(Xj - Xj_1)*norm(Xj_2 - Xj_1)+ self.eps) -np.sin(self.theta))*\
                np.array([(Xj_2 - Xj_1)[1] - ((Xj - Xj_1)[0])/(norm(Xj-Xj_1)**2 + self.eps)*np.cross((Xj - Xj_1),(Xj_2 -Xj_1)),-(Xj_2 - Xj_1)[0] - ((Xj - Xj_1)[1])/(norm(Xj-Xj_1)**2 + self.eps)*np.cross((Xj - Xj_1),(Xj_2 -Xj_1))])
        
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
                    (np.array([(Xj_1 - Xj)[1] - ((Xj1-Xj)[0])/(norm(Xj1 - Xj)**2 + self.eps)*np.cross((Xj1 - Xj),(Xj_1 -Xj)),-(Xj_1 - Xj)[0] - ((Xj1-Xj)[1])/(norm(Xj1 - Xj)**2)*np.cross((Xj1 - Xj),(Xj_1 -Xj))]) +\
                        np.array([-(Xj1 -Xj)[1] + ((Xj_1 -Xj)[0])/(norm(Xj_1 -Xj)**2 + self.eps)*np.cross((Xj1 - Xj),(Xj_1 - Xj)),(Xj1 -Xj)[0] + ((Xj_1 -Xj)[1])/(norm(Xj_1 -Xj)**2)*np.cross((Xj1 - Xj),(Xj_1 - Xj))]))

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
                (np.array([(Xj_1 - Xj)[1] - ((Xj1-Xj)[0])/(norm(Xj1 - Xj)**2 + self.eps)*np.cross((Xj1 - Xj),(Xj_1 -Xj)),-(Xj_1 - Xj)[0] - ((Xj1-Xj)[1])/(norm(Xj1 - Xj)**2)*np.cross((Xj1 - Xj),(Xj_1 -Xj))]) +\
                    np.array([-(Xj1 -Xj)[1] + ((Xj_1 -Xj)[0])/(norm(Xj_1 -Xj)**2 + self.eps)*np.cross((Xj1 - Xj),(Xj_1 - Xj)),(Xj1 -Xj)[0] + ((Xj_1 -Xj)[1])/(norm(Xj_1 -Xj)**2)*np.cross((Xj1 - Xj),(Xj_1 - Xj))]))

            # term in j,j+1,j+2 with sin
            V -= self.kt_bot/(norm(Xj2 - Xj1)*norm(Xj - Xj1) + self.eps)*((np.cross((Xj2 - Xj1),(Xj -Xj1)))/(norm(Xj2 - Xj1)*norm(Xj - Xj1)+ self.eps) - np.sin(self.theta))*\
                np.array([-(Xj2-Xj1)[1]+(((Xj-Xj1)[0])/(norm(Xj-Xj1)**2 + self.eps))*np.cross((Xj-Xj1),(Xj2-Xj1)),(Xj2 - Xj1)[0] + ((Xj - Xj1)[1])/(norm(Xj -Xj1)**2 + self.eps)*np.cross((Xj-Xj1),(Xj2-Xj1))])
        
        return V

    def Euler_explicit(self,dt):
        """Calculates the position of all the Disks by using 
        an Euler explicit integrator"""

        for i in range(self.p_i):
            self.Disks[i].X = self.Disks[i].X + dt*self.Disks[i].V
        

if __name__=="__main__":
    print("test")