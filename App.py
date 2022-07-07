from os import environ
from turtle import Turtle, width
environ['PYGAME_HIDE_SUPPORT_PROMPT'] = '1'
import pygame
import numpy  as np
import disk
import bacterium as bact
import model
from numpy.linalg import norm
import subprocess
import os
import datetime
import random

class Application(model.Model):
    """Class for handleling the graphical aspect of the simulations
    It inherits of the Tk class from the tkinter library
    
    model components :
    - list of bacteria 
    
    Application component :
    - height of window
    - width of window
    - value of a graduation
    - pixel length of a graduation
    - axis origin (y-axis is reverted)"""

    def __init__(self):
        """Constructor for the display/tkinter class"""
        #Creation of the window/tk object
        pygame.init()
        model.Model.__init__(self)
        
        """Model parameters"""
        # self.generate_bacterium()
        for i in range(0,1):
            self.generate_random_bacteria(1)

        """Display parameters"""
        #Size of the window 
        self.height = 500
        self.width = 700
        self.window = pygame.display.set_mode([self.width,self.height])
        pygame.display.set_caption("Bacteria micro-colonies simulator")
        self.black =  (0,0,0)

        #Running state 
        self.running = True

        #Micrometer to pixels conversion
        self.graduation = 1     # graduation in micrometers
        self.convert = 40 # pixel length of a graduation

        #axis offset from the side of the window
        self.axis_origin = 30

        #axis origin (to handle the movements)   
        self.x_origin = -1
        self.y_origin = -1

        #Drawing the axis
        self.axis_state = True
        self.zoom_state = 1 #zoom and dezoom state
        self.draw_axis()

        #Drawing the bacteria
        self.drawing_bacteria_state=1
        self.draw_bacteria()

    def mainloop(self):
        """Main loop of the application/simulation"""

        while self.running:            
            # Events
            self.event()

            # Fill the background with white
            # self.window.fill((255, 255, 255))
            self.window.fill((220,220,220))

            #Move the bacteria
            self.bacteria_processes()

            # Draw the bacteria
            self.draw_bacteria()

            # Draw the axises
            if(self.axis_state):
                self.draw_axis()

            #Draw the informations
            self.draw_informations()

            # Increasing the time
            self.increment_time()

            # Updating the screen
            pygame.display.flip()
            pygame.time.delay(50)
    
    ### ------------------ Drawing methods ----------------------------------###

    def draw_axis(self):
        "Draw the axises"

        #axis graduation width
        grad_width = 5
        grad_text_offset = 0
        grad_text_xoffset = -3 

        #Font loading
        self.font = pygame.font.SysFont("msreferencesansserif",11)

        #x-axis
        pygame.draw.line(self.window,(0,0,0),start_pos=[self.axis_origin,self.height-self.axis_origin],end_pos=[((self.width-self.axis_origin)//self.convert)*self.convert + self.axis_origin ,
                                self.height-self.axis_origin],width=1)

        #y-axis
        pygame.draw.line(self.window,(0,0,0),start_pos=[self.axis_origin,self.height-self.axis_origin],end_pos=[self.axis_origin,
                                    (self.height-self.axis_origin)-((self.height-self.axis_origin)//self.convert)*self.convert],width=1)

        #drawing the graduation on the x-axis
        xpos = self.axis_origin
        ypos = self.height - self.axis_origin
        i=0
        while(xpos <= self.width):
            pygame.draw.line(self.window,(0,0,0),start_pos=[xpos,self.height-self.axis_origin],end_pos=[
                                    xpos,self.height-self.axis_origin+grad_width],width=1)
            if(self.graduation<1):
                text = self.font.render("{:.2f}".format(self.x_origin + i*self.graduation),True,self.black)
                self.window.blit(text,(xpos + grad_text_xoffset,self.height-self.axis_origin+grad_width + grad_text_offset))

            else:
                text = self.font.render(f"{int(self.x_origin + i*self.graduation)}",True,self.black)
                number_offset = len(f"{self.x_origin + i*self.graduation}")
                if(number_offset==1):
                    number_offset=0
                self.window.blit(text,(xpos + grad_text_xoffset - number_offset*2,self.height-self.axis_origin+grad_width + grad_text_offset))
            i+=1
            xpos += self.convert
        
        i=0
        grad_text_offset = 17
        grad_text_xoffset = -8 
        while(ypos > 0):
            pygame.draw.line(self.window,(0,0,0),start_pos=[self.axis_origin-grad_width,ypos],end_pos=[self.axis_origin,
                                    ypos],width=1)
            if(self.graduation<1):
                text = self.font.render("{:.2f}".format(self.y_origin + i*self.graduation),True,self.black)
                self.window.blit(text,(self.axis_origin-grad_width - grad_text_offset,ypos + grad_text_xoffset))
            else:
                text = self.font.render(f"{int(self.y_origin + i*self.graduation)}",True,self.black)
                number_offset = len(f"{int(self.y_origin + i*self.graduation)}")
                if(number_offset==1):
                    number_offset=0
                self.window.blit(text,(self.axis_origin-grad_width - grad_text_offset - number_offset*2,ypos + grad_text_xoffset))
            i+=1
            ypos -= self.convert
    
    def draw_bacterium(self,bact : bact.Bacterium):
        """Draw a single bacterium"""

        #Drawing the cells
        for i in range(0,bact.p_i):
            cdisk = bact.Disks[i] #Current disk
            
            #converting the positions from micrometers to pixels
            p_x = (cdisk.X[0]-self.x_origin)*self.convert/self.graduation + self.axis_origin
            p_y = self.height - (cdisk.X[1]-self.y_origin)*self.convert/self.graduation - self.axis_origin
            p_ray = cdisk.radius*self.convert/self.graduation


            #drawing the lines if it note the last
            if(i<bact.p_i-1):
                ndisk = bact.Disks[i+1] #Current cell

                #converting the positions from micrometers to pixels for the next cell
                p_x_next = (ndisk.X[0]-self.x_origin)*self.convert/self.graduation + self.axis_origin
                p_y_next = self.height - (ndisk.X[1]-self.y_origin)*self.convert/self.graduation - self.axis_origin
                p_ray_next = ndisk.radius*self.convert/self.graduation

                #drawing the lines
                pygame.draw.line(self.window,(0,0,255),[p_x,p_y],[p_x_next,p_y_next],width=1)
            
            #Drawing the actual cell
            if(i==0):
                pygame.draw.circle(self.window,self.black,(p_x,p_y),p_ray,width=1)
            else:
                pygame.draw.circle(self.window,bact.color,(p_x,p_y),p_ray,width=1)

    def draw_bacteria(self):
        """Draw all the bacteria of the simulation"""

        for bact in self.bacteria:
            if(self.drawing_bacteria_state==1):
                self.draw_bacterium(bact)
            elif(self.drawing_bacteria_state==2):
                self.draw_bacterium(bact)
                self.draw_bacterium_hull(bact)
            else:
                self.draw_bacterium_full(bact)

    def draw_informations(self):
        """Draw the simulations informations"""
    
        #Font loading
        self.font = pygame.font.SysFont("msreferencesansserif",11)

        #Informations positions
        x = self.width - 150
        y = 5
        distance = 20

        # Number of bacteria
        text = self.font.render(f"Number of bacteria : {self.N}",True,self.black)
        self.window.blit(text,(x,y))
        y+= distance

        # elapsed time
        text = self.font.render(f"Elapsed time : {self.time}",True,self.black)
        self.window.blit(text,(x,y))
        y+= distance

        # delta t
        text = self.font.render(f"dt : {self.dt}",True,self.black)
        self.window.blit(text,(x,y))

    def draw_bacterium_hull(self,bact : bact.Bacterium):
        """Draw the hull of the bacteria : spherocylinder"""
        # Points du dessus
        x = np.zeros((bact.p_i,2))
        x_next = np.zeros((bact.p_i,2))

        # Points du dessous
        y = np.zeros((bact.p_i,2))
        y_next = np.zeros((bact.p_i,2))

        #ray offset
        ray_offset = 0

        #iteration on the disks
        for i in range(bact.p_i-1):
            cdisk = bact.Disks[i] #Current disk
            ndisk = bact.Disks[i+1]

            #converting the positions from micrometers to pixels
            p_x = (cdisk.X[0]-self.x_origin)*self.convert/self.graduation + self.axis_origin
            p_y = self.height - (cdisk.X[1]-self.y_origin)*self.convert/self.graduation - self.axis_origin
            p_ray = cdisk.radius*self.convert/self.graduation + ray_offset
            p_x_next = (ndisk.X[0]-self.x_origin)*self.convert/self.graduation + self.axis_origin
            p_y_next = self.height - (ndisk.X[1]-self.y_origin)*self.convert/self.graduation - self.axis_origin
            p_ray_next = ndisk.radius*self.convert/self.graduation + ray_offset
            
            #Calculating the angle's vectors
            vec1 = np.array([p_x_next - p_x,p_y_next-p_y])
            vec2 = np.array([p_x,0])
            
            alpha = 0
            if(p_y_next - p_y <=0):
                alpha = 2*np.pi - np.arccos(np.dot(vec1,vec2)/(norm(vec1)*norm(vec2)))
            else:
                alpha = np.arccos(np.dot(vec1,vec2)/(norm(vec1)*norm(vec2)))

            #angle 
            alpha1 = np.pi/2 + alpha
            alpha2 = alpha - np.pi/2

            #Calculation of the new points for  the current cell
            x1 = np.array([p_x + p_ray*np.cos(alpha1),p_y + p_ray*np.sin(alpha1)])
            x2 = np.array([p_x + p_ray*np.cos(alpha2),p_y + p_ray*np.sin(alpha2)])

            # #Calculation of the new points for the next cell
            alpha = alpha - np.pi
            alpha1 = alpha + np.pi/2
            alpha2 = alpha - np.pi/2
    
            x3 = [p_x_next + p_ray*np.cos(alpha1),p_y_next + p_ray*np.sin(alpha1)]
            x4 = [p_x_next + p_ray*np.cos(alpha2),p_y_next + p_ray*np.sin(alpha2)]
            
            #drawing the lines
            pygame.draw.line(self.window,bact.color,x1,x4)
            pygame.draw.line(self.window,bact.color,x2,x3)

    def draw_bacterium_full(self,bact : bact.Bacterium):
        """Draw the hull of the bacteria : spherocylinder and colors  the inside"""
        # Points du dessus
        x = np.zeros((bact.p_i,2))
        x_next = np.zeros((bact.p_i,2))

        # Points du dessous
        y = np.zeros((bact.p_i,2))
        y_next = np.zeros((bact.p_i,2))

        #ray offset
        ray_offset = 0

        #iteration on the disks
        for i in range(bact.p_i-1):
            cdisk = bact.Disks[i] #Current disk
            ndisk = bact.Disks[i+1]

            #converting the positions from micrometers to pixels
            p_x = (cdisk.X[0]-self.x_origin)*self.convert/self.graduation + self.axis_origin
            p_y = self.height - (cdisk.X[1]-self.y_origin)*self.convert/self.graduation - self.axis_origin
            p_ray = cdisk.radius*self.convert/self.graduation + ray_offset
            p_x_next = (ndisk.X[0]-self.x_origin)*self.convert/self.graduation + self.axis_origin
            p_y_next = self.height - (ndisk.X[1]-self.y_origin)*self.convert/self.graduation - self.axis_origin
            p_ray_next = ndisk.radius*self.convert/self.graduation + ray_offset
            
            #Calculating the angle's vectors
            vec1 = np.array([p_x_next - p_x,p_y_next-p_y])
            vec2 = np.array([p_x,0])
            
            alpha = 0
            if(p_y_next - p_y <=0):
                alpha = 2*np.pi - np.arccos(np.dot(vec1,vec2)/(norm(vec1)*norm(vec2)))
            else:
                alpha = np.arccos(np.dot(vec1,vec2)/(norm(vec1)*norm(vec2)))

            #angle 
            alpha1 = np.pi/2 + alpha
            alpha2 = alpha - np.pi/2

            #Calculation of the new points for  the current cell
            x1 = np.array([p_x + p_ray*np.cos(alpha1),p_y + p_ray*np.sin(alpha1)])
            x2 = np.array([p_x + p_ray*np.cos(alpha2),p_y + p_ray*np.sin(alpha2)])

            # #Calculation of the new points for the next cell
            alpha = alpha - np.pi
            alpha1 = alpha + np.pi/2
            alpha2 = alpha - np.pi/2

            x3 = [p_x_next + p_ray*np.cos(alpha1),p_y_next + p_ray*np.sin(alpha1)]
            x4 = [p_x_next + p_ray*np.cos(alpha2),p_y_next + p_ray*np.sin(alpha2)]
            
            #drawing the lines
            pygame.draw.circle(self.window,bact.color,(p_x,p_y),p_ray)
            pygame.draw.polygon(self.window,bact.color,[x1,x4,x3,x2])

        cdisk = bact.Disks[bact.p_i -1] #Current disk
        #converting the positions from micrometers to pixels
        p_x = (cdisk.X[0]-self.x_origin)*self.convert/self.graduation + self.axis_origin
        p_y = self.height - (cdisk.X[1]-self.y_origin)*self.convert/self.graduation - self.axis_origin
        p_ray = cdisk.radius*self.convert/self.graduation + ray_offset
        pygame.draw.circle(self.window,bact.color,(p_x,p_y),p_ray)

    ### ----------------------- Events methods ---------------------------------------

    def event(self):
        """Manage the events of the simulation"""
        
        # Event loop
        for event in pygame.event.get():
                # Supporting the quiting button
                if event.type == pygame.QUIT:
                    self.running = False
                
                # Supporting a key press
                elif event.type == pygame.KEYDOWN:
                    # ESC to close the simulation
                    if event.key == pygame.K_ESCAPE:
                        self.stop()
                    
                    # K to dezoom
                    if event.key == pygame.K_k:
                        self.dezoom()
                    
                    # L to zoom
                    if event.key == pygame.K_l:
                        self.zoom()
                    
                    # A to hide/show the axises
                    if event.key == pygame.K_a:
                        self.axis_state = not(self.axis_state)
                    
                    # Z to add a new bacterium
                    if event.key == pygame.K_z:
                        self.generate_random_bacteria(1)
                        self.N_bacteria()
                    
                    # T to add a new disk into the first bacterium
                    if event.key == pygame.K_t:
                        self.bacteria[0].add_disk1()
                    
                    # R to change the drawing method
                    if event.key == pygame.K_r:
                        if(self.drawing_bacteria_state==1):
                            self.drawing_bacteria_state += 1
                        elif(self.drawing_bacteria_state==2):
                            self.drawing_bacteria_state += 1
                        elif(self.drawing_bacteria_state>2):
                            self.drawing_bacteria_state = 1
                    
                    # RIGTH ARROW to go right
                    if event.key == pygame.K_RIGHT:
                        self.x_origin += self.graduation
                    
                    # LEFT ARROW to go left
                    if event.key == pygame.K_LEFT:
                        self.x_origin -=self.graduation
                    
                     # UP ARROW to go up
                    if event.key == pygame.K_UP:
                        self.y_origin +=self.graduation
                    
                    # DOWN ARROW to go down
                    if event.key == pygame.K_DOWN:
                        self.y_origin -=self.graduation
                    
                    # P to take a screenshot
                    if event.key == pygame.K_p:
                        e = datetime.datetime.now()
                        pygame.image.save(self.window,"./screenshots/screenshot_" + "%s_%s_%s_" % (e.day, e.month, e.year) + "%s_%s_%s" % (e.hour, e.minute, e.second) + ".png")
                        print("> screenshot saved")

    def zoom(self):
        """Do a zoom"""
        if(1<self.zoom_state <=5):
            if(self.zoom_state==1):
                self.graduation = 1
                self.convert = 40
            elif(self.zoom_state==5):
                self.graduation = 1
                self.convert = 25
            else:
                self.convert += 5
            self.zoom_state -= 1

        elif(5<self.zoom_state<=10):
            if(self.zoom_state==10):
                self.graduation = 10
                self.convert = 50
            else:
                self.convert += 5
            self.zoom_state -=1

        elif(10<=self.zoom_state<=15):
            self.convert += 5
            self.zoom_state -=1

    def dezoom(self):
        """Do a zoom"""

        if(self.zoom_state <5):
            self.convert -= 6
            self.zoom_state += 1

        elif(5<=self.zoom_state<10):
            if(self.zoom_state==5):
                self.graduation = 10
                self.convert = 50
            else:
                self.convert -= 5
            self.zoom_state +=1

        elif(10<=self.zoom_state<15):
            if(self.zoom_state==10):
                self.graduation = 100
                self.convert = 50
            else:
                self.convert -= 5
            self.zoom_state +=1
    
    def stop(self):
        """Close the window"""
        self.running=False

        
if __name__=="__main__":
    app = Application()
    app.mainloop()
    pygame.quit()
