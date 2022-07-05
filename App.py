from os import environ
from turtle import Turtle
import black
environ['PYGAME_HIDE_SUPPORT_PROMPT'] = '1'
import pygame
import numpy  as np
import disk
import bacterium as bact
import model


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
        self.generate_bacterium()

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

        #Drawing the axis
        self.axis_state = True
        self.zoom_state = 1 #zoom and dezoom state
        self.draw_axis()

        #Drawing the bacteria
        self.draw_bacteria()


    def mainloop(self):
        """Main loop of the application/simulation"""

        while self.running:
            # Increasing the time
            self.increment_time()
            
            # Events
            self.event()

            # Fill the background with white
            # self.window.fill((255, 255, 255))
            self.window.fill((220,220,220))

            # Draw the axises
            if(self.axis_state):
                self.draw_axis()

            #Draw the informations
            self.draw_informations()
            
            #Move the bacteria
            self.Move_bacteria()
            # print(self.bacteria[1])

            # Draw the bacteria
            self.draw_bacteria()

            # Updating the screen
            pygame.display.flip()
            # pygame.time.delay(50)
    
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
                text = self.font.render("{:.2f}".format(i*self.graduation),True,self.black)
                self.window.blit(text,(xpos + grad_text_xoffset,self.height-self.axis_origin+grad_width + grad_text_offset))
                # self.text_axis.append(self.canva.create_text(xpos,self.height-self.axis_origin+grad_width + grad_text_offset,text="{:.2f}".format(i*self.graduation),font=('sans-serif 10')))
            else:
                text = self.font.render(f"{int(i*self.graduation)}",True,self.black)
                self.window.blit(text,(xpos + grad_text_xoffset,self.height-self.axis_origin+grad_width + grad_text_offset))
                # self.text_axis.append(self.canva.create_text(xpos,self.height-self.axis_origin+grad_width + grad_text_offset,text=f"{int(i*self.graduation)}",font=('sans-serif 10')))
            i+=1
            xpos += self.convert
        
        i=0
        grad_text_offset = 17
        grad_text_xoffset = -8 
        while(ypos > 0):
            pygame.draw.line(self.window,(0,0,0),start_pos=[self.axis_origin-grad_width,ypos],end_pos=[self.axis_origin,
                                    ypos],width=1)
            if(self.graduation<1):
                text = self.font.render("{:.2f}".format(i*self.graduation),True,self.black)
                self.window.blit(text,(self.axis_origin-grad_width - grad_text_offset,ypos + grad_text_xoffset))
                # self.canva.create_text(self.axis_origin-grad_width - grad_text_offset,ypos,text="{:.2f}".format(i*self.graduation),font=('sans-serif 10')))
            else:
                text = self.font.render(f"{int(i*self.graduation)}",True,self.black)
                self.window.blit(text,(self.axis_origin-grad_width - grad_text_offset,ypos + grad_text_xoffset))
                # self.text_axis.append(self.canva.create_text(self.axis_origin-grad_width - grad_text_offset,ypos,text=f"{int(i*self.graduation)}",font=('sans-serif 10')))
            i+=1
            ypos -= self.convert
    
    def draw_bacterium(self,bact : bact.Bacterium):
        """Draw a single bacterium"""

        #Drawing the cells
        for i in range(0,bact.p_i):
            ccell = bact.Disks[i] #Current cell
            
            #converting the positions from micrometers to pixels
            p_x = ccell.X[0]*self.convert/self.graduation + self.axis_origin
            p_y = self.height - ccell.X[1]*self.convert/self.graduation - self.axis_origin
            p_ray = ccell.ray*self.convert/self.graduation
            # print(p_ray)

            #Drawing the actual cell
            pygame.draw.circle(self.window,bact.color,(p_x,p_y),p_ray,width=1)

            #drawing the lines if it note the last
            if(i<bact.p_i-1):
                ncell = bact.Disks[i+1] #Current cell

                #converting the positions from micrometers to pixels for the next cell
                p_x_next = ncell.X[0]*self.convert/self.graduation + self.axis_origin
                p_y_next = self.height - ncell.X[1]*self.convert/self.graduation - self.axis_origin
                p_ray_next = ncell.ray*self.convert/self.graduation

                #drawing the lines
                pygame.draw.line(self.window,(0,0,255),[p_x,p_y],[p_x_next,p_y_next],width=1)


    def draw_bacteria(self):
        """Draw all the bacteria of the simulation"""

        for bact in self.bacteria:
            self.draw_bacterium(bact)

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



    ### ----------------------- Events methods ---------------------------------------

    def event(self):
        """Manage the events of the simulation"""
        
        # Event loop
        for event in pygame.event.get():
                if event.type == pygame.QUIT:
                    self.running = False
                elif event.type == pygame.KEYDOWN:
                    if event.key == pygame.K_ESCAPE:
                        self.stop()
                    if event.key == pygame.K_k:
                        self.dezoom()
                    if event.key == pygame.K_l:
                        self.zoom()
                    if event.key == pygame.K_a:
                        self.axis_state = not(self.axis_state)

    def zoom(self):
        """Do a zoom"""
        if(1<self.zoom_state <=5):
            if(self.zoom_state==5):
                self.graduation = 1
                self.convert = 50
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
            # if(self.zoom==15):
            #     self.graduation = 100
            #     self.convert = 50
            # else:
            self.convert += 5
            self.zoom_state -=1

    
    def dezoom(self):
        """Do a zoom"""

        if(self.zoom_state <5):
            self.convert -= 5
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
