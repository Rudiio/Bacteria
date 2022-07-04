from lib2to3.pytree import convert
import tkinter as tk
from turtle import color
from bleach import clean
import numpy  as np
import cell
import bacterium as bact
import model



class Application(tk.Tk,model.Model):
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
        tk.Tk.__init__(self)
        model.Model.__init__(self)
        
        """Model parameters"""
        self.generate_bacterium()

        """Display parameters"""
        #Size of the window 
        self.height = 500
        self.width = 700
        self.geometry(str(self.width) + "x" + str(self.height))
        self.delay = 50 #delay of the callback in ms

        #Micrometer to pixels conversion
        self.graduation = 1     # graduation in micrometers
        self.convert = 40 # pixel length of a graduation

        #axis offset from the side of the window
        self.axis_origin = 30   

        #Canvas creation
        self.canva = tk.Canvas(self,bg="white",width=self.width,height=self.height)
        self.canva.pack()

        #Drawing the axis
        self.zoom_state = 1 #zoom and dezoom state
        self.axis = []
        self.text_axis = []
        self.draw_axis()

        #Drawing the bacteria
        self.list_disp_bacteria = [] #list of tuples of list (ovals and lines ) (for the display)
        self.draw_bacteria()

        # self.canva.create_line(self.width/2,self.height/2,self.width/2,self.height/2)

        #Action binding
        self.bind("<Escape>",self.stop)
        self.bind("<KeyPress-l>",self.zoom)  #zoom (modification of graduation and length)
        self.bind("<KeyPress-k>",self.dezoom)  #zoom (modification of graduation and length)
        
        #Launching the simulation
        self.simulation()

    def simulation(self):
        """Execution of all the model processes"""

        #Calculating the new positions of the bacteria
        self.Move_bacteria()
        print(self.bacteria[0])
        self.clean_bacteria()
        self.draw_bacteria()
        self.after(self.delay,self.simulation)

    def draw_axis(self):
        "Draw the axises"

        #axis graduation width
        grad_width = 5
        grad_text_offset = 9
        #x-axis
        self.axis.append(self.canva.create_line(self.axis_origin,self.height-self.axis_origin,((self.width-self.axis_origin)//self.convert)*self.convert + self.axis_origin ,
                                self.height-self.axis_origin,smooth=1))

        #y-axis
        self.axis.append(self.canva.create_line(self.axis_origin,self.height-self.axis_origin,self.axis_origin,
                                    (self.height-self.axis_origin)-((self.height-self.axis_origin)//self.convert)*self.convert,smooth=1))

        #drawing the graduation on the x-axis
        xpos = self.axis_origin
        ypos = self.height - self.axis_origin
        i=0
        while(xpos <= self.width):
            self.axis.append(self.canva.create_line(xpos,self.height-self.axis_origin,xpos,
                                    self.height-self.axis_origin+grad_width,smooth=1))
            if(self.graduation<1):
                self.text_axis.append(self.canva.create_text(xpos,self.height-self.axis_origin+grad_width + grad_text_offset,text="{:.2f}".format(i*self.graduation),font=('sans-serif 10')))
            else:
                self.text_axis.append(self.canva.create_text(xpos,self.height-self.axis_origin+grad_width + grad_text_offset,text=f"{int(i*self.graduation)}",font=('sans-serif 10')))
            i+=1
            xpos += self.convert
        
        i=0
        while(ypos > 0):
            self.axis.append(self.canva.create_line(self.axis_origin-grad_width+1,ypos,self.axis_origin+1,
                                    ypos,smooth=1))
            if(self.graduation<1):
                self.text_axis.append(self.canva.create_text(self.axis_origin-grad_width - grad_text_offset,ypos,text="{:.2f}".format(i*self.graduation),font=('sans-serif 10')))
            else:
                self.text_axis.append(self.canva.create_text(self.axis_origin-grad_width - grad_text_offset,ypos,text=f"{int(i*self.graduation)}",font=('sans-serif 10')))
            i+=1
            ypos -= self.convert
    
    def draw_bacterium(self,bact : bact.Bacterium):
        """Draw a single bacterium"""

        #Creating the lists
        ovals = []
        lines = []

        #Drawing the cells
        for i in range(0,bact.N):
            ccell = bact.Cells[i] #Current cell
            
            #converting the positions from micrometers to pixels
            p_x = ccell.X[0]*self.convert/self.graduation + self.axis_origin
            p_y = self.height - ccell.X[1]*self.convert/self.graduation - self.axis_origin
            p_ray = ccell.ray*self.convert/self.graduation

            #Drawing the actual cell
            ovals.append(self.canva.create_oval(p_x-p_ray,p_y-p_ray,p_x+p_ray,p_y+p_ray,width=1))

            #drawing the lines if it note the last
            if(i<bact.N-1):
                ncell = bact.Cells[i+1] #Current cell

                #converting the positions from micrometers to pixels for the next cell
                p_x_next = ncell.X[0]*self.convert/self.graduation + self.axis_origin
                p_y_next = self.height - ncell.X[1]*self.convert/self.graduation - self.axis_origin
                p_ray_next = ncell.ray*self.convert/self.graduation

                #drawing the lines
                lines.append(self.canva.create_line(p_x,p_y,p_x_next,p_y_next,width=1,fill="blue",smooth=1))
                # lines.append(self.canva.create_line(p_x,p_y+p_ray,p_x_next,p_y_next+p_ray_next,width=1))
                # lines.append(self.canva.create_line(p_x,p_y-p_ray,p_x_next,p_y_next-p_ray_next,width=1))

        return (ovals,lines)

    def draw_bacteria(self):
        """Draw all the bacteria of the simulation"""

        for bact in self.bacteria:
            self.list_disp_bacteria.append(self.draw_bacterium(bact))
    
    def clean_axis(self):
        """Delete the axises on the window"""
        
        #cleaning the lines
        for i in self.axis:
            self.canva.delete(i)
        self.axis.clear()
        
        #cleaning the texts
        for i in self.text_axis:
            self.canva.delete(i)
        self.text_axis.clear()

    def clean_bacteria(self):
        """Delete the bacteria on the display"""
        
        #iterating on the bacteria
        for ids in self.list_disp_bacteria:
            for j in ids[0]:
                self.canva.delete(j)
            for j in ids[1]:
                self.canva.delete(j)
            ids[0].clear()
            ids[1].clear()

    def zoom(self,k):
        """Do a zoom"""
        if(1<self.zoom_state<=5):
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

        #deleting the axises and the bacterias
        self.clean_axis()
        self.clean_bacteria()

        #re drawing the axises and the bacterias
        self.draw_axis()
        self.draw_bacteria()
    
    def dezoom(self,l):
        """Do a zoom"""

        if(self.zoom_state<5):
            self.convert -= 5
            self.zoom_state += 1

        elif(5<=self.zoom_state<10):
            if(self.zoom_state==5):
                self.graduation = 10
                self.convert = 100
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

        #deleting the axises and the bacterias
        self.clean_axis()
        self.clean_bacteria()

        #re drawing the axises and the bacterias
        self.draw_axis()
        self.draw_bacteria()

    def stop(self,esc):
        """Close the simulation"""
        self.quit()

        

if __name__=="__main__":
    app = Application()
    app.title("Bacteria micro-colonies simulator")
    app.mainloop()
