from lib2to3.pytree import convert
import tkinter as tk
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
        self.width = 800
        self.geometry(str(self.width) + "x" + str(self.height))

        #Micrometer to pixels conversion
        self.graduation = 1     # graduation in micrometers
        self.convert = 50 # pixel length of a graduation

        #axis offset from the side of the window
        self.axis_origin = 40   

        #Canvas creation
        self.canva = tk.Canvas(self,bg="light gray",width=self.width,height=self.height)
        self.canva.pack()

        #Drawing the axis
        self.axis = []
        self.text_axis = []
        self.draw_axis()

        #Drawing the bacteria
        self.list_disp_bacteria = [] #list of tuples of list (ovals and lines ) (for the display)
        self.draw_bacteria()

    def draw_axis(self):
        "Draw the axis"

        #axis graduation width
        grad_width = 5
        grad_text_offset = 9
        #x-axis
        self.axis.append(self.canva.create_line(self.axis_origin,self.height-self.axis_origin,((self.width-self.axis_origin)//self.convert)*self.convert + self.axis_origin ,
                                self.height-self.axis_origin))

        #y-axis
        self.axis.append(self.canva.create_line(self.axis_origin,self.height-self.axis_origin,self.axis_origin,
                                    (self.height-self.axis_origin)-((self.height-self.axis_origin)//self.convert)*self.convert))

        #drawing the graduation on the x-axis
        xpos = self.axis_origin
        ypos = self.height - self.axis_origin
        i=0
        while(xpos <= self.width):
            self.axis.append(self.canva.create_line(xpos,self.height-self.axis_origin-grad_width,xpos,
                                    self.height-self.axis_origin+grad_width))
            self.text_axis.append(self.canva.create_text(xpos,self.height-self.axis_origin+grad_width + grad_text_offset,text=str(i*self.graduation),font=('sans-serif 10')))
            i+=1
            xpos += self.convert
        
        i=0
        while(ypos > 0):
            self.axis.append(self.canva.create_line(self.axis_origin-grad_width,ypos,self.axis_origin+grad_width,
                                    ypos))
            self.text_axis.append(self.canva.create_text(self.axis_origin-grad_width - grad_text_offset,ypos,text=str(i*self.graduation),font=('sans-serif 10')))
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

            #drawing the line if it note the last

        return (ovals,lines)

    def draw_bacteria(self):
        """Draw all the bacteria of the simulation"""

        for bact in self.bacteria:
            self.list_disp_bacteria.append(self.draw_bacterium(bact))

        

if __name__=="__main__":
    app = Application()
    app.title("Bacteria micro-colonies simulator")
    app.mainloop()
