import tkinter as tk
import numpy  as np
import cell
import bacterium as bact



class Application(tk.Tk):
    """Class for handleling the graphical aspect of the simulations
    It inherits of the Tk class from the tkinter library"""

    def __init__(self):
        """Constructor for the display/tkinter class"""

        #Creation of the window/tk object
        tk.Tk.__init__(self)

        #Size of the window 
        self.height = 500
        self.width = 800
        self.geometry(str(self.width) + "x" + str(self.height))

        #Creation of the list of bacteria
        Cell1 = cell.Cell(X=np.array([[self.width/2,self.height/2]]),mass=0.1,ray=0.5)
        # print(Cell1)
        self.bacteria = [bact.Bacterium(N=1,Cells=[Cell1],Torques=[])]
        # print(self.bacteria[0])

        #Micrometer to pixels conversion
        self.convert = 50 #1 micrometer = 50 pixels

        #axis offset from the side of the window 
        self.axis_offset = 40 

        #Canvas creation
        self.canva = tk.Canvas(self,bg="light gray",width=self.width,height=self.height)
        self.canva.pack()

        #Drawing the axis
        self.draw_axis()
    
    def draw_axis(self):
        "Draw the axis"

        #axis graduation width
        grad = 5

        #x-axis
        self.canva.create_line(self.axis_offset,self.height-self.axis_offset,((self.width-self.axis_offset)//self.convert)*self.convert + self.axis_offset ,
                                self.height-self.axis_offset)

        #y-axis
        self.canva.create_line(self.axis_offset,self.height-self.axis_offset,self.axis_offset,
                                    (self.height-self.axis_offset)-((self.height-self.axis_offset)//self.convert)*self.convert)

        #drawing the graduation on the x-axis
        xpos = self.axis_offset
        ypos = self.height - self.axis_offset

        while(xpos <= self.width):
            self.canva.create_line(xpos,self.height-self.axis_offset-grad,xpos,
                                    self.height-self.axis_offset+grad)
            xpos += self.convert
        
        while(ypos > 0):
            self.canva.create_line(self.axis_offset-grad,ypos,self.axis_offset+grad,
                                    ypos)
            ypos -= self.convert

if __name__=="__main__":
    app = Application()
    app.title("Bacteria micro-colonies simulator")
    app.mainloop()
