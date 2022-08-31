# Modelisation of the morphogenesis of bacteria micro-colonies by considering a cell-spring individual based model. 

![python](https://img.shields.io/badge/langage-Python-yellow)

## Introduction 

This project consist in modeling the morphogenesis of micro-colonies of bacteria that support bending.

## Model presentation

### The main interactions

In our model, a bacteria is constitued of a chain of disks linked by springs.
- This implies that the disks are subject to linear spring forces that controls the springs length and to torsion spring forces that control the angles of system of 3 disks. 

- Non overlapping forces are implemented as elactic collision interaction between two disks that are near enough.

- Exponential growth of the bacteria are also implemented.

- When the bacteria lengths reach a threshold, the bacteria divide themselves into two daugthters.

We use an energy based approach to calculate the velocity of each disks and then integrate Newton's equation of motion.

![model](https://github.com/Rudiio/Images-factory/blob/main/model.png)

## Install

The program relies on the pygame API for the graphical interface, on Numpy for the arrays and vector computations, Scipy for stats and Joblib.

To install pygame (with pip package manager)

``` pip install pygame```

To install Numpy (with pip package manager)

``` pip install numpy ```

To install Scipy (with pip package manager)

``` pip install Scipy ```

To install Joblib (with pip package manager)

``` pip install Joblib ```

To install Jupyter Notebook

``` pip install jupyter ```


## Files and scripts explanations

### **model.py bacterium.py disk.py**

They are the main scripts of the model. **disk.py** contains the class to represent a disk, **bacterium.py** contains the bacterium class and **model.py** contains the model class.

To launch a simulation and write its data into a text file, execute **model.py** with python. 
The parameters of the simulation can be modified from the constructor of the model class. 

The data as text files can be exploited with plot_model.py and the quantifiers.ipynb.

It is important to note that the increment.npy file goes on par with bacterium.py as it contains the increment distribution for E.coli and is used to calculate the increment of each bacterium.

### Simulations and Quantifiers

Simulations can be launched through the model.py file. They are saved into the simualtions folder.
These simulations are exploited to calculates the quantifiers of the model :

1. The aspect ratio **$\alpha_r$** : quantify the shape of the colony, more precisely its elongation. It needs two lengths, $l$ that can be interpreted as a Width and $L$ that is similar to a Length. In our case, we calculate the ellipse that fit the convex hull of the colony : $l$ is the semi-minor axis and $L$ the semi-major axis.
    Then $\alpha_r= \frac{l}{L}$.
    
2. The density **$\delta$**: it is calculated with image analysis tool. It is the ration between the surface occupied by the bacteria and the surface of the hull of the colony.
    
3. The length quantifier **$s$** : it represents the global distance of the linear springs from their rest length/equilibrium.
    
    $s = \frac{1}{N} \sum^N_{i=1} \frac{1}{p_i - 1} \sum^{p_i-1}_{j=1} ||X_{j+1}^i - X_j^i| - l|$
    
4. The bending energy **$E_t$** : it quantifies the bending of the bacteria. 0 represents the perfect bending equilibrium. It is calculated from the second and third terms of the equation (1).

Thy can be calculated by the quantifiers jupyter notebook and ploted by the quantifiers_plot notebook.


### App.py (with model.py, bacterium.py and disk.py)

This script can run live simulations with graphical interface (pygame). It uses the previous file to simulate and contains an App class that handle the graphical part.

The graphical parameters (ex : window's width or heigth ...) can be modified in the class's constructor.

Here are the main control for the graphical interface:

- **A** : show/hide the axises
- **K** : zoom
- **L** : dezoom
- **E** : Recenter
- **P** : take a screenshot that will be saved in the screenshots folder
- **Z** : generate a random bacterium
- **R** : change the display mode of the bacteria
- **← → ↑ ↓** : movements (change the origin)

**Note:** the interface is optimized for windows and my computer so it can produce glitch on another computer or OS. The solution is to modify the parameters to fit your computer.

### Plot_model.py

This scripts displays a colony from a text file produced by **model.py** with the same interface as **App.py**. The parameters can be modified in the constructor of the Plot class.
The time can be chosen when creating an object but it is by default the final time.
This file also allows to have an animated evolution of the colony or even store the screen for each time.

### Plot_model.py

It does the same things as **Plot_model.py** but for the text file from the experimental data produced by **data_covert.ipynb**.

Technically, it should support also the other text files but it was not tested yet. The reverse is strictly not possible.

### data_convert.ipynb

Converts the experimental data from Duvernoy into a compatible data for our files and model.
The results are text files in the same format that with **model.py**.

### data.ipynb

It extracts the data from the experimental data from Duvernoy. The results are plots of the distributions of the different important features. Everything is explained into the file.

## Simulations

The simulations I did are contained in the simulations folder with the file to exploit them.
The converted data from the experimental data are contained in the Data_Duvernoy_converted folder as well as the scripts to exploit them.

![simu](https://github.com/Rudiio/Images-factory/blob/main/simu11.png)