# Modelisation of the morphogenesis of bacteria micro-colonies by considering a cell-spring individual based model. 

![python](https://img.shields.io/badge/langage-Python-yellow)

## Introduction 

This project consist in modeling the morphogenesis of micro-colonies of bacteria that support bending.

## Model presentation

In our model, a bacteria is constitued of a chain of disks linked by springs.
- This implies that the disks are subject to linear spring forces that controls the springs length and to torsion spring forces that control the angles of system of 3 disks. 

- Non overlapping forces are implemented as elactic collision interaction between two disks that are near enough.

- Exponential growth of the bacteria are also implemented.

- When the bacteria lengths reach a threshold, the bacteria divide themselves into two daugthters.

We use an energy based approach to calculate the velocity of each disks and then integrate Newton's equation of motion.

![model](https://github.com/Rudiio/Images-factory/blob/main/model.png)

## Simulation parameters

The data notebook extracts the numerical parameters from datasets of experimental bacteria (Coli or Pseudomonas).
Distributions of the parameters are made to understand the simulations parameters and to compare the results.

Here are the studied parameters:
- **cname** : the name of a bacteria
- **frame** : the frame during wich the bacteria was observed.
- **time** : the time of the observation
- **cellno** : the number of the cell on the current frame
- **schnitzno** : the "id" of a bacterium
- **birth** : the birth date
- **lifespan** : lifetime of a bacterium before division
- **area** : area of the bacterium
- **lentgh** :
- **angle**:
- **grate** : growth rate of the bacterium
- **width** 

## Simulations and Quantifiers

Simulations can be launched through the model.py file. They are saved into the simualtions folder.
These simulations are exploited to calculates the quantifiers of the model :

1. The aspect ratio **$\alpha_r$** : quantify the shape of the colony, more precisely its elongation. It needs two lengths, $l$ that can be interpreted as a Width and $L$ that is similar to a Length. In our case, we calculate the ellipse that fit the convex hull of the colony : $l$ is the semi-minor axis and $L$ the semi-major axis.
    Then $\alpha_r= \frac{l}{L}$.
    
2. The density **$\delta$**: it is calculated with image analysis tool. It is the ration between the surface occupied by the bacteria and the surface of the hull of the colony.
    
3. The length quantifier **$s$** : it represents the global distance of the linear springs from their rest length/equilibrium.
    
    $s = \frac{1}{N} \sum^N_{i=1} \frac{1}{p_i - 1} \sum^{p_i-1}_{j=1} ||X_{j+1}^i - X_j^i| - l|$
    
4. The bending energy **$E_t$** : it quantifies the bending of the bacteria. 0 represents the perfect bending equilibrium. It is calculated from the second and third terms of the equation (1).

Thy can be calculated by the quantifiers jupyter notebook and ploted by the quantifiers_plot notebook.

## Commands

- **A** : show/hide the axises
- **K** : zoom
- **L** : dezoom
- **E** : Recenter
- **P** : take a screenshot
- **Z** : generate a random bacterium
- **R** : change the display mode of the bacteria
- **← → ↑ ↓** : movements (change the origin)

## TODO

- [X] drawing the case of the bacteria (needs improvements)
- [X] compute spring/internal interactions
- [X] external forces 
- [X] graphical improvement (auto scaling)
- [X] collision detection and interbacteria interactions
- [ ] Calculation improvements
- [X] growth
- [X] bacteria division
- [ ] In depth parameters study


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

## Simulations

![simu](https://github.com/Rudiio/Images-factory/blob/main/simu11.png)