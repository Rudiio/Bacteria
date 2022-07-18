# Modelisation of the morphogenesis of bacteria micro-colonies by considering a cell-spring individual based model. 

![python](https://img.shields.io/badge/langage-Python-yellow)

## Introduction 

This project consist in modeling the morphogenesis of micro-colonies of bacteria that support bending.

## Model presentation
In our model, a bacteria is constitued of a chain of disks linked by springs.
- This implies that the disks are subject to linear spring forces that controls the springs length and to torsion spring forces that control the angles of system of 3 disks. 

- Exponential growth of the bacteria are also implemented.

![model](https://github.com/Rudiio/Images-factory/blob/main/model.png)

## Simulation dimensions

The dimension for a bacteria is around 10 micrometers. On the display, it is possible to change the scale micrometer/pixel, but it is natively fixed to 1 um = 40 p


## Commands
- **A** : show/hide the axises
- **K** : zoom
- **L** : dezoom
- **p** : take a screenshot
- **z** : generate a random bacterium
- **r** : change the display mode of the bacteria
- **← → ↑ ↓** : movements (change the origin)

## TODO

- [X] drawing the case of the bacteria (needs improvements)
- [X] compute spring/internal interactions
- [X] external forces 
- [ ] graphical improvement (auto scaling)
- [ ] collision detection
- [ ] growth
- [ ] bacteria division
- [ ] external interactions
- [ ] bacteria movement

## Dependencies

The program relies on the pygame API for the graphical interface and on Numpy for the arrays and vector computations. There is also a version using Tkinter but it not complete.

To install pygame (with pip package manager)

``` pip install pygame```

To install Numpy (with pip package manager)

``` pip install numpy ```

## Simulations

![simu](https://github.com/Rudiio/Images-factory/blob/main/bacteria2.png)