# Coupled opinion-environmnetal dynamics in polarized and prejudiced populations

# Abstract

## Reproducibility instructions
The python scripts provided can be used to reproduce all of the figures from the manuscript. The file and function names correspond to the section of the figure they reproduce. In general, each function takes in a given number of parameters corresponding to the constants set in that given part of the figure. The function also takes in a parameter which is set to 'show' if you want to display the plot with captions or 'export' if you want to export the plot in the exact size it is used in the figure. For example to display the bifurcation diagrams from Figure 4a you would run the script 'figure4a-c.py' and then run `gen_fig4ac(q=0, form='show')`. Further details are provided in each of the scripts.
## Required Libraries
 - **numpy**: 2.3.1
 - **matplotlib**: 3.10.3
 - **mpltern**: 1.0.4
## Repo Organization
``` bash
Code
  |-- figure1b.py
  |-- figure2.py
  |-- figure3a-c.py
  |-- figure4a-c.py
  |-- figure4e.py
  |-- figure5.py
  |-- figureD1.py
  |-- figureD2.py
Figures
  |-- Figure1.pdf
  |-- Figure2.pdf
  |-- Figure3.pdf
  |-- Figure4.pdf
  |-- Figure5.pdf
  |-- FigureD.1.pdf
  |-- FigureD.2.pdf
|-- environment.yml
|-- README.md
```
