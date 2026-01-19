# Coupled opinion-environmnetal dynamics in polarized and prejudiced populations

# Abstract
Public opinion on environmental issues remains polarized in many countries, posing a significant barrier to the implementation of effective policies. Behind this polarization, empirical studies have identified social susceptibility, personal prejudice, and personal experience as dominant factors in opinion formation on environmental issues. However, current coupled human-environment models have not yet incorporated all three factors in polarized populations. We developed a stylized coupled human-environment model to investigate how social susceptibility, personal prejudice, and personal experience shape opinion formation and the environment in polarized populations. Using analytical and numerical methods, we characterized the conditions under which polarization, consensus, opinion changes, and cyclic dynamics emerge depending on the costs of mitigation, environmental damage, and the factors influencing opinion formation. Our results demonstrate that the ratio between the strength of prejudice and objectivity, rather than social susceptibility, determines long-term opinion-environmental dynamics. Additionally, we find that prejudice is the key driver of persistent polarization, with even slightly prejudiced populations maintaining indefinite polarization. We predict that polarization can be reduced by decreasing the role of prejudice or increasing the willingness to consider opposing opinions. Finally, our model shows that cost reduction methods are less effective at reducing environmental impact in populations that are more prejudiced and less objective. Our model generates thresholds for when reducing costs or emissions is more useful depending on the factors which influence the population’s opinion formation. Overall, our model provides a framework for investigating the importance of cognitive and social structures in determining human-environment dynamics.

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
