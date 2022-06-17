# The effect of toxicant exposure on infected Flying Foxes
## A mathematical model that explores infection dynamics in toxicant-contaminated landscapes
Joukje Kloosterman & Job Maathuis  
17-06-2022

[[./images/complete_figure.jpg | width=100px]]

The model gives an insight in how the flying fox population size, infection prevalence and zoonotic spillover depend on the contamination rate of the landscape and the effect of toxicant exposure on infection, movement and survival. The model was obtained from https://royalsocietypublishing.org/doi/full/10.1098/rsbl.2020.0559 . The model was written in R and (version 4.1.3) and the DeSolve package (version 1.32) is imported in order to model the differential equations.

### Project description
A report was written, this was divided into two parts. First the research by Sánchez et al. is being reproduced in which the effect of different transmission rates of the virus in the toxicant-contaminated habitat is being modulated. Then the effect of the cost of toxicants to dispersal of the animal is being modulated.

Important parameters from this model are:
- f: This is the fraction of a landscape that is contaminated by toxicants.
- c_σ : This is the cost of toxicants to dispersal of the animal. That is, the ability of an animal of switching from a toxic contaminated habitat to a pristine habitat.
- ß_T and ß_P: These are the transmission rates in a toxic contaminated habitat (ß_T) and the transmission rate in a pristine habitat (ß_P). Transmission in this case is the transferring of a virus from host to host. 

### Run or reproduce project
To reproduce the research, the differential equations have to be put in a function in order to use the ODE method from the DeSolve package. Then different simulations can be run with the parameter values given in figure above. The simulations are run for different f values, ranging from 0.01 to 0.99. The ODE function is run for each f and a timeframe of 50 years. The initial variables need to be calculated for each f aswell, therefor the population size and how many are infected with a virus are needed. They can be calculated as follows:

- SP = (population size - infected) * (1 - f)  
- IP = infected * (1 - f)  
- ST = (population size - infected) * f  
- IT = infected * f  

To plot the population size, infection prevalence and zoonotic spillover, another function can be created or the function from the EndAssignment.Rmd can be used. If a new function is written a few extra equations have to be executed. They calculate the population size, infection prevalence and zoonotic spillover:

- Population size = SP + IP + ST + IT
- Infection prevalance = IT / population size
- Spillover risk = IT / f

Another option is to use the functions written in the EndAssignment.Rmd file, in order to reproduce the research. The first function is the function needed for the DeSolve package. The second function runs the model with the parameter values. The next obtains the data for different c_σ values, the different ß_T values are already in the function. The last function plots the different scenarios, three plots are generated: one for the population size, one for infection prevelance and the last for the spill-over risk.  
To run the simulation, run the third function with a c_σ value and then use the plot function to plot the different scenarios. 

#### Support
Email us if you have any questions:  
j.kloosterman@st.hanze.nl or j.maathuis@st.hanze.nl
