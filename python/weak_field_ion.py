import numpy as np
import matplotlib.pyplot as plt

#define the ions of interest
fe13a=np.array([1.5,10746.8,6.22,55.845*1.672621E-27,1])
fe13b=np.array([1.5,10746.8,6.22,55.845*1.672621E-27,1])
si10=np.array([1.5,10746.8,6.22,55.845*1.672621E-27,1])
si09=np.array([1.5,10746.8,6.22,55.845*1.672621E-27,1])

# define variables of common units need
 
l_speed=2.9979E+8             ## speed of light [m s^-1]  
kb=1.3806488E-23        ## Boltzman constant [m^2 kg s^-2 K^-1] 
e_mass=9.10938356e-31         ## Electron mass SI [Kg]
e_charge=1.602176634e-19      ## Electron charge SI [C]

# Ion specific parameters
lande_g=1.5					  ## effective lande g factor.							
line_w=10746.8*1e-10          ## Ion referential wavelength [A] ;converted to [m]
ion_temp=6.22			  	  ## log Ion temperature SI [K]
ion_mass=55.845*1.672621E-27  ## Ion mass SI [Kg] ## for ion XIII this needs to be computed for all 
xi=1000                       ## microturbulent velocity SI [m s^-1]

def weakfield(lande_g,line_w,ion_temp,ion_mass,xi):
    wf=4*3.141592654 * e_mass*l_speed/(lande_g*line_w*e_charge)*np.sqrt((2*kb*(10**ion_temp)/ion_mass) + xi**2)
    print('The field should be under: %2.3e Gauss' %wf)
