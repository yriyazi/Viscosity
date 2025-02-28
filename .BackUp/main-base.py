"""Glycerol_calculator.py"""
## This python code is based on the MatLab code orginaly provided by Chris Westbrook
## http://www.met.reading.ac.uk/~sws04cdw/viscosity_calc.html
__author__  = "Matthew Partridge"
__license__ = "GPL"
__version__ = "1.0"
__credits__ = "Chris Westbrook"


#Required packages ----------------

import numpy
import math
import matplotlib.pyplot    as plt

#Variables ----------------
waterVol        = 3     #volume of water required (ml)
glycerolVol     = 0.75	#volume of Glycerol used (ml)

plt.figure(figsize=[15,15])
for i in range(1,70,1):
    glycerolVol     = i
    waterVol        = 99 - i

    T               = 21    #temperature (degrees Celcius)
    #Densities ----------------
    glycerolDen = (1273.3-0.6121*T)/1000                #Density of Glycerol    (g/cm3)
    waterDen    = (1-math.pow(((abs(T-4))/622),1.7)) 	#Density of water       (g/cm3)


    #Fraction cacluator ----------------
    glycerolMass    =   glycerolDen*glycerolVol
    waterMass       =   waterDen*waterVol
    totalMass       =   glycerolMass+waterMass
    mass_fraction   =   glycerolMass/totalMass
    vol_fraction    =   glycerolVol/(glycerolVol+waterVol)

    # print ("Mass fraction of mixture =", round(mass_fraction,5))
    # print ("Volume fraction of mixture =", round(vol_fraction,5))


    #Viscosity calcualtor ----------------
    glycerolVisc    = 0.001*12100*numpy.exp((-1233+T)*T/(9900+70*T))
    waterVisc       = 0.001*1.790*numpy.exp((-1230-T)*T/(36100+360*T))

    a       =   0.705-0.0017*T
    b       =   (4.9+0.036*T)*numpy.power(a,2.5)
    alpha   =   1-mass_fraction+(a*b*mass_fraction*(1-mass_fraction))/(a*mass_fraction+b*(1-mass_fraction))
    A       =   numpy.log(waterVisc/glycerolVisc)

    viscosity_mix=glycerolVisc*numpy.exp(A*alpha)

    print (f"Viscosity = {round(viscosity_mix,5):05f} Ns/m2, Mass fraction = {round(mass_fraction,5):05f}")
    plt.scatter(round(mass_fraction,5),round(viscosity_mix,5),color="red")

plt.show()


