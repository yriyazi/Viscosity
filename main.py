import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt
from datetime import datetime

def get_time_stamp():
    return datetime.now().strftime("%m-%d-%H")


# Set global font to Times New Roman
matplotlib.rcParams['font.family'] = "Times New Roman"

# Variables ----------------
num_segments        = 12
viscosity_values    = []
mass_fractions      = []

for i in range(300, 770, 1):
    i           = i/10
    glycerolVol = i
    waterVol    = 99 - i
    T           = 21  # temperature (degrees Celsius)

    # Densities ----------------
    glycerolDen = (1273.3 - 0.6121 * T) / 1000  # Density of Glycerol (g/cm3)
    waterDen = (1 - math.pow(((abs(T - 4)) / 622), 1.7))  # Density of water (g/cm3)

    # Fraction calculator ----------------
    glycerolMass = glycerolDen * glycerolVol
    waterMass = waterDen * waterVol
    totalMass = glycerolMass + waterMass
    mass_fraction = glycerolMass / totalMass

    # Viscosity calculator ----------------
    glycerolVisc    = 0.001 * 12100 * np.exp((-1233 + T) * T / (9900    + 70    * T))
    waterVisc       = 0.001 * 1.790 * np.exp((-1230 - T) * T / (36100   + 360   * T))

    a = 0.705 - 0.0017 * T
    b = (4.9 + 0.036 * T) * np.power(a, 2.5)
    alpha = 1 - mass_fraction + (a * b * mass_fraction * (1 - mass_fraction)) / (a * mass_fraction + b * (1 - mass_fraction))
    A = np.log(waterVisc / glycerolVisc)

    viscosity_mix = glycerolVisc * np.exp(A * alpha)

    viscosity_values.append(viscosity_mix)
    mass_fractions.append(mass_fraction)

# Segmenting viscosity values into 12 even bins
min_viscosity   = min(viscosity_values)
max_viscosity   = max(viscosity_values)
segment_bounds  = np.linspace(min_viscosity, max_viscosity, num_segments + 1)


plt.figure(figsize=[15, 10],dpi=600)
j = 0
for i in range(len(viscosity_values)):
    if segment_bounds[j] <= viscosity_values[i]:
        j += 1
        plt.scatter(mass_fractions[i], viscosity_values[i], color='r')

        # # Adding annotation for each point
        # plt.annotate(f'({mass_fractions[i]:.2f}, {viscosity_values[i]:.4f})',
        #              (mass_fractions[i], viscosity_values[i]),
        #              textcoords="offset points", xytext=(5,5), ha='left', fontsize=8)
        plt.annotate(f'({mass_fractions[i]:.3f}, {viscosity_values[i]:.6f})',
                     (mass_fractions[i], viscosity_values[i]+0.001),
                     textcoords="offset points", xytext=(0, 10), ha='center', va='bottom', fontsize=9,
                     bbox=dict(boxstyle="round,pad=0.3", edgecolor="black", facecolor="yellow"),
                     arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=0", color='black', relpos=(0.5, 0.5)))


plt.xlabel("Mass Fraction of Glycerol")
plt.ylabel(r"Viscosity ($\frac{Ns}{m^{2}}$)")
plt.title(f"Segmented Viscosity Plot \n {get_time_stamp()}")
plt.grid(visible=True, which='major', color='black', linestyle='-')
plt.grid(visible=True, which='minor', color='black', linestyle='--')
plt.ylim(0,0.07)
plt.savefig("viscosity_plot.pdf", format="pdf", dpi=300, bbox_inches='tight')
plt.show()
