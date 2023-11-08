########################################################################
#This script generates Figure 6 in the following manuscript:

#Bower, S.J., Shobe, C.M., Maxwell, A.E., and Campforts, B. (2023) The 
#uncertain future of mountaintop-removal-mined landscapes 2: Modeling the influence of
#topography and vegetation. Geomorphology.

#Please cite the paper if you use this code in any way.

#Brief description: this script uses plots mean annual precipitation from 2010-2100
#derived from NASA BioClim data:

#https://www.nccs.nasa.gov/services/data-collections/land-based-products/bioclim

#Pearson, R.G., Stanton. J.C., Shoemaker, K.T., Aiello-Lammens, M.E., Ersts, P.J., 
#Horning, N., Fordham, D.A., Raxworthy, C.J., Ryu, H.Y., McNees, J., & Akçakaya, 
#H.R. Life history and spatial traits predict extinction risk due to climate change. 
#Nature Climate Change 4:217-221.

#We used BioClim's "Policy" scenario, which assumes CO2 stabilization at 450 ppm.
#We extracted the mean MAP for each year over each of our five study watersheds
#using raster zonal statistics.

########################################################################

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import mplcursors

df = pd.read_csv('climate_data.csv', delimiter = ',').to_numpy()

SF = df[0,2:]
WO = df[1,2:]
LC = df[2,2:]
BC = df[3,2:]
MR = df[4,2:]

fig, ax = plt.subplots(figsize = (8,4))

time =np.arange(0,91,1)
model_domain_time = np.arange(0,10000,1)
model_domain_precip = np.zeros(len(model_domain_time))
model_domain_precip[0:91] = SF
model_domain_precip[91:] = np.max(SF)

import matplotlib
from matplotlib import cm
from mycolorpy import colorlist as mcp

cmap = mcp.gen_color(cmap="viridis",n=5)

ax.plot(time,SF, label = 'Spruce Fork', color = cmap[0], lw = 4.5)#, ls = '-.')
ax.plot(time,MR, label = 'Mud River', color = cmap[1], lw = 4.5)#, ls = '-.')
ax.plot(time,LC, label = 'Laurel Creek', color = cmap[2], lw = 4.5)#, ls = '-.')
ax.plot(time,WO, label = 'White Oak', color = cmap[3], lw = 4.5)#, ls = '-.')
ax.plot(time,BC, label = 'Ben Creek', color = cmap[4], lw = 4.5)#, ls = '-.')
ax.set_xlim(0,90)
ax.tick_params(axis='both', which='major', labelsize=14)
ax.set_ylabel('precipitation (mm/yr)', fontsize = 14)
ax.set_xlabel('time (yr)', fontsize = 14)
ax.legend(title = 'Watershed')

plt.tight_layout()
mplcursors.cursor()

fig.savefig('fig6.png', dpi = 1000)