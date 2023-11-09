########################################################################
#This script computes data for, and generates, Figure 5 in the following manuscript:

#Bower, S.J., Shobe, C.M., Maxwell, A.E., and Campforts, B. (2023) The 
#uncertain future of mountaintop-removal-mined landscapes 2: Modeling the influence of
#topography and vegetation. Geomorphology.

#Please cite the paper if you use this code in any way.

#Brief description: this script uses the erodibility values derived from gully mapping
#(see `fig4.py`) and shifts those values such that the minimum erodibility is that
#established for the region in prior work (Gallen, 2018) while the ratio of maximum to 
#minimum erodibility we observed is preserved. See the manuscript for further details.

#This script creates the time series of K values used in our numerical simulations

########################################################################

import pandas as pd
import numpy as np
import matplotlib
from matplotlib import cm
from mycolorpy import colorlist as mcp

# Load gullies; calculate K
gullies = pd.read_csv("final_gully_data.csv", delimiter = ",")
gullies['age'] = 2021 - gullies['min_year'] #measurements were made in 2021
gullies['E'] = None
gullies['E'] = gullies['avg_depth']/gullies['age']
m = 0.5 #stream power area exponent
n = 1 #stream power slope exponent
gullies['K'] = gullies['E'] / ( (gullies['DA'] ** m) * (gullies['slope'] ** n) )


# Remove data with bad drainage area measurements (i.e., that don't fall on a mapped flow path)
DA_threshold = 100
DA_filter = np.where(gullies.DA <= DA_threshold, 1,0)
gullies['DA_filter'] = DA_filter
gullies_cleaned = gullies.loc[gullies.DA_filter == 0]
gullies['logDA'] = np.log10(gullies.DA)
gullies['logslope'] = np.log10(gullies.slope)

#calculate rank statistic
import scipy.stats as sp
slope_DA_spearmanr = sp.spearmanr(gullies_cleaned.DA,gullies_cleaned.slope)

# Calculate Kmin and Kmax as the absolute minimum of the cleaned data, and the 96th perentile respectively
K_all = gullies_cleaned['K'].to_numpy()
K_all = K_all[~np.isnan(K_all)]
Kmax = np.percentile(K_all, 50) #we assume the median K value represents "Kmax" (anything greater is an outlier)
Kmin = np.min(K_all)

#calculate the ratio of max to min measured erodibility
k_prop = Kmax / Kmin
Kmin = 1.3e-6 #Kmin is reassigned to be the value from Gallen 2018
Kmax = Kmin * k_prop #Kmax is calculated to preserve ratio Kmax/Kmin

per_rec_scale = np.arange(0,1.1,0.1)
Kmin_star = np.zeros(len(per_rec_scale))
Kmin_star = Kmin_star /100
for i in range(0,len(per_rec_scale)):
    Kmin_star[i] = Kmax -((Kmax - Kmin) * per_rec_scale[i])

t = np.arange(0,200,1)
K = np.zeros(len(t))

K_scenarios = np.zeros((len(t),len(per_rec_scale)))

for i in range(0,len(per_rec_scale)):
    m = (Kmax-Kmin_star[i])/(200**0.25)
    K[:] = Kmax - (m *t ** 0.25)
    K_scenarios[:,i] = K[:]

#############create and export the time series of K for use in numerical simulations
K_ts = np.zeros((10000,11))

for i in range(0,len(per_rec_scale)):
    K_ts[:200,i] = K_scenarios[:,i]
    K_ts[200:,i] = np.min(K_scenarios[:,i])
    
np.savetxt("K_timeseries_scaled.txt", K_ts)

#############create figure showing time series of K
import matplotlib
from matplotlib import cm
from mycolorpy import colorlist as mcp

cmap = mcp.gen_color(cmap="plasma",n=11)

fig5 = plt.figure(figsize = (7.5,4), tight_layout = True)

spec2 = fig5.add_gridspec(1, 2, wspace = 0.1, hspace= 0.0)

axl = fig5.add_subplot(spec2[:,0], zorder = 2)
axr = fig5.add_subplot(spec2[:,1])


axl.set_xlim(-10, 110)
axl.hlines(Kmax*1e5, -10, 110, ls = '--', color = 'k',clip_on=False)
axl.hlines(Kmin*1e5, -10, 110, ls = '--', color = 'k',clip_on=False)

axl.scatter(per_rec_scale[0]*100, Kmin_star[0]*1e5,zorder = 2, color = cmap[0], alpha = .75, edgecolor = 'k', s = 225, label = '0 %')
axl.scatter(per_rec_scale[5]*100, Kmin_star[5]*1e5,zorder = 2, color = cmap[5], alpha = .75, edgecolor = 'k', s = 225, label = '50 %')
axl.scatter(per_rec_scale[10]*100, Kmin_star[10]*1e5,zorder = 2, color = cmap[8], alpha = .75, edgecolor = 'k', s = 225, label = '100 %')

axr.plot(t,K_scenarios[:,0], color = cmap[0], lw = 7, alpha = .75, zorder = 4, label = '0 %')
axr.plot(t,K_scenarios[:,5], color = cmap[5], lw = 7, alpha = .75, zorder = 2, label = '50 %')
axr.plot(t,K_scenarios[:,10], color = cmap[8], lw = 7, alpha = .75, zorder = 0, label = '100 %')
       
axr.set_yticklabels([])
axr.set_ylabel(None)
axl.set_ylabel(r'$K$ ($yr^{-1}$) $\times 10^{-5}$', fontsize = 14)
axl.ticklabel_format(axis = 'y',style = 'sci',scilimits = (0,0))
axr.set_xlabel('time (yr)', fontsize = 14)
axl.set_xlabel('percent recovery', fontsize = 14)
axr.set_xlim(-10,200)
axr.hlines(Kmax, -10, 225, ls = '--', color = 'k')
axr.hlines(Kmin, -10, 225, ls = '--', color = 'k')

axl.text(0,3e-1,'$K_{min}$', fontsize = 16)
axl.text(85,3.10,'$K_{max}$', fontsize = 16)

legend = axl.legend(loc = 'center right', title = 'percent\nrecovery',markerscale = 0.5)
frame = legend.get_frame()
frame.set_edgecolor('k')

fig5.savefig('fig5.png', dpi=1000)
