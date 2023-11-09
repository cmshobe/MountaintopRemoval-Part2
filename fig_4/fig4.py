########################################################################
#This script computes data for, and generates, Figure 4 in the following manuscript:

#Bower, S.J., Shobe, C.M., Maxwell, A.E., and Campforts, B. (2023) The 
#uncertain future of mountaintop-removal-mined landscapes 2: Modeling the influence of
#topography and vegetation. Geomorphology.

#Please cite the paper if you use this code in any way.

#Brief description: this script uses the results of our gully mapping
#(mapping shapefiles can be found in the folder `gully_shapefiles`) to extract the 
#erodibility (K) values we use in our model simulations of mined and unmined landscapes.
########################################################################

import pandas as pd
import numpy as np

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

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatch
import seaborn as sbn
import seaborn as sbn
import matplotlib as mpl
from matplotlib import cm
from mycolorpy import colorlist as mcp

cmapplot = mpl.cm.viridis
norm = mpl.colors.Normalize(vmin=np.min(gullies.E), vmax=np.max(gullies.E))

gullies['color'] = None
gullies.loc[gullies.shed == 'bencreek', 'color'] = '#fde725'
gullies.loc[gullies.shed == 'sprucefork', 'color'] = '#5ec962'
gullies.loc[gullies.shed == 'laurelcreek', 'color'] = '#21918c'
gullies.loc[gullies.shed == 'whiteoak', 'color'] = '#3b528b'
gullies.loc[gullies.shed == 'mudriver', 'color'] = '#440154'

####################make figure

fig2 = plt.figure(figsize = (7.5,5), tight_layout = False)

spec = fig2.add_gridspec(5, 10, wspace = 0.0, hspace= 0.0)

ax1 = fig2.add_subplot(spec[2:5,:4])
ax2 = fig2.add_subplot(spec[1:2,:4])
ax3 = fig2.add_subplot(spec[2:,4:5])
ax4 = fig2.add_subplot(spec[1:,5:])

ax1.scatter(gullies.logslope,gullies.logDA,c = gullies.E, edgecolor = 'k', alpha = 0.75, s = 50,cmap = 'viridis')
ax1.text(-.9,5.25,r"$\rho$ = "+str(np.round(slope_DA_spearmanr[0],2))+"\np = "+str(np.round(slope_DA_spearmanr[1],2)),fontsize = 12)
ax1.tick_params(axis='both', which='major', labelsize=12)
ax1.set_xlabel(r'$log_{10}slope$', fontsize = 14)
ax1.set_ylabel(r'$log_{10}DA$', fontsize = 14)
ax1.scatter(gullies.logslope[gullies.DA_filter == 1],
            gullies.logDA[gullies.DA_filter == 1],
            edgecolor = 'red',
            s = 50, facecolor = 'none')

sbn.kdeplot(gullies.logslope,ax = ax2, c = 'k', lw = 2)
ax2.set_yticks([])
ax2.set_xticks([])
ax2.set_xlabel(None)
ax2.set_ylabel(None)
ax2.patch.set_alpha(0)
ax2.spines['left'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['bottom'].set_visible(False)

sbn.kdeplot(data = gullies,ax = ax3,y = 'logDA', c = 'k', lw = 2)
ax3.set_yticks([])
ax3.set_xticks([])
ax3.set_xlabel(None)
ax3.set_ylabel(None)
ax3.patch.set_alpha(0)
ax3.spines['left'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax3.spines['top'].set_visible(False)
ax3.spines['bottom'].set_visible(False)


sbn.kdeplot(data = gullies_cleaned['K'], color = 'grey', lw = 3, fill = True, cut = 0, ax = ax4)
ax4.ticklabel_format(axis = 'x',style = 'sci',scilimits = (0,0))
ax4.yaxis.tick_right()
ax4.yaxis.set_label_position("right")
ax4.vlines(x = Kmax, ymin = 0, ymax = 90, ls = '--', lw = 2, color = 'k')
ax4.vlines(x = Kmin, ymin = 0, ymax = 75, ls = '--', lw = 2, color = 'k')
ax4.set_yticks([])
ax4.set_ylabel(None)
ax4.set_xlabel('$K$ ($yr^{-1}$)', fontsize = 12)
ax4.spines['right'].set_visible(False)
ax4.spines['top'].set_visible(False)
ax4.spines['left'].set_visible(False)
x_ticks = [Kmin,Kmax]
x_labels = [r'$1.3*10^{-4}$',r'$1.4*10^{-2}$']
ax4.set_xticks(x_ticks)
ax4.set_xticklabels(x_labels)
ax4.tick_params(axis='both', which='major', labelsize=12)
ax4.set_xscale('log')
ax4.text(10**-4,78,'$K_{min}$',fontsize = 14)
ax4.text(.0035,90,'$K_{max}$',fontsize = 14)

cbar_ax = ax4.inset_axes([0.03, 15, .01, 50], transform = ax4.transData)
cb1 = mpl.colorbar.ColorbarBase(cbar_ax, cmap=cmapplot,
                                norm=norm,
                                orientation='vertical')
cb1.set_label(r'incision rate $myr^{-1}$')

fig2.savefig('fig4.png', dpi = 1000)