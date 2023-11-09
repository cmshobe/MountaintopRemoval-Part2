########################################################################
#This script generates Figure 13 in the following manuscript:

#Bower, S.J., Shobe, C.M., Maxwell, A.E., and Campforts, B. (2023) The 
#uncertain future of mountaintop-removal-mined landscapes 2: Modeling the influence of
#topography and vegetation. Geomorphology.

#Please cite the paper if you use this code in any way.

#Brief description: this script plots transects taken across different elevation surfaces
#(pre-mining, post-mining, and post-10,000 years fo erosion for each model scenario).
#Transect locations are shown in Figure 12.

########################################################################

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from mycolorpy import colorlist as mcp

a = pd.read_csv('transects/transect_a.csv', delimiter = ',')
b = pd.read_csv('transects/transect_b.csv', delimiter = ',')
c = pd.read_csv('transects/transect_c.csv', delimiter = ',')
d = pd.read_csv('transects/transect_d.csv', delimiter = ',')

cmap = mcp.gen_color(cmap="plasma",n=11)
cmapplot = mpl.cm.plasma
norm = mpl.colors.Normalize(vmin=0, vmax=1)

fig, (ax1,ax2,ax3,ax4) = plt.subplots(nrows = 4, figsize = (6,8),tight_layout = True)

ax1.plot(a['10'],a['Unnamed: 9'], color = cmap[8],label = '100 %')
ax1.plot(a['5'],a['Unnamed: 11'], color = cmap[5],label = '50 %')
ax1.plot(a['0'],a['Unnamed: 13'], color = cmap[0],label = '0 %')
ax1.plot(a['pre'],a['Unnamed: 1'], color = 'grey',label = 'unmined IC')
ax1.plot(a['post'],a['Unnamed: 3'], color = 'k',label = 'mined IC')
ax1.plot(a['unmined control'],a['Unnamed: 5'], color = 'grey', ls = '-.',label = 'control\nunmined')
ax1.plot(a['mined control'],a['Unnamed: 7'], color = 'k', ls = '--',label = 'control\nmined')

ax2.plot(b['10'],b['Unnamed: 9'], color = cmap[8])
ax2.plot(b['5'],b['Unnamed: 11'], color = cmap[5])
ax2.plot(b['0'],b['Unnamed: 13'], color = cmap[0])
ax2.plot(b['pre'],b['Unnamed: 1'], color = 'grey')
ax2.plot(b['post'],b['Unnamed: 3'], color = 'k')
ax2.plot(b['unmined control'],b['Unnamed: 5'], color = 'grey', ls = '-.')
ax2.plot(b['mined control'],b['Unnamed: 7'], color = 'k', ls = '--')


ax3.plot(c['10'],c['Unnamed: 9'], color = cmap[8])
ax3.plot(c['5'],c['Unnamed: 11'], color = cmap[5])
ax3.plot(c['0'],c['Unnamed: 13'], color = cmap[0])
ax3.plot(c['pre'],c['Unnamed: 1'], color = 'grey')
ax3.plot(c['post'],c['Unnamed: 3'], color = 'k')
ax3.plot(c['unmined control'],c['Unnamed: 5'], color = 'grey', ls = '-.')
ax3.plot(c['mined control'],c['Unnamed: 7'], color = 'k', ls = '--')

ax4.plot(d['10'][:-7],d['Unnamed: 9'][5:-2], color = cmap[8])
ax4.plot(d['5'][:-7],d['Unnamed: 11'][5:-2], color = cmap[5])
ax4.plot(d['0'][:-7],d['Unnamed: 13'][5:-2], color = cmap[0])
ax4.plot(d['pre'][:-7],d['Unnamed: 1'][5:-2], color = 'grey')
ax4.plot(d['post'][:-7],d['Unnamed: 3'][5:-2], color = 'k')
ax4.plot(d['unmined control'][:-7],d['Unnamed: 5'][5:-2], color = 'grey', ls = '-.')
ax4.plot(d['mined control'][:-7],d['Unnamed: 7'][5:-2], color = 'k', ls = '--')

ax4.set_xlabel('distance (m)', fontsize = 12)
ax1.set_ylabel('elevation (m)', fontsize = 12)
ax2.set_ylabel('elevation (m)', fontsize = 12)
ax3.set_ylabel('elevation (m)', fontsize = 12)
ax4.set_ylabel('elevation (m)', fontsize = 12)

ax1.legend(loc = 'upper center')
ax1.legend(bbox_to_anchor=(0.062, 1.05),
          ncol=4, fancybox=True, shadow=True)

fig.savefig('fig13.png', dpi = 1000)