########################################################################
#This script generates Figure 10 in the following manuscript:

#Bower, S.J., Shobe, C.M., Maxwell, A.E., and Campforts, B. (2023) The 
#uncertain future of mountaintop-removal-mined landscapes 2: Modeling the influence of
#topography and vegetation. Geomorphology.

#Please cite the paper if you use this code in any way.

#Brief description: this script uses the results of our numerical simulations
#(found in the folder `model_simulations`) to plot the distribution of run-averaged
#erosion rates for each model scenario in each watershed.

########################################################################

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats
from scipy.stats import norm
import landlab
from landlab.plot import imshow_grid
from landlab import RasterModelGrid
from landlab.io.esri_ascii import read_esri_ascii
import seaborn as sbn
import matplotlib as mpl
from matplotlib import cm
from mycolorpy import colorlist as mcp

################import results of model runs: post-run topography for each scenario and watershed

#ben creek
mg_pre_control, bc_z_pre_control = read_esri_ascii('../model_simulations/bencreek/bencreek_pre_10m.asc', name='topographic__elevation')
np.all(mg_pre_control.at_node['topographic__elevation'] == bc_z_pre_control)
bc_z_pre_control[bc_z_pre_control < 1] = -9999

mg_mod, bc_z_mod_control = read_esri_ascii('../model_simulations/bencreek/control_unmined_elev_01.asc', name='topographic__elevation')
np.all(mg_mod.at_node['topographic__elevation'] == bc_z_mod_control)
bc_z_mod_control[bc_z_mod_control < 1] = -9999
bc_diff_control = bc_z_mod_control - bc_z_pre_control
bc_diff_control = bc_diff_control[bc_diff_control > -9999]

mg_post, bc_z_post = read_esri_ascii('../model_simulations/bencreek/bencreek_post_10m.asc', name='topographic__elevation')
np.all(mg_post.at_node['topographic__elevation'] == bc_z_post)
bc_z_post[bc_z_post < 1] = -9999

mg_mod, bc_z_mod_0 = read_esri_ascii('../model_simulations/bencreek/scenario_0_elev_01.asc', name='topographic__elevation')
np.all(mg_mod.at_node['topographic__elevation'] == bc_z_mod_0)
bc_z_mod_0[bc_z_mod_0 < 1] = -9999
bc_diff_0 = bc_z_mod_0 - bc_z_post
bc_diff_0 = bc_diff_0[bc_diff_0 > -9999]


mg_mod, bc_z_mod_5 = read_esri_ascii('../model_simulations/bencreek/scenario_5_elev_01.asc', name='topographic__elevation')
np.all(mg_mod.at_node['topographic__elevation'] == bc_z_mod_5)
bc_z_mod_5[bc_z_mod_5 < 1] = -9999
bc_diff_5 = bc_z_mod_5 - bc_z_post
bc_diff_5 = bc_diff_5[bc_diff_5 > -9999]


mg_mod, bc_z_mod_10 = read_esri_ascii('../model_simulations/bencreek/scenario_10_elev_01.asc', name='topographic__elevation')
np.all(mg_mod.at_node['topographic__elevation'] == bc_z_mod_10)
bc_z_mod_10[bc_z_mod_10 < 1] = -9999
bc_diff_10 = bc_z_mod_10 - bc_z_post
bc_diff_10 = bc_diff_10[bc_diff_10 > -9999]

mg_mod, bc_z_mod_control_mined = read_esri_ascii('../model_simulations/bencreek/control_mined_elev_01.asc', name='topographic__elevation')
np.all(mg_mod.at_node['topographic__elevation'] == bc_z_mod_control_mined)
bc_z_mod_control_mined[bc_z_mod_control_mined < 1] = -9999
bc_diff_control_mined = bc_z_mod_control_mined - bc_z_post
bc_diff_control_mined = bc_diff_control_mined[bc_diff_control_mined > -9999]


#laurel creek
mg_pre_control, lc_z_pre_control = read_esri_ascii('../model_simulations/laurelcreek/laurelcreek_pre_10m.asc', name='topographic__elevation')
np.all(mg_pre_control.at_node['topographic__elevation'] == lc_z_pre_control)
lc_z_pre_control[lc_z_pre_control < 1] = -9999

mg_mod, lc_z_mod_control = read_esri_ascii('../model_simulations/laurelcreek/control_unmined_elev_01.asc', name='topographic__elevation')
np.all(mg_mod.at_node['topographic__elevation'] == lc_z_mod_control)
lc_z_mod_control[lc_z_mod_control < 1] = -9999
lc_diff_control = lc_z_mod_control - lc_z_pre_control
lc_diff_control = lc_diff_control[lc_diff_control > -9999]

mg_post, lc_z_post = read_esri_ascii('../model_simulations/laurelcreek/laurelcreek_post_10m.asc', name='topographic__elevation')
np.all(mg_post.at_node['topographic__elevation'] == lc_z_post)
lc_z_post[lc_z_post < 1] = -9999

mg_mod, lc_z_mod_0 = read_esri_ascii('../model_simulations/laurelcreek/scenario_0_elev_01.asc', name='topographic__elevation')
np.all(mg_mod.at_node['topographic__elevation'] == lc_z_mod_0)
lc_z_mod_0[lc_z_mod_0 < 1] = -9999
lc_diff_0 = lc_z_mod_0 - lc_z_post
lc_diff_0 = lc_diff_0[lc_diff_0 > -9999]

mg_mod, lc_z_mod_5 = read_esri_ascii('../model_simulations/laurelcreek/scenario_5_elev_01.asc', name='topographic__elevation')
np.all(mg_mod.at_node['topographic__elevation'] == lc_z_mod_5)
lc_z_mod_5[lc_z_mod_5 < 1] = -9999
lc_diff_5 = lc_z_mod_5 - lc_z_post
lc_diff_5 = lc_diff_5[lc_diff_5 > -9999]

mg_mod, lc_z_mod_10 = read_esri_ascii('../model_simulations/laurelcreek/scenario_10_elev_01.asc', name='topographic__elevation')
np.all(mg_mod.at_node['topographic__elevation'] == lc_z_mod_10)
lc_z_mod_10[lc_z_mod_10 < 1] = -9999
lc_diff_10 = lc_z_mod_10 - lc_z_post
lc_diff_10 = lc_diff_10[lc_diff_10 > -9999]

mg_mod, lc_z_mod_control_mined = read_esri_ascii('../model_simulations/laurelcreek/control_mined_elev_01.asc', name='topographic__elevation')
np.all(mg_mod.at_node['topographic__elevation'] == lc_z_mod_control_mined)
lc_z_mod_control_mined[lc_z_mod_control_mined < 1] = -9999
lc_diff_control_mined = lc_z_mod_control_mined - lc_z_post
lc_diff_control_mined = lc_diff_control_mined[lc_diff_control_mined > -9999]


#mud river
mg_pre_control, mr_z_pre_control = read_esri_ascii('../model_simulations/mudriver/mudriver_pre_10m.asc', name='topographic__elevation')
np.all(mg_pre_control.at_node['topographic__elevation'] == mr_z_pre_control)
mr_z_pre_control[mr_z_pre_control < 1] = -9999

mg_mod, mr_z_mod_control = read_esri_ascii('../model_simulations/mudriver/control_unmined_elev_01.asc', name='topographic__elevation')
np.all(mg_mod.at_node['topographic__elevation'] == mr_z_mod_control)
mr_z_mod_control[mr_z_mod_control < 1] = -9999
mr_diff_control = mr_z_mod_control - mr_z_pre_control
mr_diff_control = mr_diff_control[mr_diff_control > -9999]

mg_post, mr_z_post = read_esri_ascii('../model_simulations/mudriver/mudriver_post_10m.asc', name='topographic__elevation')
np.all(mg_post.at_node['topographic__elevation'] == mr_z_post)
mr_z_post[mr_z_post < 1] = -9999

mg_mod, mr_z_mod_0 = read_esri_ascii('../model_simulations/mudriver/scenario_0_elev_01.asc', name='topographic__elevation')
np.all(mg_mod.at_node['topographic__elevation'] == mr_z_mod_0)
mr_z_mod_0[mr_z_mod_0 < 1] = -9999
mr_diff_0 = mr_z_mod_0 - mr_z_post
mr_diff_0 = mr_diff_0[mr_diff_0 > -9999]

mg_mod, mr_z_mod_5 = read_esri_ascii('../model_simulations/mudriver/scenario_5_elev_01.asc', name='topographic__elevation')
np.all(mg_mod.at_node['topographic__elevation'] == mr_z_mod_5)
mr_z_mod_5[mr_z_mod_5 < 1] = -9999
mr_diff_5 = mr_z_mod_5 - mr_z_post
mr_diff_5 = mr_diff_5[mr_diff_5 > -9999]

mg_mod, mr_z_mod_10 = read_esri_ascii('../model_simulations/mudriver/scenario_10_elev_01.asc', name='topographic__elevation')
np.all(mg_mod.at_node['topographic__elevation'] == mr_z_mod_10)
mr_z_mod_10[mr_z_mod_10 < 1] = -9999
mr_diff_10 = mr_z_mod_10 - mr_z_post
mr_diff_10 = mr_diff_10[mr_diff_10 > -9999]

mg_mod, mr_z_mod_control_mined = read_esri_ascii('../model_simulations/mudriver/control_mined_elev_01.asc', name='topographic__elevation')
np.all(mg_mod.at_node['topographic__elevation'] == mr_z_mod_control_mined)
mr_z_mod_control_mined[mr_z_mod_control_mined < 1] = -9999
mr_diff_control_mined = mr_z_mod_control_mined - mr_z_post
mr_diff_control_mined = mr_diff_control_mined[mr_diff_control_mined > -9999]


#spruce fork
mg_pre_control, sf_z_pre_control = read_esri_ascii('../model_simulations/sprucefork/sprucefork_pre_10m.asc', name='topographic__elevation')
np.all(mg_pre_control.at_node['topographic__elevation'] == sf_z_pre_control)
sf_z_pre_control[sf_z_pre_control < 1] = -9999

mg_mod, sf_z_mod_control = read_esri_ascii('../model_simulations/sprucefork/control_unmined_elev_01.asc', name='topographic__elevation')
np.all(mg_mod.at_node['topographic__elevation'] == sf_z_mod_control)
sf_z_mod_control[sf_z_mod_control < 1] = -9999
sf_diff_control = sf_z_mod_control - sf_z_pre_control
sf_diff_control = sf_diff_control[sf_diff_control > -9999]

mg_post, sf_z_post = read_esri_ascii('../model_simulations/sprucefork/sprucefork_post_10m.asc', name='topographic__elevation')
np.all(mg_post.at_node['topographic__elevation'] == sf_z_post)
sf_z_post[sf_z_post < 1] = -9999

mg_mod, sf_z_mod_0 = read_esri_ascii('../model_simulations/sprucefork/scenario_0_elev_01.asc', name='topographic__elevation')
np.all(mg_mod.at_node['topographic__elevation'] == sf_z_mod_0)
sf_z_mod_0[sf_z_mod_0 < 1] = -9999
sf_diff_0 = sf_z_mod_0 - sf_z_post
sf_diff_0 = sf_diff_0[sf_diff_0 > -9999]

mg_mod, sf_z_mod_5 = read_esri_ascii('../model_simulations/sprucefork/scenario_5_elev_01.asc', name='topographic__elevation')
np.all(mg_mod.at_node['topographic__elevation'] == sf_z_mod_5)
sf_z_mod_5[sf_z_mod_5 < 1] = -9999
sf_diff_5 = sf_z_mod_5 - sf_z_post
sf_diff_5 = sf_diff_5[sf_diff_5 > -9999]

mg_mod, sf_z_mod_10 = read_esri_ascii('../model_simulations/sprucefork/scenario_10_elev_01.asc', name='topographic__elevation')
np.all(mg_mod.at_node['topographic__elevation'] == sf_z_mod_10)
sf_z_mod_10[sf_z_mod_10 < 1] = -9999
sf_diff_10 = sf_z_mod_10 - sf_z_post
sf_diff_10 = sf_diff_10[sf_diff_10 > -9999]

mg_mod, sf_z_mod_control_mined = read_esri_ascii('../model_simulations/sprucefork/control_mined_elev_01.asc', name='topographic__elevation')
np.all(mg_mod.at_node['topographic__elevation'] == sf_z_mod_control_mined)
sf_z_mod_control_mined[sf_z_mod_control_mined < 1] = -9999
sf_diff_control_mined = sf_z_mod_control_mined - sf_z_post
sf_diff_control_mined = sf_diff_control_mined[sf_diff_control_mined > -9999]


#white oak creek
mg_pre_control, wo_z_pre_control = read_esri_ascii('../model_simulations/whiteoak/whiteoak_pre_10m.asc', name='topographic__elevation')
np.all(mg_pre_control.at_node['topographic__elevation'] == wo_z_pre_control)
wo_z_pre_control[wo_z_pre_control < 1] = -9999

mg_mod, wo_z_mod_control = read_esri_ascii('../model_simulations/whiteoak/control_unmined_elev_01.asc', name='topographic__elevation')
np.all(mg_mod.at_node['topographic__elevation'] == wo_z_mod_control)
wo_z_mod_control[wo_z_mod_control < 1] = -9999
wo_diff_control = wo_z_mod_control - wo_z_pre_control
wo_diff_control = wo_diff_control[wo_diff_control > -9999]

mg_post, wo_z_post = read_esri_ascii('../model_simulations/whiteoak/whiteoak_post_10m.asc', name='topographic__elevation')
np.all(mg_post.at_node['topographic__elevation'] == wo_z_post)
wo_z_post[wo_z_post < 1] = -9999

mg_mod, wo_z_mod_0 = read_esri_ascii('../model_simulations/whiteoak/scenario_0_elev_01.asc', name='topographic__elevation')
np.all(mg_mod.at_node['topographic__elevation'] == wo_z_mod_0)
wo_z_mod_0[wo_z_mod_0 < 1] = -9999
wo_diff_0 = wo_z_mod_0 - wo_z_post
wo_diff_0 = wo_diff_0[wo_diff_0 > -9999]

mg_mod, wo_z_mod_5 = read_esri_ascii('../model_simulations/whiteoak/scenario_5_elev_01.asc', name='topographic__elevation')
np.all(mg_mod.at_node['topographic__elevation'] == wo_z_mod_5)
wo_z_mod_5[wo_z_mod_5 < 1] = -9999
wo_diff_5 = wo_z_mod_5 - wo_z_post
wo_diff_5 = wo_diff_5[wo_diff_5 > -9999]

mg_mod, wo_z_mod_10 = read_esri_ascii('../model_simulations/whiteoak/scenario_10_elev_01.asc', name='topographic__elevation')
np.all(mg_mod.at_node['topographic__elevation'] == wo_z_mod_10)
wo_z_mod_10[wo_z_mod_10 < 1] = -9999
wo_diff_10 = wo_z_mod_10 - wo_z_post
wo_diff_10 = wo_diff_10[wo_diff_10 > -9999]

mg_mod, wo_z_mod_control_mined = read_esri_ascii('../model_simulations/whiteoak/control_mined_elev_01.asc', name='topographic__elevation')
np.all(mg_mod.at_node['topographic__elevation'] == wo_z_mod_control_mined)
wo_z_mod_control_mined[wo_z_mod_control_mined < 1] = -9999
wo_diff_control_mined = wo_z_mod_control_mined - wo_z_post
wo_diff_control_mined = wo_diff_control_mined[wo_diff_control_mined > -9999]



##########################produce figure 10
cmap = mcp.gen_color(cmap="plasma",n=11)
cmapplot = mpl.cm.plasma
norm = mpl.colors.Normalize(vmin=0, vmax=1)

fig = plt.figure(figsize = (8,6),
                  tight_layout = False)

spec = fig.add_gridspec(10,
                         10,
                         wspace = 0.25,
                         hspace= 0.1)


fig, (ax,ax2,ax4,ax6,ax8) = plt.subplots(nrows = 5, figsize = (6,8), tight_layout = True)

plot_range = (-0.5, 3.5) 

ax.hist(mr_diff_0*-1 /10, range = plot_range, bins = 50, color = cmap[0],histtype = u'step',label = '0 %')
ax.hist(mr_diff_5*-1 /10, range = plot_range, bins = 50, color = cmap[5],histtype = u'step',label = '50 %')
ax.hist(mr_diff_10*-1 /10, range = plot_range, bins = 50, color = cmap[8],histtype = u'step',label = '100 %')
ax.hist(mr_diff_control*-1 /10, range = plot_range, bins = 50,histtype = u'step', ls = '-.',color='grey', label = 'control\nunmined')
ax.hist(mr_diff_control_mined*-1 /10, range = plot_range, bins = 50,histtype = u'step', ls = '--',color='k', label = 'control\nmined')
ax.set_yticklabels([])
ax.set_xticklabels([])

ax2.hist(wo_diff_0*-1 /10, range = plot_range, bins = 50, color = cmap[0],histtype = u'step')
ax2.hist(wo_diff_5*-1 /10, range = plot_range, bins = 50, color = cmap[5],histtype = u'step')
ax2.hist(wo_diff_10*-1 /10, range = plot_range, bins = 50, color = cmap[8],histtype = u'step')
ax2.hist(wo_diff_control*-1 /10, range = plot_range, bins = 50,histtype = u'step', ls = '-.',color='grey')
ax2.hist(wo_diff_control_mined*-1 /10, range = plot_range, bins = 50,histtype = u'step', ls = '--',color='k')
ax2.set_xticklabels([])

ax4.hist(bc_diff_0*-1 /10, range = plot_range, bins = 50, color = cmap[0],histtype = u'step')
ax4.hist(bc_diff_5*-1 /10, range = plot_range, bins = 50, color = cmap[5],histtype = u'step')
ax4.hist(bc_diff_10*-1 /10, range = plot_range, bins = 50, color = cmap[8],histtype = u'step')
ax4.hist(bc_diff_control*-1 /10, range = plot_range, bins = 50,histtype = u'step', ls = '-.',color='grey')
ax4.hist(bc_diff_control_mined*-1 /10, range = plot_range, bins = 50,histtype = u'step', ls = '--',color='k')
ax4.set_xticklabels([])

ax6.hist(lc_diff_0*-1 /10, range = plot_range, bins = 50, color = cmap[0],histtype = u'step')
ax6.hist(lc_diff_5*-1 /10, range = plot_range, bins = 50, color = cmap[5],histtype = u'step')
ax6.hist(lc_diff_10*-1 /10, range = plot_range, bins = 50, color = cmap[8],histtype = u'step')
ax6.hist(lc_diff_control*-1 /10, range = plot_range, bins = 50,histtype = u'step', ls = '-.',color='grey')
ax6.hist(lc_diff_control_mined*-1 /10, range = plot_range, bins = 50,histtype = u'step', ls = '--',color='k')
ax6.set_xticklabels([])

ax8.hist(sf_diff_0*-1 /10, range = plot_range, bins = 50, color = cmap[0],histtype = u'step')
ax8.hist(sf_diff_5*-1 /10, range = plot_range, bins = 50, color = cmap[5],histtype = u'step')
ax8.hist(sf_diff_10*-1 /10, range = plot_range, bins = 50, color = cmap[8],histtype = u'step')
ax8.hist(sf_diff_control*-1 /10, range = plot_range, bins = 50,histtype = u'step', ls = '-.',color='grey')
ax8.hist(sf_diff_control_mined*-1 /10, range = plot_range, bins = 50,histtype = u'step', ls = '--',color='k')


ax.set_yticklabels([])
ax2.set_yticklabels([])
ax4.set_yticklabels([])
ax6.set_yticklabels([])
ax8.set_yticklabels([])
ax8.set_xlabel('erosion rate (mm yr$^{-1}$)', fontsize = 12)
ax.set_ylabel('count', fontsize = 10)
ax2.set_ylabel('count', fontsize = 10)
ax4.set_ylabel('count', fontsize = 10)
ax6.set_ylabel('count', fontsize = 10)
ax8.set_ylabel('count', fontsize = 10)

ax.set_yscale('log')
ax2.set_yscale('log')
ax4.set_yscale('log')
ax6.set_yscale('log')
ax8.set_yscale('log')

ax.set_ylim(1e2,10e5)
ax2.set_ylim(1e2,10e5)
ax4.set_ylim(1e2,10e5)
ax6.set_ylim(1e2,10e5)
ax8.set_ylim(1e2,10e5)

ax1 = ax.twinx()
ax1.set_ylabel('Mud River')
ax1.set_yticks([])

ax3 = ax2.twinx()
ax3.set_ylabel('White Oak')
ax3.set_yticks([])

ax5 = ax4.twinx()
ax5.set_ylabel('Ben Creek')
ax5.set_yticks([])

ax7 = ax6.twinx()
ax7.set_ylabel('Laurel Creek')
ax7.set_yticks([])

ax9 = ax8.twinx()
ax9.set_ylabel('Spruce Fork')
ax9.set_yticks([])

ax.legend(loc = 'upper right')
ax.legend(ncol=3, fancybox=True, shadow=True)

fig.savefig('fig10.png', dpi = 1000)