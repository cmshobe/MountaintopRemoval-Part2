########################################################################
#This script generates Figure 7 in the following manuscript:

#Bower, S.J., Shobe, C.M., Maxwell, A.E., and Campforts, B. (2023) The 
#uncertain future of mountaintop-removal-mined landscapes 2: Modeling the influence of
#topography and vegetation. Geomorphology.

#Please cite the paper if you use this code in any way.

#Brief description: this script uses the results of our numerical simulations
#(found in the folder `model_simulations`) to plot the total sediment export
#after 10,000 years from each study watershed in each scenario we considered.

########################################################################

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

cm = matplotlib.colormaps['plasma']

###############import pre and post data

#mud river
mr_mined_pre = np.loadtxt('../model_simulations/mud_river/inputs/mudriver_post_10m.asc', skiprows=6)
mr0_post = np.loadtxt('../model_simulations/mud_river/outputs/scenario_0_elev_01.asc', skiprows=6)
mr5_post = np.loadtxt('../model_simulations/mud_river/outputs/scenario_5_elev_01.asc', skiprows=6)
mr10_post = np.loadtxt('../model_simulations/mud_river/outputs/scenario_10_elev_01.asc', skiprows=6)
mr_ctrl_mined_post = np.loadtxt('../model_simulations/mud_river/outputs/control_mined_elev_01.asc', skiprows=6)
mr_unmined_pre = np.loadtxt('../model_simulations/mud_river/inputs/mudriver_pre_10m.asc', skiprows=6)
mr_ctrl_unmined_post = np.loadtxt('../model_simulations/mud_river/outputs/control_unmined_elev_01.asc', skiprows=6)

#white oak creek
wo_mined_pre = np.loadtxt('../model_simulations/white_oak/inputs/whiteoak_post_10m.asc', skiprows=6)
wo0_post = np.loadtxt('../model_simulations/white_oak/outputs/scenario_0_elev_01.asc', skiprows=6)
wo5_post = np.loadtxt('../model_simulations/white_oak/outputs/scenario_5_elev_01.asc', skiprows=6)
wo10_post = np.loadtxt('../model_simulations/white_oak/outputs/scenario_10_elev_01.asc', skiprows=6)
wo_ctrl_mined_post = np.loadtxt('../model_simulations/white_oak/outputs/control_mined_elev_01.asc', skiprows=6)
wo_unmined_pre = np.loadtxt('../model_simulations/white_oak/inputs/whiteoak_pre_10m.asc', skiprows=6)
wo_ctrl_unmined_post = np.loadtxt('../model_simulations/white_oak/outputs/control_unmined_elev_01.asc', skiprows=6)

#ben creek
bc_mined_pre = np.loadtxt('../model_simulations/ben_creek/inputs/bencreek_post_10m.asc', skiprows=6)
bc0_post = np.loadtxt('../model_simulations/ben_creek/outputs/scenario_0_elev_01.asc', skiprows=6)
bc5_post = np.loadtxt('../model_simulations/ben_creek/outputs/scenario_5_elev_01.asc', skiprows=6)
bc10_post = np.loadtxt('../model_simulations/ben_creek/outputs/scenario_10_elev_01.asc', skiprows=6)
bc_ctrl_mined_post = np.loadtxt('../model_simulations/ben_creek/outputs/control_mined_elev_01.asc', skiprows=6)
bc_unmined_pre = np.loadtxt('../model_simulations/ben_creek/inputs/bencreek_pre_10m.asc', skiprows=6)
bc_ctrl_unmined_post = np.loadtxt('../model_simulations/ben_creek/outputs/control_unmined_elev_01.asc', skiprows=6)

#laurel creek
lc_mined_pre = np.loadtxt('../model_simulations/laurel_creek/inputs/laurelcreek_post_10m.asc', skiprows=6)
lc0_post = np.loadtxt('../model_simulations/laurel_creek/outputs/scenario_0_elev_01.asc', skiprows=6)
lc5_post = np.loadtxt('../model_simulations/laurel_creek/outputs/scenario_5_elev_01.asc', skiprows=6)
lc10_post = np.loadtxt('../model_simulations/laurel_creek/outputs/scenario_10_elev_01.asc', skiprows=6)
lc_ctrl_mined_post = np.loadtxt('../model_simulations/laurel_creek/outputs/control_mined_elev_01.asc', skiprows=6)
lc_unmined_pre = np.loadtxt('../model_simulations/laurel_creek/inputs/laurelcreek_pre_10m.asc', skiprows=6)
lc_ctrl_unmined_post = np.loadtxt('../model_simulations/laurel_creek/outputs/control_unmined_elev_01.asc', skiprows=6)

#spruce fork
sf_mined_pre = np.loadtxt('../model_simulations/spruce_fork/inputs/sprucefork_fixed_post_10m.asc', skiprows=6)
sf0_post = np.loadtxt('../model_simulations/spruce_fork/outputs/scenario_0_elev_01.asc', skiprows=6)
sf5_post = np.loadtxt('../model_simulations/spruce_fork/outputs/scenario_5_elev_01.asc', skiprows=6)
sf10_post = np.loadtxt('../model_simulations/spruce_fork/outputs/scenario_10_elev_01.asc', skiprows=6)
sf_ctrl_mined_post = np.loadtxt('../model_simulations/spruce_fork/outputs/control_mined_elev_01.asc', skiprows=6)
sf_unmined_pre = np.loadtxt('../model_simulations/spruce_fork/inputs/sprucefork_fixed_pre_10m.asc', skiprows=6)
sf_ctrl_unmined_post = np.loadtxt('../model_simulations/spruce_fork/outputs/control_unmined_elev_01.asc', skiprows=6)


#take the volumetric topographic difference between beginning and end of each run

cellsize = 100 #m^2

mr0_diff = np.sum(mr_mined_pre[1:] - mr0_post) * cellsize #pre for some reason has an extra line of NODATA
mr5_diff = np.sum(mr_mined_pre[1:] - mr5_post) * cellsize
mr10_diff = np.sum(mr_mined_pre[1:] - mr10_post) * cellsize
mr_ctrl_mined_diff = np.sum(mr_mined_pre[1:] - mr_ctrl_mined_post) * cellsize
mr_ctrl_unmined_diff = np.sum(mr_unmined_pre[1:] - mr_ctrl_unmined_post) * cellsize

wo0_diff = np.sum(wo_mined_pre[1:] - wo0_post) * cellsize #pre for some reason has an extra line of NODATA
wo5_diff = np.sum(wo_mined_pre[1:] - wo5_post) * cellsize
wo10_diff = np.sum(wo_mined_pre[1:] - wo10_post) * cellsize
wo_ctrl_mined_diff = np.sum(wo_mined_pre[1:] - wo_ctrl_mined_post) * cellsize
wo_ctrl_unmined_diff = np.sum(wo_unmined_pre[1:] - wo_ctrl_unmined_post) * cellsize

bc0_diff = np.sum(bc_mined_pre[1:] - bc0_post) * cellsize #pre for some reason has an extra line of NODATA
bc5_diff = np.sum(bc_mined_pre[1:] - bc5_post) * cellsize
bc10_diff = np.sum(bc_mined_pre[1:] - bc10_post) * cellsize
bc_ctrl_mined_diff = np.sum(bc_mined_pre[1:] - bc_ctrl_mined_post) * cellsize
bc_ctrl_unmined_diff = np.sum(bc_unmined_pre[1:] - bc_ctrl_unmined_post) * cellsize

lc0_diff = np.sum(lc_mined_pre[1:] - lc0_post) * cellsize #pre for some reason has an extra line of NODATA
lc5_diff = np.sum(lc_mined_pre[1:] - lc5_post) * cellsize
lc10_diff = np.sum(lc_mined_pre[1:] - lc10_post) * cellsize
lc_ctrl_mined_diff = np.sum(lc_mined_pre[1:] - lc_ctrl_mined_post) * cellsize
lc_ctrl_unmined_diff = np.sum(lc_unmined_pre[1:] - lc_ctrl_unmined_post) * cellsize

sf0_diff = np.sum(sf_mined_pre[1:] - sf0_post) * cellsize #pre for some reason has an extra line of NODATA
sf5_diff = np.sum(sf_mined_pre[1:] - sf5_post) * cellsize
sf10_diff = np.sum(sf_mined_pre[1:] - sf10_post) * cellsize
sf_ctrl_mined_diff = np.sum(sf_mined_pre[1:] - sf_ctrl_mined_post) * cellsize
sf_ctrl_unmined_diff = np.sum(sf_unmined_pre[1:] - sf_ctrl_unmined_post) * cellsize


#prepare x labels
xvals = np.arange(5)
mr_arr = np.array([mr_ctrl_unmined_diff, mr_ctrl_mined_diff, mr10_diff, mr5_diff, mr0_diff])
wo_arr = np.array([wo_ctrl_unmined_diff, wo_ctrl_mined_diff, wo10_diff, wo5_diff, wo0_diff])
bc_arr = np.array([bc_ctrl_unmined_diff, bc_ctrl_mined_diff, bc10_diff, bc5_diff, bc0_diff])
lc_arr = np.array([lc_ctrl_unmined_diff, lc_ctrl_mined_diff, lc10_diff, lc5_diff, lc0_diff])
sf_arr = np.array([sf_ctrl_unmined_diff, sf_ctrl_mined_diff, sf10_diff, sf5_diff, sf0_diff])

fig2,ax = plt.subplots(figsize=(4, 4))
colors = ['grey', 'k', cm(0.8), cm(0.5), cm(0)]
ax.scatter(xvals, mr_arr, s = 100, marker = 's', c = colors, edgecolor = 'k',alpha = 0.75, label = 'Mud River')
ax.scatter(xvals, wo_arr, s = 450, marker = '.', c = colors, edgecolor = 'k',alpha = 0.75, label = 'White Oak')
ax.scatter(xvals, bc_arr, s = 140, marker = 'p', c = colors, edgecolor = 'k',alpha = 0.75, label = 'Ben Creek')
ax.scatter(xvals, lc_arr, s = 140, marker = '^', c = colors, edgecolor = 'k',alpha = 0.75, label = 'Laurel Creek')
ax.scatter(xvals, sf_arr, s = 160, marker = '*', c = colors, edgecolor = 'k',alpha = 0.75, label = 'Spruce Fork')

ax.set_yscale('log')
ax.tick_params(axis='x', which='both', length=0)
plt.setp(ax.get_xticklabels(), visible=False)
ax.text(0, -0.05, 'unmined', transform=ax.transAxes)
ax.text(0.15, -0.1, 'mined control', transform=ax.transAxes)
ax.set_xlim(-0.5, 4.5)

ax.xaxis.set_label_coords(.73, -0.01)
ax.set_ylabel('total sediment export (m$^3$)')

handles, labels = ax.get_legend_handles_labels()
legend1 = ax.legend(handles[::-1], labels[::-1], framealpha = 1)

#second legend to do reforestation progress colors
import matplotlib.lines as mlines

blue_star = mlines.Line2D([], [], color=cm(0.8), marker='s', linestyle='None',markeredgecolor = 'k',alpha = 0.75,
                          markersize=10, label='100%')
red_square = mlines.Line2D([], [], color=cm(0.5), marker='s', linestyle='None',markeredgecolor = 'k',alpha = 0.75,
                          markersize=10, label='50%')
purple_triangle = mlines.Line2D([], [], color=cm(0), marker='s', linestyle='None',markeredgecolor = 'k',alpha = 0.75,
                          markersize=10, label='0%')

legend2 = ax.legend(title = 'Percent' '\n' 'recovery', handles=[blue_star, red_square, purple_triangle], bbox_to_anchor=(0.6,0.45))

frame = legend1.get_frame()
frame.set_edgecolor('k')
frame2 = legend2.get_frame()
frame2.set_edgecolor('k')
ax.add_artist(legend1)

ax.annotate('altered K scenarios', xy=(0.7, 0.03), xytext=(0.7, -0.08), xycoords='axes fraction', 
            ha='center', va='bottom',
            bbox=None,
            arrowprops=dict(arrowstyle='-[, widthB=6.0, lengthB=2', lw=1.0, color='k'))
            
fig2.savefig('fig7.png', dpi=1000, bbox_inches='tight')