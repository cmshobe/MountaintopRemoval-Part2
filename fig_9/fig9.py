########################################################################
#This script generates Figure 9 in the following manuscript:

#Bower, S.J., Shobe, C.M., Maxwell, A.E., and Campforts, B. (2023) The 
#uncertain future of mountaintop-removal-mined landscapes 2: Modeling the influence of
#topography and vegetation. Geomorphology.

#Please cite the paper if you use this code in any way.

#Brief description: this script uses the results of our numerical simulations
#(found in the folder `model_simulations`) to plot the percent change in sediment export
#between year 200 and year 10,000 from each study watershed in each scenario we considered.
#This is an indicator of whether the basin is on a trajectory of increasing or decreasing 
#sediment flux over time.

########################################################################

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

cm = matplotlib.colormaps['plasma']

#############import sediment flux record
#NOTE: sediment flux record comes out of the model in units of L/T, so need
#to multiply by grid cell size to get to L3/T

cellsize = 100 #m^3

#ben creek
bc0 = np.genfromtxt('../model_simulations/ben_creek/outputs/scenario_0_sedflux.txt') * cellsize
bc5 = np.genfromtxt('../model_simulations/ben_creek/outputs/scenario_5_sedflux.txt') * cellsize
bc10 = np.genfromtxt('../model_simulations/ben_creek/outputs/scenario_10_sedflux.txt') * cellsize
bc_control_mined = np.genfromtxt('../model_simulations/ben_creek/outputs/control_mined_sedflux.txt') * cellsize
bc_control_unmined = np.genfromtxt('../model_simulations/ben_creek/outputs/control_unmined_sedflux.txt') * cellsize

#mud river
mr0 = np.genfromtxt('../model_simulations/mud_river/outputs/scenario_0_sedflux.txt') * cellsize
mr5 = np.genfromtxt('../model_simulations/mud_river/outputs/scenario_5_sedflux.txt') * cellsize
mr10 = np.genfromtxt('../model_simulations/mud_river/outputs/scenario_10_sedflux.txt') * cellsize
mr_control_mined = np.genfromtxt('../model_simulations/mud_river/outputs/control_mined_sedflux.txt') * cellsize
mr_control_unmined = np.genfromtxt('../model_simulations/mud_river/outputs/control_unmined_sedflux.txt') * cellsize

#laurel creek
lc0 = np.genfromtxt('../model_simulations/laurel_creek/outputs/scenario_0_sedflux.txt') * cellsize
lc5 = np.genfromtxt('../model_simulations/laurel_creek/outputs/scenario_5_sedflux.txt') * cellsize
lc10 = np.genfromtxt('../model_simulations/laurel_creek/outputs/scenario_10_sedflux.txt') * cellsize
lc_control_mined = np.genfromtxt('../model_simulations/laurel_creek/outputs/control_mined_sedflux.txt') * cellsize
lc_control_unmined = np.genfromtxt('../model_simulations/laurel_creek/outputs/control_unmined_sedflux.txt') * cellsize

#spruce fork
sf0 = np.genfromtxt('../model_simulations/spruce_fork/outputs/scenario_0_sedflux.txt') * cellsize
sf5 = np.genfromtxt('../model_simulations/spruce_fork/outputs/scenario_5_sedflux.txt') * cellsize
sf10 = np.genfromtxt('../model_simulations/spruce_fork/outputs/scenario_10_sedflux.txt') * cellsize
sf_control_mined = np.genfromtxt('../model_simulations/spruce_fork/outputs/control_mined_sedflux.txt') * cellsize
sf_control_unmined = np.genfromtxt('../model_simulations/spruce_fork/outputs/control_unmined_sedflux.txt') * cellsize

#white oak creek
wo0 = np.genfromtxt('../model_simulations/white_oak/outputs/scenario_0_sedflux.txt') * cellsize
wo5 = np.genfromtxt('../model_simulations/white_oak/outputs/scenario_5_sedflux.txt') * cellsize
wo10 = np.genfromtxt('../model_simulations/white_oak/outputs/scenario_10_sedflux.txt') * cellsize
wo_control_mined = np.genfromtxt('../model_simulations/white_oak/outputs/control_mined_sedflux.txt') * cellsize
wo_control_unmined = np.genfromtxt('../model_simulations/white_oak/outputs/control_unmined_sedflux.txt') * cellsize

####################make figure

xvals = np.arange(5)

#calculate percent change between year 200 (timestep 400 because we used dt = 0.5 for 
#the first 200 years) and year 10,000 (timestep -1 in Python indexing)
#order: unmined, mined_ctrl, 0, 5, 10
mr_arr_flux = np.array([(mr_control_unmined[-1] - mr_control_unmined[400]) / mr_control_unmined[400],
                        (mr_control_mined[-1] - mr_control_mined[400]) / mr_control_mined[400],
                        (mr10[-1] - mr10[400]) / mr10[400],
                        (mr5[-1] - mr5[400]) / mr5[400],
                        (mr0[-1] - mr0[400]) / mr0[400]
                       ])*100

wo_arr_flux = np.array([(wo_control_unmined[-1] - wo_control_unmined[400]) / wo_control_unmined[400],
                        (wo_control_mined[-1] - wo_control_mined[400]) / wo_control_mined[400],
                        (wo10[-1] - wo10[400]) / wo10[400],
                        (wo5[-1] - wo5[400]) / wo5[400],
                        (wo0[-1] - wo0[400]) / wo0[400]
                       ])*100

bc_arr_flux = np.array([(bc_control_unmined[-1] - bc_control_unmined[400]) / bc_control_unmined[400],
                        (bc_control_mined[-1] - bc_control_mined[400]) / bc_control_mined[400],
                        (bc10[-1] - bc10[400]) / bc10[400],
                        (bc5[-1] - bc5[400]) / bc5[400],
                        (bc0[-1] - bc0[400]) / bc0[400]
                       ])*100

lc_arr_flux = np.array([(lc_control_unmined[-1] - lc_control_unmined[400]) / lc_control_unmined[400],
                        (lc_control_mined[-1] - lc_control_mined[400]) / lc_control_mined[400],
                        (lc10[-1] - lc10[400]) / lc10[400],
                        (lc5[-1] - lc5[400]) / lc5[400],
                        (lc0[-1] - lc0[400]) / lc0[400]
                       ])*100

sf_arr_flux = np.array([(sf_control_unmined[-1] - sf_control_unmined[400]) / sf_control_unmined[400],
                        (sf_control_mined[-1] - sf_control_mined[400]) / sf_control_mined[400],
                        (sf10[-1] - sf10[400]) / sf10[400],
                        (sf5[-1] - sf5[400]) / sf5[400],
                        (sf0[-1] - sf0[400]) / sf0[400]
                       ])*100


fig3,ax = plt.subplots(figsize=(4, 4))
colors = ['grey', 'k', cm(0.8), cm(0.5), cm(0)]
ax.scatter(xvals, mr_arr_flux, s = 100, marker = 's', c = colors, edgecolor = 'k',alpha = 0.75, label = 'Mud River')
ax.scatter(xvals, wo_arr_flux, s = 450, marker = '.', c = colors, edgecolor = 'k',alpha = 0.75, label = 'White Oak')
ax.scatter(xvals, bc_arr_flux, s = 140, marker = 'p', c = colors, edgecolor = 'k',alpha = 0.75, label = 'Ben Creek')
ax.scatter(xvals, lc_arr_flux, s = 140, marker = '^', c = colors, edgecolor = 'k',alpha = 0.75, label = 'Laurel Creek')
ax.scatter(xvals, sf_arr_flux, s = 160, marker = '*', c = colors, edgecolor = 'k',alpha = 0.75, label = 'Spruce Fork')


ax.tick_params(axis='x', which='both', length=0)
plt.setp(ax.get_xticklabels(), visible=False)
ax.text(0, -0.05, 'unmined', transform=ax.transAxes)
ax.text(0.15, -0.1, 'mined control', transform=ax.transAxes)
ax.set_xlim(-0.5, 4.5)

ax.xaxis.set_label_coords(.73, -0.01)
ax.set_ylabel('% change in sediment flux')

ax.axhline(0, -1, 5, linestyle = '-', color = 'grey')

handles, labels = ax.get_legend_handles_labels()
legend1 = ax.legend(handles[::-1], labels[::-1], framealpha = 1, bbox_to_anchor = (0.48,0.22))

#second legend to do reforestation progress colors
import matplotlib.lines as mlines

blue_star = mlines.Line2D([], [], color=cm(0.8), marker='s', linestyle='None',markeredgecolor = 'k',alpha = 0.75,
                          markersize=10, label='100%')
red_square = mlines.Line2D([], [], color=cm(0.5), marker='s', linestyle='None',markeredgecolor = 'k',alpha = 0.75,
                          markersize=10, label='50%')
purple_triangle = mlines.Line2D([], [], color=cm(0), marker='s', linestyle='None',markeredgecolor = 'k',alpha = 0.75,
                          markersize=10, label='0%')

legend2 = ax.legend(title = 'Percent' '\n' 'recovery', handles=[blue_star, red_square, purple_triangle], bbox_to_anchor = (0.4,1))

frame = legend1.get_frame()
frame.set_edgecolor('k')
frame2 = legend2.get_frame()
frame2.set_edgecolor('k')
ax.add_artist(legend1)

ax.text(2.2, 1, 'inc. flux over time', color = 'grey')
ax.text(2.2, -2, 'dec. flux over time', color = 'grey')

ax.annotate('altered K scenarios', xy=(0.7, 0.03), xytext=(0.7, -0.08), xycoords='axes fraction', 
            ha='center', va='bottom',
            bbox=None,
            arrowprops=dict(arrowstyle='-[, widthB=6.0, lengthB=2', lw=1.0, color='k'))
            
fig3.savefig('fig9.png', dpi=1000, bbox_inches='tight', facecolor = 'white')