########################################################################
#This script generates Figure 8 in the following manuscript:

#Bower, S.J., Shobe, C.M., Maxwell, A.E., and Campforts, B. (2023) The 
#uncertain future of mountaintop-removal-mined landscapes 2: Modeling the influence of
#topography and vegetation. Geomorphology.

#Please cite the paper if you use this code in any way.

#Brief description: this script uses the results of our numerical simulations
#(found in the folder `model_simulations`) to plot the cumulative sediment export
#over 10,000 years from each study watershed in each scenario we considered.

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

#convert from Qs (vol/time) to total sed export (volume) by multiplying each entry by dt:
#0-200 get multiplied by 0.5, 200 onward get multiplied by 1

bc0[:400] *= 0.5
bc5[:400] *= 0.5
bc10[:400] *= 0.5
bc_control_mined[:400] *= 0.5
bc_control_unmined[:400] *= 0.5

lc0[:400] *= 0.5
lc5[:400] *= 0.5
lc10[:400] *= 0.5
lc_control_mined[:400] *= 0.5
lc_control_unmined[:400] *= 0.5

mr0[:400] *= 0.5
mr5[:400] *= 0.5
mr10[:400] *= 0.5
mr_control_mined[:400] *= 0.5
mr_control_unmined[:400] *= 0.5

sf0[:400] *= 0.5
sf5[:400] *= 0.5
sf10[:400] *= 0.5
sf_control_mined[:400] *= 0.5
sf_control_unmined[:400] *= 0.5

wo0[:400] *= 0.5
wo5[:400] *= 0.5
wo10[:400] *= 0.5
wo_control_mined[:400] *= 0.5
wo_control_unmined[:400] *= 0.5

###################make figure
from matplotlib import gridspec

fig = plt.figure(figsize=(8, 9))
gs = gridspec.GridSpec(5, 2)

#200 yr plots
mr_recov = fig.add_subplot(gs[0, 0])
wo_recov = fig.add_subplot(gs[1, 0])
bc_recov = fig.add_subplot(gs[2, 0])
lc_recov = fig.add_subplot(gs[3, 0])
sf_recov = fig.add_subplot(gs[4, 0])

#10,000 yr plots
mr_long = fig.add_subplot(gs[0, 1])
wo_long = fig.add_subplot(gs[1, 1])
bc_long = fig.add_subplot(gs[2, 1])
lc_long = fig.add_subplot(gs[3, 1])
sf_long = fig.add_subplot(gs[4, 1])

#mr short and long figure

mr_recov.plot(x_vals, np.cumsum(mr0), label = '0%', color = cm(0))
mr_recov.plot(x_vals, np.cumsum(mr5), label = '50%', color = cm(0.5))
mr_recov.plot(x_vals, np.cumsum(mr10), label = '100%', color = cm(0.8))
mr_recov.plot(x_vals, np.cumsum(mr_control_mined), label = 'Post-MTR control', color = 'k', linestyle = '--')
mr_recov.plot(x_vals, np.cumsum(mr_control_unmined), label = 'Pre-MTR control', color = 'grey', linestyle = '-.')
mr_recov.set_xlim(0, 200)
mr_recov.set_ylim(1e3,2e6)
mr_recov.set_yscale('log')


mr_long.plot(x_vals/1000, np.cumsum(mr0), label = '0%', color = cm(0))
mr_long.plot(x_vals/1000, np.cumsum(mr5), label = '50%', color = cm(0.5))
mr_long.plot(x_vals/1000, np.cumsum(mr10), label = '100%', color = cm(0.8))
mr_long.plot(x_vals/1000, np.cumsum(mr_control_mined), label = 'Post-MTR control', color = 'k', linestyle = '--')
mr_long.plot(x_vals/1000, np.cumsum(mr_control_unmined), label = 'Pre-MTR control', color = 'grey', linestyle = '-.')
mr_long.set_xlim(0, 10)
mr_long.set_yscale('log')
mr_long.legend(ncol = 5, bbox_to_anchor=(-0.1, 1.15), loc='center', framealpha = 0)
mr_long.set_ylim(1e4, 1e8)

#white oak short and long figure
wo_recov.plot(x_vals, np.cumsum(wo0), label = '0%', color = cm(0))
wo_recov.plot(x_vals, np.cumsum(wo5), label = '50%', color = cm(0.5))
wo_recov.plot(x_vals, np.cumsum(wo10), label = '100%', color = cm(0.8))
wo_recov.plot(x_vals, np.cumsum(wo_control_mined), label = 'Post-MTR control ', color = 'k', linestyle = '--')
wo_recov.plot(x_vals, np.cumsum(wo_control_unmined), label = 'Pre-MTR control', color = 'grey', linestyle = '-.')
wo_recov.set_xlim(0, 200)
wo_recov.set_ylim(1e3,2e6)
wo_recov.set_yscale('log')

wo_long.plot(x_vals/1000, np.cumsum(wo0), label = '0%', color = cm(0))
wo_long.plot(x_vals/1000, np.cumsum(wo5), label = '50%', color = cm(0.5))
wo_long.plot(x_vals/1000, np.cumsum(wo10), label = '100%', color = cm(0.8))
wo_long.plot(x_vals/1000, np.cumsum(wo_control_mined), label = 'Post-MTR control', color = 'k', linestyle = '--')
wo_long.plot(x_vals/1000, np.cumsum(wo_control_unmined), label = 'Pre-MTR control', color = 'grey', linestyle = '-.')
wo_long.set_xlim(0, 10)
wo_long.set_yscale('log')
wo_long.set_ylim(1e4, 2e8)

#ben creek short and long figure
bc_recov.plot(x_vals, np.cumsum(bc0), label = '0%', color = cm(0))
bc_recov.plot(x_vals, np.cumsum(bc5), label = '50%', color = cm(0.5))
bc_recov.plot(x_vals, np.cumsum(bc10), label = '100%', color = cm(0.8))
bc_recov.plot(x_vals, np.cumsum(bc_control_mined), label = 'Post-MTR control', color = 'k', linestyle = '--')
bc_recov.plot(x_vals, np.cumsum(bc_control_unmined), label = 'Pre-MTR control', color = 'grey', linestyle = '-.')
bc_recov.set_xlim(0, 200)
bc_recov.set_ylim(1e3,2e6)
bc_recov.set_yscale('log')
bc_recov.set_ylabel('cumulative sediment export [m$^3$]')

bc_long.plot(x_vals/1000, np.cumsum(bc0), label = '0%', color = cm(0))
bc_long.plot(x_vals/1000, np.cumsum(bc5), label = '50%', color = cm(0.5))
bc_long.plot(x_vals/1000, np.cumsum(bc10), label = '100%', color = cm(0.8))
bc_long.plot(x_vals/1000, np.cumsum(bc_control_mined), label = 'Post-MTR control', color = 'k', linestyle = '--')
bc_long.plot(x_vals/1000, np.cumsum(bc_control_unmined), label = 'Pre-MTR control', color = 'grey', linestyle = '-.')
bc_long.set_xlim(0, 10)
bc_long.set_yscale('log')
bc_long.set_ylim(1e4, 2e8)

#laurel creek short and long figure
lc_recov.plot(x_vals, np.cumsum(lc0), label = '0%', color = cm(0))
lc_recov.plot(x_vals, np.cumsum(lc5), label = '50%', color = cm(0.5))
lc_recov.plot(x_vals, np.cumsum(lc10), label = '100%', color = cm(0.8))
lc_recov.plot(x_vals, np.cumsum(lc_control_mined), label = 'Post-MTR control', color = 'k', linestyle = '--')
lc_recov.plot(x_vals, np.cumsum(lc_control_unmined), label = 'Pre-MTR control', color = 'grey', linestyle = '-.')
lc_recov.set_xlim(0, 200)
lc_recov.set_ylim(2e3,5e6)
lc_recov.set_yscale('log')
lc_recov.minorticks_off()

lc_long.plot(x_vals/1000, np.cumsum(lc0), label = '0%', color = cm(0))
lc_long.plot(x_vals/1000, np.cumsum(lc5), label = '50%', color = cm(0.5))
lc_long.plot(x_vals/1000, np.cumsum(lc10), label = '100%', color = cm(0.8))
lc_long.plot(x_vals/1000, np.cumsum(lc_control_mined), label = 'Post-MTR control', color = 'k', linestyle = '--')
lc_long.plot(x_vals/1000, np.cumsum(lc_control_unmined), label = 'Pre-MTR control', color = 'grey', linestyle = '-.')
lc_long.set_xlim(0, 10)
lc_long.set_ylim(1e5, 3e8)
lc_long.set_yscale('log')

#spruce fork short and long figure
sf_recov.plot(x_vals, np.cumsum(sf0), label = '0%', color = cm(0))
sf_recov.plot(x_vals, np.cumsum(sf5), label = '50%', color = cm(0.5))
sf_recov.plot(x_vals, np.cumsum(sf10), label = '100%', color = cm(0.8))
sf_recov.plot(x_vals, np.cumsum(sf_control_mined), label = 'Post-MTR control', color = 'k', linestyle = '--')
sf_recov.plot(x_vals, np.cumsum(sf_control_unmined), label = 'Pre-MTR control', color = 'grey', linestyle = '-.')
sf_recov.set_xlim(0, 200)
sf_recov.set_ylim(2e3,5e6)
sf_recov.set_yscale('log')
sf_recov.set_xlabel('time (yr)')
sf_recov.minorticks_off()

sf_long.plot(x_vals/1000, np.cumsum(sf0), label = '0%', color = cm(0))
sf_long.plot(x_vals/1000, np.cumsum(sf5), label = '50%', color = cm(0.5))
sf_long.plot(x_vals/1000, np.cumsum(sf10), label = '100%', color = cm(0.8))
sf_long.plot(x_vals/1000, np.cumsum(sf_control_mined), label = 'Post-MTR control', color = 'k', linestyle = '--')
sf_long.plot(x_vals/1000, np.cumsum(sf_control_unmined), label = 'Pre-MTR control', color = 'grey', linestyle = '-.')
sf_long.set_xlim(0, 10)
sf_long.set_ylim(1e5, 3e8)
sf_long.set_yscale('log')
sf_long.set_xlabel('time (kyr)')

#text annotations

textx = 0.03
texty = 0.05

mr_recov.text(textx, texty, 'Mud River: recovery period', transform=mr_recov.transAxes)
mr_long.text(textx, texty, 'Mud River: long term', transform=mr_long.transAxes)

wo_recov.text(textx, texty, 'White Oak: recovery period', transform=wo_recov.transAxes)
wo_long.text(textx, texty, 'White Oak: long term', transform=wo_long.transAxes)

bc_recov.text(textx, texty, 'Ben Creek: recovery period', transform=bc_recov.transAxes)
bc_long.text(textx, texty, 'Ben Creek: long term', transform=bc_long.transAxes)

lc_recov.text(textx, texty, 'Laurel Creek: recovery period', transform=lc_recov.transAxes)
lc_long.text(textx, texty, 'Laurel Creek: long term', transform=lc_long.transAxes)

sf_recov.text(textx, texty, 'Spruce Fork: recovery period', transform=sf_recov.transAxes)
sf_long.text(textx, texty, 'Spruce Fork: long term', transform=sf_long.transAxes)

fig.savefig('fig8.png', dpi=1000, bbox_inches='tight')
