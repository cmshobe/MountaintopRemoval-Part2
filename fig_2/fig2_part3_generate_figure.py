########################################################################
#This script generates Figure 2 in the following manuscript:

#Bower, S.J., Shobe, C.M., Maxwell, A.E., and Campforts, B. (2023) The 
#uncertain future of mountaintop-removal-mined landscapes 2: Modeling the influence of
#topography and vegetation. Geomorphology.

#Please cite the paper if you use this code in any way.

#Brief description: this script uses distributions of elevation, slope, and 
#area--slope product from pre- and post-mining DEMs of five study watersheds
#(calculated in the part 1 script) as well as the results of Bayesian rank
#correlations (calculated using the R scripts in this folder) to produce the figure.

########################################################################

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats
from scipy.stats import norm
import seaborn as sbn

#read in flattened morphometric data from each watershed
flat_data_input_path = '/flat_inputs/'

paths = ['bc_pre_elev','bc_post_elev','bc_pre_slope','bc_post_slope','bc_pre_SA','bc_post_SA',
        'mr_pre_elev','mr_post_elev','mr_pre_slope','mr_post_slope','mr_pre_SA','mr_post_SA',
        'wo_pre_elev','wo_post_elev','wo_pre_slope','wo_post_slope','wo_pre_SA','wo_post_SA',
        'sf_pre_elev','sf_post_elev','sf_pre_slope','sf_post_slope','sf_pre_SA','sf_post_SA',
        'lc_pre_elev','lc_post_elev','lc_pre_slope','lc_post_slope','lc_pre_SA','lc_post_SA']

for string in paths:
    globals()[f'{string}'] = np.loadtxt(flat_data_input_path+str(string)+'.txt')
    
raw_data = [bc_pre_elev,bc_post_elev,bc_pre_slope,bc_post_slope,bc_pre_SA,bc_post_SA,
        mr_pre_elev,mr_post_elev,mr_pre_slope,mr_post_slope,mr_pre_SA,mr_post_SA,
        wo_pre_elev,wo_post_elev,wo_pre_slope,wo_post_slope,wo_pre_SA,wo_post_SA,
        sf_pre_elev,sf_post_elev,sf_pre_slope,sf_post_slope,sf_pre_SA,sf_post_SA,
        lc_pre_elev,lc_post_elev,lc_pre_slope,lc_post_slope,lc_pre_SA,lc_post_SA]

#Import results of Bayesian Wilcoxon statistical analyses

bayesian_wilcoxon_path = '/R_outputs/'
bw_test = pd.read_csv(bayesian_wilcoxon_path+'sam_bc_elev_samples.csv',delimiter = ',')
bayesian_stats = pd.read_csv(bayesian_wilcoxon_path+'sam_bc_elev.csv',delimiter = ',')

paths = ['sam_bc_elev','sam_bc_slope','sam_bc_SA',
         'sam_mr_elev','sam_mr_slope','sam_mr_SA',
         'sam_wo_elev','sam_wo_slope','sam_wo_SA',
         'sam_sf_elev','sam_sf_slope','sam_sf_SA',
         'sam_lc_elev','sam_lc_slope','sam_lc_SA']
for string in paths:
    globals()[f'{string}'] = pd.read_csv(bayesian_wilcoxon_path+str(string)+'.csv',delimiter = ',')
    
paths = ['sam_bc_elev_samples','sam_bc_slope_samples','sam_bc_SA_samples',
         'sam_mr_elev_samples','sam_mr_slope_samples','sam_mr_SA_samples',
         'sam_wo_elev_samples','sam_wo_slope_samples','sam_wo_SA_samples',
         'sam_sf_elev_samples','sam_sf_slope_samples','sam_sf_SA_samples',
         'sam_lc_elev_samples','sam_lc_slope_samples','sam_lc_SA_samples']
for string in paths:
    globals()[f'{string}'] = pd.read_csv(bayesian_wilcoxon_path+str(string)+'.csv',delimiter = ',')


bayesian_wilcoxon_density = [
    sam_bc_elev_samples.x,sam_bc_slope_samples.x,sam_bc_SA_samples.x,
    sam_mr_elev_samples.x,sam_mr_slope_samples.x,sam_mr_SA_samples.x,
    sam_wo_elev_samples.x,sam_wo_slope_samples.x,sam_wo_SA_samples.x,
    sam_sf_elev_samples.x,sam_sf_slope_samples.x,sam_sf_SA_samples.x,
    sam_lc_elev_samples.x,sam_lc_slope_samples.x,sam_lc_SA_samples.x
]

bayesian_stats = [sam_bc_elev,sam_bc_slope,sam_bc_SA,
                  sam_mr_elev,sam_mr_slope,sam_mr_SA,
                  sam_wo_elev,sam_wo_slope,sam_wo_SA,
                  sam_sf_elev,sam_sf_slope,sam_sf_SA,
                  sam_lc_elev,sam_lc_slope,sam_lc_SA]
                  
#make the figure
%precision 2

fig2 = plt.figure(
    figsize = (8,9),
    tight_layout = False
)

spec = fig2.add_gridspec(
    10,
    12,
    wspace = 0.25,
    hspace= 0.1
)

y_SA_range = (0,0.15)
y_elev_range = (0,0.01)
y_slope_range = (0,3.5)
SA_range = (0,60)
elev_range = (200,700)
slope_range = (0,1.5)
bins = 100
elev_bins = 50
elev_loc = (525,0.0085)
slope_loc = (0.75,3)
SA_loc = (40,0.125)

post_color = '#440154' #pee
pre_color = '#72ed97'  #blue

# post_color = '#721f81' #pee
# pre_color = '#fcfdbf'  #purple

# post_color = '#6b3fc4' #purple
# pre_color = '#72ed97'  #green

pre_edgecolor = 'k'
post_edgecolor = 'k'

#BEN CREEK
#elevation
ax1 = fig2.add_subplot(spec[:2,:4])
ax1.hist(
    bc_post_elev,
    density = True,
    alpha = 1,
    bins = elev_bins,
    histtype = u'step',
    color = post_edgecolor,
    lw = 1,
    range = elev_range
)
ax1.hist(
    bc_pre_elev,
    density = True,
    alpha = 1,
    bins = elev_bins,
    histtype = u'step',
    color = pre_edgecolor,
    lw = 1,
    range = elev_range
)
ax1.hist(
    bc_post_elev,
    density = True,
    alpha = 0.75,
    bins = elev_bins,
    color = post_color,
    label = 'mined',
    range = elev_range
)
ax1.hist(
    bc_pre_elev,
    density = True,
    alpha = 0.75,
    bins = elev_bins,
    color = pre_color,
    label='unmined',
    range = elev_range
)
ax1.set_ylabel('Ben Creek',fontsize = 12)
#Slope
ax2 = fig2.add_subplot(spec[:2,4:8])
ax2.hist(
    bc_post_slope,
    density = True,
    alpha = 1,
    bins = bins,
    histtype = u'step',
    color = post_edgecolor,
    lw = 1,
    range = slope_range
)
ax2.hist(
    bc_pre_slope,
    density = True,
    alpha = 1,
    bins = bins,
    histtype = u'step',
    color = pre_edgecolor,
    lw = 1,
    range = slope_range
)
ax2.hist(
    bc_post_slope,
    density = True,
    alpha = 0.75,
    bins = bins,
    color = post_color,
    label = 'mined',
    range = slope_range
)
ax2.hist(
    bc_pre_slope,
    density = True,
    alpha = 0.75,
    bins = bins,
    color = pre_color,
    label='unmined',
    range = slope_range
)
#slope*area
ax3 = fig2.add_subplot(spec[:2,8:12])
ax3.hist(
    bc_post_SA,
    density = True,
    alpha = 1,
    bins = bins,
    histtype = u'step',
    color = post_edgecolor,
    lw = 1,
    range = SA_range
)
ax3.hist(
    bc_pre_SA,
    density = True,
    alpha = 1,
    bins = bins,
    histtype = u'step',
    color = pre_edgecolor,
    lw = 1,
    range = SA_range
)
ax3.hist(
    bc_post_SA,
    density = True,
    alpha = 0.75,
    bins = bins,
    color = post_color,
    label = 'mined',
    range = SA_range
)
ax3.hist(
    bc_pre_SA,
    density = True,
    alpha = 0.75,
    bins = bins,
    color = pre_color,
    label='unmined',
    range = SA_range
)
ax3.legend(fontsize = 10, loc = 'lower right')

# MUD RIVER
#elevation
ax4 = fig2.add_subplot(spec[2:4,:4])
ax4.hist(
    mr_post_elev,
    density = True,
    alpha = 1, 
    bins = elev_bins,
    histtype = u'step', 
    color = post_edgecolor, 
    lw = 1,
    range = elev_range
)
ax4.hist(
    mr_pre_elev,
    density = True,
    alpha = 1,
    bins = elev_bins,
    histtype = u'step',
    color = pre_edgecolor,
    lw = 1,
    range = elev_range
)
ax4.hist(
    mr_post_elev,
    density = True,
    alpha = 0.75,
    bins = elev_bins,
    color = post_color,
    label = 'mined',
    range = elev_range
)
ax4.hist(
    mr_pre_elev,
    density = True,
    alpha = 0.75,
    bins = elev_bins,
    color = pre_color,
    label='unmined',
    range = elev_range
)
ax4.tick_params(axis='both', which='major', labelsize=14)
ax4.set_ylabel('Mud River',fontsize = 12)
#slope
ax5 = fig2.add_subplot(spec[2:4,4:8])
ax5.hist(
    mr_post_slope,
    density = True,
    alpha = 1,
    bins = bins,
    histtype = u'step',
    color = post_edgecolor,
    lw = 1,
    range = slope_range
)
ax5.hist(
    mr_pre_slope,
    density = True,
    alpha = 1,
    bins = bins,
    histtype = u'step',
    color = pre_edgecolor,
    lw = 1,
    range = slope_range
)
ax5.hist(
    mr_post_slope,
    density = True, 
    alpha = 0.75, 
    bins = bins, 
    color = post_color, 
    label = 'mined',
    range = slope_range
)
ax5.hist(
    mr_pre_slope,
    density = True,
    alpha = 0.75,
    bins = bins, 
    color = pre_color,
    label='unmined',
    range = slope_range
)
#slope*area
ax6 = fig2.add_subplot(spec[2:4,8:12])
ax6.hist(
    mr_post_SA,
    density = True,
    alpha = 1,
    bins = bins,
    histtype = u'step',
    color = post_edgecolor,
    lw = 1,
    range = SA_range
)
ax6.hist(
    mr_pre_SA,
    density = True,
    alpha = 1,
    bins = bins,
    histtype = u'step',
    color = pre_edgecolor,
    lw = 1,
    range = SA_range
)
ax6.hist(
    mr_post_SA,
    density = True,
    alpha = 0.75,
    bins = bins,
    color = post_color,
    label = 'mined',
    range = SA_range
)
ax6.hist(
    mr_pre_SA,
    density = True,
    alpha = 0.75,
    bins = bins,
    color = pre_color,
    label='unmined',
    range = SA_range
)

# WHITE OAK
#elevation
ax7 = fig2.add_subplot(spec[4:6,:4])
ax7.hist(
    wo_post_elev,
    density = True,
    alpha = 1,
    bins = elev_bins,
    histtype = u'step', 
    color = post_edgecolor,
    lw = 1,
    range = elev_range
)
ax7.hist(
    wo_pre_elev,
    density = True,
    alpha = 1,
    bins = elev_bins,
    histtype = u'step',
    color = pre_edgecolor,
    lw = 1,
    range = elev_range
)
ax7.hist(
    wo_post_elev,
    density = True, 
    alpha = 0.75,
    bins = elev_bins,
    color = post_color,
    label = 'mined',
    range = elev_range
)
ax7.hist(
    wo_pre_elev,
    density = True, 
    alpha = 0.75,
    bins = elev_bins,
    color = pre_color,
    label='unmined',
    range = elev_range
)
ax7.tick_params(axis='both', which='major', labelsize=14)
ax7.set_ylabel('White Oak',fontsize = 12)
#slope
ax8 = fig2.add_subplot(spec[4:6,4:8])
ax8.hist(
    wo_post_slope,
    density = True,
    alpha = 1, 
    bins = bins,
    histtype = u'step', 
    color = post_edgecolor,
    lw = 1,
    range = slope_range
)
ax8.hist(
    wo_pre_slope,
    density = True,
    alpha = 1, 
    bins = bins,
    histtype = u'step',
    color = pre_edgecolor,
    lw = 1,
    range = slope_range
)
ax8.hist(
    wo_post_slope,
    density = True, 
    alpha = 0.75,
    bins = bins,
    color = post_color, 
    label = 'mined',
    range = slope_range
)
ax8.hist(
    wo_pre_slope,
    density = True,
    alpha = 0.75, 
    bins = bins, 
    color = pre_color,
    label='unmined',
    range = slope_range
)
#slope*area
ax9 = fig2.add_subplot(spec[4:6,8:12])
ax9.hist(
    wo_post_SA,
    density = True,
    alpha = 1,
    bins = bins,
    histtype = u'step',
    color = post_edgecolor,
    lw = 1,
    range = SA_range
)
ax9.hist(
    wo_pre_SA,
    density = True,
    alpha = 1,
    bins = bins,
    histtype = u'step',
    color = pre_edgecolor,
    lw = 1,
    range = SA_range
)
ax9.hist(
    wo_post_SA,
    density = True, 
    alpha = 0.75,
    bins = bins, 
    color = post_color,
    label = 'mined',
    range = SA_range
)
ax9.hist(
    wo_pre_SA,
    density = True,
    alpha = 0.75,
    bins = bins,
    color = pre_color,
    label='unmined',
    range = SA_range
)

# SPRUCE FRORK
#elevation
ax10 = fig2.add_subplot(spec[6:8,:4])
ax10.hist(
    sf_post_elev,
    density = True,
    alpha = 1,
    bins = elev_bins,
    histtype = u'step',
    color = post_edgecolor,
    lw = 1,
    range = elev_range
)
ax10.hist(
    sf_pre_elev,
    density = True,
    alpha = 1,
    bins = elev_bins,
    histtype = u'step',
    color = pre_edgecolor,
    lw = 1,
    range = elev_range
)
ax10.hist(
    sf_post_elev,
    density = True,
    alpha = 0.75,
    bins = elev_bins,
    color = post_color,
    label = 'mined',
    range = elev_range
)
ax10.hist(
    sf_pre_elev,
    density = True,
    alpha = 0.75,
    bins = elev_bins,
    color = pre_color,
    label='unmined',
    range = elev_range
)
ax10.tick_params(axis='both', which='major', labelsize=14)
ax10.set_ylabel('Spruce Fork',fontsize = 12)
#slope
ax11 = fig2.add_subplot(spec[6:8,4:8])
ax11.hist(
    sf_post_slope,
    density = True,
    alpha = 1,
    bins = bins,
    histtype = u'step',
    color = post_edgecolor,
    lw = 1,
    range = slope_range
)
ax11.hist(
    sf_pre_slope,
    density = True,
    alpha = 1,
    bins = bins,
    histtype = u'step', 
    color = pre_edgecolor,
    lw = 1,
    range = slope_range
)
ax11.hist(
    sf_post_slope,
    density = True, 
    alpha = 0.75,
    bins = bins,
    color = post_color,
    label = 'mined',
    range = slope_range
)
ax11.hist(
    sf_pre_slope,
    density = True,
    alpha = 0.75,
    bins = bins, 
    color = pre_color,
    label='unmined',
    range = slope_range
)
#slope*area
ax12 = fig2.add_subplot(spec[6:8,8:12])
ax12.hist(
    sf_post_SA,
    density = True,
    alpha = 1,
    bins = bins,
    histtype = u'step',
    color = post_edgecolor,
    lw = 1,
    range = SA_range
)
ax12.hist(
    sf_pre_SA,
    density = True,
    alpha = 1,
    bins = bins,
    histtype = u'step',
    color = pre_edgecolor,
    lw = 1,
    range = SA_range
)
ax12.hist(
    sf_post_SA,
    density = True,
    alpha = 0.75,
    bins = bins,
    color = post_color,
    label = 'mined',
    range = SA_range
)
ax12.hist(
    sf_pre_SA,
    density = True,
    alpha = 0.75,
    bins = bins,
    color = pre_color,
    label='unmined',
    range = SA_range
)

# LAUREL CREEK
#elevation
ax13 = fig2.add_subplot(spec[8:10,:4])
ax13.hist(
    lc_post_elev,
    density = True,
    alpha = 1,
    bins = elev_bins,
    histtype = u'step', 
    color = post_edgecolor,
    lw = 1,
    range = elev_range
)
ax13.hist(
    lc_pre_elev,
    density = True,
    alpha = 1,
    bins = elev_bins,
    histtype = u'step',
    color = pre_edgecolor,
    lw = 1,
    range = elev_range
)
ax13.hist(
    lc_post_elev,
    density = True,
    alpha = 0.75,
    bins = elev_bins,
    color = post_color,
    label = 'mined',
    range = elev_range
)
ax13.hist(
    lc_pre_elev,
    density = True,
    alpha = 0.75,
    bins = elev_bins,
    color = pre_color,
    label='unmined',
    range = elev_range
)
ax13.tick_params(axis='both', which='major', labelsize=12)
ax13.set_ylabel('Laurel Creek',fontsize = 12)
ax13.set_xlabel('$Elevation (m)$')
#slope
ax14 = fig2.add_subplot(spec[8:10,4:8])
ax14.hist(
    lc_post_slope,
    density = True,
    alpha = 1,
    bins = bins,
    histtype = u'step', 
    color = post_edgecolor,
    lw = 1,
    range = slope_range
)
ax14.hist(
    lc_pre_slope,
    density = True,
    alpha = 1, 
    bins = bins,
    histtype = u'step',
    color = pre_edgecolor,
    lw = 1,
    range = slope_range
)
ax14.hist(
    lc_post_slope,
    density = True,
    alpha = 0.75,
    bins = bins,
    color = post_color,
    label = 'mined',
    range = slope_range
)
ax14.hist(
    lc_pre_slope,
    density = True,
    alpha = 0.75,
    bins = bins,
    color = pre_color,
    label='unmined',
    range = slope_range
)
ax14.tick_params(axis='both', which='major', labelsize=12)
ax14.set_xlabel('$slope (\dfrac{m}{m})$')
#slope*area
ax15 = fig2.add_subplot(spec[8:10,8:12])
ax15.hist(
    lc_post_SA,
    density = True,
    alpha = 1,
    bins = bins,
    histtype = u'step',
    color = post_edgecolor,
    lw = 1,
    range = SA_range
)
ax15.hist(
    lc_pre_SA,
    density = True,
    alpha = 1,
    bins = bins,
    histtype = u'step',
    color = pre_edgecolor,
    lw = 1,
    range = SA_range
)
ax15.hist(
    lc_post_SA,
    density = True, 
    alpha = 0.75, 
    bins = bins, 
    color = post_color,
    label = 'mined',
    range = SA_range
)
ax15.hist(
    lc_pre_SA,
    density = True, 
    alpha = 0.75, 
    bins = bins, 
    color = pre_color,
    label='unmined',
    range = SA_range
)
ax15.tick_params(axis='both', which='major', labelsize=12)
ax15.set_xlabel('$\sqrt{A}S (m)$')


#Modify axes tick labels, and add statistics
axes_lst = [ax1,ax2,ax3,
            ax4,ax5,ax6,
            ax7,ax8,ax9,
            ax10,ax11,ax12,
            ax13,ax14,ax15]

for i in range(15):
    
    if i < 12:
        #remove x ticks from all grids above the bottom row
        axes_lst[i].set_xticks([]) 

    if i%3 !=0:
        # remove y ticks from all grids right of the first column
        axes_lst[i].set_yticks([])
        #indent the following line to add y-axis labels
    axes_lst[i].set_yticklabels([]) 
    

    # modify subplot seperately [2,0] due to weird shape
    if i == 3:
        #set axes limits 
        axes_lst[i].set_xlim(elev_range)
        axes_lst[i].set_ylim(y_elev_range)
        # create inset axis with data coordinates 
        SA_bw_axes = axes_lst[i].inset_axes(
            [450, 0.0060, 225, 0.003],
            transform = axes_lst[i].transData
        )
        #remove axes ticks and frame and make inset transparent
        SA_bw_axes.set_yticks([])
        SA_bw_axes.set_xticks([])
        SA_bw_axes.spines['left'].set_visible(False)
        SA_bw_axes.spines['right'].set_visible(False)
        SA_bw_axes.spines['top'].set_visible(False)
        SA_bw_axes.patch.set_alpha(0)
        #plot kernal density of bayesian wilcoxon results
        sbn.kdeplot(
            bayesian_wilcoxon_density[i],
            ax = SA_bw_axes,
            c = 'k',
            lw = 1.25,
            zorder = 0)
        # plot the bounds of the 99% confidence range as points
        # along the x-axis, and label with 'annotate'
        SA_bw_axes.scatter(
            [bayesian_stats[i].values_column[2],
             bayesian_stats[i].values_column[3]],
            [0,0],
            edgecolor = 'k',
            facecolor = 'white',
            clip_on=False,
            zorder = 4
        )
        SA_bw_axes.annotate(
            str(np.round(bayesian_stats[i].values_column[2],2)),
            (bayesian_stats[i].values_column[2],0),
            textcoords = 'offset points',
            xytext = (0,-12),
            ha = 'center'
        )
        SA_bw_axes.annotate(
            str(np.round(bayesian_stats[i].values_column[3],2)),
            (bayesian_stats[i].values_column[3],0),
            textcoords = 'offset points',
            xytext = (0,-12),
            ha = 'center'
        )


    if (i == 0 or i == 6 or i == 9 or i== 12):
        
    #set axes limits
        axes_lst[i].set_xlim(elev_range)
        axes_lst[i].set_ylim(y_elev_range)
        SA_bw_axes = axes_lst[i].inset_axes(
            [225, 0.0060, 225, 0.003],
            transform = axes_lst[i].transData
        )
        SA_bw_axes.set_yticks([])
        SA_bw_axes.set_xticks([])
        SA_bw_axes.patch.set_alpha(0)
        SA_bw_axes.spines['left'].set_visible(False)
        SA_bw_axes.spines['right'].set_visible(False)
        SA_bw_axes.spines['top'].set_visible(False)
        sbn.kdeplot(
            bayesian_wilcoxon_density[i],
            ax = SA_bw_axes,
            c = 'k',
            lw = 1.25,
            zorder = 0
        )
        SA_bw_axes.scatter(
            [bayesian_stats[i].values_column[2],
             bayesian_stats[i].values_column[3]],
            [0,0],
            edgecolor = 'k',
            facecolor = 'white',
            clip_on=False,
            zorder = 4
        )
        SA_bw_axes.annotate(
            str(np.round(bayesian_stats[i].values_column[2],2)),
            (bayesian_stats[i].values_column[2],0),
            textcoords = 'offset points',
            xytext = (0,-12),
            ha = 'center'
        )
        SA_bw_axes.annotate(
            str(np.round(bayesian_stats[i].values_column[3],2)),
            (bayesian_stats[i].values_column[3],0),
            textcoords = 'offset points',
            xytext = (0,-12),
            ha = 'center'
        )
  

    if (i == 1 or i== 4 or i == 7 or i == 10 or i== 13):
    
        axes_lst[i].set_xlim(slope_range)
        axes_lst[i].set_ylim(y_slope_range)

        SA_bw_axes = axes_lst[i].inset_axes(
            [0.75, 2.0, .68, 1.25], 
            transform = axes_lst[i].transData
        )
        SA_bw_axes.set_yticks([])
        SA_bw_axes.set_xticks([])
        SA_bw_axes.patch.set_alpha(0)
        SA_bw_axes.spines['left'].set_visible(False)
        SA_bw_axes.spines['right'].set_visible(False)
        SA_bw_axes.spines['top'].set_visible(False)
        sbn.kdeplot(
            bayesian_wilcoxon_density[i],
            ax = SA_bw_axes,
            c = 'k',
            lw = 1.25,
            zorder = 0
        )
        SA_bw_axes.scatter(
            [bayesian_stats[i].values_column[2],
             bayesian_stats[i].values_column[3]],
            [0,0],
            edgecolor = 'k',
            facecolor = 'white',
            clip_on=False,
            zorder = 4
        )
        SA_bw_axes.annotate(
            str(np.round(bayesian_stats[i].values_column[2],2)),
            (bayesian_stats[i].values_column[2],0),
            textcoords = 'offset points',
            xytext = (0,-12),
            ha = 'center'
        )

        SA_bw_axes.annotate(
            str(np.round(bayesian_stats[i].values_column[3],2)),
            (bayesian_stats[i].values_column[3],0),
            textcoords = 'offset points',
            xytext = (0,-12),
            ha = 'center'
        )
                         
    if (i == 2 or i== 5 or i == 8 or i == 11 or i== 14):

        axes_lst[i].set_xlim(SA_range)
        axes_lst[i].set_ylim(y_SA_range)
        
        SA_bw_axes = axes_lst[i].inset_axes(
            [25, 0.085, 30, 0.05],
            transform = axes_lst[i].transData
        )
        SA_bw_axes.set_yticks([])
        SA_bw_axes.set_xticks([])
        SA_bw_axes.patch.set_alpha(0)
        SA_bw_axes.spines['left'].set_visible(False)
        SA_bw_axes.spines['right'].set_visible(False)
        SA_bw_axes.spines['top'].set_visible(False)
        sbn.kdeplot(
            bayesian_wilcoxon_density[i],
            ax = SA_bw_axes,
            c = 'k',
            lw = 1.25,
            zorder = 0
        )
        SA_bw_axes.scatter(
            [bayesian_stats[i].values_column[2],
             bayesian_stats[i].values_column[3]],
            [0,0],
            edgecolor = 'k',
            facecolor = 'white',
            clip_on=False,
            zorder = 4
        )
        SA_bw_axes.annotate(
            str(np.round(bayesian_stats[i].values_column[2],2)),
            (bayesian_stats[i].values_column[2],0),
            textcoords = 'offset points',
            xytext = (0,-12),
            ha = 'center'
        )
        SA_bw_axes.annotate(
            str(np.round(bayesian_stats[i].values_column[3],2)),
            (bayesian_stats[i].values_column[3],0),
            textcoords = 'offset points',
            xytext = (0,-12),
            ha = 'center'
        )

fig2.savefig('figure_2.png', dpi = 1000)