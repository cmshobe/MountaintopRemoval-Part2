########################################################################
#This script generates data for Figure 2 in the following manuscript:

#Bower, S.J., Shobe, C.M., Maxwell, A.E., and Campforts, B. (2023) The 
#uncertain future of mountaintop-removal-mined landscapes 2: Modeling the influence of
#topography and vegetation. Geomorphology.

#Please cite the paper if you use this code in any way.

#Brief description: this script extracts distributions of elevation, slope, and 
#area--slope product from pre- and post-mining DEMs of five study watersheds.

########################################################################

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats
from scipy.stats import norm
import landlab
import seaborn as sbn
from landlab.plot import imshow_grid
from landlab import RasterModelGrid
from landlab.io.esri_ascii import read_esri_ascii
from landlab.components import (PriorityFloodFlowRouter,
                                DepressionFinderAndRouter,
                                FlowAccumulator
                                )

#For archiving simplicity, DEMs are located in the `model_simulations' folder.
input_path = '../model_simulations/'


###########calculate morphometrics for White Oak Creek
#PRE
wo_mg_pre, wo_z_pre = read_esri_ascii(input_path+'whiteoak/inputs/whiteoak_pre_10m.asc', name='topographic__elevation')
np.all(wo_mg_pre.at_node['topographic__elevation'] == wo_z_pre)
wo_z_pre[wo_z_pre < 1] = -9999
wo_mg_pre.set_closed_boundaries_at_grid_edges(True,True,True,True)
wo_mg_pre.set_nodata_nodes_to_closed(wo_z_pre, -9999) 
wo_mg_pre.set_watershed_boundary_condition(wo_z_pre, nodata_value = -9999, return_outlet_id=True)

wo_df_pre = DepressionFinderAndRouter(wo_mg_pre)
wo_fa_pre = PriorityFloodFlowRouter(wo_mg_pre,'topographic__elevation',flow_metric = 'D8')

#POST
wo_mg_post, wo_z_post = read_esri_ascii(input_path+'whiteoak/inputs/whiteoak_post_10m.asc', name='topographic__elevation')
np.all(wo_mg_post.at_node['topographic__elevation'] == wo_z_post)
wo_z_post[wo_z_post < 1] = -9999
wo_mg_post.set_closed_boundaries_at_grid_edges(True,True,True,True)
wo_mg_post.set_nodata_nodes_to_closed(wo_z_post, -9999) 
wo_mg_post.set_watershed_boundary_condition(wo_z_post, nodata_value = -9999, return_outlet_id=True)

wo_df_post = DepressionFinderAndRouter(wo_mg_post)
wo_fa_post = PriorityFloodFlowRouter(wo_mg_post,'topographic__elevation',flow_metric = 'D8')

wo_fa_pre.run_one_step()
wo_fa_post.run_one_step()
wo_mg_pre.calc_slope_at_node()
wo_mg_post.calc_slope_at_node()


wo_pre_elev = wo_mg_pre.at_node['topographic__elevation'][wo_mg_pre.core_nodes]
wo_post_elev = wo_mg_post.at_node['topographic__elevation'][wo_mg_post.core_nodes]
wo_pre_slope = wo_mg_pre.at_node['topographic__steepest_slope'][wo_mg_pre.core_nodes]
wo_post_slope = wo_mg_post.at_node['topographic__steepest_slope'][wo_mg_post.core_nodes]
wo_pre_DA = wo_mg_pre.at_node['drainage_area'][wo_mg_pre.core_nodes]
wo_post_DA = wo_mg_post.at_node['drainage_area'][wo_mg_post.core_nodes]

wo_pre_SA = wo_pre_slope * wo_pre_DA**0.5
wo_post_SA = wo_post_slope * wo_post_DA**0.5


###########calculate morphometrics for Laurel Creek
#PRE
lc_mg_pre, lc_z_pre = read_esri_ascii(input_path+'laurelcreek/inputs/laurelcreek_pre_10m.asc', name='topographic__elevation')
np.all(lc_mg_pre.at_node['topographic__elevation'] == lc_z_pre)
lc_z_pre[lc_z_pre < 1] = -9999
lc_mg_pre.set_closed_boundaries_at_grid_edges(True,True,True,True)
lc_mg_pre.set_nodata_nodes_to_closed(lc_z_pre, -9999) 
lc_mg_pre.set_watershed_boundary_condition(lc_z_pre, nodata_value = -9999, return_outlet_id=True)

lc_df_pre = DepressionFinderAndRouter(lc_mg_pre)
lc_fa_pre = PriorityFloodFlowRouter(lc_mg_pre,'topographic__elevation',flow_metric = 'D8')

#POST
lc_mg_post, lc_z_post = read_esri_ascii(input_path+'laurelcreek/inputs/laurelcreek_post_10m.asc', name='topographic__elevation')
np.all(lc_mg_post.at_node['topographic__elevation'] == lc_z_post)
lc_z_post[lc_z_post < 1] = -9999
lc_mg_post.set_closed_boundaries_at_grid_edges(True,True,True,True)
lc_mg_post.set_nodata_nodes_to_closed(lc_z_post, -9999) 
lc_mg_post.set_watershed_boundary_condition(lc_z_post, nodata_value = -9999, return_outlet_id=True)

lc_df_post = DepressionFinderAndRouter(lc_mg_post)
lc_fa_post = PriorityFloodFlowRouter(lc_mg_post,'topographic__elevation',flow_metric = 'D8')

lc_fa_pre.run_one_step()
lc_fa_post.run_one_step()
lc_mg_pre.calc_slope_at_node()
lc_mg_post.calc_slope_at_node()


lc_pre_elev = lc_mg_pre.at_node['topographic__elevation'][lc_mg_pre.core_nodes]
lc_post_elev = lc_mg_post.at_node['topographic__elevation'][lc_mg_post.core_nodes]
lc_pre_slope = lc_mg_pre.at_node['topographic__steepest_slope'][lc_mg_pre.core_nodes]
lc_post_slope = lc_mg_post.at_node['topographic__steepest_slope'][lc_mg_post.core_nodes]
lc_pre_DA = lc_mg_pre.at_node['drainage_area'][lc_mg_pre.core_nodes]
lc_post_DA = lc_mg_post.at_node['drainage_area'][lc_mg_post.core_nodes]

lc_pre_SA = lc_pre_slope * lc_pre_DA**0.5
lc_post_SA = lc_post_slope * lc_post_DA**0.5


###########calculate morphometrics for Mud River
#PRE
mr_mg_pre, mr_z_pre = read_esri_ascii(input_path+'mudriver/inputs/mudriver_pre_10m.asc', name='topographic__elevation')
np.all(mr_mg_pre.at_node['topographic__elevation'] == mr_z_pre)
mr_z_pre[mr_z_pre < 1] = -9999
mr_mg_pre.set_closed_boundaries_at_grid_edges(True,True,True,True)
mr_mg_pre.set_nodata_nodes_to_closed(mr_z_pre, -9999) 
mr_mg_pre.set_watershed_boundary_condition(mr_z_pre, nodata_value = -9999, return_outlet_id=True)

mr_df_pre = DepressionFinderAndRouter(mr_mg_pre)
mr_fa_pre = PriorityFloodFlowRouter(mr_mg_pre,'topographic__elevation',flow_metric = 'D8')

#POST
mr_mg_post, mr_z_post = read_esri_ascii(input_path+'mudriver/inputs/mudriver_post_10m.asc', name='topographic__elevation')
np.all(mr_mg_post.at_node['topographic__elevation'] == mr_z_post)
mr_z_post[mr_z_post < 1] = -9999
mr_mg_post.set_closed_boundaries_at_grid_edges(True,True,True,True)
mr_mg_post.set_nodata_nodes_to_closed(mr_z_post, -9999) 
mr_mg_post.set_watershed_boundary_condition(mr_z_post, nodata_value = -9999, return_outlet_id=True)

mr_df_post = DepressionFinderAndRouter(mr_mg_post)
mr_fa_post = PriorityFloodFlowRouter(mr_mg_post,'topographic__elevation',flow_metric = 'D8')

mr_fa_pre.run_one_step()
mr_fa_post.run_one_step()
mr_mg_pre.calc_slope_at_node()
mr_mg_post.calc_slope_at_node()

mr_pre_elev = mr_mg_pre.at_node['topographic__elevation'][mr_mg_pre.core_nodes]
mr_post_elev = mr_mg_post.at_node['topographic__elevation'][mr_mg_post.core_nodes]
mr_pre_slope = mr_mg_pre.at_node['topographic__steepest_slope'][mr_mg_pre.core_nodes]
mr_post_slope = mr_mg_post.at_node['topographic__steepest_slope'][mr_mg_post.core_nodes]
mr_pre_DA = mr_mg_pre.at_node['drainage_area'][mr_mg_pre.core_nodes]
mr_post_DA = mr_mg_post.at_node['drainage_area'][mr_mg_post.core_nodes]

mr_pre_SA = mr_pre_slope * mr_pre_DA**0.5
mr_post_SA = mr_post_slope * mr_post_DA**0.5

###########calculate morphometrics for Spruce Fork
#PRE
sf_mg_pre, sf_z_pre = read_esri_ascii(input_path+'sprucefork/inputs/sprucefork_pre_10m.asc', name='topographic__elevation')
np.all(sf_mg_pre.at_node['topographic__elevation'] == sf_z_pre)
sf_z_pre[sf_z_pre < 1] = -9999
sf_mg_pre.set_closed_boundaries_at_grid_edges(True,True,True,True)
sf_mg_pre.set_nodata_nodes_to_closed(sf_z_pre, -9999) 
sf_mg_pre.set_watershed_boundary_condition(sf_z_pre, nodata_value = -9999, return_outlet_id=True)

sf_df_pre = DepressionFinderAndRouter(sf_mg_pre)
sf_fa_pre = PriorityFloodFlowRouter(sf_mg_pre,'topographic__elevation',flow_metric = 'D8')

#POST
sf_mg_post, sf_z_post = read_esri_ascii(input_path+'sprucefork/inputs/sprucefork_post_10m.asc', name='topographic__elevation')
np.all(sf_mg_post.at_node['topographic__elevation'] == sf_z_post)
sf_z_post[sf_z_post < 1] = -9999
sf_mg_post.set_closed_boundaries_at_grid_edges(True,True,True,True)
sf_mg_post.set_nodata_nodes_to_closed(sf_z_post, -9999) 
sf_mg_post.set_watershed_boundary_condition(sf_z_post, nodata_value = -9999, return_outlet_id=True)

sf_df_post = DepressionFinderAndRouter(sf_mg_post)
sf_fa_post = PriorityFloodFlowRouter(sf_mg_post,'topographic__elevation',flow_metric = 'D8')

sf_fa_pre.run_one_step()
sf_fa_post.run_one_step()
sf_mg_pre.calc_slope_at_node()
sf_mg_post.calc_slope_at_node()


sf_pre_elev = sf_mg_pre.at_node['topographic__elevation'][sf_mg_pre.core_nodes]
sf_post_elev = sf_mg_post.at_node['topographic__elevation'][sf_mg_pre.core_nodes]
sf_pre_slope = sf_mg_pre.at_node['topographic__steepest_slope'][sf_mg_pre.core_nodes]
sf_post_slope = sf_mg_post.at_node['topographic__steepest_slope'][sf_mg_pre.core_nodes]
sf_pre_DA = sf_mg_pre.at_node['drainage_area'][sf_mg_pre.core_nodes]
sf_post_DA = sf_mg_post.at_node['drainage_area'][sf_mg_pre.core_nodes]

sf_pre_SA = sf_pre_slope * sf_pre_DA**0.5
sf_post_SA = sf_post_slope * sf_post_DA**0.5

###########calculate morphometrics for Ben Creek
#PRE
bc_mg_pre, bc_z_pre = read_esri_ascii(input_path+'bencreek/inputs/bencreek_pre_10m.asc', name='topographic__elevation')
np.all(bc_mg_pre.at_node['topographic__elevation'] == bc_z_pre)
bc_z_pre[bc_z_pre < 1] = -9999
bc_mg_pre.set_closed_boundaries_at_grid_edges(True,True,True,True)
bc_mg_pre.set_nodata_nodes_to_closed(bc_z_pre, -9999) 
bc_mg_pre.set_watershed_boundary_condition(bc_z_pre, nodata_value = -9999, return_outlet_id=True)

bc_df_pre = DepressionFinderAndRouter(bc_mg_pre)
bc_fa_pre = PriorityFloodFlowRouter(bc_mg_pre,'topographic__elevation',flow_metric = 'D8')

#POST
bc_mg_post, bc_z_post = read_esri_ascii(input_path+'bencreek/inputs/bencreek_post_10m.asc', name='topographic__elevation')
np.all(bc_mg_post.at_node['topographic__elevation'] == bc_z_post)
bc_z_post[bc_z_post < 1] = -9999
bc_mg_post.set_closed_boundaries_at_grid_edges(True,True,True,True)
bc_mg_post.set_nodata_nodes_to_closed(bc_z_post, -9999) 
bc_mg_post.set_watershed_boundary_condition(bc_z_post, nodata_value = -9999, return_outlet_id=True)

bc_df_post = DepressionFinderAndRouter(bc_mg_post)
bc_fa_post = PriorityFloodFlowRouter(bc_mg_post,'topographic__elevation',flow_metric = 'D8')

bc_fa_pre.run_one_step()
bc_fa_post.run_one_step()
bc_mg_pre.calc_slope_at_node()
bc_mg_post.calc_slope_at_node()

bc_pre_elev = bc_mg_pre.at_node['topographic__elevation'][bc_mg_pre.core_nodes]
bc_post_elev = bc_mg_post.at_node['topographic__elevation'][bc_mg_post.core_nodes]
bc_pre_slope = bc_mg_pre.at_node['topographic__steepest_slope'][bc_mg_pre.core_nodes]
bc_post_slope = bc_mg_post.at_node['topographic__steepest_slope'][bc_mg_post.core_nodes]
bc_pre_DA = bc_mg_pre.at_node['drainage_area'][bc_mg_pre.core_nodes]
bc_post_DA = bc_mg_post.at_node['drainage_area'][bc_mg_post.core_nodes]

bc_pre_SA = bc_pre_slope * bc_pre_DA**0.5
bc_post_SA = bc_post_slope * bc_post_DA**0.5


#output flatted distributions of elevation, slope, and area--slope product
flat_data_input_path = '/flat_inputs/'
datas = [bc_pre_elev,bc_post_elev,bc_pre_slope,bc_post_slope,bc_pre_SA,bc_post_SA,
        mr_pre_elev,mr_post_elev,mr_pre_slope,mr_post_slope,mr_pre_SA,mr_post_SA,
        wo_pre_elev,wo_post_elev,wo_pre_slope,wo_post_slope,wo_pre_SA,wo_post_SA,
        sf_pre_elev,sf_post_elev,sf_pre_slope,sf_post_slope,sf_pre_SA,sf_post_SA,
        lc_pre_elev,lc_post_elev,lc_pre_slope,lc_post_slope,lc_pre_SA,lc_post_SA]

paths = ['bc_pre_elev','bc_post_elev','bc_pre_slope','bc_post_slope','bc_pre_SA','bc_post_SA',
        'mr_pre_elev','mr_post_elev','mr_pre_slope','mr_post_slope','mr_pre_SA','mr_post_SA',
        'wo_pre_elev','wo_post_elev','wo_pre_slope','wo_post_slope','wo_pre_SA','wo_post_SA',
        'sf_pre_elev','sf_post_elev','sf_pre_slope','sf_post_slope','sf_pre_SA','sf_post_SA',
        'lc_pre_elev','lc_post_elev','lc_pre_slope','lc_post_slope','lc_pre_SA','lc_post_SA']

for i in range(0,len(datas)):
    np.savetxt(flat_data_input_path+paths[i]+'.txt',datas[i])