######################################

#This model simulates the evolution of mountaintop-removal-mined landscapes 
#that are subject to a vegetation recovery trajectory

#Written by Sam Bower and Charlie Shobe, 2022-2023

#Please cite Bower et al., in review

######################################


def landscape_evolution(scenario_id):
    
    import sys
    import time
    import copy
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import landlab
    from landlab.plot import imshow_grid, imshowhs_grid
    from landlab import RasterModelGrid
    from landlab.io.esri_ascii import read_esri_ascii
    from landlab.components import (
        PriorityFloodFlowRouter,
        LinearDiffuser,
        SpaceLargeScaleEroder,
        ChannelProfiler)
    from scipy.signal import convolve2d

    path = 'inputs/'
    output_path = '/scratch/cs00048/MTR/wsheds/white_oak/outputs/'
    scenario_str = 'control_unmined_'
    #load DEM
    mg, z = read_esri_ascii(path+'whiteoak_pre_10m.asc', name='topographic__elevation')
    np.all(mg.at_node['topographic__elevation'] == z)
    mg.set_closed_boundaries_at_grid_edges(True,True,True,True)
    mg.set_nodata_nodes_to_closed(z, -99999)
    outlet_node = mg.set_watershed_boundary_condition(z,nodata_value = -99999,return_outlet_id=True)
    print('outlet: ' + str(outlet_node))

    _ = mg.add_zeros("soil__depth", at="node")
    mg.at_node["soil__depth"] += 500#z
    _ = mg.add_zeros("bedrock__elevation", at = 'node')
    mg.at_node["bedrock__elevation"] = mg.at_node["topographic__elevation"] - mg.at_node["soil__depth"]

    mg.at_node['topographic__elevation'][outlet_node[0]]

    total_time = 10000
    recovery_period = 200
    recovery_dt = 0.5 
    fast_dt = 1
    D = 3e-3
    K_br = 0
    m_sp = 0.5
    n_sp = 1
    F_f = 0
    phi = 0.3
    H_star = 1.
    v_s = 0.01
    sp_crit_sed = 0
    sp_crit_br = 0
    BL_lowering_rate = 2.7e-5
    window_size = 9 #this is the total size of the moving window for K convolution


    scenario = scenario_id #minimum vegetation recovery
    K_measured = np.loadtxt(path+'K_timeseries_scaled_CMS.txt')

    #CLIMATE TIMESERIES same length as K TIMESERIES
    climate_data = pd.read_csv(path+'yearly_data.csv').to_numpy()
    climate_data = np.round((climate_data[1,2:].astype('float')/1000),2)
    
    #charlie modified for changing recovery dt:
    precip_forcing = np.zeros(len(K_measured)+ 1)
    
    precip_forcing[0:len(climate_data)], precip_forcing[len(climate_data):] = climate_data, climate_data[-1]
    K_modeled = K_measured[:,scenario]* (v_s/np.min(precip_forcing) + 1)
    K_unmined = np.min(K_measured)* (v_s/np.min(precip_forcing) + 1)
    
    #charlie modified for changing recovery dt:
    K_max = K_modeled[int(recovery_period/1)]
        
    count = 0

    FA = PriorityFloodFlowRouter(mg, flow_metric="D8", suppress_out=True)
    FA.run_one_step()
    LD = LinearDiffuser(mg, linear_diffusivity = D) 
    mg.at_node["soil__depth"] = mg.at_node['topographic__elevation']-mg.at_node["bedrock__elevation"]
    mg.at_node["surface_water__discharge"] = mg.at_node["drainage_area"] * precip_forcing[0]


    #Make arrays for plotting sediment flux, main channel profile, and difference DEM
    #deepcopy of topographic elevation for differencing
    z1 = copy.deepcopy(z)

    #instantiate array to hold sediment flux outlet node
    sed_flux_via_diff = np.zeros(int(((recovery_period / recovery_dt)+((total_time - recovery_period)/fast_dt))))

    timer = 0
    counter = 0
    print('starting SPACE control mined')

    SPACE = SpaceLargeScaleEroder(
        mg,
        K_sed = K_unmined,
        K_br = 0,
        F_f = F_f,
        phi = phi,
        H_star = H_star,
        v_s = v_s,
        m_sp = m_sp,
        n_sp = n_sp,
        sp_crit_sed = sp_crit_sed,
        sp_crit_br = sp_crit_sed,
        erode_flooded_nodes = False
    )
        
    while timer < recovery_period:
    
        pre_dt_topo = copy.deepcopy(mg.at_node['topographic__elevation'])

        FA.run_one_step()
        FA.remove_depressions()
        mg.at_node["flood_status_code"] = np.where(mg.at_node["depression_free_elevation"] == mg.at_node["topographic__elevation"],0,3)
        mg.at_node["surface_water__discharge"] = mg.at_node["drainage_area"] * precip_forcing[int(np.floor(timer))]
        SPACE.run_one_step(dt = recovery_dt)
        LD.run_one_step(dt = recovery_dt)
            
        #we need to account for the fact that LD only modifies topographic__elevation
        #by recalculating soil depth as f(topo and br), which means it's
        #easiest to lower TOPO at baselevel which should have the same
        #effect as what Sam originally did, which was to lower SOIL
        mg.at_node['topographic__elevation'][outlet_node] -= BL_lowering_rate * recovery_dt
        mg.at_node["soil__depth"] = mg.at_node["topographic__elevation"] - mg.at_node["bedrock__elevation"] #CMS change
            
            
        #save sediment flux to array
        sed_flux_via_diff[counter] = np.sum(pre_dt_topo[mg.core_nodes] - mg.at_node['topographic__elevation'][mg.core_nodes]) / recovery_dt

        timer += recovery_dt
            
        #cleaning up machine precision decimal errors in timekeeping
        timer = np.round(timer, decimals = 2)
        counter += 1
        print('White Oak ',scenario, 'recovery: ', timer)


    while timer < total_time:#else:
    
        pre_dt_topo = copy.deepcopy(mg.at_node['topographic__elevation'])

        FA.run_one_step()
        FA.remove_depressions()
        mg.at_node["flood_status_code"] = np.where(mg.at_node["depression_free_elevation"] == mg.at_node["topographic__elevation"],0,3)
        mg.at_node["surface_water__discharge"] = mg.at_node["drainage_area"] * precip_forcing[int(timer)]
        SPACE.run_one_step(dt = fast_dt)
        LD.run_one_step(dt = fast_dt)
            
        #we need to account for the fact that LD only modifies topographic__elevation
        #by recalculating soil depth as f(topo and br), which means it's
        #easiest to lower TOPO at baselevel which should have the same
        #effect as what Sam originally did, which was to lower SOIL
        mg.at_node['topographic__elevation'][outlet_node] -= BL_lowering_rate * fast_dt
        mg.at_node["soil__depth"] = mg.at_node["topographic__elevation"] - mg.at_node["bedrock__elevation"] #CMS change

            
        sed_flux_via_diff[counter] = np.sum(pre_dt_topo[mg.core_nodes] - mg.at_node['topographic__elevation'][mg.core_nodes]) / fast_dt
        timer += fast_dt
        counter += 1
        print('White Oak ',scenario, 'fast dt: ', timer)

	#save modeled topography
    mg.save(output_path+scenario_str+'elev_01.asc',names = 'topographic__elevation')
    mg.save(output_path+scenario_str+'DA_01.asc',names = 'drainage_area')
    mg.save(output_path+scenario_str+'slope_01.asc',names = 'topographic__steepest_slope')

    #save sediment flux
    np.savetxt(output_path+scenario_str+'sedflux.txt', sed_flux_via_diff)

#now run it
import numpy as np

#you can run it with a list of scenario numbers, or with a single scenario.
#For the former uncomment the loop below and pass in the list_of_scenarios.
list_of_scenarios = np.arange(0,11,1)

#for id in list_of_scenarios:
landscape_evolution(0)
    #print("FINISHED WHITE OAK SCENARIO ID ",id)