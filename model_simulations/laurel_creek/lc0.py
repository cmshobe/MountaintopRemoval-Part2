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
    output_path = '/scratch/cs00048/MTR/wsheds/laurel_creek/outputs/'
    scenario_str = 'scenario_'+str(scenario_id)+'_'
    #load DEM
    mg, z = read_esri_ascii(path+'laurelcreek_post_10m.asc', name='topographic__elevation')
    np.all(mg.at_node['topographic__elevation'] == z)
    mg.set_closed_boundaries_at_grid_edges(True,True,True,True)
    mg.set_nodata_nodes_to_closed(z, -99999)
    outlet_node = mg.set_watershed_boundary_condition(z,nodata_value = -99999,return_outlet_id=True)
    print('outlet: ' + str(outlet_node))
    #Load rasterized Skytruth data
    k_array, ksp = read_esri_ascii(path+'LC_K_mask_edited.asc', name='k')
    k_array.set_nodata_nodes_to_closed(ksp, -9999)


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
    climate_data = np.round((climate_data[2,2:].astype('float')/1000),2)
    
    #charlie modified for changing recovery dt:
    precip_forcing = np.zeros(len(K_measured)+ 1)
    
    precip_forcing[0:len(climate_data)], precip_forcing[len(climate_data):] = climate_data, climate_data[-1]
    K_modeled = K_measured[:,scenario]* (v_s/np.min(precip_forcing) + 1)
    K_unmined = np.min(K_measured)* (v_s/np.min(precip_forcing) + 1)
    
    #charlie modified for changing recovery dt:
    K_max = K_modeled[int(recovery_period/1)]
    shape = ksp
    scape = np.zeros(len(shape),dtype = float)
    
    #charlie modified for changing recovery dt:
    K_grids = np.zeros((len(scape),int(recovery_period/1)),dtype = float)
    K_grid_recovered = np.zeros(len(scape),dtype = 'float')

    scape[shape.all() < 1 and shape.all() > -9999] = K_unmined
    scape[shape == 1] = K_max
    K_grid_recovered = scape
        
    count = 0
    #charlie modified for changing recovery dt:
    for i in np.arange(0,int(recovery_period/1),1):
        scape[shape.all() < 1 and shape.all() > -9999] = K_unmined #set K in unmined areas
        scape[shape == 1] = K_modeled[count] #set K in mined areas
        K_grids[:,i] = scape
        count += 1
        
        #Instantiate components for post-mining model with K recovery curve
    print(
        "K IC no-recovery:\nmined = %s,\nunmined = %s,\nFC no-recovery:\nmined = %s,\nunmined = %s"%(
            np.max(K_grids[:,0]),
            np.min(K_grids[:,0]),
            np.max(K_grids[:,-1]),
            np.min(K_grids[:,-1])
        )
    )

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
    print('starting SPACE on full recovery')

        
    while timer < recovery_period:
    
        pre_dt_topo = copy.deepcopy(mg.at_node['topographic__elevation'])
            
        #to allow timesteps smaller than one, I need to access the same
        #element of the K grid 1/dt times (i.e. twice if dt = 0.5yr)
        mg.at_node['erodibility'] = K_grids[:, int(np.floor(timer))]
        print('no recovery K is: ', np.max(K_grids[:, int(np.floor(timer))]))
            
        ########K CONVOLUTION
            
        #first, assign a K value to the closed nodes
        #this is any node outside the basin
        #that K value will be 10^ (the average of the log10 of Kmin and Kmax)
        outside_K_value = 10 ** (np.mean(np.array([np.log10(np.min(K_grids[:, int(np.floor(timer))])), np.log10(np.max(K_grids[:, int(np.floor(timer))]))])))

        #set closed nodes to have outside K val
        mg.at_node['erodibility'][mg.closed_boundary_nodes] = outside_K_value
            
        #define the grid to be worked on
        A = mg.at_node['erodibility'].reshape((mg.number_of_node_rows, mg.number_of_node_columns))

        #define the convolution kernel (simple moving-window average in this case
        B = np.ones((window_size, window_size)) / np.power(window_size, 2)

        #do the convolution
        logged_window_average = 10 ** convolve2d(np.log10(A), B, mode = 'same', boundary = 'fill', fillvalue = np.log10(outside_K_value))
        mg.at_node['erodibility'][:] = logged_window_average.flatten()
            
        ########END K CONVOLUTION

        #instantiate SPACE
        SPACE = SpaceLargeScaleEroder(
            mg,
            K_sed = mg.at_node['erodibility'],
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
        print('Laurel Creek ',scenario, 'recovery: ', timer)

    #K grid deletion, K convolution, and SPACE instantiation
    #moved out of the post-recovery loop for speed
        
    try:
        del(K_grids)
        print("k_grids deleted")
    except NameError:
        pass
        
    #instantiate SPACE
    SPACE = SpaceLargeScaleEroder(
        mg,
        K_sed = mg.at_node['erodibility'],
        K_br = K_br,
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
        print('Laurel Creek ',scenario, 'fast dt: ', timer)

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
    #print("FINISHED Laurel Creek SCENARIO ID ",id)