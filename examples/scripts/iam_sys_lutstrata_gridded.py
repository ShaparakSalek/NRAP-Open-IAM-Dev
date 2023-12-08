"""
Example demonstrating the use of the LookupTableStratigraphy component. Here, 
the component is used to interpolate unit thicknesses and depths at many 
locations. The x and y coordinates are arranged in a 2-dimensional grid. The 
LookupTableStratigraphy component can take the x and y grids and produce unit 
thicknesses and depths at each grid node.

Example of run:
$ python iam_sys_lutstrata_gridded.py
"""

import sys
import os
import logging
import numpy as np
from numpy import matlib
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, os.sep.join(['..', '..', 'source']))

from openiam import (SystemModel, StratigraphyDataInterpolator,
                     LookupTableStratigraphy)


if __name__ == '__main__':
    logging.basicConfig(level=logging.WARNING)

    min_val = 0
    max_val = 20000
    interval = 1000
    
    def make_x_y_grid(min_val, max_val, interval, one_dim=True):
        """
        Returns x and y values for the grid
        """
        coord_vals =  np.arange(min_val, max_val + interval, interval)
        y_val = 0
        
        if one_dim:
            coord_vals =  np.arange(min_val, max_val + interval, interval)
            
            x_vals = matlib.repmat(coord_vals, 1, len(coord_vals))
            x_vals = x_vals.reshape(-1, 1)[:, 0]
            
            y_vals = []
            for i in range(len(coord_vals)):
                for j in range(len(coord_vals)):
                     y_vals.append(y_val)
                
                y_val += interval
            
            y_vals = np.array(y_vals)
            
        else:
            x_vals = matlib.repmat(coord_vals, len(coord_vals), 1)
            
            y_vals = np.zeros((len(coord_vals), len(coord_vals)))
            for i in range(len(coord_vals)):
                for j in range(len(coord_vals)):
                     y_vals[i, j] = y_val
                
                y_val += interval
        
        return x_vals, y_vals
    
    x_vals, y_vals = make_x_y_grid(min_val, max_val, interval)
    
    numberOfShaleLayers = 3
    
    reservoirThickness = np.zeros(len(x_vals))
    shale1Thickness = np.zeros(len(x_vals))
    aquifer1Thickness = np.zeros(len(x_vals))
    shale2Thickness = np.zeros(len(x_vals))
    aquifer2Thickness =np.zeros(len(x_vals))
    shale3Thickness = np.zeros(len(x_vals))
    for i in range(len(x_vals)):
        reservoirThickness[i] = 25 + (0.00025 * x_vals[i]) + (0.0001 * y_vals[i])
        shale1Thickness[i] = 320 + (0.00033 * x_vals[i]) + (0.00015 * y_vals[i])
        aquifer1Thickness[i] = 30 + (0.0005 * x_vals[i]) + (0.0002 * y_vals[i])
        shale2Thickness[i] = 295 + (0.00025 * x_vals[i]) + (0.00033 * y_vals[i])
        aquifer2Thickness[i] = 45 + (0.0001 * x_vals[i]) + (0.0005 * y_vals[i])
        shale3Thickness[i] = 175 + (0.00005 * x_vals[i]) + (0.00075 * y_vals[i])
    
    # A .csv file will be written to this directory
    file_directory = os.sep.join([os.getcwd(), '..', '..', 'examples', 'scripts', 
                                  'lutstrata_gridded_data'])
    
    if not os.path.exists(file_directory):
        os.mkdir(file_directory)
    
    file_name = 'stratigraphy.csv'
    file_path = os.path.join(file_directory, file_name)
    
    data = {'x': x_vals, 
            'y': y_vals, 
            'reservoirThickness': reservoirThickness, 
            'shale1Thickness': shale1Thickness, 
            'aquifer1Thickness': aquifer1Thickness, 
            'shale2Thickness': shale2Thickness, 
            'aquifer2Thickness': aquifer2Thickness, 
            'shale3Thickness': shale3Thickness}
    
    strat_data = pd.DataFrame(data=data)
    
    strat_data.to_csv(file_path, index=False)
    
    # Now select x and y values within the larger domain defined by that file
    min_val = 2500
    max_val = 17500
    interval = 2500
    
    locX_grid, locY_grid = make_x_y_grid(min_val, max_val, interval, one_dim=False)

    # Define keyword arguments of the system model
    num_years = 10
    time_array = 365.25*np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array} # time is given in days
    
    numberOfShaleLayers = 3
    
    locX = 2000
    locY = 2000
    
    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)
    
    intpr = sm.add_interpolator(StratigraphyDataInterpolator(
        name='int1', parent=sm,
        header_file_dir=file_directory,
        data_file='stratigraphy.csv'), 
        intr_family='stratigraphy')

    # Add stratigraphy component
    luts = sm.add_component_model_object(LookupTableStratigraphy(
        name='luts', parent=sm, intr_family='stratigraphy', 
        locX=locX_grid, locY=locY_grid))
    
    luts.add_par('numberOfShaleLayers', value=numberOfShaleLayers, vary=False)
    
    thickness_obs = luts.get_thickness_obs_names()
    
    for ob_nm in thickness_obs:
        luts.add_obs(ob_nm, index=[0])
    
    depth_obs = luts.get_depth_obs_names()
    
    for ob_nm in depth_obs:
        luts.add_obs(ob_nm, index=[0])
    
    # Run system model using current values of its parameters
    sm.forward()  # system model is run deterministically
    
    print('Easting and Northing distances (m) for the stratigraphy information:')
    print('    Easting, x:', luts.locX)
    print('    Northing, y:', luts.locY)
    print('')
    
    # Collect thickness observations
    reservoirThickness = sm.collect_observations_as_time_series(
        luts, 'reservoirThickness', indices=[0])
    shale1Thickness = sm.collect_observations_as_time_series(
        luts, 'shale1Thickness', indices=[0])
    aquifer1Thickness = sm.collect_observations_as_time_series(
        luts, 'aquifer1Thickness', indices=[0])
    shale2Thickness = sm.collect_observations_as_time_series(
        luts, 'shale2Thickness', indices=[0])
    aquifer2Thickness = sm.collect_observations_as_time_series(
        luts, 'aquifer2Thickness', indices=[0])
    shale3Thickness = sm.collect_observations_as_time_series(
        luts, 'shale3Thickness', indices=[0])
    
    print('Unit thicknesses (m):')
    print('    reservoirThickness: ', reservoirThickness)
    print('    shale1Thickness: ', shale1Thickness)
    print('    aquifer1Thickness: ', aquifer1Thickness)
    print('    shale2Thickness: ', shale2Thickness)
    print('    aquifer2Thickness: ', aquifer2Thickness)
    print('    shale3Thickness: ', shale3Thickness)
    print('')
    
    # Collect depth observations
    reservoirTopDepth = sm.collect_observations_as_time_series(
        luts, 'reservoirTopDepth', indices=[0])
    shale1TopDepth = sm.collect_observations_as_time_series(
        luts, 'shale1TopDepth', indices=[0])
    aquifer1TopDepth = sm.collect_observations_as_time_series(
        luts, 'aquifer1TopDepth', indices=[0])
    shale2TopDepth = sm.collect_observations_as_time_series(
        luts, 'shale2TopDepth', indices=[0])
    aquifer2TopDepth = sm.collect_observations_as_time_series(
        luts, 'aquifer2TopDepth', indices=[0])
    shale3TopDepth = sm.collect_observations_as_time_series(
        luts, 'shale3TopDepth', indices=[0])
    
    # Bottom of shale 1, which is also the top of the reservoir
    shale1Depth = sm.collect_observations_as_time_series(
        luts, 'shale1Depth', indices=[0])
    
    print('Depths (m) to the top of each unit:')
    print('    reservoirTopDepth: ',  reservoirTopDepth)
    print('    shale1TopDepth: ', shale1TopDepth)
    print('    aquifer1TopDepth: ', aquifer1TopDepth)
    print('    shale2TopDepth: ', shale2TopDepth)
    print('    aquifer2TopDepth: ', aquifer2TopDepth)
    print('    shale3TopDepth: ', shale3TopDepth)
    
    # Make x and y grids for the 3D plot
    locX_grid_plots, locY_grid_plots = make_x_y_grid(min_val, max_val, interval, 
                                                     one_dim=False)
    
    # Reformat the shale1Depth array into a 2D grid
    dims = locY_grid_plots.shape
    shale1Depth_grid = np.zeros((dims[0], dims[1]))
    # These indices are used for the rows and columns of shale1Depth_grid
    i = -1
    j = 0
    for ind, val in enumerate(shale1Depth[0]):
        if (i + 1) >= dims[0]:
            i = 0
            j += 1
        else:
            i += 1
        
        shale1Depth_grid[i,j] = val
    
    shale1Depth_grid_min = np.ones((dims[0], dims[1])) * np.min(shale1Depth[0])
    
    def fmt(x, pos):
        a, b = '{:.2e}'.format(x).split('e')
        b = int(b)
        return r'${} \times 10^{{{}}}$'.format(a, b)
    
    # Make a 3-dimensional figure
    plt3d = plt.figure(1, figsize=(10, 8), dpi=200)
    plt3d_ax = plt3d.add_subplot(projection='3d')
    
    plt3d_ax.plot_surface(locX_grid_plots / 1000, locY_grid_plots / 1000, 
                          shale1Depth_grid, cmap='YlOrRd', alpha=0.67)
    plt3d_ax.plot_surface(locX_grid_plots / 1000, locY_grid_plots / 1000, 
                          shale1Depth_grid_min, color=[0.33, 0.33, 0.33], 
                          alpha=0.5, shade=False)
    
    # These vertical lines are used to highlight the perspective in the figure
    plt3d_ax.plot([locX_grid_plots[0, 0] / 1000, locX_grid_plots[0, 0] / 1000],
                  [locY_grid_plots[0, 0] / 1000, locY_grid_plots[0, 0] / 1000],
                  [np.min(shale1Depth[0]), shale1Depth_grid[0, 0]], color='k', linewidth=1)

    plt3d_ax.plot([locX_grid_plots[-1, -1] / 1000, locX_grid_plots[-1, -1] / 1000],
                  [locY_grid_plots[0, 0] / 1000, locY_grid_plots[0, 0] / 1000],
                  [np.min(shale1Depth[0]), shale1Depth_grid[0, -1]], color='k', linewidth=1)

    plt3d_ax.plot([locX_grid_plots[-1, -1] / 1000, locX_grid_plots[-1, -1] / 1000],
                  [locY_grid_plots[-1, -1] / 1000, locY_grid_plots[-1, -1] / 1000],
                  [np.min(shale1Depth[0]), shale1Depth_grid[-1, -1]], color='k', linewidth=1)

    plt3d_ax.plot([locX_grid_plots[0, 0] / 1000, locX_grid_plots[0, 0] / 1000],
                  [locY_grid_plots[-1, -1] / 1000, locY_grid_plots[-1, -1] / 1000],
                  [np.min(shale1Depth[0]), shale1Depth_grid[-1, 0]], color='k', linewidth=1)
    
    plt.title('Depth to the Top of the Reservoir', fontsize=12, fontweight='bold')
    
    axislabel_pad = 10
    plt.xlabel('Easting (km)', fontsize=12,
               fontweight='bold', labelpad=axislabel_pad)
    plt.ylabel('Northing (km)', fontsize=12,
               fontweight='bold', labelpad=axislabel_pad)
    plt3d_ax.set_zlabel('Depth (m)', fontsize=12,
                        fontweight='bold', labelpad=axislabel_pad)
    
