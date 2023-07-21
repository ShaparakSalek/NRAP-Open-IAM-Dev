"""
This example demonstrates the use of the LocationGenerator and WellDepthRiskConfigurer 
components. The x, y, and z values (easting, northing, and depth) for wellbores 
are randomly created by the LocationGenerator. The well depths are then provided 
to the WellDepthRiskConfigurer - if the well does not reach the top of the 
reservoir, then the WellDepthRiskConfigurer disables the well by setting the 
run_frequency to 0 (so the wellbore component never runs).

This simulation uses the Latin Hypercube Sampling (lhs) analysis type. The 

Note that randomly assigning coordinates with the LocationGenerator is different 
from first creating random x, y, and z values and then directly providing those 
to a reservoir component. The difference lies in the fact that the LocationGenerator 
creates the x, y, and z values during the actual simulation - if you run a 
simulation with Latin Hypercube Sampling, each realization will have different 
x, y, and z values because they are generated during the realization. Having 
random well locations generated in each realization can be important if you know 
that well records are incomplete - this approach allows one to assess the potential 
significance of unknown wellbores. Some of the unknown wells may not penetrate 
the storage reservoir and may therefore be excluded as potential leakage pathways.
"""
# @author: Nate Mitchell
# Nathaniel.Mitchell@netl.doe.gov

import sys
import os
import numpy as np
import datetime
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import ticker
from matplotlib.lines import Line2D

sys.path.insert(0, os.sep.join(['..', '..', 'source']))

from openiam import (SystemModel, Stratigraphy, SimpleReservoir, 
                     MultisegmentedWellbore, LocationGenerator, WellDepthRiskConfigurer)
from matk.sampleset import percentile, mean

if __name__ == "__main__":
    # Simulation input
    num_wells = 25

    # These are the realizations that will have map-view images of the site. 
    # The randomly generated wells will be different in each realization.
    selected_realizations = [0, 1, 2]

    # This is used as the seed input for the LocationGenerator
    LocSeed = 100
    minLocSeed = 1
    maxLocSeed = 1000

    # This is the seed used when running the LHS simulation
    lhsSeedVal = 589

    num_years = 50

    shale1Thickness = 250
    shale2Thickness = 250
    shale3Thickness = 250
    aquifer1Thickness = 100
    aquifer2Thickness = 100
    reservoirThickness = 50

    # Define boundaries of random wells domains
    x_min = -5000.0
    x_max = 5000.0

    y_min = -5000.0
    y_max = 5000.0

    injection_x_m = 0
    injection_y_m = 0

    # Boundaries defining well depths
    # Here, the z values can range from the bottom of aquifer 1 (top of shale 2) to 
    # the bottom of the reservoir. To represent depths beneath the surface, the z 
    # values need to be negative.

    # Bottom of aquifer 1
    z_max = -(shale3Thickness + aquifer2Thickness + shale2Thickness 
              + aquifer1Thickness)

    # Bottom of the reservoir - the z values need to be negative, so this is the 
    # minimum value. The generator will swap the min and max if they are mixed up though.
    z_min = -(shale3Thickness + aquifer2Thickness + shale2Thickness 
              + aquifer1Thickness + shale1Thickness  + reservoirThickness)

    # The reservoirDepth parameter is taken as positive, however. This is the depth 
    # to the top of the reservoir
    reservoirDepth = (shale3Thickness + aquifer2Thickness + shale2Thickness 
                      + aquifer1Thickness + shale1Thickness)

    # Figure formatting input
    figsize_unused_wells = (14, 6)
    figsize_used_wells = (16, 6)
    
    dpi_ref = 300
    
    genfontsize = 8
    axislabelfontsize = 10
    titlefontsize = 10
    lgndfontsize = 8

    labelfontweight = 'bold'
    colormap = 'RdYlBu_r'

    edgeWidth = 1
    edgeWidthSite = 2

    unusedWellColor = [0.33, 0.33, 0.33]
    unusedWellEdgeColor = [unusedWellColor[0] * 0.5, unusedWellColor[1] * 0.5, 
                           unusedWellColor[2] * 0.5]

    unusedWellMarkerSize = 4
    usedWellMarkerSize = 7
    injectionSiteMarkerSize = 5

    injectionSiteMarker = 's'
    wellMarker = 'o'
    
    reals_color = 'darkslateblue'
    reals_alpha = 0.8
    reals_linewidth = 2.5
    
    output_folder = 'output_locgen_depthconfig_lhs_{date_time_stamp}'.format(
        date_time_stamp=datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
    
    IAM_output_dir = os.path.join(os.getcwd(), '..', '..', 'output')
    output_dir = os.path.join(IAM_output_dir, output_folder)
    
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)


    def get_distance_from_site_km(sm, compRef, injection_x_m, injection_y_m, 
                                  locX, locY):
        """
        Returns the current component's distance from the injection site.
        """
        distance_from_site_km = (
            (((injection_x_m - locX) ** 2) + ((injection_y_m - locY) ** 2)) 
            ** 0.5) / 1000
        
        return distance_from_site_km

    # Define keyword arguments of the system model
    time_array = 365.25 * np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array} # time is given in days

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Add stratigraphy component
    strata = sm.add_component_model_object(Stratigraphy(name='strata', parent=sm))

    # Add parameters of stratigraphy component model
    strata.add_par('numberOfShaleLayers', value=3, vary=False)
    strata.add_par('shale1Thickness', value=shale1Thickness, vary=False)
    strata.add_par('shale2Thickness', value=shale2Thickness, vary=False)
    strata.add_par('shale3Thickness', value=shale3Thickness, vary=False)
    strata.add_par('aquifer1Thickness', value=aquifer1Thickness, vary=False)
    strata.add_par('aquifer2Thickness', value=aquifer2Thickness, vary=False)
    strata.add_par('reservoirThickness', value=reservoirThickness, vary=False)

    # Add generator component
    gen = sm.add_component_model_object(LocationGenerator(
        name='gen', parent=sm, x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max, 
        num_locations=num_wells, reproducible=True, z_min=z_min, z_max=z_max))

    # Add parameters of generator component
    gen.add_par('seed', value=LocSeed, min=minLocSeed, max=maxLocSeed)

    gen.add_obs_to_be_linked('locX', obs_type='grid')
    gen.add_obs_to_be_linked('locY', obs_type='grid')
    gen.add_obs_to_be_linked('locZ', obs_type='grid')

    # Locations can be obtained as observations of the generator component
    # Since the observations considered to be of gridded type the outputs
    # can only be obtained by adding them using add_obs method and
    # then reading the files with data after simulation is complete

    sres = []
    msw = []
    configs = []

    for compRef in range(num_wells):
        # This is needed for the WellDepthRiskConfigurer component
        resName = 'sres{}'.format(compRef)
        wellName = 'msw{}'.format(compRef)
        
        gen.add_local_obs('locX{}'.format(compRef), grid_obs_name='locX',
                          constr_type='array', loc_ind=compRef, index=[0])
        gen.add_local_obs('locY{}'.format(compRef), grid_obs_name='locY',
                          constr_type='array', loc_ind=compRef, index=[0])
        gen.add_local_obs('locZ{}'.format(compRef), grid_obs_name='locZ',
                          constr_type='array', loc_ind=compRef, index=[0])
        
        configs.append(sm.add_component_model_object(WellDepthRiskConfigurer(
            name='config{}'.format(compRef), parent=sm, cmpnt_nms=[resName, wellName])))
        
        configs[-1].add_par('reservoirDepth', value=reservoirDepth, vary=False)
        
        configs[-1].add_kwarg_linked_to_obs(
            'wellDepth', gen.linkobs['locZ'], obs_type='grid', 
            constr_type='array', loc_ind=[compRef])
        
        # Add simple reservoir component
        # We don't specify the location for the reservoir component, although
        # there is a default reservoir setup with locX=100, locY=100
        sres.append(sm.add_component_model_object(SimpleReservoir(
            name=resName, parent=sm, 
            injX=injection_x_m, injY=injection_y_m)))

        # Add parameters of reservoir component model
        # All parameters of the reservoir component are deterministic so all 
        # uncertainty in the simulation comes from the uncertainty of the well location
        sres[-1].add_par('numberOfShaleLayers', value=3, vary=False)
        sres[-1].add_par('injRate', value=0.5, vary=False)
        sres[-1].add_par('shale1Thickness', value=shale1Thickness, vary=False)
        sres[-1].add_par('shale2Thickness', value=shale2Thickness, vary=False)
        sres[-1].add_par('shale3Thickness', value=shale3Thickness, vary=False)
        sres[-1].add_par('aquifer1Thickness', value=aquifer1Thickness, vary=False)
        sres[-1].add_par('aquifer2Thickness', value=aquifer2Thickness, vary=False)
        sres[-1].add_par('reservoirThickness', value=reservoirThickness, vary=False)
        sres[-1].add_par('logResPerm', value=-12, vary=False) # min=-12.5, max=-11.5)

        # Simple reservoir component has keyword arguments of the model method:
        # locX and locY which we would link to the output of the generator component
        # Generator component outputs 5 random locations according to the setup above.
        # We can link reservoir component locX and locY to any of these produced locations
        # using arguments constr_type='array' and loc_ind=[compRef].
        sres[-1].add_kwarg_linked_to_obs('locX', gen.linkobs['locX'], obs_type='grid', 
                                          constr_type='array', loc_ind=[compRef])
        sres[-1].add_kwarg_linked_to_obs('locY', gen.linkobs['locY'], obs_type='grid', 
                                          constr_type='array', loc_ind=[compRef])

        # Add observations of reservoir component model
        sres[-1].add_obs('pressure')
        sres[-1].add_obs('CO2saturation')
        
        sres[-1].add_obs_to_be_linked('pressure')
        sres[-1].add_obs_to_be_linked('CO2saturation')
        
        # Add multisegmented wellbore component
        msw.append(sm.add_component_model_object(MultisegmentedWellbore(
            name='msw{}'.format(compRef), parent=sm)))

        # Add parameters of multisegmented wellbore component
        msw[-1].add_par('wellRadius', value=0.015, vary=False)
        msw[-1].add_par('numberOfShaleLayers', value=3, vary=False)
        msw[-1].add_par('shale1Thickness', value=shale1Thickness, vary=False)
        msw[-1].add_par('shale2Thickness', value=shale2Thickness, vary=False)
        msw[-1].add_par('shale3Thickness', value=shale3Thickness, vary=False)
        msw[-1].add_par('aquifer1Thickness', value=aquifer1Thickness, vary=False)
        msw[-1].add_par('aquifer2Thickness', value=aquifer2Thickness, vary=False)
        msw[-1].add_par('reservoirThickness', value=reservoirThickness, vary=False)
        msw[-1].add_par('logWellPerm', value=-13.0, vary=False) # min=-13.5, max=-12.5)
        msw[-1].add_par('logAquPerm', value=-11.0, vary=False)

        # Add keyword arguments linked to the output provided by reservoir model
        msw[-1].add_kwarg_linked_to_obs('pressure', sres[-1].linkobs['pressure'])
        msw[-1].add_kwarg_linked_to_obs('CO2saturation', sres[-1].linkobs['CO2saturation'])

        # Add observations of multisegmented wellbore component model
        msw[-1].add_obs('CO2_aquifer1')
        msw[-1].add_obs('CO2_aquifer2')
        msw[-1].add_obs('brine_aquifer1')
        msw[-1].add_obs('brine_aquifer2')


    num_samples = 30
    ncpus = 2
    # Draw Latin hypercube samples of parameter values
    s = sm.lhs(siz=num_samples, seed=lhsSeedVal)

    # Run model using values in samples for parameter values
    results = s.run(cpus=ncpus, verbose=False)


    locX = np.ones((num_samples, 1))
    locY = np.ones((num_samples, 1))
    locZ = np.ones((num_samples, 1))

    distance_range_km = [0, 0]
    
    # First, get the range of x and y values for the selected realization. The 
    # range is used for the colorbar below.
    for compRef in range(num_wells):
        locX[:, 0] = s.recarray['gen.locX{}_0'.format(compRef)]
        locY[:, 0] = s.recarray['gen.locY{}_0'.format(compRef)]
        locZ[:, 0] = s.recarray['gen.locZ{}_0'.format(compRef)]
        
        # Go through all realizations
        for real in range(num_samples):
            # Extract the x, y, and z values for the realization chosen for Figure 1
            locX_current = locX[real, 0]
            locY_current = locY[real, 0]
            locZ_current = locZ[real, 0]
            
            # Making them both positive just for conceptual clarity
            if -locZ_current >= reservoirDepth:
                distance_from_site_km = get_distance_from_site_km(
                    sm, compRef, injection_x_m, injection_y_m, locX_current, locY_current)
                
                # If it's the first time, set the value
                if distance_range_km[0] == 0:
                    distance_range_km[0] = distance_from_site_km
                
                if distance_range_km[1] == 0:
                    distance_range_km[1] = distance_from_site_km
                
                if distance_from_site_km < distance_range_km[0]:
                    distance_range_km[0] = distance_from_site_km
                
                if distance_from_site_km > distance_range_km[1]:
                    distance_range_km[1] = distance_from_site_km

    min_level = distance_range_km[0]
    max_level = distance_range_km[1]

    if min_level == 0 and max_level == 0:
        max_level = 1
    elif min_level == max_level:
        max_level *= 1.1

    interval = (max_level - min_level) / 100

    distance_levels_km = np.arange(min_level, max_level + interval, interval)


    cmap = plt.cm.get_cmap(colormap)

    font = {'family': 'Arial',
            'weight': 'normal',
            'size': genfontsize}
    plt.rc('font', **font)
    
    plt.rcParams['figure.dpi'] = dpi_ref
    plt.rcParams['image.cmap'] = 'jet'


    locX = np.ones((num_samples, 1))
    locY = np.ones((num_samples, 1))
    locZ = np.ones((num_samples, 1))

    # First, plot the map-view images for each of the selected realizations
    for real in selected_realizations:
        plt.figure(100 + real, dpi=dpi_ref)
        ax = plt.gca()
        
        plt.plot(injection_x_m / 1000, injection_y_m / 1000, color='k', 
                  linestyle='none', marker=injectionSiteMarker, markerfacecolor='none', 
                  markersize=injectionSiteMarkerSize, markeredgewidth=edgeWidthSite, zorder=3)
        
        for compRef in range(num_wells):
            locX[:, 0] = s.recarray['gen.locX{}_0'.format(compRef)]
            locY[:, 0] = s.recarray['gen.locY{}_0'.format(compRef)]
            locZ[:, 0] = s.recarray['gen.locZ{}_0'.format(compRef)]
            
            locX_current = locX[real, 0]
            locY_current = locY[real, 0]
            locZ_current = locZ[real, 0]
            
            if -locZ_current < reservoirDepth:
                # Well does not penetrate the reservoir
                ax.plot(locX_current / 1000, locY_current / 1000, 
                        linestyle='none', markeredgewidth=edgeWidth, marker=wellMarker, 
                        markerfacecolor=unusedWellColor, markersize=unusedWellMarkerSize, 
                        color=unusedWellEdgeColor, zorder=1)
            else:
                # Well penetrates the reservoir
                distance_from_site_km = get_distance_from_site_km(
                    sm, compRef, injection_x_m, injection_y_m, locX_current, locY_current)
                
                rgba = cmap(
                    (distance_from_site_km  - np.min(distance_range_km))
                    / (np.max(distance_range_km) - np.min(distance_range_km)))
                
                markerfacecolor = rgba[0:3]
                markeredgecolor = (np.array(rgba[0:3]) * 0.5).tolist()
                
                ax.plot(locX_current / 1000, locY_current / 1000, 
                        linestyle='none', markeredgewidth=edgeWidth, marker=wellMarker, 
                        markerfacecolor=markerfacecolor, markersize=usedWellMarkerSize, 
                        color=markeredgecolor, zorder=2)
    
    
    pressure = np.ones((num_samples, len(time_array)))
    brine_aquifer1 = np.ones((num_samples, len(time_array)))
    
    # Now plot the output for all realizations and all components
    for compRef in range(num_wells):
        locX[:, 0] = s.recarray['gen.locX{}_0'.format(compRef)]
        locY[:, 0] = s.recarray['gen.locY{}_0'.format(compRef)]
        locZ[:, 0] = s.recarray['gen.locZ{}_0'.format(compRef)]
        
        for ind in range(len(time_array)):
            pressure[:, ind] = s.recarray['sres{}.pressure_{}'.format(compRef, ind)]
            brine_aquifer1[:, ind] = s.recarray['msw{}.brine_aquifer1_{}'.format(compRef, ind)]
        
        for real in range(num_samples):
            locX_current = locX[real, 0]
            locY_current = locY[real, 0]
            locZ_current = locZ[real, 0]
            
            if -locZ_current < reservoirDepth:
                # Well does not penetrate the reservoir
                plt.figure(2, figsize=figsize_unused_wells, dpi=dpi_ref)
                
                color = unusedWellColor
            else:
                # Well penetrates the reservoir
                plt.figure(3, figsize=figsize_used_wells, dpi=dpi_ref)
                
                distance_from_site_km = get_distance_from_site_km(
                    sm, compRef, injection_x_m, injection_y_m, locX_current, locY_current)
                
                rgba = cmap(
                    (distance_from_site_km  - np.min(distance_range_km))
                    / (np.max(distance_range_km) - np.min(distance_range_km)))
                
                color = rgba[0:3]
                color = (np.array(rgba[0:3]) * 0.5).tolist()
            
            ax = plt.subplot(1, 2, 1)
            ax.plot(time_array / 365.25, pressure[real, :] / 1.0e+6, color=color, 
                    linestyle='-', linewidth=reals_linewidth, alpha=reals_alpha, 
                    zorder=100)
            
            ax = plt.subplot(1, 2, 2)
            ax.plot(time_array / 365.25, brine_aquifer1[real, :], color=color, 
                    linestyle='-', linewidth=reals_linewidth, alpha=reals_alpha, 
                    zorder=100)


    # Manually create the legend items for the map-view figures
    mapview_fig_handle_list = []
    mapview_fig_label_list = []

    injectionSite = Line2D(
        [0], [0], color='k', markeredgewidth=edgeWidthSite, linestyle='none', 
        marker=injectionSiteMarker, markerfacecolor='none', 
        markersize=injectionSiteMarkerSize)

    mapview_fig_handle_list.append(injectionSite)
    mapview_fig_label_list.append('Injection Site')

    unusedWells = Line2D(
        [0], [0], linestyle='none', marker=wellMarker, markerfacecolor=unusedWellColor, 
        markersize=unusedWellMarkerSize, color=unusedWellEdgeColor, 
        markeredgewidth=edgeWidth)

    mapview_fig_handle_list.append(unusedWells)
    mapview_fig_label_list.append('Well, No Connection to Reservoir')

    rgba = cmap(1)

    usedWells = Line2D(
        [0], [0], linestyle='none', marker=wellMarker, markerfacecolor=rgba[0:3], 
        markersize=usedWellMarkerSize, color=(np.array(rgba[0:3]) * 0.5).tolist(), 
        markeredgewidth=edgeWidth)

    mapview_fig_handle_list.append(usedWells)
    mapview_fig_label_list.append('Well Penetrating the Reservoir')


    # Format and save plots
    for real in selected_realizations:
        fig = plt.figure(100 + real)

        plt.xlabel('x (km)', fontsize=axislabelfontsize, fontweight=labelfontweight)
        plt.ylabel('y (km)', fontsize=axislabelfontsize, fontweight=labelfontweight)
        plt.title('Map-View Image of the Site,\nRealization {}'.format(real), 
                  fontsize=titlefontsize, fontweight=labelfontweight)

        # Make the colorbar
        cbar = fig.colorbar(cm.ScalarMappable(cmap=cmap), ax=plt.gca(),
                            values=distance_levels_km)
        cbar.set_label('Distance From Injection Site (km)', rotation=90, fontsize=axislabelfontsize,
                        fontweight=labelfontweight)
        tick_locator = ticker.MaxNLocator(nbins=5)
        cbar.locator = tick_locator
        cbar.update_ticks()

        plt.xlim(x_min / 1000, x_max / 1000)
        plt.ylim(y_min / 1000, y_max / 1000)

        plt.subplots_adjust(left=0.15, bottom=0.15, right=0.85,
                            top=0.85, wspace=0.1, hspace=0.1)
        
        # Add Legend
        plt.gca().legend(mapview_fig_handle_list, mapview_fig_label_list, fancybox=False, 
                         fontsize=lgndfontsize, ncol=1, 
                         edgecolor=[0, 0, 0], loc='upper center', 
                         bbox_to_anchor=(0.5, -0.15), framealpha=0.67)
        
        plt.savefig(os.path.join(output_dir, 'Wells_Map_View_Realizaiton{}.png'.format(
            real)), dpi=dpi_ref)
    
    
    # All Realizations, Only Wells Without a Connection to the Reservoir
    fig = plt.figure(2)
    
    ax = plt.subplot(1, 2, 1)

    plt.xlabel('Time (years)', fontsize=axislabelfontsize, fontweight=labelfontweight)
    plt.ylabel('Pressure (MPa)', fontsize=axislabelfontsize, fontweight=labelfontweight)
    
    ax = plt.subplot(1, 2, 2)

    plt.xlabel('Time (years)', fontsize=axislabelfontsize, fontweight=labelfontweight)
    plt.ylabel('Brine Leakage Rate (kg/s)', fontsize=axislabelfontsize, 
               fontweight=labelfontweight)
    
    plt.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0,0), 
                         useMathText = True)
    
    plt.suptitle('Results Over Time, Only Wells Without a Connection to the Reservoir', 
                 fontsize=titlefontsize, fontweight=labelfontweight)
    
    plt.subplots_adjust(left=0.25, bottom=0.2, right=0.75,
                        top=0.9, wspace=0.3, hspace=0.2)
    
    plt.savefig(os.path.join(output_dir, 'Results_Over_Time_Wells_No_Connection.png'), 
                dpi=dpi_ref)
    
    
    # All Realizations, Only Wells Penetrating the Reservoir
    fig = plt.figure(3)
    
    ax = plt.subplot(1, 2, 1)

    plt.xlabel('Time (years)', fontsize=axislabelfontsize, fontweight=labelfontweight)
    plt.ylabel('Pressure (MPa)', fontsize=axislabelfontsize, fontweight=labelfontweight)
    
    # Make the colorbar
    cbar = fig.colorbar(cm.ScalarMappable(cmap=cmap), ax=plt.gca(),
                        values=distance_levels_km)
    cbar.set_label('Distance From Injection Site (km)', rotation=90, fontsize=axislabelfontsize,
                    fontweight=labelfontweight)
    tick_locator = ticker.MaxNLocator(nbins=5)
    cbar.locator = tick_locator
    cbar.update_ticks()
    
    ax = plt.subplot(1, 2, 2)

    plt.xlabel('Time (years)', fontsize=axislabelfontsize, fontweight=labelfontweight)
    plt.ylabel('Brine Leakage Rate (kg/s)', fontsize=axislabelfontsize, 
               fontweight=labelfontweight)
    
    plt.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0,0), 
                         useMathText = True)
    
    # Make the colorbar
    cbar = fig.colorbar(cm.ScalarMappable(cmap=cmap), ax=plt.gca(),
                        values=distance_levels_km)
    cbar.set_label('Distance From Injection Site (km)', rotation=90, fontsize=axislabelfontsize,
                    fontweight=labelfontweight)
    tick_locator = ticker.MaxNLocator(nbins=5)
    cbar.locator = tick_locator
    cbar.update_ticks()
    
    plt.suptitle('Results Over Time, Only Wells Penetrating the Reservoir', fontsize=titlefontsize, 
                 fontweight=labelfontweight)
    
    plt.subplots_adjust(left=0.25, bottom=0.2, right=0.75,
                        top=0.9, wspace=0.3, hspace=0.2)
    
    plt.savefig(os.path.join(output_dir, 'Results_Over_Time_Wells_Penetrating_Reservoir.png'), 
                dpi=dpi_ref)

