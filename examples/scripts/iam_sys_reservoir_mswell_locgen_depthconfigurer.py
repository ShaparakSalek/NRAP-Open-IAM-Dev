"""
This example demonstrates the use of the LocationGenerator and WellDepthRiskConfigurer 
components. The x, y, and z values (easting, northing, and depth) for wellbores 
are randomly created by the LocationGenerator. The well depths are then provided 
to the WellDepthRiskConfigurer - if the well does not reach the top of the 
reservoir, then the WellDepthRiskConfigurer disables the well by setting the 
run_frequency to 0 (so the wellbore component never runs).

This simulation uses a forward analysis.

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
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import ticker
from matplotlib.lines import Line2D

sys.path.insert(0, os.sep.join(['..', '..', 'source']))

from openiam import (SystemModel, Stratigraphy, SimpleReservoir, 
                     MultisegmentedWellbore, LocationGenerator, WellDepthRiskConfigurer)


# Simulation input
num_wells = 50

seedVal = 11

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
genfontsize = 8
axislabelfontsize = 10
titlefontsize = 10
lgndfontsize = 6

labelfontweight = 'bold'
colormap = 'RdYlBu_r'

lineWidth = 1.5

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


def get_distance_from_site_km(sm, compRef, injection_x_m, injection_y_m):
    """
    Returns the current component's distance from the injection site.
    """
    x_val_m = sm.obs['gen.locX{}_0'.format(compRef)].sim
    y_val_m = sm.obs['gen.locY{}_0'.format(compRef)].sim
    
    distance_from_site_km = (
        (((injection_x_m - x_val_m) ** 2) + ((injection_y_m - y_val_m) ** 2)) 
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
gen.add_par('seed', value=seedVal, vary=False)

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
    wellName = 'msw{}'.format(compRef)
    
    gen.add_local_obs('locX{}'.format(compRef), grid_obs_name='locX',
                      constr_type='array', loc_ind=compRef, index=[0])
    gen.add_local_obs('locY{}'.format(compRef), grid_obs_name='locY',
                      constr_type='array', loc_ind=compRef, index=[0])
    gen.add_local_obs('locZ{}'.format(compRef), grid_obs_name='locZ',
                      constr_type='array', loc_ind=compRef, index=[0])
    
    configs.append(sm.add_component_model_object(WellDepthRiskConfigurer(
        name='config{}'.format(compRef), parent=sm, cmpnt_nms=[wellName])))
    
    configs[-1].add_par('reservoirDepth', value=reservoirDepth, vary=False)
    
    configs[-1].add_kwarg_linked_to_obs(
        'wellDepth', gen.linkobs['locZ'], obs_type='grid', 
        constr_type='array', loc_ind=[compRef])
    
    # Add simple reservoir component
    # We don't specify the location for the reservoir component, although
    # there is a default reservoir setup with locX=100, locY=100
    sres.append(sm.add_component_model_object(SimpleReservoir(
        name='sres{}'.format(compRef), parent=sm, 
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
    msw[-1].add_par('logWellPerm', value=-13.0, vary=False)
    msw[-1].add_par('logAquPerm', value=-11.0, vary=False)

    # Add keyword arguments linked to the output provided by reservoir model
    msw[-1].add_kwarg_linked_to_obs('pressure', sres[-1].linkobs['pressure'])
    msw[-1].add_kwarg_linked_to_obs('CO2saturation', sres[-1].linkobs['CO2saturation'])

    # Add observations of multisegmented wellbore component model
    msw[-1].add_obs('CO2_aquifer1')
    msw[-1].add_obs('CO2_aquifer2')
    msw[-1].add_obs('brine_aquifer1')
    msw[-1].add_obs('brine_aquifer2')

# Run forward simulation
sm.forward()


no_results_check = True
print('__________________Results of forward simulation___________________')
print('The simulation used {} wells with randomly generated x, y, and z.'.format(num_wells))
print('Wells penetrating the reservoir: ')
for compRef in range(num_wells):
    # Making them both positive just for conceptual clarity
    if -sm.obs['gen.locZ{}_0'.format(compRef)].sim >= reservoirDepth:
        if no_results_check:
            no_results_check = False
        
        print('compRef: ', compRef)
        print('x = {:.2f} m, y = {:.2f} m'.format(sm.obs['gen.locX{}_0'.format(compRef)].sim, 
                                          sm.obs['gen.locY{}_0'.format(compRef)].sim))
        print('z = {:.2f} m'.format(sm.obs['gen.locZ{}_0'.format(compRef)].sim))
        print('    Pressure [Pa]: ', 
              sm.collect_observations_as_time_series(sres[compRef], 'pressure'), sep='\n')
        print('    Brine Leakage Rate [kg/s]: ', 
              sm.collect_observations_as_time_series(msw[compRef], 'brine_aquifer1'), sep='\n')
        print(' ')

if no_results_check:
    print('None')


distance_range_km = [0, 0]
# First, get the range of x and y values. The range is used for the colorbar below.
for compRef in range(num_wells):
    # Making them both positive just for conceptual clarity
    if -sm.obs['gen.locZ{}_0'.format(compRef)].sim >= reservoirDepth:
        distance_from_site_km = get_distance_from_site_km(
            sm, compRef, injection_x_m, injection_y_m)
        
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
    
else:
    if min_level == max_level:
        max_level *= 1.1

interval = (max_level - min_level) / 100

distance_levels_km = np.arange(min_level, max_level + interval, interval)


cmap = plt.cm.get_cmap(colormap)

font = {'family': 'Arial',
        'weight': 'normal',
        'size': genfontsize}
plt.rc('font', **font)
dpi_ref = 300
plt.rcParams['figure.dpi'] = dpi_ref
plt.rcParams['image.cmap'] = 'jet'


for compRef in range(num_wells):
    # Making them both positive just for conceptual clarity
    if -sm.obs['gen.locZ{}_0'.format(compRef)].sim < reservoirDepth:
        # Well does not penetrate the reservoir
        plt.figure(1)
        plt.plot(sm.obs['gen.locX{}_0'.format(compRef)].sim / 1000, 
                 sm.obs['gen.locY{}_0'.format(compRef)].sim / 1000, 
                 linestyle='none', markeredgewidth=edgeWidth, marker=wellMarker, 
                 markerfacecolor=unusedWellColor, markersize=unusedWellMarkerSize, 
                 color=unusedWellEdgeColor, zorder=1)
        
    else:
        # Well penetrates the reservoir
        distance_from_site_km = get_distance_from_site_km(
            sm, compRef, injection_x_m, injection_y_m)
        
        rgba = cmap(
            (distance_from_site_km  - np.min(distance_range_km))
            / (np.max(distance_range_km) - np.min(distance_range_km)))
        
        markerfacecolor = rgba[0:3]
        markeredgecolor = (np.array(rgba[0:3]) * 0.5).tolist()
        
        plt.figure(1)
        plt.plot(sm.obs['gen.locX{}_0'.format(compRef)].sim / 1000, 
                 sm.obs['gen.locY{}_0'.format(compRef)].sim / 1000, 
                 linestyle='none', markeredgewidth=edgeWidth, marker=wellMarker, 
                 markerfacecolor=markerfacecolor, markersize=usedWellMarkerSize, 
                 color=markeredgecolor, zorder=2)
        
        pressure = sm.collect_observations_as_time_series(
            sres[compRef], 'pressure')
        
        plt.figure(2)
        plt.plot(time_array / 365.25, pressure / 1.0e+6, color=markerfacecolor, 
                 linewidth=lineWidth)
        
        brine_aquifer1 = sm.collect_observations_as_time_series(
            msw[compRef], 'brine_aquifer1')
        
        plt.figure(3)
        plt.plot(time_array / 365.25, brine_aquifer1, color=markerfacecolor, 
                 linewidth=lineWidth)

fig = plt.figure(1)

plt.plot(injection_x_m / 1000, injection_y_m / 1000, color='k', 
         linestyle='none', marker=injectionSiteMarker, markerfacecolor='none', 
         markersize=injectionSiteMarkerSize, markeredgewidth=edgeWidthSite, zorder=3)


# Manually create the legend items for figure 1
fig1_handle_list = []
fig1_label_list = []

injectionSite = Line2D(
    [0], [0], color='k', markeredgewidth=edgeWidthSite, linestyle='none', 
    marker=injectionSiteMarker, markerfacecolor='none', 
    markersize=injectionSiteMarkerSize)

fig1_handle_list.append(injectionSite)
fig1_label_list.append('Injection Site')

unusedWells = Line2D(
    [0], [0], linestyle='none', marker=wellMarker, markerfacecolor=unusedWellColor, 
    markersize=unusedWellMarkerSize, color=unusedWellEdgeColor, 
    markeredgewidth=edgeWidth)

fig1_handle_list.append(unusedWells)
fig1_label_list.append('Well, No Connection to Reservoir')

rgba = cmap(1)

usedWells = Line2D(
    [0], [0], linestyle='none', marker=wellMarker, markerfacecolor=rgba[0:3], 
    markersize=usedWellMarkerSize, color=(np.array(rgba[0:3]) * 0.5).tolist(), 
    markeredgewidth=edgeWidth)

fig1_handle_list.append(usedWells)
fig1_label_list.append('Well Penetrating the Reservoir')


# Format plots
fig = plt.figure(1)

plt.xlabel('x (km)', fontsize=axislabelfontsize, fontweight=labelfontweight)
plt.ylabel('y (km)', fontsize=axislabelfontsize, fontweight=labelfontweight)
plt.title('Map-View Image of the Site', fontsize=titlefontsize, fontweight=labelfontweight)

# Make the colorbar
cbar = fig.colorbar(cm.ScalarMappable(cmap=cmap), ax=plt.gca(),
                    values=distance_levels_km)
cbar.set_label('Distance From Injection Site (km)', rotation=90, fontsize=axislabelfontsize,
               fontweight=labelfontweight)
tick_locator = ticker.MaxNLocator(nbins=5)
cbar.locator = tick_locator
cbar.update_ticks()

plt.xlim(x_min / 1000, x_max / 1000)
plt.xlim(y_min / 1000, y_max / 1000)

plt.subplots_adjust(left=0.2, bottom=0.2, right=0.8,
                    top=0.8, wspace=0.1, hspace=0.1)

# Add Legend
plt.gca().legend(fig1_handle_list, fig1_label_list, fancybox=False, 
                 fontsize=lgndfontsize, ncol=1, 
                 edgecolor=[0, 0, 0], loc='upper center', 
                 bbox_to_anchor=(0.5, -0.2), framealpha=0.67)

fig = plt.figure(2)

plt.xlabel('Time (years)', fontsize=axislabelfontsize, fontweight=labelfontweight)
plt.ylabel('Pressure (MPa)', fontsize=axislabelfontsize, fontweight=labelfontweight)
plt.title('Reservoir Pressure Over Time', fontsize=titlefontsize, fontweight=labelfontweight)

# Make the colorbar
cbar = fig.colorbar(cm.ScalarMappable(cmap=cmap), ax=plt.gca(),
                    values=distance_levels_km)
cbar.set_label('Distance From Injection Site (km)', rotation=90, fontsize=axislabelfontsize,
               fontweight=labelfontweight)
tick_locator = ticker.MaxNLocator(nbins=5)
cbar.locator = tick_locator
cbar.update_ticks()

fig = plt.figure(3)

plt.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0,0), useMathText = True)

plt.xlabel('Time (years)', fontsize=axislabelfontsize, fontweight=labelfontweight)
plt.ylabel('Brine Leakage Rate (kg/s)', fontsize=axislabelfontsize, fontweight=labelfontweight)
plt.title('Brine Leakage Over Time', fontsize=titlefontsize, fontweight=labelfontweight)

# Make the colorbar
cbar = fig.colorbar(cm.ScalarMappable(cmap=cmap), ax=plt.gca(),
                    values=distance_levels_km)
cbar.set_label('Distance From Injection Site (km)', rotation=90, fontsize=axislabelfontsize,
               fontweight=labelfontweight)
tick_locator = ticker.MaxNLocator(nbins=5)
cbar.locator = tick_locator
cbar.update_ticks()
