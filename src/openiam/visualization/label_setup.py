# -*- coding: utf-8 -*-
"""
Module containing dictionaries defining default titles, y-axis and legend labels
for time series plots created.
"""
# These additional keys are used for the '_SUBPLOT' versions of the dictionaries.
# The 'SUBPLOT' versions are used only when subplots are used in time_series.py.
# These dictionaries specify the aquifer name in the axis labels - this approach 
# can be problematic when plotting all results in one plot. In this case, for 
# example, the axis label might show 'aquifer 2' while the results included 
# actually reflect aquifers 1 and 2.
additional_keys = [
    'brine_aquifer{}'.format(unitRef) for unitRef in range(1, 30)
    ] + ['CO2_aquifer{}'.format(unitRef) for unitRef in range(1, 30)
         ] + ['mass_brine_aquifer{}'.format(unitRef) for unitRef in range(1, 30)
              ] + ['mass_CO2_aquifer{}'.format(unitRef) for unitRef in range(1, 30)]

additional_ylabels = [
    'Leakage rate of brine to aquifer {} (kg/s)'.format(unitRef) for unitRef in range(1, 30)
    ] + ['Leakage rate of CO$_2$ to aquifer {} (kg/s)'.format(unitRef) for unitRef in range(1, 30)
         ] + ['Mass of brine leaked to aquifer {} (kg)'.format(unitRef) for unitRef in range(1, 30)
              ] + ['Mass of CO$_2$ leaked to aquifer {} (kg)'.format(unitRef) for unitRef in range(1, 30)]

additional_ylabels_2rows = [
    'Leakage rate of\nbrine to aquifer {} (kg/s)'.format(unitRef) for unitRef in range(1, 30)
    ] + ['Leakage rate of\nO$_2$ to aquifer {} (kg/s)'.format(unitRef) for unitRef in range(1, 30)
         ] + ['Mass of brine\nleaked to aquifer {} (kg)'.format(unitRef) for unitRef in range(1, 30)
              ] + ['Mass of CO$_2$\nleaked to aquifer {} (kg)'.format(unitRef) for unitRef in range(1, 30)]

additional_titles = [
    'Brine Leakage Rate to Aquifer {}'.format(unitRef) for unitRef in range(1, 30)
    ] + ['CO$_2$ Leakage Rate to Aquifer {}'.format(unitRef) for unitRef in range(1, 30)
         ] + ['Brine Mass Leaked to Aquifer {}'.format(unitRef) for unitRef in range(1, 30)
              ] + ['CO$_2$ Mass Leaked to Aquifer {}'.format(unitRef) for unitRef in range(1, 30)]

# Dictionary with the y-label names for different component output. Because single
# plots can display results for multiple aquifers, I have removed the aquifer numbers here.
Y_LABEL_DICT = {
    # Reservoir components
    'pressure': 'Reservoir pressure at location (Pa)',
    'CO2saturation': 'Reservoir CO$_2$ saturation at location [-]',
    'mass_CO2_reservoir': 'Mass of CO$_2$ in the reservoir (kg)',
    # Wellbore and adapter components
    'CO2_aquifer': 'Leakage rate of CO$_2$ to aquifer (kg/s)',
    'CO2_atm': 'Leakage rate of CO$_2$ to atmosphere (kg/s)',
    'brine_aquifer': 'Leakage rate of brine to aquifer (kg/s)',
    'brine_atm': 'Leakage rate of brine to atmosphere (kg/s)',
    'mass_CO2_aquifer': 'Mass of CO$_2$ leaked to aquifer (kg)',
    'mass_brine_aquifer': 'Mass of brine leaked to aquifer (kg)',
    'mass_CO2_atm': 'Mass of CO$_2$ leaked to atmosphere (kg)',
    'mass_brine_atm': 'Mass of brine leaked to atmosphere (kg)',
    'mass_CO2_gas_aquifer': 'Mass of CO$_2$ gas leaked to aquifer (kg)',
    'mass_methane_oil_aquifer': 'Mass of CH$_4$ oil leaked to aquifer (kg)',
    'mass_methane_gas_aquifer': 'Mass of CH$_4$ gas leaked to aquifer (kg)',
    'mass_oil_aquifer': 'Total mass of oil leaked to aquifer (kg)',
    'mass_gas_aquifer': 'Total mass of gas leaked to aquifer (kg)',
    # Fault flow and Seal horizon components
    'CO2_aquifer_total': 'Cumulative leakage rate of CO$_2$ to aquifer (kg/s)',
    'brine_aquifer_total': 'Cumulative leakage rate of brine to aquifer (kg/s)',
    'mass_CO2_aquifer_total': 'Cumulative mass of CO$_2$ leaked to aquifer (kg)',
    'mass_brine_aquifer_total': 'Cumulative mass of brine leaked to aquifer (kg)',
    # Aquifer components
    'Benzene_volume': 'Volume of plume above benzene threshold (m$^3$)',
    'Naphthalene_volume': 'Volume of plume above naphthalene threshold (m$^3$)',
    'Phenol_volume': 'Volume of plume above phenol threshold (m$^3$)',
    'As_volume': 'Volume of plume above arsenic threshold (m$^3$)',
    'Pb_volume': 'Volume of plume above lead threshold (m$^3$)',
    'Cd_volume': 'Volume of plume above cadmium threshold (m$^3$)',
    'Ba_volume': 'Volume of plume above barium threshold (m$^3$)',
    'Flux': 'CO$_2$ leakage rate to atmosphere (kg/s)',
    # FutureGen2 aquifer components
    'TDS_volume': 'Volume of plume above TDS threshold (m$^3$)',
    'TDS_dx': 'Length of plume above TDS threshold (m)',
    'TDS_dy': 'Width of plume above TDS threshold (m)',
    'TDS_dz': 'Height of plume above TDS threshold (m)',
    'pH_volume': 'Volume of plume below pH threshold (m$^3$)',
    'pH_dx': 'Length of plume below pH threshold (m)',
    'pH_dy': 'Width of plume below pH threshold (m)',
    'pH_dz': 'Height of plume below pH threshold (m)',
    'Pressure_volume': 'Volume of plume above baseline pressure change (m$^3$)',
    'Pressure_dx': 'Length of plume above baseline pressure change (m)',
    'Pressure_dy': 'Width of plume above baseline pressure change (m)',
    'Pressure_dz': 'Height of plume above baseline pressure change (m)',
    'Dissolved_CO2_volume': 'Volume of plume above baseline dissolved CO$_2$ (m$^3$)',
    'Dissolved_CO2_dx': 'Length of plume above baseline dissolved CO$_2$ (m)',
    'Dissolved_CO2_dy': 'Width of plume above baseline dissolved CO$_2$ (m)',
    'Dissolved_CO2_dz': 'Height of plume above baseline dissolved CO$_2$ (m)',
    'Temperature_volume': 'Volume of plume above baseline temperature change (m$^3$)',
    'Temperature_dx': 'Length of plume above baseline temperature change (m)',
    'Temperature_dy': 'Width of plume above baseline temperature change (m)',
    'Temperature_dz': 'Height of plume above baseline temperature change (m)',
    # Generic aquifer component
    'Dissolved_salt_volume': 'Volume of plume above baseline salt mass fraction (m$^3$)',
    'Dissolved_salt_dr': 'Radius of plume above baseline salt mass fraction (m)',
    'Dissolved_salt_dz': 'Height of plume above baseline salt mass fraction (m)',
    'Dissolved_salt_mass_fraction': 'Mass fraction of salt in aquifer [-]',
    'Dissolved_CO2_dr': 'Radius of plume above baseline dissolved CO$_2$ (m)',
    'Dissolved_CO2_dz': 'Height of plume above baseline dissolved CO$_2$ (m)',
    'Dissolved_CO2_mass_fraction': 'Mass fraction of CO$_2$ in aquifer [-]',
    # Plume stability component
    'pressure_areas': 'Pressure plume area (m$^2$)',
    'pressure_areas_dt': 'Change in pressure plume area (m$^2$/year)',
    'pressure_mobility': 'Velocity of pressure plume centroid (m/year)',
    'pressure_mobility_angles': 'Direction of pressure plume centroid [-]',
    'pressure_spreading': 'Dispersion of pressure plume (m$^2$/year)',
    'pressure_spreading_angles': 'Direction of pressure plume dispersion [-]',
    'CO2saturation_areas': 'CO$_2$ saturation plume area (m$^2$)',
    'CO2saturation_areas_dt': 'Change in CO$_2$ saturation plume area (m$^2$/year)',
    'CO2saturation_mobility': 'Velocity of CO$_2$ saturation plume centroid (m/year)',
    'CO2saturation_mobility_angles': 'Direction of CO$_2$ saturation plume centroid [-]',
    'CO2saturation_spreading': 'Dispersion of CO$_2$ saturation plume (m$^2$/year)',
    'CO2saturation_spreading_angles': 'Direction of CO$_2$ saturation plume dispersion [-]',
    }

Y_LABEL_SUBPLOT_DICT = Y_LABEL_DICT.copy()
for ind, (key, val) in enumerate(zip(additional_keys, additional_ylabels)):
    Y_LABEL_SUBPLOT_DICT[key] = val

Y_LABEL_SPLIT_DICT = {
    'CO': Y_LABEL_DICT['CO2_aquifer'],
    'brine_aquifer': Y_LABEL_DICT['brine_aquifer'],
    'brine_aquifer_cell_': Y_LABEL_DICT['brine_aquifer'],
    'brine_aquifer_segm_': Y_LABEL_DICT['brine_aquifer'],
    'mass_CO': Y_LABEL_DICT['mass_CO2_aquifer'],
    'mass_brine_aquifer': Y_LABEL_DICT['mass_brine_aquifer'],
    'mass_brine_aquifer_cell_': Y_LABEL_DICT['mass_brine_aquifer'],
    'mass_brine_aquifer_segm_': Y_LABEL_DICT['mass_brine_aquifer'],
    'Dissolved_CO': Y_LABEL_DICT['Dissolved_CO2_mass_fraction'],
    'Dissolved_salt_mass_fraction_coord_': Y_LABEL_DICT['Dissolved_salt_mass_fraction']}

# These ylabels are used if the one row case is too long on the figure. Because single
# plots can display results for multiple aquifers, I have removed the aquifer numbers here.
Y_LABEL_2ROWS_DICT = {
    # Reservoir components
    'pressure': 'Reservoir pressure\nat location (Pa)',
    'CO2saturation': 'Reservoir CO$_2$\nsaturation at location [-]',
    'mass_CO2_reservoir': 'Mass of CO$_2$\nin the reservoir (kg)',
    # Wellbore and adapter components
    'CO2_aquifer': 'Leakage rate of\nCO$_2$ to aquifer (kg/s)',
    'CO2_atm': 'Leakage rate of\nCO$_2$ to atmosphere (kg/s)',
    'brine_aquifer': 'Leakage rate of\nbrine to aquifer (kg/s)',
    'brine_atm': 'Leakage rate of\nbrine to atmosphere (kg/s)',
    'mass_CO2_aquifer': 'Mass of CO$_2$\nleaked to aquifer (kg)',
    'mass_brine_aquifer': 'Mass of brine\nleaked to aquifer (kg)',
    'mass_CO2_atm': 'Mass of CO$_2$\nleaked to atmosphere (kg)',
    'mass_brine_atm': 'Mass of brine\nleaked to atmosphere (kg)',
    'mass_CO2_gas_aquifer': 'Mass of CO$_2$ gas\nleaked to aquifer (kg)',
    'mass_methane_oil_aquifer': 'Mass of CH$_4$ oil\nleaked to aquifer (kg)',
    'mass_methane_gas_aquifer': 'Mass of CH$_4$ gas\nleaked to aquifer (kg)',
    'mass_oil_aquifer': 'Total mass of oil\nleaked to aquifer (kg)',
    'mass_gas_aquifer': 'Total mass of gas\nleaked to aquifer (kg)',
    # Fault flow and Seal horizon components
    'CO2_aquifer_total': 'Cumulative leakage rate\nof CO$_2$ to aquifer (kg/s)',
    'brine_aquifer_total': 'Cumulative leakage rate\nof brine to aquifer (kg/s)',
    'mass_CO2_aquifer_total': 'Cumulative mass of\nCO$_2$ leaked to aquifer (kg)',
    'mass_brine_aquifer_total': 'Cumulative mass of\nbrine leaked to aquifer (kg)',
    # Aquifer components
    'Benzene_volume': 'Volume of plume\nabove benzene threshold (m$^3$)',
    'Naphthalene_volume': 'Volume of plume\nabove naphthalene threshold (m$^3$)',
    'Phenol_volume': 'Volume of plume\nabove phenol threshold (m$^3$)',
    'As_volume': 'Volume of plume\nabove arsenic threshold (m$^3$)',
    'Pb_volume': 'Volume of plume\nabove lead threshold (m$^3$)',
    'Cd_volume': 'Volume of plume\nabove cadmium threshold (m$^3$)',
    'Ba_volume': 'Volume of plume\nabove barium threshold (m$^3$)',
    'Flux': 'CO$_2$ leakage rate\nto atmosphere (kg/s)',
    # FutureGen2 aquifer components
    'TDS_volume': 'Volume of plume\nabove TDS threshold (m$^3$)',
    'TDS_dx': 'Length of plume\nabove TDS threshold (m)',
    'TDS_dy': 'Width of plume\nabove TDS threshold (m)',
    'TDS_dz': 'Height of plume\nabove TDS threshold (m)',
    'pH_volume': 'Volume of plume\nnbelow pH threshold (m$^3$)',
    'pH_dx': 'Length of plume\nbelow pH threshold (m)',
    'pH_dy': 'Width of plume\nbelow pH threshold (m)',
    'pH_dz': 'Height of plume\nbelow pH threshold (m)',
    'Pressure_volume': 'Volume of plume above\nbaseline pressure change (m$^3$)',
    'Pressure_dx': 'Length of plume above\nbaseline pressure change (m)',
    'Pressure_dy': 'Width of plume above\nbaseline pressure change (m)',
    'Pressure_dz': 'Height of plume above\nbaseline pressure change (m)',
    'Dissolved_CO2_volume': 'Volume of plume above\nbaseline dissolved CO$_2$ (m$^3$)',
    'Dissolved_CO2_dx': 'Length of plume above\nbaseline dissolved CO$_2$ (m)',
    'Dissolved_CO2_dy': 'Width of plume above\nbaseline dissolved CO$_2$ (m)',
    'Dissolved_CO2_dz': 'Height of plume above\nbaseline dissolved CO$_2$ (m)',
    'Temperature_volume': 'Volume of plume above\nbaseline temperature change (m$^3$)',
    'Temperature_dx': 'Length of plume above\nbaseline temperature change (m)',
    'Temperature_dy': 'Width of plume above\nbaseline temperature change (m)',
    'Temperature_dz': 'Height of plume above\nbaseline temperature change (m)',
    # Generic aquifer component
    'Dissolved_salt_volume': 'Volume of plume above\nbaseline salt mass fraction (m$^3$)',
    'Dissolved_salt_dr': 'Radius of plume above\nbaseline salt mass fraction (m)',
    'Dissolved_salt_dz': 'Height of plume above\nbaseline salt mass fraction (m)',
    'Dissolved_salt_mass_fraction': 'Mass fraction of\nsalt in aquifer [-]',
    'Dissolved_CO2_dr': 'Radius of plume above\nbaseline dissolved CO$_2$ (m)',
    'Dissolved_CO2_dz': 'Height of plume above\nbaseline dissolved CO$_2$ (m)',
    'Dissolved_CO2_mass_fraction': 'Mass fraction of\nCO$_2$ in aquifer [-]',
    # Plume stability component
    'pressure_areas': 'Pressure\nplume area (m$^2$)',
    'pressure_areas_dt': 'Change in pressure\nplume area (m$^2$/year)',
    'pressure_mobility': 'Velocity of pressure\nplume centroid (m/year)',
    'pressure_mobility_angles': 'Direction of pressure\nplume centroid [-]',
    'pressure_spreading': 'Dispersion of\npressure plume (m$^2$/year)',
    'pressure_spreading_angles': 'Direction of pressure\nplume dispersion [-]',
    'CO2saturation_areas': 'CO$_2$ saturation\nplume area (m$^2$)',
    'CO2saturation_areas_dt': 'Change in CO$_2$\nsaturation plume area (m$^2$/year)',
    'CO2saturation_mobility': 'Velocity of CO$_2$\nsaturation plume centroid (m/year)',
    'CO2saturation_mobility_angles': 'Direction of CO$_2$\nsaturation plume centroid [-]',
    'CO2saturation_spreading': 'Dispersion of CO$_2$\nsaturation plume (m$^2$/year)',
    'CO2saturation_spreading_angles': 'Direction of CO$_2$\nsaturation plume dispersion [-]',
    }

Y_LABEL_2ROWS_SUBPLOT_DICT = Y_LABEL_2ROWS_DICT.copy()
for ind, (key, val) in enumerate(zip(additional_keys, additional_ylabels_2rows)):
    Y_LABEL_2ROWS_SUBPLOT_DICT[key] = val

Y_LABEL_2ROWS_SPLIT_DICT = {
    'CO': Y_LABEL_2ROWS_DICT['CO2_aquifer'],
    'brine_aquifer': Y_LABEL_2ROWS_DICT['brine_aquifer'],
    'brine_aquifer_cell_': Y_LABEL_2ROWS_DICT['brine_aquifer'],
    'brine_aquifer_segm_': Y_LABEL_2ROWS_DICT['brine_aquifer'],
    'mass_CO': Y_LABEL_2ROWS_DICT['mass_CO2_aquifer'],
    'mass_brine_aquifer': Y_LABEL_2ROWS_DICT['mass_brine_aquifer'],
    'mass_brine_aquifer_cell_': Y_LABEL_2ROWS_DICT['mass_brine_aquifer'],
    'mass_brine_aquifer_segm_': Y_LABEL_2ROWS_DICT['mass_brine_aquifer'],
    'Dissolved_CO': Y_LABEL_2ROWS_DICT['Dissolved_CO2_mass_fraction'],
    'Dissolved_salt_mass_fraction_coord_': Y_LABEL_2ROWS_DICT['Dissolved_salt_mass_fraction']}

# Dictionary with the title names for different component output. Because single
# plots can display results for multiple aquifers, I have removed the aquifer numbers here.
TITLE_DICT = {
    # Reservoir components
    'pressure': 'Reservoir Pressure',
    'CO2saturation': 'Reservoir CO$_2$ Saturation',
    'mass_CO2_reservoir': 'Mass of CO$_2$ in the Reservoir',
    # Wellbore and adapter components
    'CO2_aquifer': 'CO$_2$ Leakage Rate to Aquifer',
    'CO2_atm': 'CO$_2$ Leakage Rate to Atmosphere',
    'brine_aquifer': 'Brine Leakage Rate to Aquifer',
    'brine_atm': 'Brine Leakage Rate to Atmosphere',
    'mass_CO2_aquifer': 'CO$_2$ Mass Leaked to Aquifer',
    'mass_brine_aquifer': 'Brine Mass Leaked to Aquifer',
    'mass_CO2_atm': 'CO$_2$ Mass Leaked to Atmosphere',
    'mass_brine_atm': 'Brine Mass Leaked to Atmosphere',
    'mass_CO2_gas_aquifer': 'CO$_2$ Gas Leaked to Aquifer',
    'mass_methane_oil_aquifer': 'CH$_4$ Oil Leaked to Aquifer',
    'mass_methane_gas_aquifer': 'CH$_4$ Gas Leaked to Aquifer',
    'mass_oil_aquifer': 'Oil Leaked to Aquifer',
    'mass_gas_aquifer': 'Gas Leaked to Aquifer',
    # Fault flow and Seal horizon components
    'CO2_aquifer_total': 'Cumulative CO$_2$ Leakage Rate to Aquifer',
    'brine_aquifer_total': 'Cumulative Brine Leakage Rate to Aquifer',
    'mass_CO2_aquifer_total': 'Cumulative CO$_2$ Mass Leaked to Aquifer',
    'mass_brine_aquifer_total': 'Cumulative Brine Mass Leaked to Aquifer',
    # Aquifer components
    'Benzene_volume': 'Benzene Plume Volume',
    'Naphthalene_volume': 'Naphthalene Plume Volume',
    'Phenol_volume': 'Phenol Plume Volume',
    'As_volume': 'Arsenic Plume Volume',
    'Pb_volume': 'Lead Plume Volume',
    'Cd_volume': 'Cadmium Plume Volume',
    'Ba_volume': 'Barium Plume Volume',
    'Flux': 'CO$_2$ Leakage Rate to Atmosphere',
    'Dissolved_salt_volume': 'Dissolved Salt Volume',
    'Dissolved_CO2_mass_fraction': 'Dissolved CO$_2$ Mass Fraction',
    'Dissolved_salt_mass_fraction': 'Dissolved Salt Mass Fraction',
    # FutureGen2 aquifer components
    'TDS_volume': 'TDS Plume Volume',
    'TDS_dx': 'TDS Plume Length',
    'TDS_dy': 'TDS Plume Width',
    'TDS_dz': 'TDS Plume Height',
    'pH_volume': 'Low-pH Plume Volume',
    'pH_dx': 'Low-pH Plume Length',
    'pH_dy': 'Low-pH Plume Width',
    'pH_dz': 'Low-pH Plume Height',
    'Pressure_volume': 'High-Pressure Plume Volume',
    'Pressure_dx': 'High-Pressure Plume Length',
    'Pressure_dy': 'High-Pressure Plume Width',
    'Pressure_dz': 'High-Pressure Plume Height',
    'Dissolved_CO2_volume': 'Dissolved-CO$_2$ Plume Volume',
    'Dissolved_CO2_dx': 'Dissolved-CO$_2$ Plume Length',
    'Dissolved_CO2_dy': 'Dissolved-CO$_2$ Plume Width',
    'Dissolved_CO2_dz': 'Dissolved-CO$_2$ Plume Height',
    'Temperature_volume': 'Temperature-Change Plume Volume',
    'Temperature_dx': 'Temperature-Change Plume Length',
    'Temperature_dy': 'Temperature-Change Plume Width',
    'Temperature_dz': 'Temperature-Change Plume Height',
    # Generic aquifer component
    # key 'Dissolved_salt_volume' is already present for aquifer components above
    'Dissolved_salt_dr': 'Dissolved Salt Plume Radius',
    'Dissolved_salt_dz': 'Dissolved Salt Plume Height',
    # key 'Dissolved_salt_mass_fraction' is already present for aquifer components above
    'Dissolved_CO2_dr': 'Dissolved CO$_2$ Plume Radius',
    'Dissolved_CO2_dz': 'Dissolved CO$_2$ Plume Height',
    # key 'Dissolved_CO2_mass_fraction' is already present for aquifer components above
    # Plume stability component
    'pressure_areas': 'Pressure Plume Area',
    'pressure_areas_dt': 'Change in Pressure Plume Area',
    'pressure_mobility': 'Velocity of Pressure Plume',
    'pressure_mobility_angles': 'Direction of Pressure Plume',
    'pressure_spreading': 'Dispersion of Pressure Plume',
    'pressure_spreading_angles': 'Direction of Pressure Plume Dispersion',
    'CO2saturation_areas': 'CO$_2$ Saturation Plume Area',
    'CO2saturation_areas_dt': 'Change in CO$_2$ Saturation Plume Area',
    'CO2saturation_mobility': 'Velocity of CO$_2$ Saturation Plume',
    'CO2saturation_mobility_angles': 'Direction of CO$_2$ Saturation Plume',
    'CO2saturation_spreading': 'Dispersion of CO$_2$ Saturation Plume',
    'CO2saturation_spreading_angles': 'Direction of CO$_2$ Saturation Plume Dispersion',
    }

TITLE_SUBPLOT_DICT = TITLE_DICT.copy()
for ind, (key, val) in enumerate(zip(additional_keys, additional_titles)):
    TITLE_SUBPLOT_DICT[key] = val

TITLE_SPLIT_DICT = {
    'CO': TITLE_DICT['CO2_aquifer'],
    'brine_aquifer': TITLE_DICT['brine_aquifer'],
    'brine_aquifer_cell_': TITLE_DICT['brine_aquifer'],
    'brine_aquifer_segm_': TITLE_DICT['brine_aquifer'],
    'mass_CO': TITLE_DICT['mass_CO2_aquifer'],
    'mass_brine_aquifer': TITLE_DICT['mass_brine_aquifer'],
    'mass_brine_aquifer_cell_': TITLE_DICT['mass_brine_aquifer'],
    'mass_brine_aquifer_segm_': TITLE_DICT['mass_brine_aquifer'],
    'Dissolved_CO': TITLE_DICT['Dissolved_CO2_mass_fraction'],
    'Dissolved_salt_mass_fraction_coord_': TITLE_DICT['Dissolved_salt_mass_fraction']}

# Dictionary with the partial legend names for different component output.
# Some graphs can display results for multiple aquifers, and this dictionary is
# used to differentiate between such results. Otherwise, the LEGEND_DICT entries are empty.
LEGEND_DICT = {
    # Reservoir components
    'pressure': '',
    'CO2saturation': '',
    'mass_CO2_reservoir': '',
    # Wellbore and adapter components
    'CO2_aquifer': '',
    'CO2_atm': '',
    'brine_aquifer': '',
    'brine_atm': '',
    'mass_CO2_aquifer': '',
    'mass_brine_aquifer': '',
    'mass_CO2_atm': '',
    'mass_brine_atm': '',
    'mass_CO2_gas_aquifer': '',
    'mass_methane_oil_aquifer': '',
    'mass_methane_gas_aquifer': '',
    'mass_oil_aquifer': '',
    'mass_gas_aquifer': '',
    # Fault flow and Seal horizon components
    'CO2_aquifer_total': '',
    'brine_aquifer_total': '',
    'mass_CO2_aquifer_total': '',
    'mass_brine_aquifer_total': '',
    # Aquifer components
    'Benzene_volume': '',
    'Naphthalene_volume': '',
    'Phenol_volume': '',
    'As_volume': '',
    'Pb_volume': '',
    'Cd_volume': '',
    'Ba_volume': '',
    'Flux': '',
    # FutureGen2 aquifer components
    'TDS_volume': '',
    'TDS_dx': '',
    'TDS_dy': '',
    'TDS_dz': '',
    'pH_volume': '',
    'pH_dx': '',
    'pH_dy': '',
    'pH_dz': '',
    'Pressure_volume': '',
    'Pressure_dx': '',
    'Pressure_dy': '',
    'Pressure_dz': '',
    'Dissolved_CO2_volume': '',
    'Dissolved_CO2_dx': '',
    'Dissolved_CO2_dy': '',
    'Dissolved_CO2_dz': '',
    'Temperature_volume': '',
    'Temperature_dx': '',
    'Temperature_dy': '',
    'Temperature_dz': '',
    # Generic aquifer component
    'Dissolved_salt_volume': '',
    'Dissolved_salt_dr': '',
    'Dissolved_salt_dz': '',
    'Dissolved_salt_mass_fraction': '',
    'Dissolved_CO2_dr': '',
    'Dissolved_CO2_dz': '',
    'Dissolved_CO2_mass_fraction': '',
    # Plume stability component
    'pressure_areas': '',
    'pressure_areas_dt': '',
    'pressure_mobility': '',
    'pressure_mobility_angles': '',
    'pressure_spreading': '',
    'pressure_spreading_angles': '',
    'CO2saturation_areas': '',
    'CO2saturation_areas_dt': '',
    'CO2saturation_mobility': '',
    'CO2saturation_mobility_angles': '',
    'CO2saturation_spreading': '',
    'CO2saturation_spreading_angles': '',
    }
