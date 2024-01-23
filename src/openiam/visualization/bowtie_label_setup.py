# -*- coding: utf-8 -*-
"""
Module containing dictionaries defining the units for each type of metric. This 
dictionary is used in bowtie.py.
"""
additional_keys = [
    'brine_aquifer{}'.format(unitRef) for unitRef in range(1, 30)
    ] + ['CO2_aquifer{}'.format(unitRef) for unitRef in range(1, 30)
         ] + ['mass_brine_aquifer{}'.format(unitRef) for unitRef in range(1, 30)
              ] + ['mass_CO2_aquifer{}'.format(unitRef) for unitRef in range(1, 30)]

additional_units = [
    'kg/s' for unitRef in range(1, 30)] + [
        'kg/s' for unitRef in range(1, 30)
        ] + ['kg' for unitRef in range(1, 30)
             ] + ['kg' for unitRef in range(1, 30)]

additional_labels = [
    'Brine Leakage Rate to Aquifer' for unitRef in range(1, 30)] + [
        'CO$_2$ Leakage Rate to Aquifer' for unitRef in range(1, 30)
        ] + ['Brine Mass Leaked to Aquifer' for unitRef in range(1, 30)
             ] + ['CO$_2$ Mass Leaked to Aquifer' for unitRef in range(1, 30)]

UNIT_DICT = {
    # Reservoir components
    'pressure': 'Pa',
    'CO2saturation': '[-]',
    'mass_CO2_reservoir': 'kg',
    # Wellbore and adapter components
    'CO2_aquifer': 'kg/s',
    'CO2_atm': 'kg/s',
    'brine_aquifer': 'kg/s',
    'brine_atm': 'kg/s',
    'mass_CO2_aquifer': 'kg',
    'mass_brine_aquifer': 'kg',
    'mass_CO2_atm': 'kg',
    'mass_brine_atm': 'kg',
    'mass_CO2_gas_aquifer': 'kg',
    'mass_methane_oil_aquifer': 'kg',
    'mass_methane_gas_aquifer': 'kg',
    'mass_oil_aquifer': 'kg',
    'mass_gas_aquifer': 'kg',
    # Fault flow and Seal horizon components
    'CO2_aquifer_total': 'kg/s',
    'brine_aquifer_total': 'kg/s',
    'mass_CO2_aquifer_total': 'kg',
    'mass_brine_aquifer_total': 'kg',
    # Aquifer components
    'Benzene_volume': 'm$^3$',
    'Naphthalene_volume': 'm$^3$',
    'Phenol_volume': 'm$^3$',
    'As_volume': 'm$^3$',
    'Pb_volume': 'm$^3$',
    'Cd_volume': 'm$^3$',
    'Ba_volume': 'm$^3$',
    'Flux': 'kg/s',
    # FutureGen2 aquifer components
    'TDS_volume': 'm$^3$',
    'TDS_dx': 'm',
    'TDS_dy': 'm',
    'TDS_dz': 'm',
    'pH_volume': 'm$^3$',
    'pH_dx': 'm',
    'pH_dy': 'm',
    'pH_dz': 'm',
    'Pressure_volume': 'm$^3$',
    'Pressure_dx': 'm',
    'Pressure_dy': 'm',
    'Pressure_dz': 'm',
    'Dissolved_CO2_volume': 'm$^3$',
    'Dissolved_CO2_dx': 'm',
    'Dissolved_CO2_dy': 'm',
    'Dissolved_CO2_dz': 'm',
    'Temperature_volume': 'm$^3$',
    'Temperature_dx': 'm',
    'Temperature_dy': 'm',
    'Temperature_dz': 'm',
    # Generic aquifer component
    'Dissolved_salt_volume': 'm$^3$',
    'Dissolved_salt_dr': 'm',
    'Dissolved_salt_dz': 'm',
    'Dissolved_salt_mass_fraction': '[-]',
    'Dissolved_CO2_dr': 'm',
    'Dissolved_CO2_mass_fraction': '[-]',
    # Plume stability component
    'pressure_areas': 'm$^2$',
    'pressure_areas_dt': 'm$^2$/year',
    'pressure_mobility': 'm/year',
    'pressure_mobility_angles': '[-]',
    'pressure_spreading': 'm$^2$/year',
    'pressure_spreading_angles': '[-]',
    'CO2saturation_areas': 'm$^2$',
    'CO2saturation_areas_dt': 'm$^2$/year',
    'CO2saturation_mobility': 'm/year',
    'CO2saturation_mobility_angles': '[-]',
    'CO2saturation_spreading': 'm$^2$/year',
    'CO2saturation_spreading_angles': '[-]',
    }

for ind, (key, val) in enumerate(zip(additional_keys, additional_units)):
    UNIT_DICT[key] = val

METRIC_LABEL_DICT = {
    # Reservoir components
    'pressure': 'Reservoir Pressure',
    'CO2saturation': 'CO$_2$ Saturation',
    'mass_CO2_reservoir': 'CO$_2$ Mass in Reservoir',
    # Wellbore and adapter components
    'CO2_aquifer': 'CO$_2$ Leakage Rate to Aquifer',
    'CO2_atm': 'CO$_2$ Leakage Rate to Atmosphere',
    'brine_aquifer': 'Brine Leakage Rate to Aquifer',
    'brine_atm': 'Brine Leakage Rate to Atmosphere',
    'mass_CO2_aquifer': 'CO$_2$ Mass Leaked to Aquifer',
    'mass_brine_aquifer': 'Leaked Brine Mass',
    'mass_CO2_atm': 'CO$_2$ Mass Leaked to Atmosphere',
    'mass_brine_atm': 'Brine Mass Leaked to Atmosphere',
    'mass_CO2_gas_aquifer': 'Mass of CO$_2$ Gas Leaked to Aquifer',
    'mass_methane_oil_aquifer': 'Mass of CH$_4$ Oil Leaked to Aquifer',
    'mass_methane_gas_aquifer': 'Mass of CH$_4$ Gas Leaked to Aquifer',
    'mass_oil_aquifer': 'Total Oil Mass Leaked to Aquifer',
    'mass_gas_aquifer': 'Total Gas Mass Leaked to Aquifer',
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
    # FutureGen2 aquifer components
    'TDS_volume': 'TDS Plume Volume',
    'TDS_dx': 'TDS Plume Length',
    'TDS_dy': 'TDS Plume Width',
    'TDS_dz': 'TDS Plume Height',
    'pH_volume': 'pH Plume Volume',
    'pH_dx': 'pH Plume Length',
    'pH_dy': 'pH Plume Width',
    'pH_dz': 'pH Plume Height',
    'Pressure_volume': 'Pressure Plume Volume',
    'Pressure_dx': 'Pressure Plume Length',
    'Pressure_dy': 'Pressure Plume Width',
    'Pressure_dz': 'Pressure Plume Height',
    'Dissolved_CO2_volume': 'Dissolved CO$_2$ Plume Volume',
    'Dissolved_CO2_dx': 'Dissolved CO$_2$ Plume Length',
    'Dissolved_CO2_dy': 'Dissolved CO$_2$ Plume Width',
    'Dissolved_CO2_dz': 'Dissolved CO$_2$ Plume Height',
    'Temperature_volume': 'Temperature Plume Volume',
    'Temperature_dx': 'Temperature Plume Length',
    'Temperature_dy': 'Temperature Plume Width',
    'Temperature_dz': 'Temperature Plume Height',
    # Generic aquifer component
    'Dissolved_salt_volume': 'Dissolved Salt Plume Volume',
    'Dissolved_salt_dr': 'Dissolved Salt Plume Radius',
    'Dissolved_salt_dz': 'Dissolved Salt Plume Height',
    'Dissolved_salt_mass_fraction': 'Dissolved Salt Mass Fraction',
    'Dissolved_CO2_dr': 'Dissolved CO$_2$ Plume Radius',
    'Dissolved_CO2_mass_fraction': 'Dissolved CO$_2$ Mass fraction',
    # Plume stability component
    'pressure_areas': 'Pressure Plume Area',
    'pressure_areas_dt': 'Change in Pressure Plume Area',
    'pressure_mobility': 'Pressure Plume Velocity',
    'pressure_mobility_angles': 'Pressure Plume Direction',
    'pressure_spreading': 'Pressure Plume Dispersion',
    'pressure_spreading_angles': 'Direction of Pressure Plume Dispersion',
    'CO2saturation_areas': 'CO$_2$ Saturation Plume Area',
    'CO2saturation_areas_dt': 'Change in CO$_2$ Saturation Plume Area',
    'CO2saturation_mobility': 'CO$_2$ Saturation Plume Velocity',
    'CO2saturation_mobility_angles': 'CO$_2$ Saturation Plume Direction',
    'CO2saturation_spreading': 'CO$_2$ Saturation Plume Dispersion',
    'CO2saturation_spreading_angles': 'Direction of CO$_2$ Saturation Plume Dispersion',
    }

METRIC_LABEL_CONCISE_DICT = {
    # Reservoir components
    'pressure': 'Pressure',
    'CO2saturation': 'CO$_2$ Saturation',
    'mass_CO2_reservoir': 'CO$_2$ Mass',
    # Wellbore and adapter components
    'CO2_aquifer': 'CO$_2$ Leakage Rate',
    'CO2_atm': 'CO$_2$ Leakage Rate',
    'brine_aquifer': 'Brine Leakage Rate',
    'brine_atm': 'Brine Leakage Rate',
    'mass_CO2_aquifer': 'CO$_2$ Mass',
    'mass_brine_aquifer': 'Brine Mass',
    'mass_CO2_atm': 'CO$_2$ Mass',
    'mass_brine_atm': 'Brine Mass',
    'mass_CO2_gas_aquifer': 'Mass of CO$_2$ Gas',
    'mass_methane_oil_aquifer': 'Mass of CH$_4$ Oil',
    'mass_methane_gas_aquifer': 'Mass of CH$_4$ Gas',
    'mass_oil_aquifer': 'Total Oil Mass',
    'mass_gas_aquifer': 'Total Gas Mass',
    # Fault flow and Seal horizon components
    'CO2_aquifer_total': 'CO$_2$ Leakage Rate',
    'brine_aquifer_total': 'Brine Leakage Rate',
    'mass_CO2_aquifer_total': 'CO$_2$ Mass',
    'mass_brine_aquifer_total': 'Brine Mass',
    # Aquifer components
    'Benzene_volume': 'Plume Volume',
    'Naphthalene_volume': 'Plume Volume',
    'Phenol_volume': 'Plume Volume',
    'As_volume': 'Plume Volume',
    'Pb_volume': 'Plume Volume',
    'Cd_volume': 'Plume Volume',
    'Ba_volume': 'Plume Volume',
    'Flux': 'CO$_2$ Leakage Rate',
    # FutureGen2 aquifer components
    'TDS_volume': 'Plume Volume',
    'TDS_dx': 'Plume Length',
    'TDS_dy': 'Plume Width',
    'TDS_dz': 'Plume Height',
    'pH_volume': 'Plume Volume',
    'pH_dx': 'Plume Length',
    'pH_dy': 'Plume Width',
    'pH_dz': 'Plume Height',
    'Pressure_volume': 'Plume Volume',
    'Pressure_dx': 'Plume Length',
    'Pressure_dy': 'Plume Width',
    'Pressure_dz': 'Plume Height',
    'Dissolved_CO2_volume': 'Plume Volume',
    'Dissolved_CO2_dx': 'Plume Length',
    'Dissolved_CO2_dy': 'Plume Width',
    'Dissolved_CO2_dz': 'Plume Height',
    'Temperature_volume': 'Plume Volume',
    'Temperature_dx': 'Plume Length',
    'Temperature_dy': 'Plume Width',
    'Temperature_dz': 'Plume Height',
    # Generic aquifer component
    'Dissolved_salt_volume': 'Plume Volume',
    'Dissolved_salt_dr': 'Plume Radius',
    'Dissolved_salt_dz': 'Plume Height',
    'Dissolved_salt_mass_fraction': 'Mass Fraction',
    'Dissolved_CO2_dr': 'Plume Radius',
    'Dissolved_CO2_mass_fraction': 'Mass fraction',
    # Plume stability component
    'pressure_areas': 'Plume Area',
    'pressure_areas_dt': 'Change in Plume Area',
    'pressure_mobility': 'Plume Velocity',
    'pressure_mobility_angles': 'Plume Direction',
    'pressure_spreading': 'Plume Dispersion',
    'pressure_spreading_angles': 'Plume Dispersion Direction',
    'CO2saturation_areas': 'Plume Area',
    'CO2saturation_areas_dt': 'Change in Plume Area',
    'CO2saturation_mobility': 'Plume Velocity',
    'CO2saturation_mobility_angles': 'Plume Direction',
    'CO2saturation_spreading': 'Plume Dispersion',
    'CO2saturation_spreading_angles': 'Plume Dispersion Direction',
    }
