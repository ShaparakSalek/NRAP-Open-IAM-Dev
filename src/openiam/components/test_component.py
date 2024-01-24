# -*- coding: utf-8 -*-
"""
This script simplifies execution of the tests within each component module.

Created November, 2023
"""
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from openiam.components.alluvium_aquifer_component import test_alluvium_aquifer_component
from openiam.components.alluvium_aquifer_lf_component import test_alluvium_aquifer_lf_component
from openiam.components.analytical_reservoir_component import test_analytical_reservoir_component
from openiam.components.atmRom_component import test_atmRom_component
from openiam.components.carbonate_aquifer_component import test_carbonate_aquifer_component
from openiam.components.cemented_wellbore_component import test_cemented_wellbore_component
from openiam.components.cemented_wellbore_wr_component import test_cemented_wellbore_wr_component
from openiam.components.chemical_well_sealing import test_chemical_well_sealing
from openiam.components.configurer_component import test_configurer_component
from openiam.components.deep_alluvium_aquifer_component import test_deep_alluvium_aquifer_component
from openiam.components.deep_alluvium_aquifer_ml_component import test_deep_alluvium_aquifer_ml_component
from openiam.components.fault_flow_component import test_fault_flow_component
from openiam.components.fault_leakage_component import test_fault_leakage_component
from openiam.components.futuregen2_aquifer_component import test_futuregen2_aquifer_component
from openiam.components.futuregen2_azmi_component import test_futuregen2_azmi_component
from openiam.components.generalized_flow_rate_component import test_generalized_flow_rate_component
from openiam.components.generic_aquifer_component import test_generic_aquifer_component
from openiam.components.generic_reservoir_component import test_generic_reservoir_component
from openiam.components.hydrocarbon_leakage_component import test_hydrocarbon_leakage_component
from openiam.components.kimberlina_wellbore_component import test_kimberlina_wellbore_component
from openiam.components.location_generator import test_location_generator
from openiam.components.lookup_table_reservoir_component import test_lookup_table_reservoir_component
from openiam.components.monitoring_scheduler_component import test_monitoring_scheduler_component
from openiam.components.monitoring_tool_component import test_monitoring_tool_component
from openiam.components.multisegmented_wellbore_component import test_multisegmented_wellbore_component
from openiam.components.open_wellbore_component import test_open_wellbore_component
from openiam.components.parameter_setup_component import test_parameter_setup_component
from openiam.components.plume_stability_component import test_plume_stability_component
from openiam.components.rate_to_mass_adapter import test_rate_to_mass_adapter
from openiam.components.reservoir_data_interpolator import test_reservoir_data_interpolator
from openiam.components.seal_horizon_component import test_seal_horizon_component
from openiam.components.simple_reservoir_component import test_simple_reservoir_component
from openiam.components.stratigraphy_component import test_stratigraphy_component
from openiam.components.theis_reservoir_component import test_theis_reservoir_component


TEST_FUNCTIONS = {
    'AAC': [test_alluvium_aquifer_component, None],
    'AALFC': [test_alluvium_aquifer_lf_component, None],
    'ARC': [test_analytical_reservoir_component, None],
    'ATMRC': [test_atmRom_component, None],
    'CAC': [test_carbonate_aquifer_component, None],
    'CWC': [test_cemented_wellbore_component, None],
    'CWWRC': [test_cemented_wellbore_wr_component, None],
    'CWS': [test_chemical_well_sealing, None],
    'CC': [test_configurer_component, [1, 2, 3]], # on a TODO list
    'DAAC': [test_deep_alluvium_aquifer_component, None],
    'DAAMLC': [test_deep_alluvium_aquifer_ml_component, None],
    'FFC': [test_fault_flow_component, [1, 2]],
    'FLC': [test_fault_leakage_component, None],
    'FGAQ': [test_futuregen2_aquifer_component, None],
    'FGAZ': [test_futuregen2_azmi_component, None],
    'GFRC': [test_generalized_flow_rate_component, None],
    'GAC': [test_generic_aquifer_component, None],
    'GRC': [test_generic_reservoir_component, None],
    'HLC': [test_hydrocarbon_leakage_component, None],
    'KWC': [test_kimberlina_wellbore_component, None],
    'LG': [test_location_generator, None],
    'LUTRC': [test_lookup_table_reservoir_component, [1, 2]],
    'MSC': [test_monitoring_scheduler_component, [1, 2, 3]], # on a TODO list
    'MTC': [test_monitoring_tool_component, [1, 2, 3]],
    'MSWC': [test_multisegmented_wellbore_component, None],
    'OWC': [test_open_wellbore_component, [1, 2]],
    'PASC': [test_parameter_setup_component, [1, 2, 3, 4, 5]], # on a TODO list
    'PSC': [test_plume_stability_component, [1, 2, 3]],
    'RTMA': [test_rate_to_mass_adapter, None],
    'RDI': [test_reservoir_data_interpolator, [1, 2]],
    'SHC': [test_seal_horizon_component, [1, 2]],
    'SRC': [test_simple_reservoir_component, None],
    'STRATA': [test_stratigraphy_component, None],
    'TRC': [test_theis_reservoir_component, None]
    }


def run_test(test_case='AAC', spec_val=None):
    # Get test function and argument list
    test_function = TEST_FUNCTIONS[test_case][0]
    arg_list = TEST_FUNCTIONS[test_case][1]

    # Check whether test function was provided with arguments
    if isinstance(arg_list, list):
        if spec_val in arg_list:
            arg_list = [spec_val]
        for val in arg_list:
            print('Running {} with input {}'.format(test_function.__name__, val))
            print(80*'-')
            test_function(val)
    else:
        print('Running {}...'.format(test_function.__name__))
        print(80*'-')
        test_function()


if __name__ == '__main__':

    run_test('CWC')
