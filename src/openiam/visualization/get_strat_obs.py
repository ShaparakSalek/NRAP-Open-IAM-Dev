import sys
import os
import logging
import numpy as np

SOURCE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(SOURCE_DIR)

try:
    from openiam import (IAM_DIR, SystemModel, LookupTableStratigraphy,
                         StratigraphyDataInterpolator, DippingStratigraphy)
except ImportError as err:
    print('Unable to load IAM class module: {}'.format(err))

REFERENCE_LOC_WARNING = ''.join([
    'The input given for the {} coordinate of the reference location for the ',
    'DippingStratigraphy component was not of type int or float. The input given ',
    'was {}. Check your input. The default reference {} coordinate of 0 m will ',
    'be used.'])


def get_strat_comp_obs(x_loc, y_loc, numShaleLayers, yaml_data):
    """
    This function takes x and y grids and runs a LookupTableStratigraphy component
    for each grid node. The unit depth observations from the LookupTableStratigraphy
    component are then stored in an array that is returned by the function. This
    function is used by stratigraphy_plot.py, stratigraphic_column.py, and
    ttfd.py.
    """
    try:
        x_dims = x_loc.shape
    except:
        x_loc = np.ones((1, 1)) * x_loc
        x_dims = x_loc.shape

    try:
        y_dims = y_loc.shape
    except:
        y_loc = np.ones((1, 1)) * y_loc
        y_dims = y_loc.shape

    if x_dims != y_dims:
        warning_msg = ''.join([
            'There was an error while using the function get_strat_comp_obs() ',
            'in the file get_strat_obs.py (in the folder source/openiam/visualize). ',
            'This function is used to obtain output from components like LookupTableStratigraphy ',
            'and DippingStratigraphy components for different plot types ',
            '(Stratigraphy and TTDF plots). While using the function, the dimensions of the ',
            'x array {} and y array {} provided do not match. '.format(x_dims, y_dims),
            'The function assumes that the x and y arrays provided have the ',
            'same dimensions, so the input used may results in an error.'])

        logging.error(warning_msg)

    if 'LookupTableStratigraphy' in yaml_data:
        comp_type = 'LookupTableStratigraphy'

        comp_data = yaml_data['LookupTableStratigraphy']

        # These inputs are required for a LookupTableStratigraphy component
        file_directory = comp_data['FileDirectory']
        file_directory = os.path.join(IAM_DIR, file_directory)

        data_file = comp_data['FileName']

    elif 'DippingStratigraphy' in yaml_data:
        comp_type = 'DippingStratigraphy'

        comp_data = yaml_data['DippingStratigraphy']

        dipDirection = None
        if 'Controls' in comp_data:
            if 'dipDirection' in comp_data['Controls']:
                dipDirection = comp_data['Controls']['dipDirection']

        locXRef = 0
        locYRef = 0
        locZRef = 0
        if 'ReferenceLocation' in comp_data:
            if 'coordx' in comp_data['ReferenceLocation']:
                locXRef = comp_data['ReferenceLocation']['coordx']

                if isinstance(locXRef, list):
                    if len(locXRef) == 1:
                        locXRef = locXRef[0]
                    else:
                        logging.warning(REFERENCE_LOC_WARNING.format('x', locXRef, 'x'))
                        locXRef = 0

                if isinstance(locXRef, (int, float)):
                    locXRef = comp_data['ReferenceLocation']['coordx']
                else:
                    logging.warning(REFERENCE_LOC_WARNING.format('x', locXRef, 'x'))
                    locXRef = 0

            if 'coordy' in comp_data['ReferenceLocation']:
                locYRef = comp_data['ReferenceLocation'][
                    'coordy']

                if isinstance(locYRef, list):
                    if len(locYRef) == 1:
                        locYRef = locYRef[0]
                    else:
                        logging.warning(REFERENCE_LOC_WARNING.format('y', locYRef, 'y'))
                        locYRef = 0

                if isinstance(locYRef, (int, float)):
                    locYRef = comp_data['ReferenceLocation']['coordy']
                else:
                    logging.warning(REFERENCE_LOC_WARNING.format('y', locYRef, 'y'))
                    locYRef = 0

            if 'coordz' in comp_data['ReferenceLocation']:
                locZRef = comp_data['ReferenceLocation']['coordz']

                if isinstance(locZRef, list):
                    if len(locZRef) == 1:
                        locZRef = locZRef[0]
                    else:
                        logging.warning(REFERENCE_LOC_WARNING.format('z', locZRef, 'z'))
                        locZRef = 0

                if isinstance(locZRef, (int, float)):
                    locZRef = comp_data['ReferenceLocation']['coordz']
                else:
                    logging.warning(REFERENCE_LOC_WARNING.format('z', locZRef, 'z'))
                    locZRef = 0

    if 'Parameters' in comp_data:
        pars = comp_data['Parameters']

        for key in pars.keys():
            val = pars[key]

            if not isinstance(val, dict):
                pars[key] = {}

                pars[key]['value'] = val

            if not ('min' in pars[key] and 'max' in pars[key]) or ('vary' in pars[key]):
                pars[key]['vary'] = False
    else:
        pars = {}

    # Run a fake simulation just to interpolate the stratigraphy information
    # Define keyword arguments of the system model
    num_years = 1
    time_array = 365.25 * np.arange(0.0, num_years + 1)
    sm_model_kwargs = {'time_array': time_array} # time is given in days

    # Create dummy system model
    dummy_sm = SystemModel(model_kwargs=sm_model_kwargs)

    if comp_type == 'LookupTableStratigraphy':
        _ = dummy_sm.add_interpolator(StratigraphyDataInterpolator(
            name='int1', parent=dummy_sm,
            header_file_dir=file_directory,
            data_file=data_file),
            intr_family='stratigraphy')

    strat_comps = []
    for x_ref in range(x_dims[0]):
        for y_ref in range(x_dims[1]):
            x_loc_temp = float(x_loc[x_ref, y_ref])
            y_loc_temp = float(y_loc[x_ref, y_ref])

            try:
                if comp_type == 'LookupTableStratigraphy':
                    strat_comps.append(dummy_sm.add_component_model_object(
                        LookupTableStratigraphy(name='strata_x{}_y{}'.format(x_ref, y_ref),
                                                parent=dummy_sm, intr_family='stratigraphy',
                                                locX=x_loc_temp, locY=y_loc_temp)))

                elif comp_type == 'DippingStratigraphy':
                    strat_comps.append(dummy_sm.add_component_model_object(
                        DippingStratigraphy(
                            name='strata_x{}_y{}'.format(x_ref, y_ref),
                            parent=dummy_sm, locXRef=locXRef, locYRef=locYRef, locZRef=locZRef,
                            dipDirection=dipDirection, locX=x_loc_temp, locY=y_loc_temp)))

                strat_comps[-1].add_par('numberOfShaleLayers', value=numShaleLayers,
                                        vary=False)

                thickness_obs = strat_comps[-1].get_thickness_obs_names()

                for key in pars.keys():
                    strat_comps[-1].add_par(key, **pars[key])

                for ob_nm in thickness_obs:
                    strat_comps[-1].add_obs(ob_nm, index=[0])

                depth_obs = strat_comps[-1].get_depth_obs_names()

                for ob_nm in depth_obs:
                    strat_comps[-1].add_obs(ob_nm, index=[0])
            except:
                pass

    # Run the dummy system model just to get the stratigraphy information
    try:
        dummy_sm.forward()
        successful_run = True
    except:
        # If the simulation fails, return a value of None. The simulation may
        # fail if the x and y locations are outside the domain contained in the
        # file used for a LookupTableStratigraphy component.
        successful_run = False

    if not successful_run:
        stratigraphy_by_loc = None
    else:
        # Dictionary contaiming unit thicknesses and depths at each location
        stratigraphy_by_loc = dict()

        stratigraphy_by_loc['reservoirThickness'] = np.zeros(
            (x_loc.shape[0], x_loc.shape[1])) * np.nan
        stratigraphy_by_loc['reservoirTopDepth'] = np.zeros(
            (x_loc.shape[0], x_loc.shape[1])) * np.nan
        stratigraphy_by_loc['reservoirMidDepth'] = np.zeros(
            (x_loc.shape[0], x_loc.shape[1])) * np.nan
        stratigraphy_by_loc['reservoirDepth'] = np.zeros(
            (x_loc.shape[0], x_loc.shape[1])) * np.nan

        for shaleRef in range(numShaleLayers - 1, -1, -1):
            nm_v1 = 'shale{}'.format(shaleRef + 1)
            stratigraphy_by_loc[nm_v1 + 'Thickness'] = np.zeros(
                (x_loc.shape[0], x_loc.shape[1])) * np.nan
            stratigraphy_by_loc[nm_v1 + 'Depth'] = np.zeros(
                (x_loc.shape[0], x_loc.shape[1])) * np.nan
            stratigraphy_by_loc[nm_v1 + 'MidDepth'] = np.zeros(
                (x_loc.shape[0], x_loc.shape[1])) * np.nan
            stratigraphy_by_loc[nm_v1 + 'TopDepth'] = np.zeros(
                (x_loc.shape[0], x_loc.shape[1])) * np.nan

            if (shaleRef + 1) < numShaleLayers:
                nm_v2 = 'aquifer{}'.format(shaleRef + 1)
                stratigraphy_by_loc[nm_v2 + 'Thickness'] = np.zeros(
                    (x_loc.shape[0], x_loc.shape[1])) * np.nan
                stratigraphy_by_loc[nm_v2 + 'Depth'] = np.zeros(
                    (x_loc.shape[0], x_loc.shape[1])) * np.nan
                stratigraphy_by_loc[nm_v2 + 'MidDepth'] = np.zeros(
                    (x_loc.shape[0], x_loc.shape[1])) * np.nan
                stratigraphy_by_loc[nm_v2 + 'TopDepth'] = np.zeros(
                    (x_loc.shape[0], x_loc.shape[1])) * np.nan

        for x_ref in range(x_loc.shape[0]):
            for y_ref in range(x_loc.shape[1]):
                try:
                    comp = dummy_sm.component_models['strata_x{}_y{}'.format(x_ref, y_ref)]

                    thickness_temp = dummy_sm.collect_observations_as_time_series(
                        comp, 'reservoirThickness', indices=[0])

                    stratigraphy_by_loc['reservoirThickness'][x_ref, y_ref] = thickness_temp[0]

                    top_depth_temp = dummy_sm.collect_observations_as_time_series(
                        comp, 'reservoirTopDepth', indices=[0])

                    stratigraphy_by_loc['reservoirTopDepth'][x_ref, y_ref] = top_depth_temp[0]

                    mid_depth_temp = dummy_sm.collect_observations_as_time_series(
                        comp, 'reservoirMidDepth', indices=[0])

                    stratigraphy_by_loc['reservoirMidDepth'][x_ref, y_ref] = mid_depth_temp[0]

                    depth_temp = dummy_sm.collect_observations_as_time_series(
                        comp, 'reservoirDepth', indices=[0])

                    stratigraphy_by_loc['reservoirDepth'][x_ref, y_ref] = depth_temp[0]

                    for shaleRef in range(numShaleLayers - 1, -1, -1):
                        nm_v1 = 'shale{}'.format(shaleRef + 1)

                        thickness_temp = dummy_sm.collect_observations_as_time_series(
                            comp, nm_v1+'Thickness', indices=[0])

                        depth_temp = dummy_sm.collect_observations_as_time_series(
                            comp, nm_v1+'Depth', indices=[0])

                        mid_depth_temp = dummy_sm.collect_observations_as_time_series(
                            comp, nm_v1+'MidDepth', indices=[0])

                        top_depth_temp = dummy_sm.collect_observations_as_time_series(
                            comp, nm_v1+'TopDepth', indices=[0])

                        stratigraphy_by_loc[nm_v1+'Thickness'][x_ref, y_ref] = thickness_temp[0]

                        stratigraphy_by_loc[nm_v1+'Depth'][x_ref, y_ref] = depth_temp[0]

                        stratigraphy_by_loc[nm_v1+'MidDepth'][x_ref, y_ref] = mid_depth_temp[0]

                        stratigraphy_by_loc[nm_v1+'TopDepth'][x_ref, y_ref] = top_depth_temp[0]

                        if (shaleRef + 1) < numShaleLayers:
                            nm_v2 = 'aquifer{}'.format(shaleRef + 1)

                            thickness_temp = dummy_sm.collect_observations_as_time_series(
                                comp, nm_v2+'Thickness', indices=[0])

                            depth_temp = dummy_sm.collect_observations_as_time_series(
                                comp, nm_v2+'Depth', indices=[0])

                            mid_depth_temp = dummy_sm.collect_observations_as_time_series(
                                comp, nm_v2+'MidDepth', indices=[0])

                            top_depth_temp = dummy_sm.collect_observations_as_time_series(
                                comp, nm_v2+'TopDepth', indices=[0])

                            stratigraphy_by_loc[nm_v2+'Thickness'][x_ref, y_ref] = thickness_temp[0]

                            stratigraphy_by_loc[nm_v2+'Depth'][x_ref, y_ref] = depth_temp[0]

                            stratigraphy_by_loc[nm_v2+'MidDepth'][x_ref, y_ref] = mid_depth_temp[0]

                            stratigraphy_by_loc[nm_v2+'TopDepth'][x_ref, y_ref] = top_depth_temp[0]
                except:
                    for shaleRef in range(numShaleLayers - 1, -1, -1):
                        nm_v1 = 'shale{}'.format(shaleRef + 1)

                        stratigraphy_by_loc[nm_v1+'Thickness'][x_ref, y_ref] = None

                        stratigraphy_by_loc[nm_v1+'Depth'][x_ref, y_ref] = None

                        stratigraphy_by_loc[nm_v1+'MidDepth'][x_ref, y_ref] = None

                        stratigraphy_by_loc[nm_v1+'TopDepth'][x_ref, y_ref] = None

                        if (shaleRef + 1) < numShaleLayers:
                            stratigraphy_by_loc[nm_v2+'Thickness'][x_ref, y_ref] = None

                            stratigraphy_by_loc[nm_v2+'Depth'][x_ref, y_ref] = None

                            stratigraphy_by_loc[nm_v2+'MidDepth'][x_ref, y_ref] = None

                            stratigraphy_by_loc[nm_v2+'TopDepth'][x_ref, y_ref] = None

    return stratigraphy_by_loc
