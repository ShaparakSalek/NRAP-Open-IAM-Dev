# -*- coding: utf-8 -*-
import os
import sys
import logging
import numpy as np

try:
    import openiam.components.iam_base_classes as iam_bc
except ImportError as err:
    print('Unable to load NRAP-Open-IAM base classes module: {}'.format(err))


class Stratigraphy(iam_bc.ComponentModel):
    """
    The Stratigraphy component is a component containing parameters describing
    the structure of the storage reservoir system. The stratigraphy component
    allows the number of shale (or aquitard) layers to be set, thus, setting the
    total number of layers in the system. Each shale or aquifer layer can
    take on the default thickness for that layer type or be assigned a user defined value
    with ``shale#Thickness`` or ``aquifer#Thickness`` keywords where # is replaced
    by an index of the layer (e.g., ``shale1Thickness``). Layers are numbered from
    the bottom upward: shale 1 is the layer above the storage reservoir,
    and, with N shale layers total, shale N is the surface layer.

    .. image:: ../../images/ShaleLayers.png
       :align: center
       :scale: 100%
       :alt: Layer ordering

    This component represents flat-lying stratigraphy; when using this component, 
    unit thicknesses and depths do not vary across the domain.
    
    Descriptions of the component's parameters are provided below.

    * **numberOfShaleLayers** [-] (3 to 30) - number of shale layers in the
      system (default: 3). The shale units must be separated by an aquifer.

    * **shaleThickness** [|m|] (1 to 1600) - thickness of shale layers (default
      250). Thickness of shale layer 1, for example, can be defined by
      **shale1Thickness**; otherwise, shale layers for which the thickness is not
      defined will be assigned a default thickness.

    * **aquiferThickness** [|m|] (1 to 1600) - thickness of aquifers (default: 100).
      Thickness of aquifer 1, for example, can be defined by **aquifer1Thickness**;
      otherwise, aquifers for which the thickness is not defined will be assigned
      a default thickness.

    * **reservoirThickness** [|m|] (1 to 1600) - thickness of reservoir (default: 50)

    * **datumPressure** [|Pa|] (80,000 to 300,000) - pressure at the top of the
      system (default: 101,325)

    * **depth** [|m|] (5 to 30,000) - depth to the top of reservoir (default: 950).
      The depth to the top of the reservoir can also be obtained through the 
      **reservoirTopDepth** or **shale1Depth** composite parameters (the bottom 
      of shale 1 is also the top of the reservoir).

    The depths to the bottom, middle, and top of each unit are produced as composite
    parameters (i.e., calculated as functions of the thickness parameters used):

    * **shale#Depth** [|m|] (boundaries depend on the user input) - depth
      to the base of the shale layer with index # (default value is not defined)

    * **aquifer#Depth** [|m|] (boundaries depend on the user input) - depth
      to the base of the aquifer layer with index # (default value is not defined)

    * **reservoirDepth** [|m|] (boundaries depend on the user input) - depth
      to the base of the reservoir (default value is not defined).
    
    * **shale#MidDepth** [|m|] (boundaries depend on the user input) - depth
      to the middle of the shale layer with index # (default value is not defined)

    * **aquifer#MidDepth** [|m|] (boundaries depend on the user input) - depth
      to the middle of the aquifer layer with index # (default value is not defined)

    * **reservoirMidDepth** [|m|] (boundaries depend on the user input) - depth
      to the middle of the reservoir (default value is not defined).
    
    * **shale#TopDepth** [|m|] (boundaries depend on the user input) - depth
      to the top of the shale layer with index # (default value is not defined)

    * **aquifer#TopDepth** [|m|] (boundaries depend on the user input) - depth
      to the top of the aquifer layer with index # (default value is not defined)

    * **reservoirTopDepth** [|m|] (boundaries depend on the user input) - depth
      to the top of the reservoir (default value is not defined).

    Although these composite parameters are automatically added in the control 
    file interface and graphical user interface, they must be added explicitly 
    as composite parameters in script applications. For examples of such 
    applications, see the script examples ``iam_sys_strata.py``, 
    ``iam_sys_strata_reservoir_openwell_genericaquifer.py``, or 
    ``iam_sys_strata_reservoir_openwell_genericaquifer_5locs.py``.
    """
    def __init__(self, name, parent):
        """
        Constructor method of Container Class.
        """
        super().__init__(name, parent, model='')

        # Add type attribute
        self.class_type = 'Stratigraphy'

        # Set default parameters of the component model
        self.add_default_par('numberOfShaleLayers', value=3)
        self.add_default_par('shaleThickness', value=250.0)
        self.add_default_par('aquiferThickness', value=100.0)
        self.add_default_par('reservoirThickness', value=50.0)
        self.add_default_par('datumPressure', value=101325.0)
        self.add_default_par('depth', value=950.0)   # depth to the top of the reservoir

        # Define dictionary of boundaries
        self.pars_bounds = dict()
        self.pars_bounds['numberOfShaleLayers'] = [3, 30]
        self.pars_bounds['shaleThickness'] = [1.0, 1600.0]
        self.pars_bounds['aquiferThickness'] = [1.0, 1600.0]
        self.pars_bounds['reservoirThickness'] = [1.0, 1600.0]
        self.pars_bounds['datumPressure'] = [8.0e+4, 3.0e+5]
        self.pars_bounds['depth'] = [5.0, 30000.0]

        # Indicate that the component should not be run
        self.default_run_frequency = 0
        self.run_frequency = 0

        debug_msg = 'Stratigraphy component created with name {}'.format(name)
        logging.debug(debug_msg)

    def check_input_parameters(self, p):
        """
        Check whether input parameters fall within specified boundaries.

        :param p: input parameters of component model
        :type p: dict
        """
        debug_msg = 'Input parameters of component {} are {}'.format(self.name, p)
        logging.debug(debug_msg)

        for key, val in p.items():
            if not key in self.gridded_pars:
                warn_msg = ''.join([
                    'Parameter {} of Stratigraphy component {} ',
                    'is out of boundaries.']).format(key, self.name)
                if (key[0:5] == 'shale' and key[-9:] == 'Thickness'):
                    if (val < self.pars_bounds['shaleThickness'][0]) or (
                            val > self.pars_bounds['shaleThickness'][1]):
                        logging.warning(warn_msg)
                elif key[0:7] == 'aquifer' and key[-9:] == 'Thickness':
                    if (val < self.pars_bounds['aquiferThickness'][0]) or (
                            val > self.pars_bounds['aquiferThickness'][1]):
                        logging.warning(warn_msg)
                elif key[0:5] == 'depth':
                    if (val < self.pars_bounds['depth'][0]) or (
                            val > self.pars_bounds['depth'][1]):
                        logging.warning(warn_msg)
                elif key in self.pars_bounds:
                    if (val < self.pars_bounds[key][0]) or (
                            val > self.pars_bounds[key][1]):
                        logging.warning(warn_msg)
                else:
                    warn_msg = ''.join([
                        'Parameter {} is not recognized as an input parameter ',
                        'of Stratigraphy component {}.']).format(key, self.name)
                    logging.warning(warn_msg)
    
    def get_thickness_par_names(self):
        """
        Returns a list of the thickness parameter names, given the numberOfShaleLayers.
        """
        numberOfShaleLayers = self.get_num_shale_layers()
        
        par_names = ['shale{}Thickness'.format(ind) \
                     for ind in range(1, numberOfShaleLayers + 1)] + [
                         'aquifer{}Thickness'.format(ind) \
                             for ind in range(1, numberOfShaleLayers)] + [
                                     'reservoirThickness']
        
        return par_names
    
    def get_composite_depth_names(self, cfi=False):
        """
        Returns a list of the names for all composite depth parameters, given 
        the value for the numberOfShaleLayers parameter.
        """
        numberOfShaleLayers = self.get_num_shale_layers(cfi=cfi)
        
        names = ['shale{}Depth'.format(ind) \
                 for ind in range(1, numberOfShaleLayers + 1)] + [
                         'aquifer{}Depth'.format(ind) \
                             for ind in range(1, numberOfShaleLayers)] + [
                                     'reservoirDepth']
        
        names += ['shale{}MidDepth'.format(ind) \
                  for ind in range(1, numberOfShaleLayers + 1)] + [
                          'aquifer{}MidDepth'.format(ind) \
                              for ind in range(1, numberOfShaleLayers)] + [
                                      'reservoirMidDepth']
        
        names += ['shale{}TopDepth'.format(ind) \
                  for ind in range(1, numberOfShaleLayers + 1)] + [
                          'aquifer{}TopDepth'.format(ind) \
                              for ind in range(1, numberOfShaleLayers)] + [
                                      'reservoirTopDepth']
        
        return names
    
    def get_depth_expr(self, par_nm):
        """
        Returns a string for the expression used to solve for a unit's top, middle, 
        or bottom depth. This expression can be used to set a unit depth as a 
        composite parameter.
        """
        if par_nm[:5] == 'shale':
            unit_type = 'shale'
            par_nm_end = par_nm[5:]
        elif par_nm[:7] == 'aquifer':
            unit_type = 'aquifer'
            par_nm_end = par_nm[7:]
        elif par_nm[:9] == 'reservoir':
            unit_type = 'reservoir'
            par_nm_end = par_nm[9:]
            unit_num = 0
        
        if (par_nm_end[-8:] == 'TopDepth') or (par_nm_end[-8:] == 'MidDepth'):
            if unit_type != 'reservoir':
                unit_num = int(par_nm_end[:-8])
            
            if par_nm_end[-8:-5] == 'Top':
                top_mid_bottom = 'top'
            elif par_nm_end[-8:-5] == 'Mid':
                top_mid_bottom = 'mid'
                
        elif par_nm_end[-5:] == 'Depth':
            if unit_type != 'reservoir':
                unit_num = int(par_nm_end[:-5])
            
            top_mid_bottom = 'bottom'
        
        numberOfShaleLayers = self.get_num_shale_layers()
        
        par_expr = '+'.join(
            ['{}.shale{}Thickness'.format(self.name, j) for j in range(
                unit_num+1, numberOfShaleLayers+1)]+[
                    '{}.aquifer{}Thickness'.format(self.name, j) for j in range(
                        unit_num+1, numberOfShaleLayers)])
        
        if unit_type == 'shale' and unit_num != numberOfShaleLayers:
            par_expr += '+{}.aquifer{}Thickness'.format(self.name, unit_num)
        elif unit_type == 'reservoir':
            unit_num = ''
        
        if unit_type == 'shale' and unit_num == numberOfShaleLayers:
            expr_addition_start = ''
        else:
            expr_addition_start = '+'
            
        if top_mid_bottom == 'mid':
            par_expr += '{}({}.{}{}Thickness/2)'.format(
                expr_addition_start, self.name, unit_type, unit_num)
        elif top_mid_bottom == 'bottom':
            par_expr += '{}{}.{}{}Thickness'.format(
                expr_addition_start, self.name, unit_type, unit_num)
        
        if len(par_expr) == 0:
            par_expr = '0'
        
        return par_expr
    
    def get_num_shale_layers(self, cfi=False):
        """
        Returns the value of the numberOfShaleLayers parameter.
        """
        # Check if numberOfShaleLayers is among stochastic parameters
        if 'numberOfShaleLayers' in self.pars and cfi:
            err_msg = ''.join(['Parameter numberOfShaleLayers cannot be ',
                               'stochastic for the control file interface.'])
            logging.error(err_msg)
            raise TypeError(err_msg)
            
        elif 'numberOfShaleLayers' in self.pars and not cfi:
            numberOfShaleLayers = self.pars['numberOfShaleLayers'].value
            
        elif 'numberOfShaleLayers' in self.deterministic_pars:
            numberOfShaleLayers = self.deterministic_pars['numberOfShaleLayers'].value
            
        else:
            numberOfShaleLayers = self.default_pars['numberOfShaleLayers'].value
            
            warn_msg = ''.join([
                'Parameter numberOfShaleLayers is not defined in the control ', 
                'file interface. The parameter will be assigned a default value ', 
                'of {}.'.format(numberOfShaleLayers)])
            logging.warn(warn_msg)
            
            self.add_par('numberOfShaleLayers', value=numberOfShaleLayers, vary=False)
        
        return numberOfShaleLayers
    
    def connect_with_system(self):
        """
        Code to add stratigraphy to system model for control file interface.
        """
        # Check if numberOfShaleLayers is among stochastic parameters
        _ = self.get_num_shale_layers(cfi=True)

        # Check whether all other needed parameters are defined by user
        par_names = self.get_thickness_par_names()

        for ind, par_nm in enumerate(par_names):
            if (par_nm not in self.pars) and (par_nm not in self.deterministic_pars):
                if 'shale' in par_nm:
                    default_value = self.default_pars['shaleThickness'].value
                elif 'aquifer' in par_nm:
                    default_value = self.default_pars['aquiferThickness'].value
                elif 'reservoir' in par_nm:
                    default_value = self.default_pars['reservoirThickness'].value
                warn_msg = ''.join([
                    'Parameter {} is not defined in the control file ',
                    'interface. The parameter will be assigned ',
                    'a default value of {}.']).format(par_nm, default_value)
                logging.warn(warn_msg)
                self.add_par(par_nm, value=default_value, vary=False)

        # Add composite parameters of the stratigraphy component representing 
        # the depths to the bottom, middle, and top of each shale and aquifer 
        # layer as well as the reservoir.
        composite_depth_names = self.get_composite_depth_names(cfi=True)
        
        for comp_depth in composite_depth_names:
            par_expr = self.get_depth_expr(comp_depth)
            self.add_composite_par(comp_depth, par_expr)


def test_stratigraphy_component():
    try:
        from openiam.components.analytical_reservoir_component import AnalyticalReservoir
    except ImportError as err:
        print('Unable to load IAM class module: '+str(err))

    logging.basicConfig(level=logging.WARNING)
    # Define keyword arguments of the system model
    num_years = 50
    time_array = 365.25*np.arange(0.0, num_years+1) # time is in days
    sm_model_kwargs = {'time_array': time_array}
    sm = iam_bc.SystemModel(model_kwargs=sm_model_kwargs)

    strata = sm.add_component_model_object(Stratigraphy(name='strata', parent=sm))

    # Add parameters of stratigraphy component model. Parameters that are not 
    # entered will resort to default values.
    strata.add_par('numberOfShaleLayers', value=3, vary=False)
    strata.add_par('shale1Thickness', min=30.0, max=50., value=40.0)
    strata.add_par('shale2Thickness', min=40.0, max=60., value=50.0)

    # Add reservoir component
    ares = sm.add_component_model_object(AnalyticalReservoir(name='ares', parent=sm))

    # Add parameters of reservoir component model
    ares.add_par('injRate', min=0.4, max=0.6, value=0.5)
    ares.add_par('reservoirRadius', value=10000, vary=False)
    ares.add_par_linked_to_par('numberOfShaleLayers',
                               strata.deterministic_pars['numberOfShaleLayers'])
    ares.add_par_linked_to_par('shale1Thickness', strata.pars['shale1Thickness'])
    ares.add_par_linked_to_par('shale2Thickness', strata.pars['shale2Thickness'])
    ares.add_par_linked_to_par('shale3Thickness',
                               strata.default_pars['shaleThickness'])

    ares.add_par_linked_to_par('aquifer1Thickness',
                               strata.default_pars['aquiferThickness'])
    ares.add_par_linked_to_par('aquifer2Thickness',
                               strata.default_pars['aquiferThickness'])
    ares.add_par_linked_to_par('aquifer3Thickness',
                               strata.default_pars['aquiferThickness'])

    ares.add_par_linked_to_par('reservoirThickness',
                               strata.default_pars['reservoirThickness'])

    ares.add_par_linked_to_par('datumPressure',
                               strata.default_pars['datumPressure'])

    ares.add_obs('pressure')

    # Run system model using current values of its parameters
    sm.forward()

    print('------------------------------------------------------------------')
    print('                  Forward method illustration ')
    print('------------------------------------------------------------------')

    # Print pressure
    print('Pressure', sm.collect_observations_as_time_series(ares, 'pressure'),
          sep='\n')
