# -*- coding: utf-8 -*-
import os
import sys
import logging
import numpy as np
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from openiam import SystemModel, ComponentModel
except ImportError as err:
    print('Unable to load IAM class module: {}'.format(err))

DIP_DIRECTION_DEGREE_OPTIONS = [0, 90, 180, 270, 360]
DIP_DIRECTION_OPTIONS = ['N', 'S', 'E', 'W', 'NE', 'SE', 'SW', 'NW']
DEFAULT_DIP_DIRECTION = 'NE'

CHECK_CONDITIONS_MSG = ''.join([
    'Check your input. Alternatively, calculating unit thicknesses with these ',
    'strike and dip values may not be appropriate for these locations and/or ',
    'conditions.'])

REFERENCE_LOC_WARNING = ''.join([
    'The input given for the {} coordinate of the reference location for the ', 
    'DippingStratigraphy component was not of type int or float. The input given ', 
    'was {}. Check your input. The reference {} coordinate will instead be set to ', 
    '{} m.'])

DIP_DIRECTION_WARNING = ''.join([
    'The dipDirection ({}) provided to the  DippingStratigraphy component was ', 
    'not one of the available options ({}). The dipDirection will be set to ', 
    '90 degrees clockwise from the direction of strike(when viewed in a ', 
    'map-view figure).'])

class DippingStratigraphy(ComponentModel):
    """
    The Dipping Stratigraphy component is similar to to the Stratigraphy component, 
    but it portrays dipping units. This component takes strike and dip values 
    as well as unit thicknesses for a reference location. It then uses that input 
    to produce unit thicknesses and depths at for the given output location(s). 
    This component 

    The arrangement of geologic units for this component is the same as the that 
    for the ``Stratigraphy`` component. The reservoir is the lowest unit. Alternating 
    shale and aquifer layers overlie the reservoir. These units each have an index, 
    with the index number increasing closer to the surface. Shale 1 is above the 
    reservoir, aquifer 1 is above shale 1, and shale 2 is above aquifer 1. This 
    pattern continues until shale ``N``, where ``N`` is the number of shale units 
    used (the *numberOfShaleLayers* parameter). There are ``N - 1`` aquifers. The 
    *numberOfShaleLayers* parameter cannot vary across the domain.
    
    The ``dipDirection`` keyword can be used to specify the direction of dip. The 
    ``dipDirection`` can be given as a cardinal direction (``N``, ``E``, ``S``, 
    ``W``, ``NE``, ``SE``, ``SW``, and ``NW``). In the control file interface, 
    the ``dipDirection`` can be specified under the ``Controls`` entry, which 
    is entered directly under ``DippingStratigraphy``. For example, the *strike* 
    and *dip* parameters could be given as 45 |deg| and 3 |deg|, respectively. 
    Because the strike is to the northwest, the ``dipDirection`` could be to the 
    southeast or northwest (``SE`` or ``NW``). In this example, providing a 
    ``dipDirection`` to the south or north (``S`` or ``N``) will also be taken 
    as a dip to the southeast or northeast, respectively (i.e., a southwards 
    component of dip vs. a northwards component). If the ``dipDirection`` keyword 
    is not specified, then the ``dipDirection`` will be taken as 90 |deg| clockwise 
    from the strike direction, if viewed in a map-view image (i.e., the 
    right-hand rule in geology; when facing in the direction of strike, 
    the dip direction is to your right).
    
    In the NRAP-Open-IAM control file interface, the type name for the Dipping 
    Stratigraphy component is ``DippingStratigraphy``.
    
    Descriptions of the component's parameters are provided below.
    
    * **strike** [-] (0 to 360) - strike of the units in degrees (default: 315). 
      In a map-view image, this parameter is taken as increasing from zero with 
      rotation clockwise from north. For example, *strike* values of 0, 90, 180, 
      and 270 would represent units striking to the north, east, south, and west, 
      repsectively.

    * **dip** [-] (0.1 to 89.9) - dip of the units in degrees (default: 5). 
      If the dip causes any unit thicknesses in the domain to fall outside the 
      range of 1 |m| to 1600 |m|, that unit thickness will be set to the closest 
      limit.
    
    * **numberOfShaleLayers** [-] (3 to 30) - number of shale layers in the
      system (default: 3). The shale units must be separated by an aquifer.

    * **shaleThickness** [|m|] (1 to 1600) - thickness of shale layers at the 
      reference location (default: 250). Thickness of shale layer 1, for example, 
      can be defined by **shale1Thickness**; otherwise, shale layers for which 
      the thickness is not defined will be assigned a default thickness.

    * **aquiferThickness** [|m|] (1 to 1600) - thickness of aquifers at the reference 
      location (default: 100). Thickness of aquifer 1, for example, can be defined 
      by **aquifer1Thickness**; otherwise, aquifers for which the thickness is 
      not defined will be assigned a default thickness.

    * **reservoirThickness** [|m|] (1 to 1600) - thickness of the reservoir at 
      the reference location (default: 50)

    * **datumPressure** [|Pa|] (80,000 to 300,000) - pressure at the top of the
      system (default: 101,325)

    Unit thicknesses and the depths to the bottom, middle, and top of each unit 
    at specific locations (``x`` and ``y`` values) are produced as observations:
    
    The observations from the Lookup Table Stratigraphy component are:
    
    * **shale#Thickness** [|m|] - thickness of shale ``#`` at the output location(s), 
      where ``#`` is an index ranging from one to *numberOfShaleLayers* .

    * **aquifer#Thickness** [|m|] - thickness of aquifer ``#`` at the output 
      location(s), where ``#`` is an index ranging from one to 
      (*numberOfShaleLayers* - 1).

    * **reservoirThickness** [|m|] - reservoir thickness at the output location(s).
    
    * **shale#Depth** [|m|] - depth to the bottom of the shale unit with index 
      ``#`` at the output location(s).

    * **aquifer#Depth** [|m|] - depth to the bottom of aquifer unit with index 
      ``#`` at the output location(s).

    * **reservoirDepth** [|m|] - depth to the base of the reservoir at the 
      output location(s).
    
    * **shale#MidDepth** [|m|] - depth to the middle of the shale layer with 
      index ``#`` at the output location(s).

    * **aquifer#MidDepth** [|m|] - depth to the middle of the aquifer layer with 
      index ``#`` at the output location(s).

    * **reservoirMidDepth** [|m|] - depth to the middle of the reservoir at the 
      output location(s).
    
    * **shale#TopDepth** [|m|] - depth to the top of the shale unit with index 
      ``#`` at the output location(s).

    * **aquifer#TopDepth** [|m|] - depth to the top of the aquifer unit with 
      index ``#`` at the output location(s).

    * **reservoirTopDepth** [|m|] - depth to the top of the reservoir at the 
      output location(s).

    If the component produces a unit thickness for a shale, aquifer, or reservoir 
    that falls outside the range of 1 |m| to 1600 |m|, then that thickness will 
    be set to the lower or upper limit (whichever is closer). For example, if the 
    unit thickness would be 1700 |m|, based on the parameters and output locations 
    used, then that thickness observation will instead be set to 1600 |m|.
    
    For control file examples using the ``DippingStratigraphy`` component, see 
    *ControlFile_ex33a.yaml* to *ControlFile_ex35.yaml*, and *ControlFile_ex39b.yaml*. 
    For script examples, see *iam_sys_dippingstrata.py*, *iam_sys_dippingstrata_gridded.py*, 
    *iam_sys_reservoir_mswell_stratplot_dipping_strata.py*, and 
    *iam_sys_reservoir_mswell_futuregen_ttfdplot_dipping_strata.py*
    """
    def __init__(self, name, parent, locXRef=0, locYRef=0, locZRef=0, 
                 dipDirection=None, locX=None, locY=None):
        """
        Constructor method of DippingStratigraphy Class.
        
        :param name: name of component model
        :type name: str

        :param parent: the SystemModel object that the component model
            belongs to
        :type parent: SystemModel object
        
        :param locXRef: x-coordinate of the reference location (default: 0 |m|). 
            The unit thicknesses at the reference location (assigned as parameters) 
            are used to calculate unit thicknesses and depths at different 
            locations.
        :type locXRef: float

        :param locYRef: y-coordinate of the reference location (default: 0 |m|). 
            The unit thicknesses at the reference location (assigned as parameters) 
            are used to calculate unit thicknesses and depths at different 
            locations.
        :type locYRef: float
        
        :param locX: x-coordinate of the location at which stratigraphy information 
            is to to be calculated.
            By default, it is None, which means that the whole grid data is requested
        :type locX: float or array-like of floats

        :param locY: y-coordinate of the location at which stratigraphy information 
            is to to be calculated.
            By default, it is None, which means that the whole grid data is requested
        :type locY: float or array-like of floats
        """
        # Set up keyword arguments of the 'model' method provided by the system model
        model_kwargs = {'time_point': 365.25}  # default value of 365.25 days
        
        super().__init__(name, parent, model=self.simulation_model,
                         model_kwargs=model_kwargs)
        
        # The model only needs to be run once
        self.run_frequency = 1
        self.default_run_frequency = 1

        # Add type attribute
        self.class_type = 'DippingStratigraphy'
        
        if dipDirection:
            if not dipDirection in DIP_DIRECTION_OPTIONS:
                logging.warning(DIP_DIRECTION_WARNING.format(
                    dipDirection, DIP_DIRECTION_OPTIONS))
                self.dipDirection = None
                
            else:
                self.dipDirection = dipDirection
        else:
            self.dipDirection = dipDirection
        
        self.locXRef = locXRef
        self.locYRef = locYRef
        self.locZRef = locZRef
        
        # Set up location attributes
        if locX is not None and locY is not None:
            self.locX = np.asarray(locX).flatten()
            self.locY = np.asarray(locY).flatten()
            if len(self.locX) != len(self.locY):
                err_msg = 'Length of locX and locY arrays are not the same.'
                logging.error(err_msg)
                raise ValueError(err_msg)

            if len(self.locX) > 1:
                self.grid_obs_keys = []
                self.grid_obs_requested = True
            else:
                self.grid_obs_requested = False
            
            self.num_points = len(self.locX)
        else:
            self.grid_obs_keys = []
            self.locX = np.array([None])
            self.locY = np.array([None])
            # Add flag indicating whether whole grid data is needed
            self.grid_obs_requested = True
            self.num_points = None

        # Set default parameters of the component model
        self.add_default_par('numberOfShaleLayers', value=3)
        self.add_default_par('shaleThickness', value=250.0)
        self.add_default_par('aquiferThickness', value=100.0)
        self.add_default_par('reservoirThickness', value=50.0)
        self.add_default_par('datumPressure', value=101325.0)
        self.add_default_par('strike', value=315)
        self.add_default_par('dip', value=5)

        # Define dictionary of boundaries
        self.pars_bounds = dict()
        self.pars_bounds['numberOfShaleLayers'] = [3, 30]
        self.pars_bounds['shaleThickness'] = [1.0, 1600.0]
        self.pars_bounds['aquiferThickness'] = [1.0, 1600.0]
        self.pars_bounds['reservoirThickness'] = [1.0, 1600.0]
        self.pars_bounds['depth'] = [5.0, 30000.0]
        self.pars_bounds['datumPressure'] = [8.0e+4, 3.0e+5]
        self.pars_bounds['strike'] = [0, 360]
        self.pars_bounds['dip'] = [0.1, 89.9]

        # Indicate that the component should not be run
        self.default_run_frequency = 1
        self.run_frequency = 1

        debug_msg = 'DippingStratigraphy component created with name {}'.format(name)
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
    
    def parameter_assignment_warning_msg(self):
        """
        Generates a warning message. This warning message is meant to be used 
        in the connect_with_system() methods of components in the control file 
        interface. Specifically, it is used when the parameter of the component 
        is assigned as a stochastic or deterministic parameter when it should 
        instead be linked to an observation from the DippingStratigraphy 
        component.
        """
        warning_msg = ''.join([
            'When using a DippingStratigraphy component, the {} parameter ', 
            'cannot be specified as a stochastic or deterministic ', 
            'parameter. The parameter needs to be linked to the ', 
            '{} observation of the DippingStratigraphy component. ', 
            'The input provided for that parameter will not be used. ', 
            'In the control file interface, the parameter linkage is performed ', 
            'automatically. In a script application, the parameter can be linked ', 
            'with the add_par_linked_to_obs() method of the component with the ', 
            'parameter being linked to a DippingStratigraphy observation. ', 
            'Before linking the DippingStratigraphy observation to the ', 
            'parameter, use the add_obs() and add_obs_to_be_linked() methods ', 
            'of the DippingStratigraphy component. The add_obs() method ', 
            'can include the obsertvation name and index=[0]. For example, ', 
            'DipStratComponentName.add_obs(obs_name, index=[0]).'])
        
        return warning_msg
    
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
    
    def get_unit_thicknesses_ref_point(self):
        """
        Returns lists of shale and aqufier thicknesses at the reference location 
        as well as the reservoir thickness. The shale and aquifer thickness lists 
        are arranged with the lower units at the beginning of the list. For example, 
        shale1Thickness is stored in index 0 of the list, while shale3Thickness 
        is stored in index 2.
        """
        numberOfShaleLayers = self.get_num_shale_layers()
        
        shaleThicknessList = []
        aquiferThicknessList = []
        
        warn_msg = ''.join([
            'The DippingStratigraphy parameter {} is not defined. The parameter ', 
            'will be assigned a default value of {} m. Note that this ', 
            'parameter is representative of the reference location for the ', 
            'DippingStratigraphy component. The thicknesses and depths at the ', 
            'output location are produced as observations by the ', 
            'DippingStratigraphy component.'])
        
        for unitRef in range(numberOfShaleLayers):
            unit_name = 'shale{}Thickness'.format(unitRef + 1)
            
            if unit_name in self.pars:
                shaleThicknessList.append(self.pars[unit_name].value)
                
            elif unit_name in self.deterministic_pars:
                shaleThicknessList.append(self.deterministic_pars[unit_name].value)
                
            else:
                logging.warn(warn_msg.format(
                    unit_name, self.default_pars['shaleThickness'].value))
                shaleThicknessList.append(self.default_pars['shaleThickness'].value)
            
            if (unitRef + 1) < numberOfShaleLayers:
                unit_name = 'aquifer{}Thickness'.format(unitRef + 1)
                
                if unit_name in self.pars:
                    aquiferThicknessList.append(self.pars[unit_name].value)
                    
                elif unit_name in self.deterministic_pars:
                    aquiferThicknessList.append(self.deterministic_pars[unit_name].value)
                    
                else:
                    logging.warn(warn_msg.format(
                        unit_name, self.default_pars['aquiferThickness'].value))
                    aquiferThicknessList.append(self.default_pars['aquiferThickness'].value)
        
        unit_name = 'reservoirThickness'
        
        if unit_name in self.pars:
            reservoirThickness = self.pars[unit_name].value
            
        elif unit_name in self.deterministic_pars:
            reservoirThickness = self.deterministic_pars[unit_name].value
            
        else:
            logging.warn(warn_msg.format(
                unit_name, self.default_pars[unit_name].value))
            reservoirThickness = self.default_pars[unit_name].value
        
        return shaleThicknessList, aquiferThicknessList, reservoirThickness
    
    def get_par_values(self, par_names):
        """
        Returns the values of the given parameter names.
        """
        par_values = []
        for par in par_names:
            if par in self.pars:
                par_values.append(self.pars[par].value)
                
            elif par in self.deterministic_pars:
                par_values.append(self.deterministic_pars[par].value)
                
            else:
                par_values.append(self.default_pars[par].value)
        
        return par_values
    
    def layer_thickness_bound_debug_message(self, unit, bound, locx, locy, msg_option=1):
        """
        Returns string delivering debug message regarding setting unit thickness to
        either minimum or maximum value.

        unit: string, 'shale' or 'aquifer'
        index: integer, index of the unit: e.g., 2 for aquifer 2
        bound: string, 'minimum' or 'maximum'
        locx: float, x-cordinate of location of interest
        locy: float, y-cordinate of location of interest

        """
        side = {'minimum': 'below', 'maximum': 'above'}

        if msg_option == 1:
            msg = ''.join([
                'The thickness of {} went {} the {} value allowed at ',
                'x = {} m and y = {} m. Setting the parameter to the {} value ',
                'and accomplishing the required depth changes by altering ',
                'the thicknesses of the underlying units.']).format(
                    unit, side[bound], bound, locx, locy, bound)
        else:
            msg = ''.join([
                'The thickness of {} went {} the {} value allowed at ',
                'x = {} m and y = {} m. Setting the parameter to the {} value.']).format(
                    unit, side[bound], bound, locx, locy, bound)

        return msg
    
    def make_xyz_points_from_strike_and_dip(self, L1=100,L2=100, L3=100, L4=100,
                                            point0_xyz=None):
        """
        This function calculates a set of x, y, and z values based on a reference
        location, and dip direction. Here, x and y are horizontal distances and z
        is elevation relative to the point0_xy (negative means downwards).
        The point0_xy is assumed to be at 0 m above the surface.

        The dx and dy values used to scale each distance are calculated based on
        the dip direction (in degrees) and the input values for the corresponding
        length scale (L0, L1, L2, L3, L4). Only three points are needed to create
        a three-dimensional plane in the depth_change() function (p0, p1, and p2)
        but all five points are needed in the stratigraphy_plot() function in
        stratigraphy_plot.py.

        Point (p0) is at the location in point0_xyz. The 2nd (p1) is at a distance of 
        dx and dy from p0 in the direction of dip. The 3rd and 4th points (p2 and p3)
        are at distances dx and dy away from p0 in the directions of strike (one
        in each direction). The 5th point (p4) is at a distances of dx and dy from
        p0 in the direction of dip, and is used for the location of the number in
        the dip symbol. The dx and dy values are calculated with the length scales
        (L1, L2, L3, and L4) in a way that depends on dip direction.

        Note that L4 is doubled for certain values of dipDirectionDegrees because
        in Stratigraphy plots, those cases the dip number would often look like 
        it is right on top of the strike and dip symbol (for the default perspectives 
        set by view_elev and view_azimuth).

        :param L1: Length scale (m) used to calculate the position of point 1,
            which is described by x_p1, y_p1, and z_p1. Point 1 is in the direction
            of dip.
        :type L1: int or float

        :param L2: Length scale (m) used to calculate the position of point 2,
            which is described by x_p2, y_p2, and z_p2. Point 2 is in the direction
            of strike, 90 degrees counter clockwise from the dip direction.
        :type L2: int or float

        :param L3: Length scale (m) used to calculate the position of point 3,
            which is described by x_p3, y_p3, and z_p3. Point 3 is in the direction
            of strike, 90 degrees clockwise from the dip direction.
        :type L3: int or float

        :param L4: Length scale (m) used to calculate the position of point 4,
            which is described by x_p4, y_p4, and z_p4. Point 4 is in the direction
            of dip, but is meant to be farther down dip than point 1. This point is
            used for the placement of dip numbers in the stratigraphy_plot()
            function.
        :type L3: int or float

        :param point0_xyz: Array of length 2 containing the x, y, and z values for
            point 0 (x_p0, y_p0, and z_p0). In the stratigraphy _plot() function,
            point 0 is taken as the location for the strike and dip symbol
            (SandD_Location).
        :type point0_xy: array

        :returns: x_points, y_points, and z_points
        """
        par_values = self.get_par_values(['dip'])
        dip = par_values[0]
        
        dipDirectionDegrees = self.obtain_dip_direction_degrees()
        
        if point0_xyz is None:
            x_p0, y_p0, z_p0 = 0, 0, 0
        else:
            x_p0 = point0_xyz[0]
            y_p0 = point0_xyz[1]
            z_p0 = point0_xyz[2]

        if dipDirectionDegrees not in DIP_DIRECTION_DEGREE_OPTIONS:
            if 0 < dipDirectionDegrees < 90:
                dx = L1 * np.sin(np.radians(dipDirectionDegrees))
                dy = L1 * np.cos(np.radians(dipDirectionDegrees))
                x_p1 = x_p0 + dx
                y_p1 = y_p0 + dy
                z_p1 = (-L1 * np.tan(np.radians(dip))) + z_p0

                dx = L2 * np.cos(np.radians(dipDirectionDegrees))
                dy = L2 * np.sin(np.radians(dipDirectionDegrees))
                x_p2 = x_p0 - dx
                y_p2 = y_p0 + dy
                z_p2 = z_p0

                dx = L3 * np.cos(np.radians(dipDirectionDegrees))
                dy = L3 * np.sin(np.radians(dipDirectionDegrees))
                x_p3 = x_p0 + dx
                y_p3 = y_p0 - dy
                z_p3 = z_p0

                dx = L4 * np.sin(np.radians(dipDirectionDegrees))
                dy = L4 * np.cos(np.radians(dipDirectionDegrees))
                x_p4 = x_p0 + dx
                y_p4 = y_p0 + dy
                z_p4 = (-L4 * np.tan(np.radians(dip))) + z_p0

            elif 90 < dipDirectionDegrees < 180:
                dx = L1 * np.cos(np.radians(dipDirectionDegrees - 90))
                dy = L1 * np.sin(np.radians(dipDirectionDegrees - 90))
                x_p1 = x_p0 + dx
                y_p1 = y_p0 - dy
                z_p1 = (-L1 * np.tan(np.radians(dip))) + z_p0

                dx = L2 * np.cos(np.radians(180 - dipDirectionDegrees))
                dy = L2 * np.sin(np.radians(180 - dipDirectionDegrees))
                x_p2 = x_p0 + dx
                y_p2 = y_p0 + dy
                z_p2 = z_p0

                dx = L3 * np.cos(np.radians(180 - dipDirectionDegrees))
                dy = L3 * np.sin(np.radians(180 - dipDirectionDegrees))
                x_p3 = x_p0 - dx
                y_p3 = y_p0 - dy
                z_p3 = z_p0

                dx = L4 * 2 * np.cos(np.radians(dipDirectionDegrees - 90))
                dy = L4 * 2 * np.sin(np.radians(dipDirectionDegrees - 90))
                x_p4 = x_p0 + dx
                y_p4 = y_p0 - dy
                z_p4 = (-L4 * 2 * np.tan(np.radians(dip))) + z_p0

            elif 180 < dipDirectionDegrees < 270:
                dx = L1 * np.sin(np.radians(dipDirectionDegrees - 180))
                dy = L1 * np.cos(np.radians(dipDirectionDegrees - 180))
                x_p1 = x_p0 - dx
                y_p1 = y_p0 - dy
                z_p1 = (-L1 * np.tan(np.radians(dip))) + z_p0

                dx = L2 * np.cos(np.radians(dipDirectionDegrees - 180))
                dy = L2 * np.sin(np.radians(dipDirectionDegrees - 180))
                x_p2 = x_p0 + dx
                y_p2 = y_p0 - dy
                z_p2 = z_p0

                dx = L3 * np.cos(np.radians(dipDirectionDegrees - 180))
                dy = L3 * np.sin(np.radians(dipDirectionDegrees - 180))
                x_p3 = x_p0 - dx
                y_p3 = y_p0 + dy
                z_p3 = z_p0

                dx = L4 * np.sin(np.radians(dipDirectionDegrees - 180))
                dy = L4 * np.cos(np.radians(dipDirectionDegrees - 180))
                x_p4 = x_p0 - dx
                y_p4 = y_p0 - dy
                z_p4 = (-L4 * np.tan(np.radians(dip))) + z_p0

            elif 270 < dipDirectionDegrees < 360:
                dx = L1 * np.cos(np.radians(dipDirectionDegrees - 270))
                dy = L1 * np.sin(np.radians(dipDirectionDegrees - 270))
                x_p1 = x_p0 - dx
                y_p1 = y_p0 + dy
                z_p1 = (-L1 * np.tan(np.radians(dip))) + z_p0

                dx = L2 * np.cos(np.radians(360 - dipDirectionDegrees))
                dy = L2 * np.sin(np.radians(360 - dipDirectionDegrees))
                x_p2 = x_p0 - dx
                y_p2 = y_p0 - dy
                z_p2 = z_p0

                dx = L3 * np.cos(np.radians(360 - dipDirectionDegrees))
                dy = L3 * np.sin(np.radians(360 - dipDirectionDegrees))
                x_p3 = x_p0 + dx
                y_p3 = y_p0 + dy
                z_p3 = z_p0

                dx = L4 * 2 * np.cos(np.radians(dipDirectionDegrees - 270))
                dy = L4 * 2 * np.sin(np.radians(dipDirectionDegrees - 270))
                x_p4 = x_p0 - dx
                y_p4 = y_p0 + dy
                z_p4 = (-L4 * 2 * np.tan(np.radians(dip))) + z_p0

        elif dipDirectionDegrees in [0, 360]:
            dx = 0
            dy = L1
            x_p1 = x_p0 + dx
            y_p1 = y_p0 + dy
            z_p1 = (-np.tan(np.radians(dip)) * L1) + z_p0

            dx = L2
            dy = 0
            x_p2 = x_p0 + dx
            y_p2 = y_p0 + dy
            z_p2 = z_p0

            dx = L3
            dy = 0
            x_p3 = x_p0 - dx
            y_p3 = y_p0 - dy
            z_p3 = z_p0

            dx = 0
            dy = L4
            x_p4 = x_p0 + dx
            y_p4 = y_p0 + dy
            z_p4 = (-np.tan(np.radians(dip)) * L4) + z_p0

        elif dipDirectionDegrees == 180:
            dx = 0
            dy = L1
            x_p1 = x_p0 + dx
            y_p1 = y_p0 - dy
            z_p1 = (-np.tan(np.radians(dip)) * L1) + z_p0

            dx = L2
            dy = 0
            x_p2 = x_p0 + dx
            y_p2 = y_p0 + dy
            z_p2 = z_p0

            dx = L3
            dy = 0
            x_p3 = x_p0 - dx
            y_p3 = y_p0 - dy
            z_p3 = z_p0

            dx = 0
            dy = L4 * 2
            x_p4 = x_p0 + dx
            y_p4 = y_p0 - dy
            z_p4 = (-np.tan(np.radians(dip)) * L4 * 2) + z_p0

        elif dipDirectionDegrees == 90:
            dx = L1
            dy = 0
            x_p1 = x_p0 + dx
            y_p1 = y_p0 + dy
            z_p1 = (-np.tan(np.radians(dip)) * L1) + z_p0

            dx = 0
            dy = L2
            x_p2 = x_p0 + dx
            y_p2 = y_p0 + dy
            z_p2 = z_p0

            dx = 0
            dy = L3
            x_p3 = x_p0 - dx
            y_p3 = y_p0 - dy
            z_p3 = z_p0

            dx = L4
            dy = 0
            x_p4 = x_p0 + dx
            y_p4 = y_p0 + dy
            z_p4 = (-np.tan(np.radians(dip)) * L4) + z_p0

        elif dipDirectionDegrees == 270:
            dx = L1
            dy = 0
            x_p1 = x_p0 - dx
            y_p1 = y_p0 + dy
            z_p1 = (-np.tan(np.radians(dip)) * L1) + z_p0

            dx = 0
            dy = L2
            x_p2 = x_p0 + dx
            y_p2 = y_p0 + dy
            z_p2 = z_p0

            dx = 0
            dy = L3
            x_p3 = x_p0 - dx
            y_p3 = y_p0 - dy
            z_p3 = z_p0

            dx = L4
            dy = 0
            x_p4 = x_p0 - dx
            y_p4 = y_p0 + dy
            z_p4 = (-np.tan(np.radians(dip)) * L4) + z_p0

        x_points = [x_p0, x_p1, x_p2, x_p3, x_p4]
        y_points = [y_p0, y_p1, y_p2, y_p3, y_p4]
        z_points = [z_p0, z_p1, z_p2, z_p3, z_p4]

        return x_points, y_points, z_points
    
    def obtain_dip_direction_degrees(self):
        """
        Function that provides the dip direction in degrees clockwise from north,
        so that 90 is east, 180 is south, and 270 is west.

        :param strike: strike of the units in degrees clockwise from north, where
            0 is north, 90 is east, 180 is south, and 270 is west.
        :type strike: int or float

        :param dipDirection: Direction of dip, provided in a cardinal direction -
            N, E, S, W, NE, SE, SW, or NW.
        :type dipDirection: str

        :returns: dipDirectionDegrees
        """
        par_values = self.get_par_values(['strike'])
        strike = par_values[0]
        
        if not self.dipDirection:
            dipDirectionDegrees = strike + 90
            
            if dipDirectionDegrees >= 360:
                dipDirectionDegrees -= 360
            
        else:
            if self.dipDirection not in DIP_DIRECTION_OPTIONS:
                err_msg = ''.join([
                    'Dip direction provided to DippingStratigraphy component ', 
                    '{} does not match any of the available options ({}).'.format(
                        self.dipDirection, DIP_DIRECTION_OPTIONS)])
                raise KeyError(err_msg)

            dipDirectionDegrees = None
            if 0 < strike < 90:
                if self.dipDirection in ['S', 'E', 'SE']:
                    dipDirectionDegrees = strike + 90
                elif self.dipDirection in ['N', 'W', 'NW']:
                    dipDirectionDegrees = strike + 270
            elif 90 < strike < 180:
                if self.dipDirection in ['N', 'E', 'NE']:
                    dipDirectionDegrees = strike - 90
                elif self.dipDirection in ['S', 'W', 'SW']:
                    dipDirectionDegrees = strike + 90
            elif 180 < strike < 270:
                if self.dipDirection in ['N', 'W', 'NW']:
                    dipDirectionDegrees = strike + 90
                elif self.dipDirection in ['S', 'E', 'SE']:
                    dipDirectionDegrees = strike - 90
            elif 270 < strike < 360:
                if self.dipDirection in ['N', 'E', 'NE']:
                    dipDirectionDegrees = (strike + 90) - 360
                elif self.dipDirection in ['S', 'W', 'SW']:
                    dipDirectionDegrees = strike - 90
            elif strike in [0, 180, 360]:
                if self.dipDirection == 'E':
                    dipDirectionDegrees = 90
                elif self.dipDirection == 'W':
                    dipDirectionDegrees = 270
            elif strike in [90, 270]:
                if self.dipDirection == 'N':
                    dipDirectionDegrees = 0
                elif self.dipDirection == 'S':
                    dipDirectionDegrees = 180

            if dipDirectionDegrees is None:
                err_msg = ''.join(['The dip direction provided does not work ',
                                   'with the strike provided. Dip must be ',
                                   'orthogonal to strike. For example, the strike ',
                                   'and dip cannot both be to the NE. A unit '
                                   'striking to the NE can have dips of E, S, ',
                                   'or SE (all for dip to the SE) or N, W, or NW',
                                   '(all for dip to the NW). A unit striking ',
                                   'N or S (0 or 180 degrees) can only dip E or W. ',
                                   'A unit striking E or W (90 or 270 degrees) ',
                                   'can only dip N or S.'])
                raise ValueError(err_msg)

        return dipDirectionDegrees
    
    def depth_increases_from_strike_dip(self, x_locations, y_locations,
                                        output_type='single_point'):
        """
        This function calculates the changes in unit depths due to a prescribed
        strike and dip.

        :param x_locations: x value or a list of x values for the location(s) at which
            the function estimates changes in depth due to strike and dip. Note
            that x is assumed to increase to the east.
        :type x_locations: int, float, or list

        :param y_locations: y value or a list of y values for the location(s) at which
            the function estimates changes in depth due to strike and dip. Note
            that y is assumed to increase to the north.
        :type x_locations: int, float, or list

        :param output_type: Option specifying the type of analysis. The first
            option is 'single_point', where the function is given the x and y
            coordinates for a single location and the function then provides the
            change in unit depth at that location (due to the strike and dip data
            provided). The second option is 'point_list', where the function is
            given lists of x and y values for multiple locations. The lists are
            assumed to be aligned, so that the first x value and first y value
            represent one location (etc.). With 'point_list', the function returns
            the changes in elevation at each of the locations provided. The third
            option is 'grid'. While 'point_list' produces a 1-dimensional list of
            depth changes for a 1-dimensional list of x and y coordinates, the grid
            option produces a 2-dimensional array of depth changes for 2-D grid of
            x and y points.
        :type output_type: str

        :returns: depth_increase
        """
        if isinstance(x_locations, list) and isinstance(y_locations, list):
            if len(x_locations) != len(y_locations):
                err_msg = ''.join([
                    'The x_locations and y_locations lists provided to the ',
                    'DippingStratigraphy component method depth_increases_from_strike_dip() ',
                    'do not have equal lengths.'])
                raise ValueError(err_msg)

        if isinstance(x_locations, list) and not isinstance(y_locations, list):
            err_msg = ''.join([
                'The x_locations provided to the DippingStratigraphy component ', 
                'method depth_increases_from_strike_dip() were a list, but the ', 
                'y_locations were not a list. Check your input.'])
            raise TypeError(err_msg)

        if not isinstance(x_locations, list) and isinstance(y_locations, list):
            err_msg = ''.join([
                'The y_locations provided to the DippingStratigraphy component ', 
                'method depth_increases_from_strike_dip() were a list, but the ', 
                'x_locations were not a list. Check your input.'])
            raise TypeError(err_msg)

        # Create a 3-dimensional plane representing the changes in the depth of
        # each unit's base. To create the plane, use the change in depths for 3
        # locations. One location is the reference point (x = locXRef, y = locYRef,
        # and z = 0). The two other locations will have z values calculated with 
        # the strike and dip. Point p1 is in the dip direction, while point p2 
        # is 90 degrees counter clockwise from the dip direction (in one of the 
        # directions of strike). Points 3 and 4 aren't needed here.

        x_points, y_points, z_points = self.make_xyz_points_from_strike_and_dip(
            L1=100, L2=100, L3=100, L4=100, 
            point0_xyz=[self.locXRef, self.locYRef, self.locZRef])

        x_p0 = x_points[0]
        x_p1 = x_points[1]
        x_p2 = x_points[2]

        y_p0 = y_points[0]
        y_p1 = y_points[1]
        y_p2 = y_points[2]

        z_p0 = z_points[0]
        z_p1 = z_points[1]
        z_p2 = z_points[2]

        # Set of known x, y, and z points
        points = [[x_p0, y_p0, z_p0],
                  [x_p1, y_p1, z_p1],
                  [x_p2, y_p2, z_p2]]

        p0, p1, p2 = points
        x0, y0, z0 = p0
        x1, y1, z1 = p1
        x2, y2, z2 = p2

        ux, uy, uz = [x1 - x0, y1 - y0, z1 - z0]
        vx, vy, vz = [x2 - x0, y2 - y0, z2 - z0]

        u_cross_v = [(uy * vz) - (uz * vy),
                     (uz * vx) - (ux * vz),
                     (ux * vy) - (uy * vx)]

        point  = np.array(p0)
        normal = np.array(u_cross_v)

        d = -point.dot(normal)

        # Solve for the change in depth at the location - depth_increase is positive
        # when the depth increases and is negative when it decreases.
        if output_type == 'single_point':
            depth_increase = -(-(normal[0] * x_locations)
                            - (normal[1] * y_locations) - d) / normal[2]

        elif output_type == 'point_list':
            depth_increase = []

            for loc_ref, xval in enumerate(x_locations):
                depth_increase_val = -(-(normal[0] * xval)
                         - (normal[1] * y_locations[loc_ref]) - d) / normal[2]
                depth_increase.append(depth_increase_val)

        elif output_type == 'grid':
            depth_increase = np.zeros((x_locations.shape[0], x_locations.shape[1]))

            for x_ref in range(x_locations.shape[0]):
                for y_ref in range(x_locations.shape[1]):
                    depth_increase[x_ref, y_ref] = -(-(normal[0] * x_locations[x_ref, y_ref])
                             - (normal[1] * y_locations[x_ref, y_ref]) - d) / normal[2]

        return depth_increase
    
    def update_stratigraphy_by_strike_and_dip(self, location_x, location_y, 
                                              updated_strat=None, ind=None):
        """
        This function provides updated thicknesses for shales, aquifers, and the
        reservoir based on strike and dip information as well as the location at
        which the estimates are to be made (location_x, location_y).

        :param location_x: x value [|m|] for the location at which the function is
            estimating unit thicknesses.
        :type location_x: int, float, or numpy.ndarray (if output_in_progress)

        :param location_y: y value [|m|] for the location at which the function is
            estimating unit thicknesses.
        :type location_y: int, float, or numpy.ndarray (if output_in_progress)
        
        :param updated_strat: Output dictionary. If set to None, updated_strat 
            is made as a dictionary containing one output value per key. Otherwise, 
            for each key (e.g., 'shale1Thickness', 'aquifer2MidDepth', or 'reservoirTopDepth') 
            the dictionary is assumed to contain numpy array contaning thicknesses 
            and depths at different locations corresponding with index ind in 
            the array. If so, the dictionary is updated repeatedly with this method.
        :type  updated_strat: None or numpy.ndarray

        :returns: updated_strat, a dictionary containing the updated
            thicknesses and depths for all shales, aquifers, and the reservoir.
        """
        numberOfShaleLayers = self.get_num_shale_layers()
        
        shaleThicknessList, aquiferThicknessList, reservoirThickness = \
            self.get_unit_thicknesses_ref_point()
        
        depth_increase = self.depth_increases_from_strike_dip(
            location_x, location_y, output_type='single_point')

        reservoirThicknessUpdated = reservoirThickness

        shaleThicknessListUpdated = [None] * numberOfShaleLayers
        aquiferThicknessListUpdated = [None] * (numberOfShaleLayers - 1)

        # First, handle the highest shale
        shaleThickness = shaleThicknessList[-1] + depth_increase

        # If the depth change has made the shaleThickness fall beneath the
        # minimum value, set it to the minimum value
        if shaleThickness < self.pars_bounds['shaleThickness'][0]:
            debug_msg = self.layer_thickness_bound_debug_message(
                'shale {}'.format(numberOfShaleLayers), 'minimum', location_x, location_y)
            logging.debug(debug_msg)
            additional_depth_increase = shaleThickness - self.pars_bounds[
                'shaleThickness'][0]
            # Make the aquifer beneath this shale thinner to accomodate
            # the additional_depth_increase
            aquiferThicknessList[-1] += additional_depth_increase
            del additional_depth_increase
            shaleThickness = self.pars_bounds['shaleThickness'][0]

        # If the depth change has made the shaleThickness go above the maximum
        # value, set it to the maximum value
        if shaleThickness > self.pars_bounds['shaleThickness'][1]:
            debug_msg = self.layer_thickness_bound_debug_message(
                'shale {}'.format(numberOfShaleLayers), 'maximum', location_x, location_y)
            logging.debug(debug_msg)
            additional_depth_increase = shaleThickness - self.pars_bounds[
                'shaleThickness'][1]
            # Make the aquifer beneath this shale thicker to accomodate
            # the additional_depth_increase
            aquiferThicknessList[-1] += additional_depth_increase
            del additional_depth_increase
            shaleThickness = self.pars_bounds['shaleThickness'][1]

        shaleThicknessListUpdated[-1] = shaleThickness

        # Now, go from the highest aquifer to the lowest shale so that required
        # thickness changes can be carried down through the units
        for shaleRef in range(numberOfShaleLayers - 2, -1, -1):
            aquiferThickness = aquiferThicknessList[shaleRef]
            shaleThickness = shaleThicknessList[shaleRef]

            # If the depth change has made the aquiferThickness fall beneath the
            # minimum value, set it to the minimum value
            if aquiferThickness < self.pars_bounds['aquiferThickness'][0]:
                debug_msg = self.layer_thickness_bound_debug_message(
                    'aquifer {}'.format(shaleRef+1), 'minimum', location_x, location_y)
                logging.debug(debug_msg)
                additional_depth_increase = aquiferThickness \
                    - self.pars_bounds['aquiferThickness'][0]
                # Make the shale beneath this aquifer thinner to accomodate
                # the additional_depth_increase
                shaleThickness += additional_depth_increase
                del additional_depth_increase
                aquiferThickness = self.pars_bounds['aquiferThickness'][0]

            # If the depth change has made the aquiferThickness go above the maximum
            # value, set it to the maximum value
            if aquiferThickness > self.pars_bounds['aquiferThickness'][1]:
                debug_msg = self.layer_thickness_bound_debug_message(
                    'aquifer {}'.format(shaleRef+1), 'maximum', location_x, location_y)
                logging.debug(debug_msg)
                additional_depth_increase = aquiferThickness \
                    - self.pars_bounds['aquiferThickness'][1]
                # Make the shale beneath this aquifer thicker to accomodate
                # the additional_depth_increase
                shaleThickness += additional_depth_increase
                del additional_depth_increase
                aquiferThickness = self.pars_bounds['aquiferThickness'][1]

            aquiferThicknessListUpdated[shaleRef] = aquiferThickness

            # If the depth change has made the shaleThickness fall beneath the
            # minimum value, set it to the minimum value
            if shaleThickness < self.pars_bounds['shaleThickness'][0]:
                debug_msg = self.layer_thickness_bound_debug_message(
                    'shale {}'.format(shaleRef+1), 'minimum', location_x, location_y)
                logging.debug(debug_msg)
                additional_depth_increase = shaleThickness \
                    - self.pars_bounds['shaleThickness'][0]
                # Make the aquifer beneath this shale thinner to accomodate
                # the additional_depth_increase. If it's shale 1, make the reservoir thinner.
                if (shaleRef + 1) != 1:
                    aquiferThicknessList[shaleRef - 1] += additional_depth_increase
                elif (shaleRef + 1) == 1:
                    reservoirThicknessUpdated += additional_depth_increase
                del additional_depth_increase
                shaleThickness = self.pars_bounds['shaleThickness'][0]

            # If the depth change has made the shaleThickness go above the maximum
            # value, set it to the maximum value
            if shaleThickness > self.pars_bounds['shaleThickness'][1]:
                debug_msg = self.layer_thickness_bound_debug_message(
                    'shale {}'.format(shaleRef+1), 'maximum', location_x, location_y)
                logging.debug(debug_msg)
                additional_depth_increase = shaleThickness \
                    - self.pars_bounds['shaleThickness'][1]
                # Make the aquifer beneath this shale thicker to accomodate
                # the additional_depth_increase. If it's shale 1, make the reservoir thicker.
                if (shaleRef + 1) != 1:
                    aquiferThicknessList[shaleRef - 1] += additional_depth_increase
                elif (shaleRef + 1) == 1:
                    reservoirThicknessUpdated += additional_depth_increase
                del additional_depth_increase
                shaleThickness = self.pars_bounds['shaleThickness'][1]

            shaleThicknessListUpdated[shaleRef] = shaleThickness
        # End of loop through shales and aquifers
        del shaleThicknessList, aquiferThicknessList

        if reservoirThicknessUpdated < self.pars_bounds['reservoirThickness'][0]:
            debug_msg = self.layer_thickness_bound_debug_message(
                'reservoir', 'minimum', location_x, location_y, msg_option=2)
            logging.debug(debug_msg)
            reservoirThicknessUpdated = self.pars_bounds[
                'reservoirThickness'][0]

        if reservoirThicknessUpdated > self.pars_bounds['reservoirThickness'][1]:
            debug_msg = self.layer_thickness_bound_debug_message(
                'reservoir', 'maximum', location_x, location_y, msg_option=2)
            logging.debug(debug_msg)
            reservoirThicknessUpdated = self.pars_bounds[
                'reservoirThickness'][1]

        reservoirTopDepthUpdated = np.sum([shaleThicknessListUpdated])
        reservoirTopDepthUpdated += np.sum([aquiferThicknessListUpdated])
        
        reservoirMidDepthUpdated = reservoirTopDepthUpdated + (
            reservoirThicknessUpdated / 2)
        
        reservoirDepthUpdated = reservoirTopDepthUpdated + reservoirThicknessUpdated

        if reservoirTopDepthUpdated < self.pars_bounds['depth'][0]:
            err_msg = ''.join([
                'The depth of the reservoir top went below the minimum value ',
                'allowed at x = {} m and y = {} m. ',
                CHECK_CONDITIONS_MSG]).format(location_x, location_y)
            raise ValueError(err_msg)

        if reservoirTopDepthUpdated > self.pars_bounds['depth'][1]:
            err_msg = ''.join([
                'The depth of the reservoir top went above the maximum value ',
                'allowed at x = {} m and y = {} m. ',
                CHECK_CONDITIONS_MSG]).format(location_x, location_y)
            raise ValueError(err_msg)

        shaleTopDepthListUpdated = [None] * numberOfShaleLayers
        aquiferTopDepthListUpdated = [None] * (numberOfShaleLayers - 1)

        shaleMidDepthListUpdated = [None] * numberOfShaleLayers
        aquiferMidDepthListUpdated = [None] * (numberOfShaleLayers - 1)

        shaleDepthListUpdated = [None] * numberOfShaleLayers
        aquiferDepthListUpdated = [None] * (numberOfShaleLayers - 1)

        # The top depth of the highest shale is 0
        shaleTopDepthListUpdated[-1] = 0
        shaleMidDepthListUpdated[-1] = shaleThicknessListUpdated[-1] / 2
        shaleDepthListUpdated[-1] = shaleThicknessListUpdated[-1]

        for shaleRef in range(numberOfShaleLayers - 2, -1, -1):
            aquiferTopDepthListUpdated[shaleRef] = shaleTopDepthListUpdated[
                shaleRef + 1] + shaleThicknessListUpdated[shaleRef + 1]

            if aquiferTopDepthListUpdated[shaleRef] <= 0:
                err_msg = ''.join([
                    'The top depth of aquifer {} was <= 0 m at x = {} m and ',
                    'y = {} m. Only shale {} should reach the surface. ',
                    CHECK_CONDITIONS_MSG]).format(
                        shaleRef + 1, location_x, location_y, numberOfShaleLayers)
                raise ValueError(err_msg)

            if aquiferTopDepthListUpdated[shaleRef] <= shaleTopDepthListUpdated[shaleRef + 1]:
                err_msg = ''.join([
                    'The top depth of aquifer {} was <= the top depth of shale {} ',
                    'at x = {} m and y = {} m. This aquifer should always be ',
                    'beneath this shale. ', CHECK_CONDITIONS_MSG]).format(
                        shaleRef + 1, shaleRef + 2, location_x, location_y)
                raise ValueError(err_msg)
            
            aquiferMidDepthListUpdated[shaleRef] = aquiferTopDepthListUpdated[
                shaleRef] + (aquiferThicknessListUpdated[shaleRef] / 2)
            
            aquiferDepthListUpdated[shaleRef] = aquiferTopDepthListUpdated[
                shaleRef] + aquiferThicknessListUpdated[shaleRef]

            shaleTopDepthListUpdated[shaleRef] = aquiferTopDepthListUpdated[
                shaleRef] + aquiferThicknessListUpdated[shaleRef]

            if shaleTopDepthListUpdated[shaleRef] <= 0:
                err_msg = ''.join([
                    'The top depth of shale {} was <= 0 m at x = {} m and ',
                    'y = {} m. Only shale {} should reach the surface. ',
                    CHECK_CONDITIONS_MSG]).format(
                        shaleRef + 1, location_x, location_y, numberOfShaleLayers)
                raise ValueError(err_msg)

            if shaleTopDepthListUpdated[shaleRef] <= aquiferTopDepthListUpdated[shaleRef]:
                err_msg = ''.join([
                    'The top depth of shale {} was <= the top depth of aquifer {} ',
                    'at x = {} m and y = {} m. This shale should always be beneath this ',
                    'aquifer. ', CHECK_CONDITIONS_MSG]).format(
                        shaleRef + 1, shaleRef + 1, location_x, location_y)
                raise ValueError(err_msg)
            
            shaleMidDepthListUpdated[shaleRef] = shaleTopDepthListUpdated[
                shaleRef] + (shaleThicknessListUpdated[shaleRef] / 2)
            
            shaleDepthListUpdated[shaleRef] = shaleTopDepthListUpdated[
                shaleRef] + shaleThicknessListUpdated[shaleRef]

        if not updated_strat:
            updated_strat = dict()
        
        for shaleRef in range(numberOfShaleLayers):
            nm = 'shale{}Thickness'.format(shaleRef + 1)
            val = shaleThicknessListUpdated[shaleRef]
            
            updated_strat = self.format_output(nm, val, updated_strat, ind)
            
            nm = 'shale{}Depth'.format(shaleRef + 1)
            val = shaleDepthListUpdated[shaleRef]
            
            updated_strat = self.format_output(nm, val, updated_strat, ind)
            
            nm = 'shale{}MidDepth'.format(shaleRef + 1)
            val = shaleMidDepthListUpdated[shaleRef]
            
            updated_strat = self.format_output(nm, val, updated_strat, ind)
            
            nm = 'shale{}TopDepth'.format(shaleRef + 1)
            val = shaleTopDepthListUpdated[shaleRef]
            
            updated_strat = self.format_output(nm, val, updated_strat, ind)
            
            if (shaleRef + 1) < numberOfShaleLayers:
                nm = 'aquifer{}Thickness'.format(shaleRef + 1)
                val = aquiferThicknessListUpdated[shaleRef]
                
                updated_strat = self.format_output(nm, val, updated_strat, ind)
                
                nm = 'aquifer{}Depth'.format(shaleRef + 1)
                val = aquiferDepthListUpdated[shaleRef]
                
                updated_strat = self.format_output(nm, val, updated_strat, ind)
                
                nm = 'aquifer{}MidDepth'.format(shaleRef + 1)
                val = aquiferMidDepthListUpdated[shaleRef]
                
                updated_strat = self.format_output(nm, val, updated_strat, ind)
                
                nm = 'aquifer{}TopDepth'.format(shaleRef + 1)
                val = aquiferTopDepthListUpdated[shaleRef]
                
                updated_strat = self.format_output(nm, val, updated_strat, ind)
        
        nm = 'reservoirThickness'
        val = reservoirThicknessUpdated
        
        updated_strat = self.format_output(nm, val, updated_strat, ind)
        
        nm = 'reservoirDepth'
        val = reservoirDepthUpdated
        
        updated_strat = self.format_output(nm, val, updated_strat, ind)
        
        nm = 'reservoirMidDepth'
        val = reservoirMidDepthUpdated
        
        updated_strat = self.format_output(nm, val, updated_strat, ind)
        
        nm = 'reservoirTopDepth'
        val = reservoirTopDepthUpdated
        
        updated_strat = self.format_output(nm, val, updated_strat, ind)

        return updated_strat
    
    def format_output(self, nm, val, updated_strat, ind):
        """
        Checks if update_stratigraphy_by_strike_and_dip() is updating a dictionary 
        (updated_strat) containing numpy.ndarrays. If the index ind is not None, 
        the index of the array within updated_strat[nm] is set to val. Otherwise, 
        the output dictionary (updated_strat) for the given key nm is set to val.
        """
        if ind is None:
            updated_strat[nm] = float(val)
        else:
            updated_strat[nm][ind] = float(val)
            
        return updated_strat
    
    def get_thickness_obs_names(self):
        """
        Returns a list of the thickness observation names, given the numberOfShaleLayers.
        """
        numberOfShaleLayers = self.get_num_shale_layers()
        
        obs_names = ['shale{}Thickness'.format(ind) \
                     for ind in range(1, numberOfShaleLayers + 1)] + [
                         'aquifer{}Thickness'.format(ind) \
                             for ind in range(1, numberOfShaleLayers)] + [
                                     'reservoirThickness']
        
        return obs_names
    
    def get_depth_obs_names(self):
        """
        Returns a list of the depth observation names, given the numberOfShaleLayers.
        """
        numberOfShaleLayers = self.get_num_shale_layers()
        
        obs_names = ['shale{}Depth'.format(ind) \
                     for ind in range(1, numberOfShaleLayers + 1)] + [
                         'aquifer{}Depth'.format(ind) \
                             for ind in range(1, numberOfShaleLayers)] + [
                                     'reservoirDepth']
        
        obs_names += ['shale{}MidDepth'.format(ind) \
                     for ind in range(1, numberOfShaleLayers + 1)] + [
                         'aquifer{}MidDepth'.format(ind) \
                             for ind in range(1, numberOfShaleLayers)] + [
                                     'reservoirMidDepth']
        
        obs_names += ['shale{}TopDepth'.format(ind) \
                     for ind in range(1, numberOfShaleLayers + 1)] + [
                         'aquifer{}TopDepth'.format(ind) \
                             for ind in range(1, numberOfShaleLayers)] + [
                                     'depth', 'reservoirTopDepth']
        
        return obs_names
    
    def check_data_for_ref_loc(self, comp_data):
        """
        Checks the dictionary comp_data for input related to the reference 
        location's x, y, and z coordinates. If the input is present and formatted 
        correctly, it is given to the component.
        """
        if 'ReferenceLocation' in comp_data:
            if 'coordx' in comp_data['ReferenceLocation']:
                locXRef = comp_data['ReferenceLocation']['coordx']
                
                if isinstance(locXRef, list):
                    if len(locXRef) == 1:
                        self.locXRef = locXRef[0]
                    else:
                        logging.warning(REFERENCE_LOC_WARNING.format(
                            'x', locXRef, 'x', self.locX))
                elif isinstance(locXRef, (int, float)):
                    self.locXRef = comp_data['ReferenceLocation']['coordx']
                else:
                    logging.warning(REFERENCE_LOC_WARNING.format(
                        'x', locXRef, 'x', self.locX))
            
            if 'coordy' in comp_data['ReferenceLocation']:
                locYRef = comp_data['ReferenceLocation']['coordy']
                
                if isinstance(locYRef, list):
                    if len(locYRef) == 1:
                        self.locYRef = locYRef[0]
                    else:
                        logging.warning(REFERENCE_LOC_WARNING.format(
                            'y', locYRef, 'y', self.locYRef))
                elif isinstance(locYRef, (int, float)):
                    self.locYRef = comp_data['ReferenceLocation']['coordy']
                else:
                    logging.warning(REFERENCE_LOC_WARNING.format(
                        'y', locYRef, 'y', self.locYRef))
            
            if 'coordz' in comp_data['ReferenceLocation']:
                locZRef = comp_data['ReferenceLocation']['coordz']
                
                if isinstance(locZRef, list):
                    if len(locZRef) == 1:
                        self.locZRef = locZRef[0]
                    else:
                        logging.warning(REFERENCE_LOC_WARNING.format(
                            'z', locZRef, 'z', self.locZRef))
                elif isinstance(locZRef, (int, float)):
                    self.locZRef = comp_data['ReferenceLocation']['coordz']
                else:
                    logging.warning(REFERENCE_LOC_WARNING.format(
                        'z', locZRef, 'z', self.locZRef))
    
    def connect_with_system(self, component_data):
        """
        Code to add stratigraphy to system model for control file interface.
        """
        if 'ReferenceLocation' in component_data:
            if 'coordx' in component_data['ReferenceLocation']:
                locXRef = component_data['ReferenceLocation']['coordx']
                
                if isinstance(locXRef, list):
                    if len(locXRef) == 1:
                        locXRef = locXRef[0]
                    else:
                        logging.warning(REFERENCE_LOC_WARNING.format('x', locXRef, 'x'))
                
                if isinstance(locXRef, (int, float)):
                    self.locXRef = component_data['ReferenceLocation']['coordx']
                else:
                    logging.warning(REFERENCE_LOC_WARNING.format('x', locXRef, 'x'))
            
            if 'coordy' in component_data['ReferenceLocation']:
                locYRef = component_data['ReferenceLocation']['coordy']
                
                if isinstance(locYRef, list):
                    if len(locYRef) == 1:
                        locYRef = locYRef[0]
                    else:
                        logging.warning(REFERENCE_LOC_WARNING.format('y', locYRef, 'y'))
                
                if isinstance(locYRef, (int, float)):
                    self.locYRef = component_data['ReferenceLocation']['coordy']
                else:
                    logging.warning(REFERENCE_LOC_WARNING.format('y', locYRef, 'y'))
            
            if 'coordz' in component_data['ReferenceLocation']:
                locZRef = component_data['ReferenceLocation']['coordz']
                
                if isinstance(locZRef, list):
                    if len(locZRef) == 1:
                        locZRef = locZRef[0]
                    else:
                        logging.warning(REFERENCE_LOC_WARNING.format('z', locZRef, 'z'))
                
                if isinstance(locZRef, (int, float)):
                    self.locZRef = component_data['ReferenceLocation']['coordz']
                else:
                    logging.warning(REFERENCE_LOC_WARNING.format('z', locZRef, 'z'))
        
        if 'Controls' in component_data:
            if 'dipDirection' in component_data['Controls']:
                dipDirection = component_data['Controls']['dipDirection']
                
                if component_data['Controls']['dipDirection'] in DIP_DIRECTION_OPTIONS:
                    self.dipDirection = dipDirection
                else:
                    logging.warning(DIP_DIRECTION_WARNING.format(
                        dipDirection, DIP_DIRECTION_OPTIONS))
                    self.dipDirection = None
        
        # Run this to check is numberOfShaleLayers is stochastic, which does not 
        # work for the control file interface.
        _ = self.get_num_shale_layers(cfi=True)

        # Check whether all other needed parameters are defined by user
        par_names = self.get_thickness_obs_names()

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
        
        # For the control file interface, locX and locY are handled in the 
        # function strata.process_spatially_variable_strata()
        self.grid_obs_requested = False
        
        # Run this to check is numberOfShaleLayers is stochastic, which does not 
        # work for the control file interface.
        _ = self.get_num_shale_layers(cfi=True)
        
        thickness_obs = self.get_thickness_obs_names()
        
        for ob_nm in thickness_obs:
            self.add_obs(ob_nm, index=[0])
            self.add_obs_to_be_linked(ob_nm)
        
        depth_obs = self.get_depth_obs_names()
        
        for ob_nm in depth_obs:
            self.add_obs(ob_nm, index=[0])
            self.add_obs_to_be_linked(ob_nm)
    
    def simulation_model(self, p, time_point=365.25, locX=None, locY=None):
        """
        Return unit thicknesses and depths at the location given.
        
        :param p: input parameters of LookupTableStratigraphy model
        :type p: dict
        
        :param time_point: time point (in days). This input is not used here, 
            but it is kept for compatibility.
        :type time_point: float

        :param locX: x-coordinate of the location at which stratigraphy observations
            are to be calculated.
        :type locX: float or array-like of floats

        :param locY: y-coordinate of the location at which stratigraphy observations
            are to be calculated.
        :type locY: float or array-like of floats

        :returns: out - dictionary of observations of lookup table stratigraphy
            component model; keys (here, N is the numberOfShaleLayers): 
                ['reservoirThickness','shale1Thickness', aquifer1Thickness', 
                 'shale2Thickness', 'aquifer2Thickness', ..., 'shaleNThickness', 
                 'reservoirDepth', 'reservoirMidDepth', 'reservoirTopDepth', 
                 'shale1Depth', 'shale1MidDepth', 'shale1TopDepth', 
                 'aquifer1Depth', 'aquifer1MidDepth', 'aquifer1TopDepth', ..., 
                 'shaleNDepth', 'shaleNMidDepth', 'shaleNTopDepth']
        """
        # Check whether locations are provided as keyword arguments
        if locX is not None and locY is not None:
            # Get the provided locations reformatted
            prov_locX = np.array([locX]).flatten()
            prov_locY = np.array([locY]).flatten()

            # If single point is provided update type of observations
            if len(prov_locX) == 1:
                self.grid_obs_keys = []
            elif len(self.locX) > 1:
                self.grid_obs_keys = self.get_thickness_obs_names()

            # Check consistency of lengths
            if len(prov_locX) != len(prov_locY):
                err_msg = ''.join([
                    'Length of locX and locY arrays provided ',
                    'to the model method are not the same.'])
                logging.error(err_msg)
                raise ValueError(err_msg)

            # Check whether the previously defined locations are the same
            # as the new ones
            if (len(prov_locX) != len(self.locX)) or (
                    np.not_equal(prov_locX, self.locX).any() or np.not_equal(
                        prov_locY, self.locY).any()):
                # Update locations
                self.locX = prov_locX
                self.locY = prov_locY
                self.grid_obs_requested = False
        
        # Initialize output dictionary
        out = dict()

        if not self.grid_obs_requested:
            out = self.update_stratigraphy_by_strike_and_dip(
                self.locX, self.locY)
        else:
            thickness_obs = self.get_thickness_obs_names()
            
            for obs in thickness_obs:
                out[obs] = np.zeros(self.num_points)
            
            depth_obs = self.get_depth_obs_names()
            
            for obs in depth_obs:
                out[obs] = np.zeros(self.num_points)
            
            for ind in range(self.num_points):
                out = self.update_stratigraphy_by_strike_and_dip(
                    self.locX[ind], self.locY[ind], out, ind=ind)

        # Return dictionary of outputs
        return out


if __name__ == "__main__":
    try:
        from openiam import AnalyticalReservoir
    except ImportError as err:
        print('Unable to load IAM class module: '+str(err))

    logging.basicConfig(level=logging.WARNING)
    
    # Reference location for the DippingStratigraphy component
    locXRef = 0
    locYRef = 0
    
    # Location at which to evaluate results
    locX = 3535.35
    locY = 3535.35
    
    # Define keyword arguments of the system model
    num_years = 5
    time_array = 365.25*np.arange(0.0, num_years+1) # time is in days
    sm_model_kwargs = {'time_array': time_array}
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # dipDirection is not specified, so it will be based on the right-hand rule
    strata = sm.add_component_model_object(DippingStratigraphy(
        name='strata', parent=sm, locXRef=locXRef, locYRef=locYRef, 
        locX=locX, locY=locY))

    # Add parameters of stratigraphy component model. Parameters that are not 
    # entered will resort to default values.
    strata.add_par('numberOfShaleLayers', value=3, vary=False)
    strata.add_par('shale1Thickness', value=260.0, vary=False)
    strata.add_par('shale2Thickness', value=120.0, vary=False)
    strata.add_par('shale3Thickness', value=20.0, vary=False)
    strata.add_par('aquifer1Thickness', value=30.0, vary=False)
    strata.add_par('aquifer2Thickness', value=5.0, vary=False)
    strata.add_par('reservoirThickness', value=12.0, vary=False)
    strata.add_par('strike', value=340, vary=False)
    strata.add_par('dip', value=2, vary=False)
    
    thickness_obs = strata.get_thickness_obs_names()
    
    # The thicness and depths observations are only produced in the first time 
    # step, so use index=[0] when adding these observations.
    for ob_nm in thickness_obs:
        strata.add_obs(ob_nm, index=[0])
        strata.add_obs_to_be_linked(ob_nm)
    
    depth_obs = strata.get_depth_obs_names()
    
    for ob_nm in depth_obs:
        strata.add_obs(ob_nm, index=[0])
        strata.add_obs_to_be_linked(ob_nm)
    
    # Add reservoir component
    ares = sm.add_component_model_object(AnalyticalReservoir(name='ares', parent=sm))

    # Add parameters of reservoir component model
    ares.add_par('injRate', value=0.5, vary=False)
    ares.add_par('logResPerm', value=-14.0, vary=False)
    ares.add_par('brineResSaturation', value=0.02, vary=False)
    
    ares.add_par_linked_to_par('numberOfShaleLayers',
                               strata.deterministic_pars['numberOfShaleLayers'])
    ares.add_par_linked_to_par('datumPressure',
                               strata.default_pars['datumPressure'])
    
    ares.add_par_linked_to_obs('shale1Thickness', strata.linkobs['shale1Thickness'])
    ares.add_par_linked_to_obs('shale2Thickness', strata.linkobs['shale2Thickness'])
    ares.add_par_linked_to_obs('shale3Thickness', strata.linkobs['shale3Thickness'])
    ares.add_par_linked_to_obs('aquifer1Thickness', strata.linkobs['aquifer1Thickness'])
    ares.add_par_linked_to_obs('aquifer2Thickness', strata.linkobs['aquifer2Thickness'])
    ares.add_par_linked_to_obs('reservoirThickness', strata.linkobs['reservoirThickness'])

    ares.add_obs('pressure')
    ares.add_obs('CO2saturation')
    
    # Run system model using current values of its parameters
    sm.forward()

    print('------------------------------------------------------------------')
    print('                  Forward method illustration ')
    print('------------------------------------------------------------------')
    
    # Print results
    print('Results at x = {} m, y = {} m: '.format(locX, locY))
    print('Aquifer 1 Bottom Depth (m): ', sm.collect_observations_as_time_series(
        strata, 'aquifer1Depth', indices=[0]), sep='\n')
    print('Aquifer 2 Bottom Depth (m): ', sm.collect_observations_as_time_series(
        strata, 'aquifer2Depth', indices=[0]), sep='\n')
    print('time_array (years): ', time_array / 365.25)
    print('Pressure (Pa): ', sm.collect_observations_as_time_series(ares, 'pressure'),
          sep='\n')
    print('CO2 Saturation: ', sm.collect_observations_as_time_series(ares, 'CO2saturation'),
          sep='\n')
