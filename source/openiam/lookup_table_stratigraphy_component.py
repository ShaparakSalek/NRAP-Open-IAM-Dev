# -*- coding: utf-8 -*-
"""
Last modified: November 16th, 2023

Authors: Seth King, Veronika Vasylkivska, Nate Mitchell
"""
import os
import sys
import logging
import numpy as np

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from openiam import (SystemModel, ComponentModel,
                         StratigraphyDataInterpolator, IAM_DIR)
    from openiam.iam_gridded_observation import interp_weights
except ImportError as err:
    print('Unable to load IAM class module: {}'.format(err))


def find_weights(loc_xyz, triangulation):
    """ Find vertices of the simplex and calculate corresponding weights.

    :param loc_xyz: list of (locX, locY) where locX and locY are all array-like
    :type loc_xyz: list()

    :param triangulation: Delaunay tesselation in 2- or 3-d,
        triangulation for which weights and vertices will be found.
    :type triangulation: scipy.spatial.Delaunay

    :returns: list of num_points, tri_vertices, tri_weights
    """
    # Determine number of dimensions: 2d or 3d
    num_dims = len(loc_xyz)

    # Determine number of different points
    num_points = len(loc_xyz[0])
    tri_vertices = np.zeros((num_points, num_dims+1), dtype=np.intp)
    tri_weights = np.zeros((num_points, num_dims+1))

    # Extract coordinates common for both dimension scenarios
    locX = loc_xyz[0]
    locY = loc_xyz[1]

    # Find vertices and corresponding weights for each point
    for ind in range(num_points):
        tri_vertices[[ind]], tri_weights[[ind]] = (
            interp_weights(triangulation,
                           np.array([[locX[ind], locY[ind]]])))

    # Setup machine epsilon to replace comparison with zero
    machine_eps = np.finfo(np.float32).eps  # 1.1920929e-07

    # Check whether weights are reasonable; if they are not,
    # then it means that the point at which we need to interpolate
    # is outside the data range.
    # The first if condition checks whether any of the weights
    # are negative up to machine epsilon. The second condition
    # checks whether any of the weights exceed 1 up to machine epsilon.
    if (np.any(tri_weights < -machine_eps)) or (
            np.any(tri_weights > (1+machine_eps))):
        err_msg = ''.join(['Interpolation is requested at the ',
                           'location outside the data domain.'])
        logging.error(err_msg)
        raise ValueError(err_msg)

    return num_points, tri_vertices, tri_weights

class LookupTableStratigraphy(ComponentModel):
    """
    The Lookup Table Stratigraphy component can be used to read stratigraphy 
    data from a *.csv* or *.hdf5* file and then utilize this stratigraphy data 
    in a simulation. Unit thicknesses are specified for different coordinates 
    in the file, and spatial interpolation is used to produce unit thicknesses 
    at different locations within the domain. This approach is meant to allow 
    greater flexibility, but it also requires the user to be careful when setting 
    up the data file used. The output from the Lookup Table Stratigraphy component 
    can be checked with the ``Stratigraphy`` and ``StratigraphicColumn`` plot 
    types in the control file interface.
    
    In the NRAP-Open-IAM control file interface, the type name for the Lookup 
    Table Stratigraphy component is ``LookupTableStratigraphy``.
    
    The arrangement of geologic units for this component is the same as the that for
    the ``Stratigraphy`` component. The lowest unit considered is the reservoir. There 
    are alternating shales and aquifers overlying the reservoir. The lowest shale, 
    shale 1, is immediately above the reservoir. Aquifer 1 is above shale 1, and 
    shale 2 is above aquifer 1. This pattern continues until shale ``N``, which 
    extends to the surface, where ``N`` is the number of shale units used. There 
    are ``N - 1`` aquifers, and the number of shale layers is specified with the 
    *numberOfShaleLayers* parameter. The *numberOfShaleLayers* cannot vary across 
    the domain. If one wants a unit to effectively pinch out in certain areas, 
    however, its thickness values for that area in the file can be decreased to 
    to the minimum thickness of 1 |m|.

    In the NRAP-Open-IAM control file interface, the ``FileDirectory`` and 
    ``FileName`` keywords must be specified. The ``FileDirectory`` keyword 
    indicates the directory where the data files for the lookup tables are 
    located. The ``FileName`` keyword is the name of the *.csv* or *.hdf5* file 
    that stores the unit thicknesses. The name should be given with the extension 
    included (e.g., ``stratigraphy.csv``). The first two columns in the file should 
    be labeled as ``x`` and ``y``, with these two columns containing easting and 
    northing distances, respectively. The other columns in the file should be 
    labeled with names that correspond with shale units (``shale#Thickness``), 
    aquifer units (``aquifer#Thickness``), or the reservoir (``reservoirThickness``). 
    Here, the ``#`` character shown for shales and aquifers should be replaced 
    with the index (e.g., ``shale5Thickness`` for shale 5). The unit thicknesses 
    in each row of the file represent the stratigraphy at the corresponding ``x`` 
    and ``y`` values. Unit thicknesses can vary over space, but they cannot change 
    with time in a simulation.

    The Lookup Table Stratigraphy component produces the output using spatial 
    interpolation within the domain defined by the coordinates in the data file 
    used. An error will occur if the component is asked to produce output for 
    a location that lies outside of the domain covered by the file. If this 
    error occurs, the user can expand the domain covered by the data file by adding 
    more rows with unit thicknesses at ``x`` and ``y`` values outside the current 
    ranges.
    
    Descriptions of the component's parameters are provided below.

    * **numberOfShaleLayers** [-] (3 to 30) - number of shale layers in the
      system (default: 3). The shale units must be separated by an aquifer.

    * **datumPressure** [|Pa|] (80,000 to 300,000) - pressure at the top of the
      system (default: 101,325).
    
    The observations from the Lookup Table Stratigraphy component are:

    * **shale#Thickness** [|m|] - thickness of shale unit ``#``, where ``#`` is
      an index ranging from one to *numberOfShaleLayers*. If data for a particular 
      shale layer are not included in the file used, that layer will be assigned 
      a default thickness is 250 |m|.

    * **aquifer#Thickness** [|m|] - thickness of aquifer unit ``#``, where ``#`` 
      is an index ranging from one to (*numberOfShaleLayers* - 1). If data for a 
      particular aquifer layer are not included in the file used, that layer will 
      be assigned a default thickness of 100 |m|.
     
    * **reservoirThickness** [|m|] - thickness of the storage reservoir. If data 
      for the reservoir are not included in the file used, the default reservoir 
      thickness is 50 |m|.
    
    * **shale#Depth** [|m|] - depth to the bottom of shale unit ``#``, where 
      ``#`` is an index ranging from one to *numberOfShaleLayers*.

    * **aquifer#Depth** [|m|] - depth to the bottom of aquifer unit ``#``, where 
      ``#`` is an index ranging from one to (*numberOfShaleLayers* - 1).
      
    * **reservoirDepth** [|m|] - depth to the bottom of the storage reservoir.
    
    * **shale#MidDepth** [|m|] - depth to the middle of shale unit ``#``, where 
      ``#`` is an index ranging from one to *numberOfShaleLayers*.
    
    * **aquifer#MidDepth** [|m|] - depth to the middle of aquifer unit ``#``, where 
      ``#`` is an index ranging from one to (*numberOfShaleLayers* - 1).
    
    * **reservoirMidDepth** [|m|] - depth to the middle of the storage reservoir. 
    
    * **shale#TopDepth** [|m|] - depth to the top of shale unit ``#``, where 
      ``#`` is an index ranging from one to *numberOfShaleLayers*.
    
    * **aquifer#TopDepth** [|m|] - depth to the top of aquifer unit ``#``, where 
      ``#`` is an index ranging from one to (*numberOfShaleLayers* - 1).
    
    * **reservoirTopDepth** [|m|] - depth to the top of the storage reservoir. 
      The same value can also be produced with the observation name **depth**, 
      however. The use of the **depth** output name is included to conform with 
      the ``Stratigraphy`` component.
    
    If the component produces a unit thickness for a shale, aquifer, or reservoir 
    that falls outside the range of 1 |m| to 1600 |m|, then that thickness will 
    be set to the lower or upper limit (whichever is closer). For example, if an 
    interpolated thickness is 0.1 |m|, based on the values included in the data 
    file, then that thickness observation will be set to 1 |m| instead.
    
    Currently, this component cannot be used to produce unit thicknesses that vary 
    stochastically in a simulation. The unit thicknesses and depths produced by 
    a ``LookupTableStratigraphy`` component can vary over space within a simulation, 
    but the thicknesses and depths produced for one location cannot vary across 
    different realizations of a simulation. If one wanted to assess how altering 
    unit thicknesses would impact results in a simulation using Latin Hypercube 
    Sampling (``LHS``), for example, one should run multiple ``LHS`` simulations 
    with different files used for the ``LookupTableStratigraphy`` components 
    (with the files containing different stratigraphy inputs).
    
    For control file examples using the ``LookupTableStratigraphy`` component, 
    see *ControlFile_ex32b.yaml* to *ControlFile_ex32c.yaml*, *ControlFile_ex38a.yaml* 
    to *ControlFile_ex38c.yaml*, *ControlFile_ex39c.yaml*, *ControlFile_ex55b.yaml*, 
    and *ControlFile_ex55d.yaml*. For script example examples, see *iam_sys_lutstrata.py*, 
    *iam_sys_lutstrata_gridded.py*, *iam_sys_lutstrata_reservoir_openwell.py*, 
    *iam_sys_lutstrata_reservoir_mswell.py*, and *iam_sys_lutstrata_reservoir_cmwell.py*.
    """
    def __init__(self, name, parent, intr_family=None, locX=None, locY=None, 
                 file_directory=None, file_name=None):
        """
        Constructor method of LookupTableStratigraphy class.

        :param name: name of component model
        :type name: str

        :param parent: the SystemModel object that the component model
            belongs to
        :type parent: SystemModel object

        :param intr_family: name of interpolator family which belongs to the parent
            system model and whose data will be used for stratigraphy observations
        :type intr_family: str

        :param locX: x-coordinate of the location at which stratigraphy information 
            is to to be calculated.
            By default, it is None, which means that the whole grid data is requested
        :type locX: float or array-like of floats

        :param locY: y-coordinate of the location at which stratigraphy information 
            is to to be calculated.
            By default, it is None, which means that the whole grid data is requested
        :type locY: float or array-like of floats

        :param file_directory: location (directory) of the stratigraphy data 
            that will be used to create stratigraphy family of interpolators 
            in the case those are not created before the LookupTableStratigraphy 
            component.
        :type file_directory: str

        :param file_name: name of a *.csv or *.hdf5 file containing unit thicknesses
            at different coordinates in the domain.
        :type file_name: str

        :returns: LookupTableStratigraphy class object
        """
        # Set up keyword arguments of the 'model' method provided by the system model
        model_kwargs = {'time_point': 365.25}  # default value of 365.25 days

        super().__init__(name, parent, model=self.simulation_model,
                         model_kwargs=model_kwargs)
        
        # The model only needs to be run once
        self.run_frequency = 1
        self.default_run_frequency = 1

        # Add type attribute
        self.class_type = 'LookupTableStratigraphy'

        # Setup attributes related to interpolator family
        self.intr_family = intr_family
        self.linked_to_intr_family = False
        # Dictionary of pairs (index, name)
        self.intr_names = None

        self.file_directory = file_directory
        
        self.file_name = file_name

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
            
            self.grid_obs_requested = False

        else:
            self.grid_obs_keys = []
            self.locX = np.array([None])
            self.locY = np.array([None])
            # Add flag indicating whether whole grid data is needed
            self.grid_obs_requested = True

        # Reserve space for additional attributes
        self.num_points = None
        self.tri_vertices = None
        self.tri_weights = None
        
        # Add one of the default parameters
        self.add_default_par('numberOfShaleLayers', value=3)
        self.add_default_par('datumPressure', value=101325.0)
        
        self.pars_bounds = dict()
        self.pars_bounds['numberOfShaleLayers'] = [3, 30]
        self.pars_bounds['datumPressure'] = [8.0e+4, 3.0e+5]

        # Setup default observations of the component
        self.default_obs = {'shaleThickness': 250,
                            'aquiferThickness': 100, 
                            'reservoirThickness': 50}

        if intr_family and intr_family in self._parent.interpolators:
            self.link_to_interpolators()

        # Log creating the component
        debug_msg = 'LookupTableStratigraphy component created with name {}'.format(self.name)
        logging.debug(debug_msg)

    def connect_with_system(self, component_data, *args, **kwargs):
        """
        Code to add LookupTableStratigraphy to system model for control file interface.

        :param component_data: Dictionary of component data
        :type component_data: dict

        :returns: None
        """
        # Run this to check if numberOfShaleLayers is stochastic, which does not 
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
        
        err_msg = ''.join(['{} must be provided for ',
                           'LookupTableStratigraphy Component.'])
        
        if 'FileDirectory' in component_data:
            self.file_directory = os.path.join(IAM_DIR,
                                               component_data['FileDirectory'])
        else:
            logging.error(err_msg.format('FileDirectory'))
            raise IOError(err_msg.format('FileDirectory'))
        
        if 'FileName' in component_data:
            self.file_name = component_data['FileName']
        else:
            logging.error(err_msg.format('FileName'))
            raise IOError(err_msg.format('FileName'))

        # For the control file interface, locX and locY are handled in the 
        # function strata.process_spatially_variable_strata()
        self.grid_obs_requested = False

        # If 'stratigraphy' interpolator family was not created before,
        # we need to create it first, then link it to the stratigraphy component(s)
        if 'stratigraphy' not in self._parent.interpolators:
            # Specify the name of interpolator family to be created to make it explicit
            self.build_and_link_interpolators(intr_family='stratigraphy')

        else:
            self.link_to_interpolators(intr_family='stratigraphy')
    
    def parameter_assignment_warning_msg(self):
        """
        Generates a warning message. This warning message is meant to be used 
        in the connect_with_system() methods of components in the control file 
        interface. Specifically, it is used when the parameter of the component 
        is assigned as a stochastic or deterministic parameter when it should 
        instead be linked to an observation from the LookupTableStratigraphy 
        component.
        """
        warning_msg = ''.join([
            'When using LookupTableStratigraphy, the {} parameter ', 
            'cannot be specified as a stochastic or deterministic ', 
            'parameter. The parameter needs to be linked to the ', 
            '{} observation of the LookupTableStratigraphy component. ', 
            'The input provided for that parameter will not be used. ', 
            'In the control file interface, the parameter linkage is performed ', 
            'automatically. In a script application, the parameter can be linked ', 
            'with the add_par_linked_to_obs() method of the component with the ', 
            'parameter being linked to a LookupTableStratigraphy observation. ', 
            'Before linking the LookupTableStratigraphy observation to the ', 
            'parameter, use the add_obs() and add_obs_to_be_linked() methods ', 
            'of the LookupTableStratigraphy component. The add_obs() method ', 
            'should include the observation name and index=[0]. For example, ', 
            'LUTSComponentName.add_obs(obs_name, index=[0]).'])
        
        return warning_msg
    
    def build_and_link_interpolators(self, file_directory=None, file_name=None, 
                                     intr_family=None, default_values=None, 
                                     build_on_the_fly=False):
        """
        Builds interpolators for component model to use.

        :param file_directory: Directory that contains the time file,
            the parameters file, and all realization files.
        :type file_directory: str

        :param intr_family: name of interpolator family whose data is used
            for stratigraphy calculations
        :type intr_family: str

        :param default_values: dictionary of default values of observations
            which are not provided in the lookup tables. By default,
            it is assumed that data provided in the lookup tables is enough,
            thus, by default, the parameter is None.
            The values provided are assumed to be constant for all time steps
            and all spatial points.
        type default_values: dict()

        :param build_on_the_fly: flag variable indicating whether the data
            from the lookup tables corresponding to the linked family
            of interpolators will be read only if needed, e.g. at the first
            call of model method utilizing a particular interpolator. By default,
            all interpolators are created before the model method of the current
            component (value=False).
        :type build_on_the_fly: boolean

        :returns: None
        """
        if file_directory:
            self.file_directory = file_directory
        
        if file_name:
            self.file_name = file_name
        
        if intr_family:
            self.intr_family = intr_family
        else:
            self.intr_family = 'stratigraphy'
        
        # Only need one interpolator for LookupTableStratigraphy
        if 'stratigraphy' in self._parent.interpolators:
            if not 'int1' in self._parent.interpolators[intr_family].items():
                _ = self._parent.add_interpolator(
                    StratigraphyDataInterpolator(
                        name='int1',
                        parent=self._parent,
                        header_file_dir=self.file_directory,
                        data_file=self.file_name,
                        default_values=default_values,
                        build_on_the_fly=False),
                    intr_family=self.intr_family,
                    creator=self)
        else:
            _ = self._parent.add_interpolator(
                StratigraphyDataInterpolator(
                    name='int1',
                    parent=self._parent,
                    header_file_dir=self.file_directory,
                    data_file=self.file_name,
                    default_values=default_values,
                    build_on_the_fly=False),
                intr_family=self.intr_family,
                creator=self)

        logging.debug('Interpolator created for {}.'.format(self.name))
        # Link just created interpolators to the LookupTableStratigraphy component
        self.link_to_interpolators()

    def link_to_interpolators(self, intr_family=None):
        """
        Link the component to the interpolators to be used.

        The method links the component to the interpolators to be used, check
        interpolation weights, and sets the default parameter values.

        :param intr_family: name of interpolation family to be used.
        :type intr_family: str

        :returns: None
        """
        # Check whether interpolators family intr_family exists
        if intr_family:
            self.intr_family = intr_family
        else:
            intr_family = self.intr_family

        if intr_family in self._parent.interpolators:
            # Get the first interpolator in the family
            interpr = list(self._parent.interpolators[self.intr_family].values())[0]

            if not self.grid_obs_requested:
                self.num_points, self.tri_vertices, self.tri_weights = find_weights(
                    (self.locX, self.locY), interpr.triangulation)

            # Setup link flag
            self.linked_to_intr_family = True

        else:
            # Show useful message
            err_msg = ''.join(['Attempt to link stratigraphy component {} ',
                               'to the family of interpolators failed: ',
                               'family of interpolators {} does ',
                               'not exist.']).format(self.name, intr_family)
            logging.error(err_msg)
            raise KeyError(err_msg)
    
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
    
    def get_thickness_obs_names(self):
        """
        Returns a list of the thickness observation names, given the numberOfShaleLayers. 
        These are the observations that are expected to be in the data file.
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
        These are observations that are not expected to be in the data file, they 
        are calculated based on unit thicknesses.
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
    
    def format_output(self, out, nm, val, spat_format=False, ind=None):
        """
        Takes output from either calculate_csv_based_output() or 
        calculate_hdf5_based_output() and adjusts the out dictionary entry. 
        If spat_format is True, then the output is treated as an array, with 
        positions defined by ind. If ind is given as 'all', then it loops 
        through all positions and assigns the given value.
        """
        if spat_format:
            if not nm in out:
                out[nm] = np.zeros(self.num_points)
            
            if ind == 'all':
                for ind in range(self.num_points):
                    out[nm][ind] = val
            else:
                out[nm][ind] = val
            
        else:
            out[nm] = val
        
        return out
    
    def calculate_unit_depths(self, out):
        """
        Takes the output dictionary from simulation_model() and gives it any 
        missing unit thicknesses (default values) as well as all unit depths.
        """
        obs_names = self.get_thickness_obs_names()
        
        numberOfShaleLayers = self.get_num_shale_layers()
        if self.grid_obs_requested or self.num_points == 1:
            for ob_nm in obs_names:
                if not ob_nm in out:
                    if 'shale' and 'Thickness' in ob_nm:
                        obs_type = 'shaleThickness'
                    elif 'aquifer' and 'Thickness' in ob_nm:
                        obs_type = 'aquiferThickness'
                    else:
                        obs_type  = ob_nm
                    
                    val = self.default_obs[obs_type]
                    out = self.format_output(out, ob_nm, val)
            
            depth_temp = 0
            ob_nm = 'shale{}TopDepth'.format(numberOfShaleLayers)
            out = self.format_output(out, ob_nm, depth_temp)
            
            depth_temp = out['shale{}Thickness'.format(numberOfShaleLayers)] / 2
            
            ob_nm = 'shale{}MidDepth'.format(numberOfShaleLayers)
            out = self.format_output(out, ob_nm, depth_temp)
            
            depth_temp = out['shale{}Thickness'.format(numberOfShaleLayers)]
            
            ob_nm = 'shale{}Depth'.format(numberOfShaleLayers)
            out = self.format_output(out, ob_nm, depth_temp)
            
            for shaleRef in range(numberOfShaleLayers - 1, 0, -1):
                ob_nm = 'aquifer{}TopDepth'.format(shaleRef)
                out = self.format_output(out, ob_nm, depth_temp)
                
                depth_temp_v2 = depth_temp + (out['aquifer{}Thickness'.format(shaleRef)] / 2)
                
                ob_nm = 'aquifer{}MidDepth'.format(shaleRef)
                out = self.format_output(out, ob_nm, depth_temp_v2)
                
                depth_temp += out['aquifer{}Thickness'.format(shaleRef)]
                
                ob_nm = 'aquifer{}Depth'.format(shaleRef)
                out = self.format_output(out, ob_nm, depth_temp)
                
                ob_nm = 'shale{}TopDepth'.format(shaleRef)
                out = self.format_output(out, ob_nm, depth_temp)
                
                depth_temp_v2 = depth_temp + (out['shale{}Thickness'.format(shaleRef)] / 2)
                
                ob_nm = 'shale{}MidDepth'.format(shaleRef)
                out = self.format_output(out, ob_nm, depth_temp_v2)
                
                depth_temp += out['shale{}Thickness'.format(shaleRef)]
                
                ob_nm = 'shale{}Depth'.format(shaleRef)
                out = self.format_output(out, ob_nm, depth_temp)
                    
            # depth is also equal to shale1Depth
            ob_nm = 'depth'
            out = self.format_output(out, ob_nm, depth_temp)
            
            # Keeping the 'depth' observation to be consistent with the Stratigraphy 
            # component, but here the same value is also represented by the 
            # reservoirTopDepth observation
            ob_nm = 'reservoirTopDepth'
            out = self.format_output(out, ob_nm, depth_temp)
            
            depth_temp_v2 = depth_temp + (out['reservoirThickness'] / 2)
            
            ob_nm = 'reservoirMidDepth'
            out = self.format_output(out, ob_nm, depth_temp_v2)
            
            depth_temp += out['reservoirThickness']
            
            ob_nm = 'reservoirDepth'
            out = self.format_output(out, ob_nm, depth_temp)
            
        else:
            for ob_nm in obs_names:
                obs_type = ob_nm
                
                if not obs_type in out:
                    if 'shale' and 'Thickness' in ob_nm:
                        obs_type = 'shaleThickness'
                    elif 'aquifer' and 'Thickness' in ob_nm:
                        obs_type = 'aquiferThickness'
                    else:
                        obs_type = ob_nm
                    
                    val = self.default_obs[obs_type]
                    out = self.format_output(out, ob_nm, val, 
                                             spat_format=True, ind='all')
            
            # Interpolate over all requested points
            for ind in range(self.num_points):
                depth_temp = 0
                
                ob_nm = 'shale{}TopDepth'.format(numberOfShaleLayers)
                out = self.format_output(out, ob_nm, depth_temp, 
                                         spat_format=True, ind=ind)
                
                ob_nm = 'shale{}Thickness'.format(numberOfShaleLayers)
                depth_temp_v2 = out[ob_nm][ind] / 2
                
                ob_nm = 'shale{}MidDepth'.format(numberOfShaleLayers)
                out = self.format_output(out, ob_nm, depth_temp_v2)
                
                # The highest shale's bottom depth is just that shale's thickness
                ob_nm = 'shale{}Thickness'.format(numberOfShaleLayers)
                
                depth_temp = out[ob_nm][ind]
                
                ob_nm = 'shale{}Depth'.format(numberOfShaleLayers)
                out = self.format_output(out, ob_nm, depth_temp, 
                                         spat_format=True, ind=ind)
                
                for shaleRef in range(numberOfShaleLayers - 1, 0, -1):
                    ob_nm = 'aquifer{}TopDepth'.format(shaleRef)
                    out = self.format_output(out, ob_nm, depth_temp, 
                                             spat_format=True, ind=ind) 
                    
                    ob_nm = 'aquifer{}Thickness'.format(shaleRef)
                    
                    depth_temp_v2 = depth_temp + (out[ob_nm][ind] / 2)
                    
                    ob_nm = 'aquifer{}MidDepth'.format(shaleRef)
                    out = self.format_output(out, ob_nm, depth_temp_v2, 
                                             spat_format=True, ind=ind)
                    
                    ob_nm = 'aquifer{}Thickness'.format(shaleRef)
                    
                    depth_temp += out[ob_nm][ind]
                    
                    ob_nm = 'aquifer{}Depth'.format(shaleRef)
                    out = self.format_output(out, ob_nm, depth_temp, 
                                             spat_format=True, ind=ind) 
                    
                    ob_nm = 'shale{}TopDepth'.format(shaleRef)
                    out = self.format_output(out, ob_nm, depth_temp, 
                                             spat_format=True, ind=ind)
                    
                    ob_nm = 'shale{}Thickness'.format(shaleRef)
                    
                    depth_temp_v2 = depth_temp + (out[ob_nm][ind] / 2)
                    
                    ob_nm = 'shale{}MidDepth'.format(shaleRef)
                    out = self.format_output(out, ob_nm, depth_temp_v2, 
                                             spat_format=True, ind=ind)
                    
                    ob_nm = 'shale{}Thickness'.format(shaleRef)
                    
                    depth_temp += out[ob_nm][ind]
                    
                    ob_nm = 'shale{}Depth'.format(shaleRef)
                    out = self.format_output(out, ob_nm, depth_temp, 
                                             spat_format=True, ind=ind) 
                
                # The bottom depth of shale 1 is also the top of the reservoir
                ob_nm = 'depth'
                out = self.format_output(out, ob_nm, depth_temp, 
                                         spat_format=True, ind=ind)
                
                ob_nm = 'reservoirTopDepth'
                out = self.format_output(out, ob_nm, depth_temp, 
                                         spat_format=True, ind=ind)
                
                ob_nm = 'reservoirThickness'
                depth_temp_v2 = depth_temp + (out[ob_nm][ind] / 2)
                
                ob_nm = 'reservoirMidDepth'
                out = self.format_output(out, ob_nm, depth_temp_v2, 
                                         spat_format=True, ind=ind)
                
                ob_nm = 'reservoirThickness'
                depth_temp += out[ob_nm][ind]
                
                ob_nm = 'reservoirDepth'
                out = self.format_output(out, ob_nm, depth_temp, 
                                         spat_format=True, ind=ind)
            # End of ind loop
        
        return out

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
        loc_updated = False
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
                loc_updated = True

            if loc_updated:
                self.grid_obs_requested = False
                # Find first linked interpolator
                interpr = list(self._parent.interpolators[self.intr_family].values())[0]

                # Recalculate the weights of the simplex and corresponding weights
                self.num_points, self.tri_vertices, self.tri_weights = find_weights(
                    (self.locX, self.locY), interpr.triangulation)

        # Initialize output dictionary
        out = dict()
        
        # Only one interpolator is used in the stratigraphy family, int1
        interpr = self._parent.interpolators[self.intr_family]['int1']

        if self.grid_obs_requested:
            interpr_out = interpr()
            
            for nm in interpr_out:
                out[nm] = interpr_out[nm]
            
        else:
            if self.num_points == 1:
                interpr_out = interpr(
                    self.tri_vertices, self.tri_weights)
                
                # Create dictionary of output
                for nm in interpr_out:
                    out[nm] = interpr_out[nm][0]
                    
            else:
                # Interpolate over all requested points
                for ind in range(self.num_points):
                    interpr_out = interpr(self.tri_vertices[[ind]],
                                          self.tri_weights[[ind]])
                    
                    for nm in interpr_out:
                        if not nm in out:
                            out[nm] = np.zeros(self.num_points)
                        
                        out[nm][ind] = interpr_out[nm][0]

            debug_msg = '{} model call output {}'.format(self.name, out)
            logging.debug(debug_msg)
        
        out = self.calculate_unit_depths(out)
        
        # Return dictionary of outputs
        return out

    def reset(self):
        pass


if __name__ == "__main__":
    logging.basicConfig(level=logging.WARNING)

    file_directory = os.sep.join(['..', '..', 'examples', 'Control_Files', 
                                  'input_data', 'ex38a'])

    # Define keyword arguments of the system model
    num_years = 10
    time_array = 365.25 * np.arange(0.0, num_years + 1)
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
        name='luts', parent=sm, intr_family='stratigraphy', locX=locX, locY=locY))
    
    luts.add_par('numberOfShaleLayers', value=numberOfShaleLayers, vary=False)
    
    thickness_obs = luts.get_thickness_obs_names()
    
    for ob_nm in thickness_obs:
        luts.add_obs(ob_nm, index=[0])
    
    depth_obs = luts.get_depth_obs_names()
    
    for ob_nm in depth_obs:
        luts.add_obs(ob_nm, index=[0])
    
    # Run system model using current values of its parameters
    sm.forward()  # system model is run deterministically

    print('------------------------------------------------------------------')
    print('                  Forward method illustration ')
    print('------------------------------------------------------------------')
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
    print('    aquiferThickness: ', aquifer1Thickness)
    print('    shale2Thickness: ', shale2Thickness)
    print('    aquifer2Thickness: ', aquifer2Thickness)
    print('    shale3Thickness: ', shale3Thickness)
    print('')
    
    # Collect depth observations
    reservoirDepth = sm.collect_observations_as_time_series(
        luts, 'reservoirDepth', indices=[0])
    shale1Depth = sm.collect_observations_as_time_series(
        luts, 'shale1Depth', indices=[0])
    aquifer1Depth = sm.collect_observations_as_time_series(
        luts, 'aquifer1Depth', indices=[0])
    shale2Depth = sm.collect_observations_as_time_series(
        luts, 'shale2Depth', indices=[0])
    aquifer2Depth = sm.collect_observations_as_time_series(
        luts, 'aquifer2Depth', indices=[0])
    shale3Depth = sm.collect_observations_as_time_series(
        luts, 'shale3Depth', indices=[0])
    
    print('Depths (m) to the bottom of each unit:')
    print('    reservoirDepth: ', reservoirDepth)     
    print('    shale1Depth: ', shale1Depth)
    print('    aquiferDepth: ', aquifer1Depth)
    print('    shale2Depth: ', shale2Depth)
    print('    aquifer2Depth: ', aquifer2Depth)
    print('    shale3Depth: ', shale3Depth)
