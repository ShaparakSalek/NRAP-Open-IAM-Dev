# -*- coding: utf-8 -*-
import logging
import sys
import os
import numpy as np

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import os.path
import pandas as pd
import matplotlib.pyplot as plt

try:
    import openiam.components.iam_base_classes as iam_bc
except ImportError as err:
    print('Unable to load NRAP-Open-IAM base classes module: {}'.format(err))

try:
    import openiam.components.models.wellbore.multisegmented_ai.multisegmented_wellbore_ai_ROM as mswrom
    from openiam.components.models.wellbore.multisegmented_ai.fluid_models import FluidModels
    from openiam.components.models.wellbore.multisegmented_ai.caprock_segment_model import DLCaprockSegmentModels
    from openiam.components.models.wellbore.multisegmented_ai.aquifer_segment_model import DLAquiferSegmentModels
except ImportError:
    print('\nERROR: Unable to load ROM for AI Multisegmented Wellbore component\n')
    sys.exit()

AQUIFER_FILE_EDX_LINK = 'http://edx.netl.doe.gov/dataset/42fb2930-7bf5-48d9-821b-d801c6f29f36/resource/8b786e28-0e52-40d4-84b5-23c25dc2d338/download'
SHALE_FILE_EDX_LINK = 'http://edx.netl.doe.gov/dataset/42fb2930-7bf5-48d9-821b-d801c6f29f36/resource/42882a99-a5b0-4d68-8f4c-cad4812358c4/download'

MODEL_DATA_LOCATION = os.path.join(iam_bc.IAM_DIR, 'src', 'openiam', 'components', 
                                   'models', 'wellbore', 'multisegmented_ai')

MODEL_LOADING_ERR_MSG = ''.join([
    'While setting up the MultisegmentedWellboreAI component {}, there was an ', 
    'error loading the model files required. A zipped file must be downloaded ', 
    'from the following link: {}. Once downloaded, the zipped folder must be ', 
    'extracted and placed in the directory {}. Without the required file, the ', 
    'MultisegmentedWellboreAI component cannot be used.'])


class MultisegmentedWellboreAI(iam_bc.ComponentModel):
    """
    The Multisegmented Wellbore AI component estimates the leakage rates of 
    brine and |CO2| along wells in the presence of overlying aquifers or thief 
    zones. The component is a deep learning model based on full physics numerical 
    simulations of wellbore leakage. The model is based on the work of Nordbotten 
    et al., :cite:`N2009`.
    
    The conceptualization of this component (:cite:`N2009`, :cite:`N2011a`) is 
    similar to the Multisegmented Wellbore component (:cite:'BaekEtAl2021b'), 
    but that component is an analytical model while this component is a machine 
    learning model. This component was developed by Baek et al., :cite:'BaekEtAl2023'. 
    That study compared the analytical and machine learning Multisegmented Wellbore 
    components, and found that the machine learning model has a significantly 
    improved performance. Note that the machine learning approach is only used 
    here if the parameter **useDLmodel** is set to 1 (``True``). Otherwise, the 
    results of this component are derived with analytical modeling.
    
    The files required for this component are relatively large (about 7.2 GB), so 
    they are not included in the download of NRAP-Open-IAM. The two files must be 
    downloaded from EDX (links shown in control file example 51a). The zipped 
    files must then be unzipped, with the resulting folders placed into this 
    directory:
    
        ``src/openiam/components/models/wellbore/multisegmented_ai``
    
    The directory above is relative to the location of NRAP-Open-IAM (the ``src`` 
    folder is located in the main directory for NRAP-Open-IAM). One of the two 
    files (about 2 GB) is always required for this component, as it is used to 
    calculate fluid properties. The second file (about 5.2 GB) is only required 
    if the parameter **useDLmodel** is set to 1. That parameter controls whether 
    the component 
    
    In the NRAP-Open-IAM control file interface, the type name for the Multisegmented 
    Wellbore AI component is ``MultisegmentedWellboreAI``. Descriptions of the 
    component's parameters are provided below.
    
    * **logWellPerm** [|log10| |m^2|] (-101 to -12.3) - logarithm of well (vertical) 
      permeability passing through shale (default: -13). Logarithm of well 
      permeability passing through shale 1, for example, can be defined by 
      **logWell1Perm**. Permeability of the well along any shale layers not 
      assigned a value will be assigned a default value.

    * **logAquPerm** [|log10| |m^2|] (-16 to -12.3) - logarithm of well (vertical) 
      permeability passing through aquifer (default: -13). Logarithm of well 
      permeability passing through aquifer 1, for example, can be defined by 
      **logAqu1Perm**.  Permeability of the well along any aquifer layers not 
      assigned a value will be assigned a default value.

    * **logAquiferHPerm** [|log10| |m^2|] (-14 to -12) - logarithm of aquifer 
      horizontal permeability (default: -13). Logarithm of aquifer 1  horizontal 
      permeability, for example, can be defined by **logAquifer1HPerm**. 
      Aquifer horizontal permeability not defined by user will be assigned a 
      default value.
      
    * **logAquiferVPerm** [|log10| |m^2|] (-15 to -12.3) - logarithm of aquifer 
      vertical permeability  (default: -13). Logarithm of aquifer 1 vertical 
      permeability, for example, can be defined by **logAquifer1VPerm**. 
      Aquifer vertical permeability not defined by user will be assigned a
      default value.

    * **brineDensity** [|kg/m^3|] (900 to 1400) - density of brine phase 
      (default: 1000). Not in use but necessary for a placeholder.

    * **CO2Density** [|kg/m^3|] (20 to 1000) - density of |CO2| phase 
      (default: 479). Not in use but necessary for a placeholder.

    * **brineViscosity** [|Pa*s|] (1.0e-4 to 5.0e-3) - viscosity of brine phase 
      (default: 2.535e-3). Not in use but necessary for a placeholder.

    * **CO2Viscosity** [|Pa*s|] (1.0e-6 to 1.0e-4)  - viscosity of |CO2| phase 
      (default: 3.95e-5). Not in use but necessary for a placeholder.

    * **aquBrineResSaturation** [-] (0 to 0.99) - residual saturation of brine phase 
      in each aquifer (default: 0.0). For example, the residual brine saturation 
      of aquifer1 can be defined by **aqu1BrineResSaturation**; otherwise, aquifer 
      layers for which the residual brine saturation is not defined will be 
      assigned a default value.

    * **compressibility** [|Pa^-1|] (1.0e-13 to 1.0e-8) - total compressibility 
      (default: 5.1e-11). Necessary for analytical calculation in the initial 
      timesteps.
    
    * **wellRadius** [|m|] (0.03 to 0.042) - radius of leaking well (default: 0.03).

    * **numberOfShaleLayers** [-] (3 to 30) - number of shale layers in the 
      system (default: 3); *linked to Stratigraphy*. The shale units must be 
      separated by aquifers.

    * **reservoirThickness** [|m|] (1 to 1600) - thickness of reservoir (default: 30); 
      *linked to Stratigraphy*.

    * **shaleThickness** [|m|] (50 to 2150) - thickness of shale layers (default: 
      250); *linked to Stratigraphy*. Thickness of shale layer 1, for example, 
      can be defined by **shale1Thickness**; otherwise, shale layers for which 
      the thickness is not defined will be assigned a default thickness.

    * **aquiferThickness** [|m|] (20 to 250) - thickness of aquifers (default: 100); 
      *linked to Stratigraphy*. Thickness of aquifer 1, for example, can be defined 
      by **aquifer1Thickness**; otherwise, aquifers for which the thickness 
      is not defined will be assigned a default thickness.

    * **aquiferPorosity** [-] (0.1 to 0.5) - porosity of aquifers (default: 0.2).
      Porosity of aquifer 1, for example, can be defined by **aquifer1porosity**; 
      otherwise, aquifers for which the porosity is not defined will be assigned 
      a default thickness.
      
    * **salinity** [|kg/kg|] (1e-4 to 0.015) - initial salinity of aquifers 
      (default: 0.01). The current version applies the same salinity to all aquifers. 

    * **tempGrad** [|C/km|] (24 to 32) - initial salinity of aquifers (default: 25). 
      The current version applies the same salinity to all aquifers. 
      
    * **datumPressure** [|Pa|] (80,000 to 300,000) - pressure at the top of the 
      system (default: 101,325); *linked to Stratigraphy*.
      
    * **useDLmodel** [-] (0 or 1) - option to use the machine learning models (1) 
      or analytical models (0) for the results. The machine learning models can 
      be computationally demanding, and using multiple locations for this component 
      increases that demand. The only acceptable values for this parameter are 
      1 or 0 (``True`` or ``False``).
    
    The possible outputs from the AI Multisegmented Wellbore component are 
    leakage rates of |CO2| and brine to each of the aquifers in the system as 
    well as the atmosphere. The names of the observations are of the form:

    * **CO2_aquifer1**, **CO2_aquifer2**,..., **CO2_atm** [|kg/s|] -
      |CO2| leakage rates

    * **brine_aquifer1**, **brine_aquifer2**,..., **brine_atm** [|kg/s|] -
      brine leakage rates

    * **mass_CO2_aquifer1**, **mass_CO2_aquifer2**,..., **mass_CO2_aquiferN** [|kg|]
      - mass of the |CO2| leaked into the aquifer.

    * **mass_brine_aquifer1**, **mass_brine_aquifer2**,..., **mass_brine_aquiferN** [|kg|]
      - mass of the brine leaked into the aquifer.
      
    For control file examples using the Multisegmented Wellbore AI component, see 
    *ControlFile_ex31g*, *ControlFile_ex48b*, *ControlFile_ex48c*, *ControlFile_ex51a* 
    to *ControlFile_ex51c*, and *ControlFile_ex58e*. For script examples, see 
    *iam_sys_strata_reservoir_mswellai_genericaquifer.py* and 
    *iam_sys_lutreservoir_mswellai.py*.
    """
    def __init__(self, name, parent, useDLmodel=False):
        """
        Constructor method of MultisegmentedWellboreAI class

        :param name: name of component model
        :type name: str

        :param parent: the SystemModel object that the component model
            belongs to
        :type parent: SystemModel object

        :returns: MultisegmentedWellboreAI class object
        """
        # Set up keyword arguments of the 'model' method provided by the system model
        model_kwargs = {'time_point': 365.25, 'time_step': 365.25}   # default value of 365.25 days

        super().__init__(name, parent, model=self.simulation_model, 
                         model_kwargs=model_kwargs)

        # Add type attribute
        self.class_type = 'MultisegmentedWellboreAI'

        # Set default parameters of the component model
        self.add_default_par('numberOfShaleLayers', value=3)
        self.add_default_par('shaleThickness', value=250.0)
        self.add_default_par('aquiferThickness', value=100.0)
        self.add_default_par('reservoirThickness', value=30.0)
        self.add_default_par('logWellPerm', value=-13.0)
        self.add_default_par('logAquPerm', value=-13.0)
        self.add_default_par('datumPressure', value=101325.0)
        self.add_default_par('brineDensity', value=1000.0) # not in use
        self.add_default_par('CO2Density', value=479.0)  # not in use
        self.add_default_par('brineViscosity', value=2.535e-3)  # not in use
        self.add_default_par('CO2Viscosity', value=3.95e-5)  # not in use
        self.add_default_par('aquBrineResSaturation', value=0.0)
        self.add_default_par('compressibility', value=5.1e-11)
        self.add_default_par('wellRadius', value=0.03)
        self.add_default_par('salinity', value=0.01) # initial salinity
        self.add_default_par('tempGrad', value=25.) # C/km
        ## new
        self.add_default_par('useDLmodel', useDLmodel)
        self.add_default_par('aquiferPorosity', value=0.2) ## new
        self.add_default_par('logAquiferHPerm', value=-13.0) ## new
        self.add_default_par('logAquiferVPerm', value=-13.0) ## new

        # Define dictionary of boundaries
        self.pars_bounds = dict()
        self.pars_bounds['numberOfShaleLayers'] = [3, 30]
        self.pars_bounds['shaleThickness'] = [50, 2150.0]
        self.pars_bounds['aquiferThickness'] = [20, 250.0]
        self.pars_bounds['reservoirThickness'] = [1.0, 1600.0]
        self.pars_bounds['logWellPerm'] = [-101.0, -12.3]
        self.pars_bounds['logAquPerm'] = [-16.0, -12.3]
        self.pars_bounds['datumPressure'] = [8.0e+4, 3.0e+5]
        self.pars_bounds['brineDensity'] = [900.0, 1400.0]  # not in use
        self.pars_bounds['CO2Density'] = [20.0, 1000.0]  # not in use
        self.pars_bounds['brineViscosity'] = [1.0e-4, 5.0e-3]  # not in use
        self.pars_bounds['CO2Viscosity'] = [1.0e-6, 1.5e-4]  # not in use
        self.pars_bounds['aquBrineResSaturation'] = [0.0, 0.99]
        self.pars_bounds['compressibility'] = [1.0e-13, 1.0e-8]
        self.pars_bounds['wellRadius'] = [0.03, 0.042]
        self.pars_bounds['salinity'] = [1e-4, 0.015]
        self.pars_bounds['tempGrad'] = [24., 32.]
        
        self.pars_bounds['useDLmodel'] = [0, 1]
        self.pars_bounds['aquiferPorosity'] = [0.1, 0.5]
        self.pars_bounds['logAquiferHPerm'] = [-14, -12]
        self.pars_bounds['logAquiferVPerm'] = [-15, -12.3]
        
        # Define dictionary of parameters that are restricted to specific values
        self.acceptable_vals = dict()
        self.acceptable_vals['useDLmodel'] = [0, 1]

        # By default, the smallest number of aquifers the system can have is 2,
        # so we add two accumulators by default. Extra will be added as needed,
        # once the system knows how many there are
        self.num_accumulators = 2
        for i in range(2):
            self.add_accumulator('mass_CO2_aquifer'+str(i+1), sim=0.0)
            self.add_accumulator('mass_brine_aquifer'+str(i+1), sim=0.0)
            self.add_accumulator('CO2_saturation_aquifer'+str(i+1), sim=0.0)
        
            # Attribute num_accumulators equals to the number of aquifers
            self.add_accumulator('pressure_aquifer'+str(i+1), sim=0.0)
            self.add_accumulator('pressure_shaleBot'+str(i+1), sim=0.0)
            self.add_accumulator('pressure_shaleTop'+str(i+1), sim=0.0)
        
        self.useDLmodel = useDLmodel
        
        componentPath0 = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        
        # The extracted data might repeat the folder name. Check which approach 
        # is used and store the correct pathway.
        componentPath = os.path.join(componentPath0, 'components', 'models', 'wellbore', 
                                     'multisegmented_ai', 'shale_feb2024_models')
        componentPath_v2 = os.path.join(componentPath0, 'components', 'models', 'wellbore', 
                                     'multisegmented_ai', 'shale_feb2024_models', 
                                     'shale_feb2024_models')
        
        componentPath2 = os.path.join(componentPath0, 'components', 'models', 'wellbore', 
                                      'multisegmented_ai', 'aquifer_mar2024_models')
        componentPath2_v2 = os.path.join(componentPath0, 'components', 'models', 'wellbore', 
                                      'multisegmented_ai', 'aquifer_mar2024_models', 
                                      'aquifer_mar2024_models')
        
        # Caprock segment model files
        if os.path.exists(componentPath_v2):
            self.componentPath = componentPath_v2
        else:
            self.componentPath = componentPath
            
            # These files are required for the FluidModels object, so check if 
            # they are present. Log an error if they are not.
            if not os.path.exists(self.componentPath):
                err_msg = MODEL_LOADING_ERR_MSG.format(
                    self.name, SHALE_FILE_EDX_LINK, MODEL_DATA_LOCATION)
                
                logging.error(err_msg)
        
        # Aquifer segment model files
        if os.path.exists(componentPath2_v2):
            self.componentPath2 = componentPath2_v2
        else:
            self.componentPath2 = componentPath2
            
            # Only check these files if useDLmodel is True, they are not required 
            # for the FluidModels object.
            if not os.path.exists(self.componentPath2) and useDLmodel:
                err_msg = MODEL_LOADING_ERR_MSG.format(
                    self.name, AQUIFER_FILE_EDX_LINK, MODEL_DATA_LOCATION)
                
                logging.error(err_msg)
        
        # These are updated in simulation_model() if useDLmodel is True
        self.DLModels1 = None
        self.DLModels2 = None
        
        # This is always updated in simulation_model()
        self.FluidModels = None
        
        debug_msg = 'MultisegmentedWellboreAI component created with name {}'.format(name)
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
            
            warn_msg = ''.join([
                'Parameter {} of MultisegmentedWellboreAI component {} ',
                'is out of boundaries.']).format(key, self.name)
            
            if key.startswith('shale') and key.endswith('Thickness'):
                if (val < self.pars_bounds['shaleThickness'][0]) or (
                        val > self.pars_bounds['shaleThickness'][1]):
                    logging.warning(warn_msg)
                continue
            
            if key.startswith('aquifer') and key.endswith('Thickness'):
                if (val < self.pars_bounds['aquiferThickness'][0]) or (
                        val > self.pars_bounds['aquiferThickness'][1]):
                    logging.warning(warn_msg)
                continue
            
            if key.startswith('logWell') and key.endswith('Perm'):
                if (val < self.pars_bounds['logWellPerm'][0]) or (
                        val > self.pars_bounds['logWellPerm'][1]):
                    logging.warning(warn_msg)
                continue
            
            if key.startswith('logAqu') and key.endswith('Perm'):
                if ((val < self.pars_bounds['logAquPerm'][0]) or
                        (val > self.pars_bounds['logAquPerm'][1])):
                    logging.warning(warn_msg)
                continue
            
            if key.startswith('aqu') and key.endswith('BrineResSaturation'):
                if ((val < self.pars_bounds['aquBrineResSaturation'][0]) or
                        (val > self.pars_bounds['aquBrineResSaturation'][1])):
                    logging.warning(warn_msg)
                continue
            
            if key.startswith('aquifer') and key.endswith('Porosity'):
                if ((val < self.pars_bounds['aquiferPorosity'][0]) or
                        (val > self.pars_bounds['aquiferPorosity'][1])):
                    logging.warning(warn_msg)
                continue
            
            if key.startswith('logAquifer') and key.endswith('HPerm'):
               if ((val < self.pars_bounds['logAquiferHPerm'][0]) or
                       (val > self.pars_bounds['logAquiferHPerm'][1])):
                   logging.warning(warn_msg)
               continue
           
            if key.startswith('logAquifer') and key.endswith('VPerm'):
                if ((val < self.pars_bounds['logAquiferVPerm'][0]) or
                        (val > self.pars_bounds['logAquiferVPerm'][1])):
                    logging.warning(warn_msg)
                continue
            
            if key in self.pars_bounds:
                if ((val < self.pars_bounds[key][0]) or
                        (val > self.pars_bounds[key][1])):
                    logging.warning(warn_msg)
    
    def simulation_model(self, p, time_point=365.25, time_step=365.25, 
                         pressure=0.0, CO2saturation=0.0, salinity=0.0):
        """
        Return |CO2| and brine leakage rates corresponding to the provided input.
        Note that the names of parameters contained in the input dictionary
        coincide with those defined in this module's docstring.

        :param p: input parameters of multisegmented wellbore model
        :type p: dict

        :param pressure: pressure at the bottom of leaking well, in Pa;
            by default, its value is 0.0
        :type pressure: float

        :param CO2saturation: saturation of |CO2| phase at the bottom
            of leaking well; by default, its value is 0.0
        :type CO2saturation: float
        
        :param salinity: mass fraction of salt in water (kg/kg) at the bottom
            of leaking well; by default, its value is 0.0
        :type salinity: float

        :param time_point: time point in days at which the leakage rates are
            to be calculated; by default, its value is 365.25 (1 year in days)
        :type time_point: float

        :param time_step: difference between the current and previous
            time points in days; by default, its value is 365.25 (1 year in days)
        :type time_point: float

        :returns: out - dictionary of observations of multisegmented wellbore
            model; keys:
            ['CO2_aquifer1','CO2_aquifer2',...,'CO2_atm',
            'brine_aquifer1','brine_aquifer2',...,'brine_atm',
            'mass_CO2_aquifer1','mass_CO2_aquifer2',...,
            'mass_brine_aquifer1','mass_brine_aquifer2',...]
        """
        # Obtain the default values of the parameters from dictionary of default parameters
        actual_p = {k: v.value for k, v in self.default_pars.items()}
        # Update default values of parameters with the provided ones
        actual_p.update(p)
        inputParameters = mswrom.Parameters()

        nSL = int(actual_p['numberOfShaleLayers'])
        inputParameters.numberOfShaleLayers = nSL

        # Add extra accumulators after the number of aquifers became known
        num_extra_accumulators = nSL-1 - self.num_accumulators
        
        if num_extra_accumulators != 0: # if component has no enough accumulators
            for i in range(num_extra_accumulators):
                self.add_accumulator('mass_CO2_aquifer'+str(i+3), sim=0.0)
                self.add_accumulator('mass_brine_aquifer'+str(i+3), sim=0.0)
                self.add_accumulator('CO2_saturation_aquifer'+str(i+3), sim=0.0)
                
                self.add_accumulator('pressure_aquifer'+str(i+3), sim=0.0)
                self.add_accumulator('pressure_shaleBot'+str(i+3), sim=0.0)
                self.add_accumulator('pressure_shaleTop'+str(i+3), sim=0.0)
            
            self.add_accumulator('dl_inputarray1', sim=0.0)
            self.add_accumulator('dl_inputarray1_minus_one', sim=0.0)
            self.add_accumulator('dl_inputarray1_minus_two', sim=0.0)
            
            self.add_accumulator('dl_inputarray2', sim=0.0)
            self.add_accumulator('dl_inputarray2_minus_one', sim=0.0)
            self.add_accumulator('dl_inputarray2_minus_two', sim=0.0)
            
            self.add_accumulator('n_timestep', sim=0.0)
                
            self.num_accumulators = nSL - 1
        
        accumulator_check = ('dl_inputarray1' in self.accumulators) and (
            'dl_inputarray1_minus_one' in self.accumulators) and (
                'dl_inputarray1_minus_two' in self.accumulators) and (
                    'dl_inputarray2' in self.accumulators) and (
                        'dl_inputarray2_minus_one' in self.accumulators) and (
                            'dl_inputarray2_minus_two' in self.accumulators) and (
                                'n_timestep' in self.accumulators)
        
        if not accumulator_check:
            self.add_accumulator('dl_inputarray1', sim=0.0)
            self.add_accumulator('dl_inputarray1_minus_one', sim=0.0)
            self.add_accumulator('dl_inputarray1_minus_two', sim=0.0)
            
            self.add_accumulator('dl_inputarray2', sim=0.0)
            self.add_accumulator('dl_inputarray2_minus_one', sim=0.0)
            self.add_accumulator('dl_inputarray2_minus_two', sim=0.0)
            
            self.add_accumulator('n_timestep', sim=0.0)

        # Check whether the initial state of the system is requested
        # Then if yes, return the state without proceeding further
        if time_point == 0.0:
            # Create dictionary of leakage rates
            out = dict()
            for i in range(nSL-1):
                out['CO2_aquifer'+str(i+1)] = 0.0
                out['brine_aquifer'+str(i+1)] = 0.0
                out['mass_CO2_aquifer'+str(i+1)] = 0.0
                out['mass_brine_aquifer'+str(i+1)] = 0.0
                self.accumulators['mass_CO2_aquifer'+str(i+1)].sim = 0.0
                self.accumulators['mass_brine_aquifer'+str(i+1)].sim = 0.0
                self.accumulators['CO2_saturation_aquifer'+str(i+1)].sim = 0.0
                self.accumulators['pressure_aquifer'+str(i+1)].sim = 0.0
                self.accumulators['pressure_shaleBot'+str(i+1)].sim = 0.0
                self.accumulators['pressure_shaleTop'+str(i+1)].sim = 0.0
            out['CO2_atm'] = 0.0
            out['brine_atm'] = 0.0
            self.prevCO2Mass = 0.0
            self.prevBrineMass = 0.0
            
            self.accumulators['dl_inputarray1'].sim = dict()
            self.accumulators['dl_inputarray1_minus_one'].sim = dict()
            self.accumulators['dl_inputarray1_minus_two'].sim = dict()
            
            self.accumulators['dl_inputarray2'].sim = dict()
            self.accumulators['dl_inputarray2_minus_one'].sim = dict()
            self.accumulators['dl_inputarray2_minus_two'].sim = dict()
            self.accumulators['n_timestep'].sim = 0
            
            return out

        inputParameters.shaleThickness = actual_p['shaleThickness'] * np.ones(nSL)
        inputParameters.shalePermeability = (
            10 ** actual_p['logWellPerm'] * np.ones(nSL))
        inputParameters.aquiferThickness = actual_p['aquiferThickness'] * np.ones((nSL-1))
        inputParameters.aquiferPermeability = (
            10 ** actual_p['logAquPerm'] * np.ones((nSL-1)))
        inputParameters.aquBrineResSaturation = np.ones((nSL))

        inputParameters.aquiferPorosity = np.ones((nSL-1))
        inputParameters.logAquiferHPerm = np.ones((nSL-1))
        inputParameters.logAquiferVPerm = np.ones((nSL-1))
        
        # Set up shale, aquifer and reservoir thickness parameters
        for i in range(nSL):
            nm = 'shale{}Thickness'.format(i+1)
            if nm in p:
                inputParameters.shaleThickness[i] = p[nm]
        for i in range(nSL-1):
            nm = 'aquifer{}Thickness'.format(i+1)
            if nm in p:
                inputParameters.aquiferThickness[i] = p[nm]
        inputParameters.reservoirThickness = actual_p['reservoirThickness']

        # Set up permeability of the well segment along shale
        for i in range(nSL):
            nm = 'logWell{}Perm'.format(i+1)
            if nm in p:
                inputParameters.shalePermeability[i] = 10**p[nm]

        # Set up permeability of the wellsegment along aquifer
        for i in range(nSL-1):
            nm = 'logAqu{}Perm'.format(i+1)
            if nm in p:
                inputParameters.aquiferPermeability[i] = 10**p[nm]

        # Set up residual saturation of aquifers
        for i in range(nSL):
            nm = 'aqu{}BrineResSaturation'.format(i) # aqu0BrineResSaturation is reservoir
            if nm in p:
                inputParameters.aquBrineResSaturation[i] = p[nm]
        inputParameters.aquBrineResSaturation[0] = 0.0 # 0 is reservoir

        for i in range(nSL-1):
            nm = 'aquifer{}Porosity'.format(i+1)
            if nm in p:
                inputParameters.aquiferPorosity[i] = p[nm]
        for i in range(nSL-1):
            nm = 'logAquifer{}HPerm'.format(i+1)
            if nm in p:
                inputParameters.logAquiferHPerm[i] = p[nm]
        for i in range(nSL-1):
            nm = 'logAquifer{}HPerm'.format(i+1)
            if nm in p:
                inputParameters.logAquiferVPerm[i] = p[nm]        

        # Set up land surface pressure
        inputParameters.datumPressure = actual_p['datumPressure']
        inputParameters.tempGrad = actual_p['tempGrad']

        # Set up brine and CO2 density
        inputParameters.brineDensity = actual_p['brineDensity']
        inputParameters.CO2Density = actual_p['CO2Density']

        # Set up brine and CO2 viscosity
        inputParameters.brineViscosity = actual_p['brineViscosity']
        inputParameters.CO2Viscosity = actual_p['CO2Viscosity']

        # Set up residual saturation and compressibility
        inputParameters.compressibility = actual_p['compressibility']

        # Set up well radius parameter
        inputParameters.wellRadius = actual_p['wellRadius']
        inputParameters.flowArea = np.pi * (inputParameters.wellRadius ** 2)

        # Get parameters from keyword arguments of the 'model' method
        inputParameters.timeSpan = time_step         # in days
        inputParameters.timePoint = time_point       # in days
        inputParameters.pressure = pressure

        inputParameters.prevCO2Mass = np.zeros(nSL-1)
        inputParameters.prevBrineMass = np.zeros(nSL-1)
        inputParameters.CO2saturation = np.zeros(nSL)
        inputParameters.CO2saturation[0] = CO2saturation
        
        inputParameters.pressure_prev = np.zeros(nSL)
        inputParameters.pressure_prev[0] = pressure
        inputParameters.pressureShaleTop_prev = np.zeros(nSL)
        inputParameters.pressureShaleBot_prev = np.zeros(nSL)
        
        for i in range(nSL-1):
            inputParameters.prevCO2Mass[i] = (
                self.accumulators['mass_CO2_aquifer'+str(i+1)].sim)
            inputParameters.prevBrineMass[i] = (
                self.accumulators['mass_brine_aquifer'+str(i+1)].sim)
            inputParameters.CO2saturation[i+1] = (
                self.accumulators['CO2_saturation_aquifer{}'.format(i+1)].sim)
            inputParameters.pressure_prev[i+1] = (
                self.accumulators['pressure_aquifer{}'.format(i+1)].sim)           
            inputParameters.pressureShaleTop_prev[i] = (
                self.accumulators['pressure_shaleTop{}'.format(i+1)].sim)                   
            inputParameters.pressureShaleBot_prev[i] = (
                self.accumulators['pressure_shaleBot{}'.format(i+1)].sim)
        
        inputParameters.dl1input_minthree = self.accumulators['dl_inputarray1_minus_two'].sim
        inputParameters.dl1input_mintwo = self.accumulators['dl_inputarray1_minus_one'].sim
        inputParameters.dl1input_minone = self.accumulators['dl_inputarray1'].sim
        
        inputParameters.dl2input_minthree = self.accumulators['dl_inputarray2_minus_two'].sim
        inputParameters.dl2input_mintwo = self.accumulators['dl_inputarray2_minus_one'].sim
        inputParameters.dl2input_minone = self.accumulators['dl_inputarray2'].sim
        inputParameters.timeStep = self.accumulators['n_timestep'].sim
        
        inputParameters.pressureShaleBot_prev[0] = pressure     
        
        inputParameters.SaltMassFrac = np.zeros(nSL)
        inputParameters.SaltMassFrac[:] = actual_p['salinity'] # use initial salinity
        inputParameters.SaltMassFrac[0] = salinity # dynamic value at the bottom leaky well
        inputParameters.useDLmodel = actual_p['useDLmodel']
        
        if not inputParameters.useDLmodel in self.acceptable_vals['useDLmodel']:
            warn_msg = ''.join([
                'The {} parameter was given as {}. The only acceptable ', 
                'values for this parameter are {}, however. The default ', 
                'value of {} will be used.']).format(
                    'useDLmodel', inputParameters.useDLmodel, 
                    self.acceptable_vals['useDLmodel'], self.default_pars['useDLmodel'])
            
            logging.warning(warn_msg)
            
            inputParameters.useDLmodel = self.default_pars['useDLmodel']
        
        # Load the models, if enabled and they have not been loaded yet
        if (inputParameters.useDLmodel or self.useDLmodel) and (
                self.DLModels1 is None or self.DLModels1 is None):
            # May run out of memory while attempting to load these models
            self.DLModels1 = DLCaprockSegmentModels(self.componentPath)
            self.DLModels2 = DLAquiferSegmentModels(self.componentPath2)
        
        # This object should not be stored in the component during the __init__() 
        # function. Doing so causes errors whenever pickle is used ("TypeError: 
        # cannot pickle 'module' object"). During a stochastic simulation, pickle 
        # is used during the multiprocessing steps in matk. 
        if self.FluidModels is None:
            self.FluidModels = FluidModels(self.componentPath)
        
        # Create solution object
        sol = mswrom.Solution(inputParameters, self.FluidModels, 
                              self.DLModels1, self.DLModels2)
        
        # Find solution corresponding to the inputParameters
        sol.find()

        # Create dictionary of leakage rates
        out = dict()
        for i in range(nSL-1):
            out['CO2_aquifer'+str(i+1)] = sol.CO2LeakageRates[i]
            out['brine_aquifer'+str(i+1)] = sol.brineLeakageRates[i]
            out['mass_CO2_aquifer{}'.format(i+1)] = sol.CO2Mass[i]            
            out['mass_brine_aquifer{}'.format(i+1)] = sol.brineMass[i]      
            
            # Keep mass in cubic meters in accumulators
            self.accumulators['mass_CO2_aquifer{}'.format(i+1)].sim = sol.CO2Mass[i]
            self.accumulators['mass_brine_aquifer{}'.format(i+1)].sim = sol.brineMass[i]
            self.accumulators['CO2_saturation_aquifer{}'.format(i+1)].sim = sol.CO2SaturationAq[i]
            
            self.accumulators[
                'pressure_aquifer{}'.format(i+1)].sim = sol.pressure_prev[i]                        
            self.accumulators[
                'pressure_shaleBot{}'.format(i+1)].sim = sol.pressureShaleBot_prev[i] # saving for only effective shale (excluding top shale)
            self.accumulators[
                'pressure_shaleTop{}'.format(i+1)].sim = sol.pressureShaleTop_prev[i] # saving for only effective shale (excluding top shale)
         
        self.accumulators['dl_inputarray1_minus_two'].sim = self.accumulators[
            'dl_inputarray1_minus_one'].sim
        self.accumulators['dl_inputarray1_minus_one'].sim = self.accumulators[
            'dl_inputarray1'].sim 
        self.accumulators['dl_inputarray1'].sim = sol.inputarray1
        
        self.accumulators['dl_inputarray2_minus_two'].sim = self.accumulators[
            'dl_inputarray2_minus_one'].sim
        self.accumulators['dl_inputarray2_minus_one'].sim = self.accumulators[
            'dl_inputarray2'].sim 
        self.accumulators['dl_inputarray2'].sim = sol.inputarray2
        self.accumulators['n_timestep'].sim = sol.timeStep
         
        out['CO2_atm'] = sol.CO2LeakageRates[inputParameters.numberOfShaleLayers-1]
        out['brine_atm'] = sol.brineLeakageRates[inputParameters.numberOfShaleLayers-1]

        return out

    def reset(self):
        return

    # Attributes for system connections
    system_inputs = ['pressure',
                     'CO2saturation']
    system_params = ['numberOfShaleLayers',
                     'shale1Thickness',
                     'shale2Thickness',
                     'shale3Thickness',
                     'aquifer1Thickness',
                     'aquifer2Thickness',
                     'reservoirThickness',
                     'datumPressure']

def read_data(filename):
    """
    Routine used for reading data files.

    Read data from file and create numpy.array if file exists.
    """
    # Check whether the file with given name exists
    if os.path.isfile(filename):
        data = np.genfromtxt(filename)
        return data

    return None


if __name__ == "__main__":
    try:
        from openiam.components.analytical_reservoir_component import AnalyticalReservoir
    except ImportError as err:
        print('Unable to load IAM class module: {}'.format(err))
    
    __spec__ = None

    logging.basicConfig(level=logging.WARNING)
    
    # Define keyword arguments of the system model
    isUQAnalysisOn = 0 # 0: one forward run, 1: stochastic runs
    isPlottingOn = 1
    to_save_png = False

    num_years = 50
    time_array = 365.25*np.arange(0, num_years+1)

    sm_model_kwargs = {'time_array': time_array} # time is given in days

    # Create system model
    sm = iam_bc.SystemModel(model_kwargs=sm_model_kwargs)

    # Add reservoir component
    res = sm.add_component_model_object(AnalyticalReservoir(
        name='res', parent=sm, injX=100., injY=0., locX=0., locY=0.))

    # Add parameters of reservoir component model
    res.add_par('injRate', min=1e-5, max=1e+2, value=0.0185, vary=False)
    res.add_par('reservoirRadius', min=500.0, max=1000000, value=100000.0, vary=False)
    res.add_par('reservoirThickness', min=10.0, max=150., value=30.0, vary=False)
    res.add_par('logResPerm', min=-16.0, max=-9.0, value=-13.69897, vary=False)
    res.add_par('reservoirPorosity', min=0.01, max=0.5, value=0.15, vary=False)

    res.add_par('shaleThickness', value=970.0, vary=False)
    res.add_par('aquiferThickness', value=30.0, vary=False)

    res.add_par('brineDensity', min=900.0, max=1300., value=1045.0, vary=False)
    res.add_par('CO2Density', min=200.0, max=900., value=479.0, vary=False)
    res.add_par('brineViscosity', min=1.0e-5, max=1.0e-3, value=2.535e-4, vary=False)
    res.add_par('CO2Viscosity', min=1.0e-5, max=1.0e-4, value=3.95e-5, vary=False)
    res.add_par('brineResSaturation', min=0.0, max=0.99, value=0.0, vary=False)
    res.add_par('brineCompressibility', min=1e-9, max=1e-13, value=1.e-11, vary=False)

    res.add_par('datumPressure', value=101325.0, vary=False)

    # Add observations of reservoir component model
    res.add_obs_to_be_linked('pressure')
    res.add_obs_to_be_linked('CO2saturation')

    res.add_obs('pressure')
    res.add_obs('CO2saturation')
    res.add_obs('mass_CO2_reservoir')

    # Add multisegmented wellbore component
    ms = sm.add_component_model_object(MultisegmentedWellboreAI(
        name='ms', parent=sm, useDLmodel=True))
    
    ms.add_par('wellRadius', min=0.03, max=0.042, value=0.035)

    ms.add_par('numberOfShaleLayers', value=5, vary=False)
    ms.add_par('shale1Thickness', min=1.0, max=200, value=100.0)
    ms.add_par('shale2Thickness', min=1.0, max=200, value=100.0)
    ms.add_par('shale3Thickness', min=1.0, max=200, value=100.0)
    ms.add_par('shale4Thickness', min=1.0, max=200, value=400.0)
    ms.add_par('shale5Thickness', min=1.0, max=2500, value=2150.0)

    ms.add_par('logWell1Perm', min=-16, max=-9, value=-13)
    ms.add_par('logWell2Perm', min=-16, max=-9, value=-13)
    ms.add_par('logWell3Perm', min=-16, max=-9, value=-13)
    ms.add_par('logWell4Perm', min=-16, max=-9, value=-13)
    ms.add_par('logWell5Perm', min=-101, max=-9, value=-100)

    ms.add_par('aquifer1Thickness', min=1.0, max=200., value=30.0)
    ms.add_par('aquifer2Thickness', min=1.0, max=200., value=30.0)
    ms.add_par('aquifer3Thickness', min=1.0, max=200., value=30.0)
    ms.add_par('aquifer4Thickness', min=1.0, max=200., value=30.0)

    ms.add_par('logAqu1Perm', min=-16, max=-9, value=-13)
    ms.add_par('logAqu2Perm', min=-16, max=-9, value=-13)
    ms.add_par('logAqu3Perm', min=-16, max=-9, value=-13)
    ms.add_par('logAqu4Perm', min=-16, max=-9, value=-13)

    ms.add_par('aqu1BrineResSaturation', min=0.0, max=0.99, value=0.18)
    ms.add_par('aqu2BrineResSaturation', min=0.0, max=0.99, value=0.4)
    ms.add_par('aqu3BrineResSaturation', min=0.0, max=0.99, value=0.5)
    ms.add_par('aqu4BrineResSaturation', min=0.0, max=0.99, value=0.0)

    ms.add_par('aquifer1Porosity', min=0.0, max=0.99, value=0.4)
    ms.add_par('aquifer2Porosity', min=0.0, max=0.99, value=0.3)
    ms.add_par('aquifer3Porosity', min=0.0, max=0.99, value=0.2)
    ms.add_par('aquifer4Porosity', min=0.0, max=0.99, value=0.1)
    
    ms.add_par('logAquifer1HPerm', min=-16, max=-9, value=-12.5)
    ms.add_par('logAquifer2HPerm', min=-16, max=-9, value=-13.0)
    ms.add_par('logAquifer3HPerm', min=-16, max=-9, value=-13.5)
    ms.add_par('logAquifer4HPerm', min=-16, max=-9, value=-14.0)
    
    ms.add_par('logAquifer1VPerm', min=-16, max=-9, value=-14)
    ms.add_par('logAquifer2VPerm', min=-16, max=-9, value=-14)
    ms.add_par('logAquifer3VPerm', min=-16, max=-9, value=-14)
    ms.add_par('logAquifer4VPerm', min=-16, max=-9, value=-14)

    # Add linked parameters: common to both components
    ms.add_par_linked_to_par('reservoirThickness',
                             res.default_pars['reservoirThickness'])
    ms.add_par_linked_to_par('aqu0BrineResSaturation',
                             res.default_pars['brineResSaturation']) # reservoir
    ms.add_par_linked_to_par('datumPressure',
                             res.default_pars['datumPressure'])

    # Add keyword arguments linked to the output provided by reservoir model
    ms.add_kwarg_linked_to_obs('pressure', res.linkobs['pressure'])
    ms.add_kwarg_linked_to_obs('CO2saturation', res.linkobs['CO2saturation'])

    # Add observations of multisegmented wellbore component model
    # When add_obs method is called, system model automatically
    # (if no other options of the method is used) adds a list of observations
    # with names ms.obsnm_0, ms.obsnm_1,.... The final index is determined by the
    # number of time points in system model time_array attribute.
    # For more information, see the docstring of add_obs of ComponentModel class.
    ms.add_obs('CO2_aquifer1')
    ms.add_obs('CO2_aquifer2')
    ms.add_obs('CO2_aquifer3')
    ms.add_obs('CO2_aquifer4')
    ms.add_obs('CO2_atm')

    ms.add_obs('mass_CO2_aquifer1')
    ms.add_obs('mass_CO2_aquifer2')
    ms.add_obs('mass_CO2_aquifer3')
    ms.add_obs('mass_CO2_aquifer4')

    ms.add_obs('brine_aquifer1')
    ms.add_obs('brine_aquifer2')
    ms.add_obs('brine_aquifer3')
    ms.add_obs('brine_aquifer4')

    if isUQAnalysisOn == 0:
        # Run system model using current values of its parameters
        sm.forward()  # system model is run deterministically

        print('------------------------------------------------------------------')
        print('                  Forward method illustration ')
        print('------------------------------------------------------------------')
        # Since the observations at the particular time points are different variables,
        # method collect_observations_as_time_series creates lists of
        # values of observations belonging to a given component (e.g. cw) and having the same
        # common name (e.g. 'CO2_aquifer1', etc) but differing in indices.
        # More details are given in the docstring and documentation to the method
        # collect_observations_as_time_series of SystemModel class.

        pressure = sm.collect_observations_as_time_series(res, 'pressure')
        CO2saturation = sm.collect_observations_as_time_series(res, 'CO2saturation')
        mass_CO2_reservoir = sm.collect_observations_as_time_series(res, 'mass_CO2_reservoir')
        
        CO2_aquifer1 = sm.collect_observations_as_time_series(ms, 'CO2_aquifer1')
        CO2_aquifer2 = sm.collect_observations_as_time_series(ms, 'CO2_aquifer2')
        CO2_aquifer3 = sm.collect_observations_as_time_series(ms, 'CO2_aquifer3')
        CO2_aquifer4 = sm.collect_observations_as_time_series(ms, 'CO2_aquifer4')
        CO2_atm = sm.collect_observations_as_time_series(ms, 'CO2_atm')

        mass_CO2_aquifer1 = sm.collect_observations_as_time_series(ms, 'mass_CO2_aquifer1')
        mass_CO2_aquifer2 = sm.collect_observations_as_time_series(ms, 'mass_CO2_aquifer2')
        mass_CO2_aquifer3 = sm.collect_observations_as_time_series(ms, 'mass_CO2_aquifer3')
        mass_CO2_aquifer4 = sm.collect_observations_as_time_series(ms, 'mass_CO2_aquifer4')

        brine_aquifer1 = sm.collect_observations_as_time_series(ms, 'brine_aquifer1')
        brine_aquifer2 = sm.collect_observations_as_time_series(ms, 'brine_aquifer2')
        brine_aquifer3 = sm.collect_observations_as_time_series(ms, 'brine_aquifer3')
        brine_aquifer4 = sm.collect_observations_as_time_series(ms, 'brine_aquifer4')

        num_samples = 1
        print('Forward run is done. ')

        if isPlottingOn:
           # Plot results
           label_size = 13
           font_size = 16
           ticks_size = 12
           line_width = 1
           fig = plt.figure(figsize=(13, 12))

           ax = fig.add_subplot(331)
           plt.plot(time_array/365.25, CO2_aquifer1,
                    color='steelblue', linewidth=line_width)
           plt.xlabel('Time, years', fontsize=label_size)
           plt.ylabel('Leakage rates, kg/s', fontsize=label_size)
           plt.title(r'Leakage of CO$_2$: aquifer 1', fontsize=font_size)
           plt.tick_params(labelsize=ticks_size)
           plt.xlim([0, 50])
           ax.get_yaxis().set_label_coords(-0.13, 0.5)

           ax = fig.add_subplot(332)
           plt.plot(time_array/365.25, CO2_aquifer2,
                    color='steelblue', linewidth=line_width)
           plt.xlabel('Time, years', fontsize=label_size)
           plt.ylabel('Leakage rates, kg/s', fontsize=label_size)
           plt.title(r'Leakage of CO$_2$: aquifer 2', fontsize=font_size)
           plt.tick_params(labelsize=ticks_size)
           plt.xlim([0, 50])
           ax.get_yaxis().set_label_coords(-0.13, 0.5)
           
           ax = fig.add_subplot(333)
           plt.plot(time_array/365.25, CO2_aquifer3,
                    color='steelblue', linewidth=line_width)
           plt.xlabel('Time, years', fontsize=label_size)
           plt.ylabel('Leakage rates, kg/s', fontsize=label_size)
           plt.title(r'Leakage of CO$_2$: aquifer 3', fontsize=font_size)
           plt.tick_params(labelsize=ticks_size)
           plt.xlim([0, 50])
           ax.get_yaxis().set_label_coords(-0.13, 0.5)

           ax = fig.add_subplot(334)
           plt.plot(time_array/365.25, brine_aquifer1,
                    color='rosybrown', linewidth=line_width)
           plt.xlabel('Time, years', fontsize=label_size)
           plt.ylabel('Leakage rates, kg/s', fontsize=label_size)
           plt.title(r'Leakage of brine: aquifer 1', fontsize=font_size)
           plt.tight_layout()
           plt.tick_params(labelsize=ticks_size)
           plt.xlim([0, 50])
           ax.get_yaxis().set_label_coords(-0.14, 0.5)

           ax = fig.add_subplot(335)
           plt.plot(time_array/365.25, brine_aquifer2,
                    color='rosybrown', linewidth=line_width)
           plt.xlabel('Time, years', fontsize=label_size)
           plt.ylabel('Leakage rates, kg/s', fontsize=label_size)
           plt.title(r'Leakage of brine: aquifer 2', fontsize=font_size)
           plt.tight_layout()
           plt.tick_params(labelsize=ticks_size)
           plt.xlim([0, 50])
           ax.get_yaxis().set_label_coords(-0.14, 0.5)
           
           ax = fig.add_subplot(336)
           plt.plot(time_array/365.25, brine_aquifer3,
                    color='rosybrown', linewidth=line_width)
           plt.xlabel('Time, years', fontsize=label_size)
           plt.ylabel('Leakage rates, kg/s', fontsize=label_size)
           plt.title(r'Leakage of brine: aquifer 3', fontsize=font_size)
           plt.tight_layout()
           plt.tick_params(labelsize=ticks_size)
           plt.xlim([0, 50])
           ax.get_yaxis().set_label_coords(-0.14, 0.5)

           fig.subplots_adjust(
               left=0.1, top=0.95, right=0.95, bottom=0.05, wspace=0.22)

           # to_save_png = False
           if to_save_png:
               plt.savefig('MSWObservationsCombinedPlotVsTime.jpeg', dpi=300)
               print('"MSWObservationsCombinedPlotVsTime.jpeg" successfully saved.')
               
           fig = plt.figure(figsize=(13, 12))

           ax = fig.add_subplot(331)
           plt.plot(time_array/365.25, CO2_aquifer1-CO2_aquifer2,
                    color='steelblue', linewidth=line_width)
           plt.xlabel('Time, years', fontsize=label_size)
           plt.ylabel('Net leakage rates, kg/s', fontsize=label_size)
           plt.title(r'Net leakage of CO$_2$: aquifer 1', fontsize=font_size)
           plt.tick_params(labelsize=ticks_size)
           plt.xlim([0, 50])
           ax.get_yaxis().set_label_coords(-0.13, 0.5)

           ax = fig.add_subplot(332)
           plt.plot(time_array/365.25, CO2_aquifer2-CO2_aquifer3,
                    color='steelblue', linewidth=line_width)
           plt.xlabel('Time, years', fontsize=label_size)
           plt.ylabel('Net leakage rates, kg/s', fontsize=label_size)
           plt.title(r'Net leakage of CO$_2$: aquifer 2', fontsize=font_size)
           plt.tick_params(labelsize=ticks_size)
           plt.xlim([0, 50])
           ax.get_yaxis().set_label_coords(-0.13, 0.5)
           
           ax = fig.add_subplot(333)
           plt.plot(time_array/365.25, CO2_aquifer3,
                    color='steelblue', linewidth=line_width)
           plt.xlabel('Time, years', fontsize=label_size)
           plt.ylabel('Net leakage rates, kg/s', fontsize=label_size)
           plt.title(r'Net leakage of CO$_2$: aquifer 3', fontsize=font_size)
           plt.tick_params(labelsize=ticks_size)
           plt.xlim([0, 50])
           ax.get_yaxis().set_label_coords(-0.13, 0.5)

           ax = fig.add_subplot(334)
           plt.plot(time_array/365.25, brine_aquifer1-brine_aquifer2,
                    color='rosybrown', linewidth=line_width)
           plt.xlabel('Time, years', fontsize=label_size)
           plt.ylabel('Net leakage rates, kg/s', fontsize=label_size)
           plt.title(r'Net leakage of brine: aquifer 1', fontsize=font_size)
           plt.tight_layout()
           plt.tick_params(labelsize=ticks_size)
           plt.xlim([0, 50])
           ax.get_yaxis().set_label_coords(-0.14, 0.5)

           ax = fig.add_subplot(335)
           plt.plot(time_array/365.25, brine_aquifer2-brine_aquifer3,
                    color='rosybrown', linewidth=line_width)
           plt.xlabel('Time, years', fontsize=label_size)
           plt.ylabel('Net leakage rates, kg/s', fontsize=label_size)
           plt.title(r'Net leakage of brine: aquifer 2', fontsize=font_size)
           plt.tight_layout()
           plt.tick_params(labelsize=ticks_size)
           plt.xlim([0, 50])
           ax.get_yaxis().set_label_coords(-0.14, 0.5)
           
           ax = fig.add_subplot(336)
           plt.plot(time_array/365.25, brine_aquifer3,
                    color='rosybrown', linewidth=line_width)
           plt.xlabel('Time, years', fontsize=label_size)
           plt.ylabel('Net leakage rates, kg/s', fontsize=label_size)
           plt.title(r'Net leakage of brine: aquifer 3', fontsize=font_size)
           plt.tight_layout()
           plt.tick_params(labelsize=ticks_size)
           plt.xlim([0, 50])
           ax.get_yaxis().set_label_coords(-0.14, 0.5)

           fig.subplots_adjust(
               left=0.1, top=0.95, right=0.95, bottom=0.05, wspace=0.22)

           # to_save_png = False
           if to_save_png:
               plt.savefig('MSWObservationsCombinedPlotVsTime2.jpeg', dpi=300)
               print('"MSWObservationsCombinedPlotVsTime2.jpeg" successfully saved.')  
               
           # Pressure and Saturation - Reservoir    
           fig = plt.figure(figsize=(3.0, 6))

           ax = fig.add_subplot(211)
           plt.plot(time_array/365.25, pressure*1e-6,
                    color='steelblue', linewidth=line_width)
           plt.xlabel('Time, years', fontsize=label_size)
           plt.ylabel('Pressure, MPa', fontsize=label_size)
           plt.title(r'Reservoir Pressure at the well', fontsize=font_size)
           plt.tick_params(labelsize=ticks_size)
           plt.xlim([0, 50])
           ax.get_yaxis().set_label_coords(-0.13, 0.5)

           ax = fig.add_subplot(212)
           plt.plot(time_array/365.25, brine_aquifer1-brine_aquifer2,
                    color='rosybrown', linewidth=line_width)
           plt.xlabel('Time, years', fontsize=label_size)
           plt.ylabel('CO2 saturation', fontsize=label_size)
           plt.title(r'Reservoir CO2 Saturation at the well', fontsize=font_size)
           plt.tight_layout()
           plt.tick_params(labelsize=ticks_size)
           plt.xlim([0, 50])
           ax.get_yaxis().set_label_coords(-0.14, 0.5)

           fig.subplots_adjust(
               left=0.1, top=0.95, right=0.95, bottom=0.05, wspace=0.22)

           # to_save_png = False
           if to_save_png:
               plt.savefig('ReservoirPressureCO2saturation.jpeg', dpi=300)
               print('"ReservoirPressureCO2saturation.jpeg" successfully saved.')      

    else:

        print('------------------------------------------------------------------')
        print('                          UQ illustration ')
        print('------------------------------------------------------------------')

        import random
        num_samples = 50
        ncpus = 1
        
        # Draw Latin hypercube samples of parameter values
        seed = random.randint(500, 1100)

        s = sm.lhs(siz=num_samples, seed=seed)

        # Run model using values in samples for parameter values
        results = s.run(cpus=ncpus, verbose=False)

        # Extract results from stochastic simulations
        outputs = s.collect_observations_as_time_series()
        
        CO2_aquifer1 = outputs['ms.CO2_aquifer1']
        CO2_aquifer2 = outputs['ms.CO2_aquifer2']
        CO2_aquifer3 = outputs['ms.CO2_aquifer3']
        CO2_aquifer4 = outputs['ms.CO2_aquifer4']
        brine_aquifer1 = outputs['ms.brine_aquifer1']
        brine_aquifer2 = outputs['ms.brine_aquifer2']
        brine_aquifer3 = outputs['ms.brine_aquifer3']
        brine_aquifer4 = outputs['ms.brine_aquifer4']
        mass_CO2_aquifer1 = outputs['ms.mass_CO2_aquifer1']
        mass_CO2_aquifer2 = outputs['ms.mass_CO2_aquifer2']
        mass_CO2_aquifer3 = outputs['ms.mass_CO2_aquifer3']
        mass_CO2_aquifer4 = outputs['ms.mass_CO2_aquifer4']

        mass_CO2_reservoir = outputs['res.mass_CO2_reservoir']

        print('UQ run is done. ')

        if isPlottingOn:
            # Plot result 1 : CO2 mass rate, kg/s
            label_size = 8
            font_size = 10
            ticks_size = 6
            line_width = 1
            fig = plt.figure(1, figsize=(13, 3.5))

            ax = fig.add_subplot(141)
            for j in range(num_samples):
                plt.plot(time_array/365.25, CO2_aquifer1[j],
                         color='steelblue', linewidth=line_width)
            plt.xlabel('Time, t (years)', fontsize=label_size)
            plt.ylabel('Leakage rates, q (kg/s)', fontsize=label_size)
            plt.title(r'Leakage of CO$_2$: aquifer 1', fontsize=font_size)
            plt.tick_params(labelsize=ticks_size)
            plt.xlim([0, num_years])
            ax.get_yaxis().set_label_coords(-0.13, 0.5)

            ax = fig.add_subplot(142)
            for j in range(num_samples):
                plt.plot(time_array/365.25, CO2_aquifer2[j],
                         color='steelblue', linewidth=line_width)
            plt.xlabel('Time, t (years)', fontsize=label_size)
            plt.ylabel('Leakage rates, q (kg/s)', fontsize=label_size)
            plt.title(r'Leakage of CO$_2$: aquifer 2', fontsize=font_size)
            plt.tick_params(labelsize=ticks_size)
            plt.xlim([0, num_years])
            ax.get_yaxis().set_label_coords(-0.13, 0.5)

            ax = fig.add_subplot(143)
            for j in range(num_samples):
                plt.plot(time_array/365.25, CO2_aquifer3[j],
                         color='steelblue', linewidth=line_width)
            plt.xlabel('Time, t (years)', fontsize=label_size)
            plt.ylabel('Leakage rates, q (kg/s)', fontsize=label_size)
            plt.title(r'Leakage of CO$_2$: aquifer 3', fontsize=font_size)
            plt.tick_params(labelsize=ticks_size)
            plt.xlim([0, num_years])
            ax.get_yaxis().set_label_coords(-0.13, 0.5)

            ax = fig.add_subplot(144)
            for j in range(num_samples):
                plt.plot(time_array/365.25, CO2_aquifer4[j],
                         color='steelblue', linewidth=line_width)
            plt.xlabel('Time, t (years)', fontsize=label_size)
            plt.ylabel('Leakage rates, q (kg/s)', fontsize=label_size)
            plt.title(r'Leakage of CO$_2$: aquifer 4', fontsize=font_size)
            plt.tick_params(labelsize=ticks_size)
            plt.xlim([0, num_years])
            ax.get_yaxis().set_label_coords(-0.13, 0.5)

            # to_save_png = False
            if to_save_png:
                plt.savefig('MSWObservationsCombinedPlot1VsTime.png', dpi=300)
                print('"MSWObservationsCombinedPlot1VsTime.png" successfully saved.')

            # Plot result 2 : CO2 mass, kg
            label_size = 8
            font_size = 10
            ticks_size = 6
            line_width = 1
            fig = plt.figure(2, figsize=(13, 3.5))

            # CO2 mass
            ax = fig.add_subplot(141)
            for j in range(num_samples):
                plt.plot(time_array/365.25, mass_CO2_aquifer1[j],
                         color='steelblue', linewidth=line_width)
            plt.xlabel('Time, t (years)', fontsize=label_size)
            plt.ylabel('Leaked mass, kg', fontsize=label_size)
            plt.title(r'Leakage of CO$_2$: aquifer 1', fontsize=font_size)
            plt.tick_params(labelsize=ticks_size)
            ax.yaxis.get_offset_text().set_fontsize(ticks_size)
            plt.xlim([0, num_years])
            ax.get_yaxis().set_label_coords(-0.13, 0.5)

            ax = fig.add_subplot(142)
            for j in range(num_samples):
                plt.plot(time_array/365.25, mass_CO2_aquifer2[j],
                         color='steelblue', linewidth=line_width)
            plt.xlabel('Time, t (years)', fontsize=label_size)
            plt.ylabel('Leaked mass, kg', fontsize=label_size)
            plt.title(r'Leakage of CO$_2$: aquifer 2', fontsize=font_size)
            plt.tick_params(labelsize=ticks_size)
            ax.yaxis.get_offset_text().set_fontsize(ticks_size)
            plt.xlim([0, num_years])
            ax.get_yaxis().set_label_coords(-0.13, 0.5)

            ax = fig.add_subplot(143)
            for j in range(num_samples):
                plt.plot(time_array/365.25, mass_CO2_aquifer3[j],
                         color='steelblue', linewidth=line_width)
            plt.xlabel('Time, t (years)', fontsize=label_size)
            plt.ylabel('Leaked mass, kg', fontsize=label_size)
            plt.title(r'Leakage of CO$_2$: aquifer 3', fontsize=font_size)
            plt.tick_params(labelsize=ticks_size)
            ax.yaxis.get_offset_text().set_fontsize(ticks_size)
            plt.xlim([0, num_years])
            ax.get_yaxis().set_label_coords(-0.13, 0.5)

            ax = fig.add_subplot(144)
            for j in range(num_samples):
                plt.plot(time_array/365.25, mass_CO2_aquifer4[j],
                         color='steelblue', linewidth=line_width)
            plt.xlabel('Time, t (years)', fontsize=label_size)
            plt.ylabel('Leaked mass, kg', fontsize=label_size)
            plt.title(r'Leakage of CO$_2$: aquifer 4', fontsize=font_size)
            plt.tick_params(labelsize=ticks_size)
            ax.yaxis.get_offset_text().set_fontsize(ticks_size)
            plt.xlim([0, num_years])
            ax.get_yaxis().set_label_coords(-0.13, 0.5)

            # to_save_png = False
            if to_save_png:
                plt.savefig('MSWObservationsCombinedPlot2VsTime.png', dpi=300)
                print('"MSWObservationsCombinedPlot2VsTime.png" successfully saved.')

            # Plot result 3 : brine mass rate, kg/s
            label_size = 8
            font_size = 10
            ticks_size = 6
            line_width = 1
            fig = plt.figure(3, figsize=(13, 3.5))

            ax = fig.add_subplot(141)
            for j in range(num_samples):
                plt.plot(time_array/365.25, brine_aquifer1[j],
                         color='rosybrown', linewidth=line_width)
            plt.xlabel('Time, t (years)', fontsize=label_size)
            plt.ylabel('Leakage rates, q (kg/s)', fontsize=label_size)
            plt.title(r'Leakage of brine: aquifer 1', fontsize=font_size)
            plt.tight_layout()
            plt.tick_params(labelsize=ticks_size)
            plt.xlim([0, num_years])
            ax.get_yaxis().set_label_coords(-0.14, 0.5)

            ax = fig.add_subplot(142)
            for j in range(num_samples):
                plt.plot(time_array/365.25, brine_aquifer2[j],
                         color='rosybrown', linewidth=line_width)
            plt.xlabel('Time, t (years)', fontsize=label_size)
            plt.ylabel('Leakage rates, q (kg/s)', fontsize=label_size)
            plt.title(r'Leakage of brine: aquifer 2', fontsize=font_size)
            plt.tight_layout()
            plt.tick_params(labelsize=ticks_size)
            plt.xlim([0, num_years])
            ax.get_yaxis().set_label_coords(-0.14, 0.5)

            ax = fig.add_subplot(143)
            for j in range(num_samples):
                plt.plot(time_array/365.25, brine_aquifer3[j],
                         color='rosybrown', linewidth=line_width)
            plt.xlabel('Time, t (years)', fontsize=label_size)
            plt.ylabel('Leakage rates, q (kg/s)', fontsize=label_size)
            plt.title(r'Leakage of brine: aquifer 3', fontsize=font_size)
            plt.tight_layout()
            plt.tick_params(labelsize=ticks_size)
            plt.xlim([0, num_years])
            ax.get_yaxis().set_label_coords(-0.14, 0.5)

            ax = fig.add_subplot(144)
            for j in range(num_samples):
                plt.plot(time_array/365.25, brine_aquifer4[j],
                         color='rosybrown', linewidth=line_width)
            plt.xlabel('Time, t (years)', fontsize=label_size)
            plt.ylabel('Leakage rates, q (kg/s)', fontsize=label_size)
            plt.title(r'Leakage of brine: aquifer 4', fontsize=font_size)
            plt.tight_layout()
            plt.tick_params(labelsize=ticks_size)
            plt.xlim([0, num_years])
            ax.get_yaxis().set_label_coords(-0.14, 0.5)

            # to_save_png = False
            if to_save_png:
                plt.savefig('MSWObservationsCombinedPlot3VsTime.png', dpi=300)
                print('"MSWObservationsCombinedPlot3VsTime.png" successfully saved.')
