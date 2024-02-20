"""
The module contains the solution class for the Multisegmented Wellbore AI 
Component.

The model is based on the work of Baek et al.(:cite:'BaekEtAl2021b', :cite:'BaekEtAl2023').

Author: Seunghwan Baek, Veronika S. Vasylkivska, and Nate Mitchell
Date: 1/1/2024
"""
import logging
import os
import sys
import numpy as np
import scipy.special as scm
from scipy.interpolate import interp1d
import warnings
import pandas as pd
import time

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from multisegmented import units


class Parameters():
    """ Parameters class for multisegmented wellbore ROM."""
    def __init__(self):
        """ Create an instance of the Parameters class. """
        self.numberOfShaleLayers = None

        self.shaleThickness = None
        self.aquiferThickness = None
        self.reservoirThickness = None

        # Permeability
        self.shalePermeability = None
        self.aquiferPermeability = None
        self.aquiferPorosity = None
        self.logAquiferHPerm = None
        self.logAquiferVPerm = None

        # Land surface pressure
        self.datumPressure = None
        self.tempGrad = None

        # Input pressure, CO2 saturation and time points
        self.pressure = None
        self.CO2saturation = None
        self.timePoint = None
        self.pressure_prev = None
        self.pressureShaleTop_prev = None
        self.pressureShaleBot_prev = None
        
        # Amount of CO2 present in aquifers before simulation
        self.prevCO2Mass = None

        # Amount of brine present in aquifers before simulation
        self.prevBrineMass = None

        self.timeSpan = None
        self.timeStep = None

        # Density
        self.brineDensity = None
        self.brineDensityAquifer = None  
        self.brineDensityForm = "Constant"     # constant or functional
        self.CO2Density = None        
        self.CO2DensityShale = None  
        self.CO2DensityForm = "Constant"       # constant or functional

        # Viscosity
        self.brineViscosity = None
        self.brineViscosityAquifer = None    
        self.brineViscosityForm = "Constant"     # constant or functional
        self.CO2Viscosity = None
        self.CO2ViscosityAquifer = None  
        self.CO2ViscosityForm = "Constant"       # constant or functional

        # Residual saturation and compressibility
        self.aquBrineResSaturation = None
        self.compressibility = None

        # Well radius
        self.wellRadius = None
        self.flowArea = None
        
        self.SaltMassFrac = None
        self.useDLmodel = None

        self.inputarray1 = None
        self.dl1input_minthree = None
        self.dl1input_mintwo = None
        self.dl1input_minone = None        
        self.inputarray2 = None
        self.dlinput_minthree = None
        self.dlinput_mintwo = None
        self.dlinput_minone = None

class Solution():
    """ Solution class for multisegmented wellbore ROM."""
    def __init__(self, parameters, FluidModels=None, DLModels1=None,
                 DLModels2=None):
        """ Create an instance of the Solution class."""

        # Setup the parameters that define a solution
        self.parameters = Parameters()
        self.parameters = parameters

        # Reserve a space for the solution pieces
        self.CO2LeakageRates = None
        self.brineLeakageRates = None
        self.g = 9.8  # acceleration due to gravity

        # Reserve space for instance attributes
        self.nSL = None
        self.reservoirBottom = None
        self.depth = None

        self.interface = None
        self.thicknessH = None
        self.initialPressure = None
        self.iniTopPressure = None

        self.CO2Mass = None
        self.CO2SaturationAq = None

        self.CO2MobilityAq = None
        self.brineMobilityAq = None

        self.CO2MobilityShale = None
        self.brineMobilityShale = None
        self.effectiveMobility = None

        self.CO2Q = None
        self.brineQ = None
        self.Lambda = None
        
        self.brineDensityAquifer = None   
        self.CO2DensityAquifer = None 
        self.brineViscosityAquifer = None  
        self.CO2ViscosityAquifer = None 
  
        self.brineDensityShale = None
        self.CO2DensityShale = None
        
        self.brineDensityShale = None
        self.CO2DensityShale = None
        
        self.FluidModels = FluidModels
        self.DLModels1 = DLModels1
        self.DLModels2 = DLModels2
    
    def update_parameters(self, parameters):
        """ Update the parameter object """
        self.parameters = parameters

    def get_depth_for_pressure(self):
        """ Calculate depth of all layers for the pressure calculations."""
        # The deepest shale has the first thickness in the array
        shaleThickness = self.parameters.shaleThickness
        aquiferThickness = self.parameters.aquiferThickness
        reservoirThickness = self.parameters.reservoirThickness

        # Get the sum of shale layers thicknesses
        shaleThicknessSum = sum(shaleThickness)

        # Number of aquifer includes reservoir so
        # number of aquifers = (number of shales - 1) + 1 = number of shales
        aquiferNum = self.parameters.numberOfShaleLayers

        # Get the sum of aquifer (not including reservoir) thicknesses
        if aquiferNum == 2:
            aquiferThicknessSum = aquiferThickness
        else:
            aquiferThicknessSum = sum(aquiferThickness)

        # Get the depth of the bottom of the reservoir
        self.reservoirBottom = (shaleThicknessSum + aquiferThicknessSum +
                                reservoirThickness)
        self.reservoirTop = (shaleThicknessSum + aquiferThicknessSum)

        # We need to know the depths to the bottom of each aquifer
        # Initialize depth array for depths of the bottoms of the aquifers
        self.depth = np.zeros((aquiferNum,))

        self.depth[0] = self.reservoirBottom
        self.depth[1] = self.depth[0] - (
            shaleThickness[0] + reservoirThickness)

        if aquiferNum > 2:
            for ind in range(2, aquiferNum):
                self.depth[ind] = self.depth[ind-1] - (
                    shaleThickness[ind-1] + aquiferThickness[ind-2])

    def get_fluidProps(self, pressure):   
    
        temp_grad = self.parameters.tempGrad * 1e-3 #C/m
        SaltMassFrac = self.parameters.SaltMassFrac
        
        depth = self.depth - self.thicknessH * 0.5 
            
        Temp = depth*temp_grad + 14.7
        
        # Condi_brine = np.zeros((self.nSL,3))        
        # Condi_CO2 = np.zeros((self.nSL,2))
        
        features = pd.DataFrame(np.zeros((self.nSL,3)),columns=['P,MPa', 'T,C', 'Salinity'])
        features.loc[:,'P,MPa'] = pressure / 1e+6
        features.loc[:,'T,C'] = Temp
        features.loc[:,'Salinity'] = SaltMassFrac
        
        TempCO2 = self.FluidModels.CO2Model(features["T,C"], features["P,MPa"])
        self.CO2DensityAquifer = TempCO2[:,0]
        self.CO2ViscosityAquifer = TempCO2[:,1]

        TempBrine = self.FluidModels.BrineModel(
            features["T,C"], features["P,MPa"], features["Salinity"])
        self.brineDensityAquifer = TempBrine[:,0]
        self.brineViscosityAquifer = TempBrine[:,1] # unit taken care
        
        TempBrine2 = self.FluidModels.BrineModel(
            features["T,C"], features["P,MPa"] + 1, features["Salinity"])
        self.brinecfAquifer  = -2 / (
            self.brineDensityAquifer+TempBrine2[:,0]) * (
                self.brineDensityAquifer - TempBrine2[:,0]) / (1 * 1e6) # 1/Pa
        
    def setup_initial_conditions(self):
        """ Setup initial conditions at the abandoned well."""
        # Determine depths of the aquifers
        self.get_depth_for_pressure()

        # Setup array for the thickness of aquifers and reservoir
        self.nSL = self.parameters.numberOfShaleLayers
        self.thicknessH = np.ones(self.nSL)
        self.thicknessH[0] = self.parameters.reservoirThickness
        self.thicknessH[1:self.nSL] = self.parameters.aquiferThickness
        
        # Setup initial hydrostatic pressure
        # datumPressure is atmospheric pressure or some reference pressure
        # at the top of the upper shale layer
        pressGrad = 0.009785602 * 1e+6 # hydrostatic, Pa/m
        self.initialPressure = self.parameters.datumPressure + (
            self.depth * pressGrad)
        pressureTopReservoir = self.parameters.datumPressure + (
            (self.depth - self.thicknessH) * pressGrad)
        self.iniTopPressure  = pressureTopReservoir

        self.CO2SaturationAq = np.zeros(self.nSL-1)
        
        pressureAve = (self.initialPressure+pressureTopReservoir)/2  
        self.get_fluidProps(pressureAve) # Use initial values 
        
        # Setup initial masses of CO2 in each aquifer and reservoir
        self.CO2Mass = self.parameters.prevCO2Mass
        self.brineMass = self.parameters.prevBrineMass
  
        # Setup array for the thickness of aquifers and reservoir
        self.thicknessH = np.ones(self.nSL)
        self.thicknessH[0] = self.parameters.reservoirThickness
        self.thicknessH[1:self.nSL] = self.parameters.aquiferThickness
  
        # Setup initial interface height h at the leaking well
        self.interface = np.zeros(self.nSL)
        for j in range(self.nSL):
            self.interface[j] = self.parameters.CO2saturation[j]*self.thicknessH[j]
          
    def get_mobility_for_aquifers(self):
        """ Calculate mobility parameter for brine and CO2 in the aquifers."""

        CO2RelPermAq = np.zeros(self.nSL)
        brineRelPermAq = np.ones(self.nSL)

        Sres = self.parameters.aquBrineResSaturation

        for ind in range(self.nSL):
            CO2RelPermAq[ind] = np.minimum((
                1-Sres[ind]), self.interface[ind] / self.thicknessH[ind])
            
            if (1 - Sres[ind]) == 0:
                brineRelPermAq[ind] = 0
            else:
                brineRelPermAq[ind] = (1 / (1 - Sres[ind])) * (
                    1 - Sres[ind] - CO2RelPermAq[ind])

        self.CO2MobilityAq = CO2RelPermAq / self.CO2ViscosityAquifer
        self.brineMobilityAq = brineRelPermAq / self.brineViscosityAquifer
                
        self.effectiveMobility = ((self.interface * self.CO2MobilityAq) + (
            (self.thicknessH - self.interface) * self.brineMobilityAq)
            ) / self.thicknessH

        # For use in well functions
        for ind in range(0, self.nSL):
            if self.interface[ind] == 0:
                self.effectiveMobility[ind] = 1 / self.brineViscosityAquifer[ind]

    def find(self):
        """ Find solution of the ROM corresponding to the provided parameters."""
        timeSpan = self.parameters.timeSpan # days *units.day()
        timeSpan_sec = self.parameters.timeSpan * units.day()
        timePoint = self.parameters.timePoint
        self.parameters.timeStep += 1 # starts from zero
        self.timeStep = self.parameters.timeStep
        
        # Setup initial pressure, saturations, masses and interface height
        self.setup_initial_conditions()

        # Setup lambda and mobility
        self.find_lambda()

        self.get_mobility_for_aquifers()

        self.get_mobility_for_shales()

        # Initialize matrix for system of equations. AA*Pressure(Aquifers) = BB
        # LUT coupled: the pressure at the bottom of reservoir (or aquifer 1)
        # is known so it's not a variable
        StartingIdx = 1

        # tridiagonal (N-1)*(N-1) (N: # of aquifers including the bottom reservoir =  self.nSL)
        AA = np.zeros((self.nSL - StartingIdx, self.nSL - StartingIdx))
        BB = np.zeros(self.nSL - StartingIdx) # vector (N-1)*1

        # Calculate coefficients of Darcy equation for each aquifer for the matrix build-up
        # Refer to the document for the detailed information
        CC = self.CC_shale()
        GG = self.GG_shale()
        FF = self.FF_aquifer()
        WV = self.WV_aquifer()
        
        pressGrad = 0.009785602*1e+6 # hydrostatic, Pa/m  = 9.792 kPa/m
        self.initialAvePressure = self.parameters.datumPressure + (
            self.depth - (self.thicknessH * 0.5)) * pressGrad

        for ind in range(StartingIdx, self.nSL):

            F_below = FF[ind-1]
            if ind == StartingIdx:
                F_below = 0     # nothing below 1st shale layer

            F_zero = FF[ind]

            if ind < self.nSL-1:
                F_above = FF[ind+1]
            else:
                F_above = 0 # nothing above top shale layer

            C_below = CC[ind-1]
            C_zero = CC[ind]

            G_below = GG[ind-1]
            G_zero = GG[ind]

            # aquifer2 (bottom). Do not calculate the bottom reservoir due to LUT coupled case
            if ind == StartingIdx:
                # Top reservoir, not averaged, pressure from LUT
                ReservoirTopPressure = self.parameters.pressure

                AA[ind-StartingIdx, ind-StartingIdx] = 1 + (
                    WV[ind] * C_below) - (WV[ind] * C_zero)
                AA[ind-StartingIdx, ind+1-StartingIdx] = WV[ind] * C_zero

                BB[ind-StartingIdx] = WV[ind]*(-C_below*F_zero-G_below+C_zero*F_zero+\
                    C_zero*F_above-G_zero-C_below*ReservoirTopPressure) +self.initialAvePressure[ind]

            elif ind == (self.nSL-1): # aquifer N (top)
                DatumPressure = self.parameters.datumPressure

                AA[ind-StartingIdx, ind-1-StartingIdx] = -WV[ind] * C_below
                AA[ind-StartingIdx, ind-StartingIdx] = 1 + (
                    WV[ind] * C_below) - (WV[ind-1] * C_zero)

                BB[ind-StartingIdx] = WV[ind-1] * (
                    (-C_below * F_below) - (C_below * F_zero) - G_below + 
                    (C_zero * F_zero) - G_zero - (C_zero * DatumPressure)
                    ) + self.initialAvePressure[ind] 

            else: # aquifer3 to aquifer N-1 (intermediate)

                AA[ind-StartingIdx, ind-1-StartingIdx] = -WV[ind] * C_below
                AA[ind-StartingIdx, ind-StartingIdx] = 1 + (
                    WV[ind] * C_below) - (WV[ind] * C_zero)
                AA[ind-StartingIdx, ind+1-StartingIdx] = WV[ind] * C_zero

                BB[ind-StartingIdx] = WV[ind] * (
                    (-C_below * F_below) - (C_below * F_zero) - G_below + 
                    (C_zero * F_zero) + (C_zero * F_above) - G_zero
                    ) + self.initialAvePressure[ind] 

        # Calculate coupled pressure of each aquifer as solution of linear system
        # Vertically averaged delta pressure (compared to the initial values) at each aquifer
        Pave = np.linalg.solve(AA, BB) # Vertically averaged pressure at each aquifer
        self.pressure_prev = Pave 

        if (StartingIdx == 1): # LUT coupled 
            
            Pbot = Pave + FF[StartingIdx:] 
            Pbot = np.append(self.initialPressure[0], Pbot) # inclusion of the bottom layer of the storage reservoir pressure

            Ptop = Pave - FF[StartingIdx:] 
            Ptop = np.append(ReservoirTopPressure, Ptop) # inclusion of the top layer of the storage reservoir pressure from LUT
          
            # checked ok
            self.pressureShaleBot_prev = Ptop[0:] 
            self.pressureShaleTop_prev = self.pressureShaleBot_prev.copy()
            self.pressureShaleTop_prev[:len(Pbot[1:])] = Pbot[1:]
            self.pressureShaleTop_prev[-1] = self.parameters.datumPressure
    
            # Find change in pressure to calculate flow rate 
            deltaP = np.zeros(self.nSL)        
            deltaP[0:self.nSL-1] = Ptop[0:self.nSL-1] - Pbot[1:self.nSL]                
            deltaP[self.nSL-1]   = Ptop[self.nSL-1] - self.parameters.datumPressure

        DarcyCoefCO2      = (self.parameters.flowArea * self.parameters.shalePermeability 
                             * self.CO2MobilityShale)
        DarcyCoefCO2Shale = DarcyCoefCO2 * (1 / self.parameters.shaleThickness)
        GravityCO2Shale   = DarcyCoefCO2 * self.CO2DensityShale * self.g
        
        DarcyCoefBrine      = (self.parameters.flowArea * self.parameters.shalePermeability 
                               * self.brineMobilityShale)
        DarcyCoefBrineShale = DarcyCoefBrine * (1 / self.parameters.shaleThickness)        
        GravityBrineShale   = DarcyCoefBrine * self.brineDensityShale * self.g        
                
        # Inflow volumetric rate along shale layers = inflow volumetric rate into the aquifer above the shale layer
        self.CO2Q   = np.zeros(self.nSL)
        self.brineQ = np.zeros(self.nSL)   
        
        self.CO2Q[:]   = (DarcyCoefCO2Shale*deltaP - GravityCO2Shale) # volume rate
        self.brineQ[:] = (DarcyCoefBrineShale*deltaP - GravityBrineShale) # volume rate
              
        # Inflow mass rate along shale layers = inflow mass rate into the aquifer above the shale layer       
        self.CO2LeakageRates = np.zeros(self.nSL)
        self.brineLeakageRates = np.zeros(self.nSL)

        self.CO2LeakageRates   = self.CO2Q * self.CO2DensityShale # mass rate        
        self.brineLeakageRates = self.brineQ * self.brineDensityShale # mass rate
        
        if self.parameters.useDLmodel == 1:
            
            if self.timeStep >= 4: # due to look back = 3
                # Caprock segment
                features = self.dataprocessing4tsmodel_caproock(
                    self.parameters.dl1input_minthree, 
                    self.parameters.dl1input_mintwo, 
                    self.parameters.dl1input_minone)
                
                self.CO2LeakageRates[0] = self.DLModels1.model_rf_qC.predict(features)  
                self.brineLeakageRates[0] = self.DLModels1.model_rf_qB.predict(features)
                CO2_Saturation = self.DLModels1.model_rf_sC.predict(features)  

                # Aquifer segment
                features = self.dataprocessing4tsmodel(
                    self.parameters.dl2input_minthree, 
                    self.parameters.dl2input_mintwo, 
                    self.parameters.dl2input_minone)
           
                qB_pred = self.DLModels2.model_rf_qB.predict(features)  
                qC_pred = self.DLModels2.model_rf_qC.predict(features)  
                SC_pred = self.DLModels2.model_rf_sC.predict(features)  
            
                self.brineLeakageRates[1:-1] = qB_pred # excluding bottom aqu and atm
                self.CO2LeakageRates[1:-1] = qC_pred # excluding bottom aqu and atm
        
        for j in range(1,self.nSL):
            self.CO2LeakageRates[j] = np.minimum(self.CO2LeakageRates[j-1], 
                                                 self.CO2LeakageRates[j])
            self.CO2LeakageRates[j] = np.maximum(0, self.CO2LeakageRates[j])    
        
        # Update CO2 cumulative mass in each aquifer (not including the bottom reservoir)
        self.CO2Mass = self.CO2Mass + (       
            self.CO2LeakageRates[0:self.nSL-1] 
            - self.CO2LeakageRates[1:self.nSL]) * timeSpan_sec
            
        self.CO2Mass = np.maximum(self.CO2Mass, np.zeros(self.nSL - 1))

        # Update brine cumulative mass in each aquifer (not including the bottom reservoir)
        self.brineMass = self.brineMass + (
            self.brineLeakageRates[0:self.nSL-1] 
            - self.brineLeakageRates[1:self.nSL]) * timeSpan_sec
        
        self.brineMass = np.maximum(self.brineMass, np.zeros(self.nSL - 1))

        self.get_interface()
        for j in range(self.nSL-1):
            self.CO2SaturationAq[j] = self.interface[j+1] / self.thicknessH[j+1]
        
        if (self.parameters.useDLmodel == 1) & (self.timeStep >= 4):
            self.CO2SaturationAq[1:] = CO2_Saturation
            self.CO2SaturationAq[0] = np.minimum(CO2_Saturation, 
                                                 self.CO2SaturationAq[0])
            
        for j in range(1,self.nSL-1):
            self.CO2SaturationAq[j] = np.minimum(self.CO2SaturationAq[j-1], 
                                                 self.CO2SaturationAq[j]) 
            
        # Build inputarray for the next time step for shale segment model
        self.inputarray1 = {'TopDepth':self.depth[1],
                        'BottomDepth':self.depth[1] + self.parameters.shaleThickness[0],
                        'Sg_bot': self.interface[0] / self.thicknessH[0],
                        'WellRadius': self.parameters.wellRadius,
                        'Salinity': self.parameters.SaltMassFrac[0],
                        'PbotMPa': self.parameters.pressure * 1e-6, # LUT
                        'ShaleThickness': self.parameters.shaleThickness[0],
                        'WellPermeability':self.parameters.shalePermeability[0]}
        
        # Build inputarray for the next time step for aquifer segment model
        self.inputarray2 = {
                        'BottomDepth':np.array([self.depth[ii+1] + self.thicknessH[ii+1] 
                                                for ii in range(len(self.depth)-1)]),
                        'MidDepth':np.array([self.depth[ii+1] 
                                             for ii in range(len(self.depth) - 1)]),
                        'TopDepth':np.array([self.depth[ii+1] 
                                             - self.parameters.shaleThickness[ii+1] 
                                             for ii in range(len(self.depth) - 1)]),
                        'ShaleThickness': self.parameters.shaleThickness[1:],
         
                        'PbotMPa': self.pressureShaleTop_prev[:-1] * 1e-6, # aquifer bottom
                        'Sg_bot': np.array(self.CO2SaturationAq[:-1]), 
                        'LeakRateCO2,kg/s': self.CO2LeakageRates[:-1], # excluding release to atm
                        'LeakRateBrine,kg/s': self.brineLeakageRates[:-1], # excluding release to atm
                        'LeakAqCumMassCO2,kg':  self.CO2Mass, 
                        'LeakAqCumMassBrine,kg': self.brineMass,
                        'simtime,days': timePoint,
                        'timespan,days': timeSpan,
                        'nStep': self.timeStep,
                        'Salinity': np.array([self.parameters.SaltMassFrac[ii+1] 
                                              for ii in range(len(self.depth)-1)]), 
                        
                        'WellRadius': self.parameters.wellRadius,
                        'aquiferWellPermeability':self.parameters.aquiferPermeability,
                        'shaleWellPermeability':self.parameters.shalePermeability[1:],
                        
                        'AquPoro': self.parameters.aquiferPorosity,
                        'PermHAq': self.parameters.logAquiferHPerm,
                        'PermVAq': self.parameters.logAquiferVPerm,
                        
                        'Brine_top,kg/s': self.brineLeakageRates[1:-1],
                        'CO2_top,kg/s': self.CO2LeakageRates[1:-1],
                        'SCO_top': self.CO2SaturationAq[1:],
                        }
       
    def dataprocessing4tsmodel_caproock(self,inputarray_minthree,
                                   inputarray_mintwo,inputarray_minone):
           
                def makingFeatures(FluidModels,inputarray):
                    
                    # Feature processing for DL model prediction
                    Temp_surface     = 15 #C
                    Press_surface    = 101352.9 #Pa
                    Press_grad       = 9.785602 # kPa/m
             
                    TopDepth = inputarray['TopDepth']
                    BottomDepth  = inputarray['BottomDepth']
                    Sg_bot  = inputarray['Sg_bot']
                    WellRadius = inputarray['WellRadius']
                    Salinity = inputarray['Salinity']
                    PbotMPa = inputarray['PbotMPa']
                    ShaleThickness = inputarray['ShaleThickness']
                    WellPermeability = inputarray['WellPermeability']
            
                    features = pd.DataFrame(
                        np.zeros((1, 33)), 
                        columns=[
                            'T_bot,C', 'Pavg_bot,MPa', 'Sg_bot', 'rhoG_bot,kg/m3', 
                            'rhoL_bot,kg/m3', 'visG_bot,Pa-s', 'visL_bot,Pa-s', 
                            'Pavg_top,MPa', 'rhoG_top,kg/m3', 'rhoL_top,kg/m3', 
                            'shaleDepth,m', 'shaleThickness,m', 'shalePerm,mD', 
                            'TempGrad,C/km', 'TopDepth,m', 'Pavg,bot-top', 
                            'Pavg-Pini,bot', 'Salinity', 'WellRadius,m', 'Sb_bot', 
                            'temp_top,C', 'visG_top,Pa-s', 'visL_top,Pa-s', 
                            'gammaB', 'gammaC', 'hydCondtopB', 'hydCondbotB', 
                            'hydCondRatioB', 'hydCondtopC', 'hydCondbotC', 
                            'hydCondRatioC', 'mobRatiotop', 'mobRatiobot'])
                    
                    features["TempGrad,C/km"] = self.parameters.tempGrad
                    features['TopDepth,m'] = TopDepth
                   
                    features['Sg_bot'] = Sg_bot
                    features['Sg_bot'] = features['Sg_bot'].where(
                        features['Sg_bot'] > 0.01, 0.0)
                    features['Sb_bot'] = 1.0-features['Sg_bot'].copy()
                    
                    features['WellRadius,m'] = WellRadius # need to check currently 0.15
                
                    features["T_bot,C"]   = (
                        (BottomDepth * features["TempGrad,C/km"] * 1e-3) 
                        + Temp_surface)
                    features["Pavg_bot,MPa"] = PbotMPa
                
                    features["temp_top,C"]   = (
                        (features['TopDepth,m'] * features["TempGrad,C/km"] * 1e-3) 
                        + Temp_surface)
                    features["Pavg_top,MPa"] = (
                        (features['TopDepth,m'] * Press_grad * 1e-3) 
                        + (Press_surface * 1e-6)) # hydrostatic 
                   
                    features["Salinity"] = Salinity
                   
                    # bottom
                    BottomCO2 = FluidModels.CO2Model(
                        features["T_bot,C"], features["Pavg_bot,MPa"])
                    BottomBrine = FluidModels.BrineModel(
                        features["T_bot,C"], features["Pavg_bot,MPa"], 
                        features["Salinity"])
                   
                    # top
                    TopCO2 = FluidModels.CO2Model(
                        features["temp_top,C"], features["Pavg_top,MPa"])
                    TopBrine = FluidModels.BrineModel(
                        features["temp_top,C"], features["Pavg_top,MPa"], 
                        features["Salinity"])
                   
                    features["rhoG_bot,kg/m3"] = BottomCO2[:,0] 
                    features["visG_bot,Pa-s"] =  BottomCO2[:,1] 
                   
                    features["rhoL_bot,kg/m3"] = BottomBrine[:,0]
                    features["visL_bot,Pa-s"] = BottomBrine[:,1]
                   
                    features["rhoG_top,kg/m3"] = TopCO2[:,0]
                    features["visG_top,Pa-s"] = TopCO2[:,1]
                   
                    features["rhoL_top,kg/m3"] = TopBrine[:,0]
                    features["visL_top,Pa-s"] = TopBrine[:,1] 
                   
                    features['shaleDepth,m'] = BottomDepth
                    features["shaleThickness,m"] = ShaleThickness 
                   
                    features['Pavg,bot-top']  = (
                        features['Pavg_bot,MPa'] - features['Pavg_top,MPa'])
                    features['Pavg-Pini,bot'] = (
                        features['Pavg_bot,MPa'] - (
                            0.1013529 + BottomDepth * 0.009785602))
                   
                    shaleBrineden = (features["rhoL_bot,kg/m3"] 
                                     + features["rhoL_top,kg/m3"]) * 0.5
                    features['gammaB'] = (
                        (features['Pavg_bot,MPa'] - features['Pavg_top,MPa']) 
                        * 1e6) / (shaleBrineden * 9.785602 * features["shaleThickness,m"])
                    
                    shaleCO2den = (features["rhoG_bot,kg/m3"] 
                                   + features["rhoG_top,kg/m3"]) * 0.5
                    
                    features['gammaC']  = (features["Pavg,bot-top"] * 1e6) / (
                        shaleCO2den * 9.785602 * features["shaleThickness,m"])
                    
                    features['shalePerm,mD'] = WellPermeability * 1e15
                    features['shalePerm,mD'] = np.log10(features['shalePerm,mD'])
                   
                    shalePerm = (10.0 ** (features['shalePerm,mD'])) * 1e-15
                    
                    hydCondtopB = (shalePerm * features["rhoL_top,kg/m3"] 
                                   * (9.785602 / features["visL_top,Pa-s"]))
                    hydCondbotB = (shalePerm * features["rhoL_bot,kg/m3"] 
                                   * (9.785602 / features["visL_bot,Pa-s"]))
                    
                    features['hydCondtopB'] = np.log10(hydCondtopB) 
                    features['hydCondbotB'] = np.log10(hydCondbotB) 
                    features['hydCondRatioB'] = hydCondtopB/hydCondbotB
                
                    hydCondtopC = (shalePerm * features["rhoG_top,kg/m3"] 
                                   * (9.785602 / features["visG_top,Pa-s"]))
                    hydCondbotC = (shalePerm * features["rhoG_bot,kg/m3"] 
                                   * (9.785602 / features["visG_bot,Pa-s"]))
                    
                    features['hydCondtopC'] = np.log10(hydCondtopC) 
                    features['hydCondbotC'] = np.log10(hydCondbotC) 
                    features['hydCondRatioC'] = hydCondtopC/hydCondbotC
                
                    mobilityRatiotop = (
                        (features["rhoG_top,kg/m3"] / features["visG_top,Pa-s"]) 
                        / (features["rhoL_top,kg/m3"] / features["visL_top,Pa-s"]))
                    mobilityRatiobot = (
                        (features["rhoG_bot,kg/m3"] / features["visG_bot,Pa-s"]) 
                        / (features["rhoL_bot,kg/m3"] / features["visL_bot,Pa-s"]))
                    
                    features['mobRatiotop'] = np.log10(mobilityRatiotop)
                    features['mobRatiobot'] = np.log10(mobilityRatiobot)
                
                    return features 
                
                # For only bottom layer
                feature_minthree = makingFeatures(self.FluidModels,inputarray_minthree)
                feature_mintwo = makingFeatures(self.FluidModels,inputarray_mintwo)
                fefature_mieone = makingFeatures(self.FluidModels,inputarray_minone)
        
                feature_window_layers = np.concatenate(
                    [feature_minthree, feature_mintwo, fefature_mieone], axis=1)
                    
                return feature_window_layers
            
    def dataprocessing4tsmodel(self,inputarray_minthree,
                               inputarray_mintwo,inputarray_minone):
       
            def makingFeatures(FluidModels,inputarray,nPred):
            
                # Feature processing for DL model prediction
                Temp_surface     = 15 #C
                Press_surface    = 101352.9 #Pa
                Press_grad       = 9.785602 # kPa/m
                
                features = pd.DataFrame(np.zeros((1,68)),columns=
                ['simtime,days', 'AquBotDepth,m', 'AquThickness,m', 'AquPoro',
                'Salinity', 'PermHAq,mD', 'PermVAq,mD', 'PermWAq,mD', 'AqWellRadius,m',
                'ShaThickness,m', 'PermWSh,mD', 'ShWellRadius,m', 'TempGrad,C/km',
                'temp_bot,C', 'temp_mid,C', 'temp_top,C', 'Pavg_bot,MPa',
                'Pavg_mid,MPa', 'Pavg_top,MPa', 'rhoG_bot,kg/m3', 'rhoG_mid,kg/m3',
                'rhoG_top,kg/m3', 'rhoL_bot,kg/m3', 'rhoL_mid,kg/m3', 'rhoL_top,kg/m3',
                'visG_bot,Pa-s', 'visG_mid,Pa-s', 'visG_top,Pa-s', 'visL_bot,Pa-s',
                'visL_mid,Pa-s', 'visL_top,Pa-s', 'LeakRateCO2,kg/s',
                'LeakRateBrine,kg/s', 'CaseNum', 'Sg_bot', 'Sb_bot',
                'PermRatio_AqH_AqV', 'PermRatio_AqH_AqW', 'PermRatio_AqV_AqW',
                'PermRatio_AqH_ShW', 'PermRatio_AqV_ShW', 'PermRatio_AqW_ShW',
                'TopDepth,m', 'MidDepth,m', 'depthRatio1', 'depthRatio2',
                'Pavg,bot-top', 'Pavg-Pini,bot', 'LeakAqCumMassCO2,kg',
                'LeakAqCumMassBrine,kg', 'gammaB', 'gammaC', 'gammaB2', 'gammaC2',
                'hydCondtopB', 'hydCondtopC', 'hydCondRatio2B', 'hydCondRatio2C',
                'hydCondRatio3B', 'hydCondRatio3C', 'hydCondRatioB', 'hydCondRatioC',
                'mobRatiotop', 'mobRatiomid', 'mobRatiobot', 
                'Brine_top,kg/s','CO2_top,kg/s','SCO_top'])
                
                features['simtime,days'] = inputarray['simtime,days']
                features['AquBotDepth,m'] = inputarray['BottomDepth'][nPred]
                features['MidDepth,m'] = inputarray['MidDepth'][nPred]
                features['TopDepth,m'] = inputarray['TopDepth'][nPred]
                
                features['AquThickness,m'] = (inputarray['BottomDepth'] 
                                              - inputarray['MidDepth'])[nPred]
                features['AquPoro'] = inputarray['AquPoro'][nPred]
                features['Salinity'] = inputarray['Salinity'][nPred]
                
                features['PermHAq,mD'] = (np.log10(
                    10 ** inputarray['PermHAq'] * 1e15))[nPred]
                features['PermVAq,mD'] = (np.log10(
                    10 ** inputarray['PermVAq'] * 1e15))[nPred]
                features['PermWAq,mD'] =  (np.log10(
                    inputarray['aquiferWellPermeability'] * 1e15))[nPred]
                features['AqWellRadius,m'] = inputarray['WellRadius']
                features['ShaThickness,m'] = inputarray['ShaleThickness'][nPred]
                
                features['PermWSh,mD'] =  (np.log10(
                    inputarray['shaleWellPermeability'] * 1e15))[nPred]
                features['ShWellRadius,m'] = inputarray['WellRadius']
                features['TempGrad,C/km'] = self.parameters.tempGrad
                
                # Top
                features['temp_top,C'] = (
                    (features['TopDepth,m'] * features["TempGrad,C/km"] * 1e-3) 
                    + Temp_surface)
                features['Pavg_top,MPa'] = (
                    (features['TopDepth,m'] * Press_grad * 1e-3) 
                    + (Press_surface * 1e-6))
                # features["Pavg_bot,MPa"] = inputarray['PbotMPa']*1e-6 # dynamic
                
                Temp3 = FluidModels.CO2Model(features["temp_top,C"], 
                                             features["Pavg_top,MPa"])
                features["rhoG_top,kg/m3"] = Temp3[:,0]
                features["visG_top,Pa-s"] = Temp3[:,1]
                
                Temp2 = FluidModels.BrineModel(features["temp_top,C"], 
                                               features["Pavg_top,MPa"], 
                                               features["Salinity"])
                features["rhoL_top,kg/m3"] = Temp2[:,0]
                features["visL_top,Pa-s"] = Temp2[:,1]
                
                # Mid
                features["temp_mid,C"]   = (
                    (features['MidDepth,m'] * features["TempGrad,C/km"] * 1e-3) 
                    + Temp_surface)
                features["Pavg_mid,MPa"] = (
                    (features['MidDepth,m'] * Press_grad * 1e-3) 
                    + Press_surface * 1e-6)
                
                Temp4 = FluidModels.CO2Model(features["temp_mid,C"], 
                                             features["Pavg_mid,MPa"])
                features["rhoG_mid,kg/m3"] = Temp4[:,0]
                features["visG_mid,Pa-s"] = Temp4[:,1]
                
                Temp5 = FluidModels.BrineModel(features["temp_mid,C"], 
                                               features["Pavg_mid,MPa"], 
                                               features["Salinity"])
                features["rhoL_mid,kg/m3"] = Temp5[:,0]
                features["visL_mid,Pa-s"] = Temp5[:,1]
                
                # Bottom
                features["temp_bot,C"]   = (
                    (features['AquBotDepth,m'] * features["TempGrad,C/km"] * 1e-3) 
                    + Temp_surface)
                features["Pavg_bot,MPa"] = (
                    (features['AquBotDepth,m'] * Press_grad * 1e-3 ) 
                    + Press_surface * 1e-6)
                
                # inputarray['PbotMPa'][nPred]
                
                Temp3 = FluidModels.CO2Model(features["temp_bot,C"], 
                                             features["Pavg_bot,MPa"])
                features["rhoG_bot,kg/m3"] = Temp3[:,0]
                features["visG_bot,Pa-s"] = Temp3[:,1]
                
                Temp2 = FluidModels.BrineModel(features["temp_bot,C"], 
                                               features["Pavg_bot,MPa"], 
                                               features["Salinity"])
                features["rhoL_bot,kg/m3"] = Temp2[:,0]
                features["visL_bot,Pa-s"] = Temp2[:,1]
                
                features['PermRatio_AqH_AqV'] = (features['PermHAq,mD'] 
                                                 / features['PermVAq,mD'])
                features['PermRatio_AqH_AqW'] = np.log10(
                    (10 ** features['PermHAq,mD']) / (10 ** features['PermWAq,mD']))
                features['PermRatio_AqV_AqW'] = np.log10(
                    (10 ** features['PermVAq,mD']) / (10 ** features['PermWAq,mD']))
                features['PermRatio_AqH_ShW'] = np.log10(
                    (10 ** features['PermHAq,mD']) / (10 ** features['PermWSh,mD']))
                features['PermRatio_AqV_ShW'] = np.log10(
                    (10 ** features['PermVAq,mD']) / (10 ** features['PermWSh,mD']))
                features['PermRatio_AqW_ShW'] = np.log10(
                    (10 ** features['PermWAq,mD']) / (10 ** features['PermWSh,mD']))
                
                features['depthRatio1'] = (
                    features['AquBotDepth,m'] / features['TopDepth,m'])
                features['depthRatio2'] = (
                    features['AquBotDepth,m'] / features['MidDepth,m'])
                
                features['Pavg,bot-top'] = (
                    features['Pavg_bot,MPa'] - features['Pavg_top,MPa'])
                features['Pavg-Pini,bot'] = features['Pavg_bot,MPa'] - (
                    0.1013529 + features['AquBotDepth,m'] * 0.009785602) 
                
                features['LeakRateCO2,kg/s'] = inputarray['LeakRateCO2,kg/s'][nPred]
                features['LeakRateBrine,kg/s'] = inputarray['LeakRateBrine,kg/s'][nPred]
                features['LeakAqCumMassCO2,kg'] = inputarray['LeakAqCumMassCO2,kg'][nPred]
                features['LeakAqCumMassBrine,kg'] = inputarray['LeakAqCumMassBrine,kg'][nPred]
                
                features['CaseNum'] = 1 # will be removed in next 
                features['Sg_bot'] = np.maximum(inputarray['Sg_bot'][nPred],0.01)
                features['Sb_bot'] = 1.0-features['Sg_bot']
                
                shaleBrineden = ((features["rhoL_bot,kg/m3"] + features["rhoL_mid,kg/m3"]) 
                                 * 0.5)
                gammaB = (
                    ((features['Pavg_bot,MPa'] - features['Pavg_mid,MPa']) * 1e6) 
                    / (shaleBrineden * 9.785602 * features["AquThickness,m"]))
                
                shaleBrineden2 = (features["rhoL_bot,kg/m3"] 
                                  + features["rhoL_top,kg/m3"]) * 0.5
                gammaB2 = (
                    ((features['Pavg_bot,MPa'] - features['Pavg_top,MPa']) * 1e6) 
                    / (shaleBrineden2 * 9.785602 * (
                        features["AquThickness,m"] + features["ShaThickness,m"])))
                
                shaleCO2den = (features["rhoG_bot,kg/m3"] 
                               + features["rhoG_mid,kg/m3"]) * 0.5
                gammaC = (
                    ((features['Pavg_bot,MPa'] - features['Pavg_mid,MPa']) * 1e6) 
                    / (shaleCO2den * 9.785602 * features["AquThickness,m"]))
                
                shaleCO2den2 = (features["rhoG_bot,kg/m3"] 
                                + features["rhoG_top,kg/m3"]) * 0.5
                gammaC2 = (
                    ((features['Pavg_bot,MPa'] - features['Pavg_top,MPa']) * 1e6) 
                    / (shaleCO2den2 * 9.785602 * (
                        features["AquThickness,m"] + features["ShaThickness,m"])))
    
                features['gammaB'] = gammaB
                features['gammaC'] = gammaC
                features['gammaB2'] = gammaB2
                features['gammaC2'] = gammaC2
                
                shalePerm = (10.0 ** (features['PermWSh,mD'])) * 1e-15
                hydCondtopB = (shalePerm * features["rhoL_top,kg/m3"] 
                               * 9.785602 / features["visL_top,Pa-s"])
                hydCondmidB = (shalePerm  *features["rhoL_mid,kg/m3"] 
                               * 9.785602 / features["visL_mid,Pa-s"])
                hydCondbotB = (shalePerm * features["rhoL_bot,kg/m3"] 
                               * 9.785602 / features["visL_bot,Pa-s"])
                
                hydCondtopC = ((shalePerm * features["rhoG_top,kg/m3"] * 9.785602) 
                               / features["visG_top,Pa-s"])
                hydCondmidC = ((shalePerm * features["rhoG_mid,kg/m3"] * 9.785602) 
                               / features["visG_mid,Pa-s"])
                hydCondbotC = ((shalePerm * features["rhoG_bot,kg/m3"] * 9.785602) 
                               / features["visG_bot,Pa-s"])
    
                features['hydCondtopB'] = np.log10(hydCondtopB) 
                features['hydCondtopC'] = np.log10(hydCondtopC) 
                
                features['hydCondRatioB'] = hydCondtopB/hydCondbotB
                features['hydCondRatioC'] = hydCondtopC/hydCondbotC
                features['hydCondRatio2B'] = hydCondmidB/hydCondbotB
                features['hydCondRatio2C'] = hydCondmidC/hydCondbotC
                features['hydCondRatio3B'] = hydCondtopB/hydCondmidB
                features['hydCondRatio3C'] = hydCondtopC/hydCondmidC
                
                mobilityRatiotop = (
                    (features["rhoG_top,kg/m3"] / features["visG_top,Pa-s"]) 
                    / (features["rhoL_top,kg/m3"] / features["visL_top,Pa-s"]))
                mobilityRatiomid = (
                    (features["rhoG_mid,kg/m3"] / features["visG_mid,Pa-s"]) 
                    / (features["rhoL_mid,kg/m3"] / features["visL_mid,Pa-s"]))
                mobilityRatiobot = (
                    (features["rhoG_bot,kg/m3"] / features["visG_bot,Pa-s"]) 
                    / (features["rhoL_bot,kg/m3"] / features["visL_bot,Pa-s"]))
    
                features['mobRatiotop'] = np.log10(mobilityRatiotop)
                features['mobRatiomid'] = np.log10(mobilityRatiomid)
                features['mobRatiobot'] = np.log10(mobilityRatiobot)
                
                features['Brine_top,kg/s'] = inputarray['Brine_top,kg/s'][nPred]
                features['CO2_top,kg/s'] = inputarray['CO2_top,kg/s'][nPred] 
                features['SCO_top'] = inputarray['SCO_top'][nPred]
                
                features.drop(columns=['CaseNum'],inplace=True)
                
                return features 
            
            for nPred in range(self.nSL-2):
                feature_minthree = makingFeatures(
                    self.FluidModels, inputarray_minthree, nPred)
                feature_mintwo = makingFeatures(
                    self.FluidModels, inputarray_mintwo, nPred)
                fefature_mieone = makingFeatures(
                    self.FluidModels, inputarray_minone, nPred)
        
                feature_window_layer_i = np.concatenate([
                    feature_minthree, feature_mintwo, fefature_mieone], axis=1)
                
                if nPred == 0:
                    feature_window_layers = feature_window_layer_i.copy()
                else:
                    feature_window_layers = np.concatenate(
                        [feature_window_layers, feature_window_layer_i], axis=0)
           
            return feature_window_layers
        
    def get_mobility_for_shales(self):
        """ Calculate mobility parameter for brine and CO2 in the shale layers."""
        CO2RelPermShale = np.zeros(self.nSL)
        brineRelPermShale = np.ones(self.nSL)

        Sres = self.parameters.aquBrineResSaturation

        for ind in range(0, self.nSL):
            if self.interface[ind] > 0:
                CO2RelPermShale[ind] = min(
                    1 - Sres[ind], self.interface[ind] / self.thicknessH[ind])
                brineRelPermShale[ind] = (
                    1 / (1 - Sres[ind])) * (1 - Sres[ind] - CO2RelPermShale[ind])
            elif ind == self.nSL-1: # last shale is dummy
                CO2RelPermShale[ind] = min(
                    1-Sres[ind], self.interface[ind] / self.thicknessH[ind])

        CO2ViscosityShale = np.zeros(self.nSL)
        brineViscosityShale = np.ones(self.nSL)
        
        for ii in range(self.nSL):
            
            if ii == self.nSL-1:                
                CO2ViscosityShale[ii] = self.CO2ViscosityAquifer[ii]
                brineViscosityShale[ii] = self.brineViscosityAquifer[ii]
            else:
                CO2ViscosityShale[ii] = (
                    self.CO2ViscosityAquifer[ii] + self.CO2ViscosityAquifer[ii+1]
                    ) * 0.5
                brineViscosityShale[ii] = (
                    self.brineViscosityAquifer[ii] + self.brineViscosityAquifer[ii+1]
                    ) * 0.5
                                                
        self.CO2MobilityShale = CO2RelPermShale / CO2ViscosityShale
        self.brineMobilityShale = brineRelPermShale / brineViscosityShale

    def find_lambda(self):
        """ Calculate lambda constant."""
        # Lambda is a constant in the formula for the outer egde of the plume
        # It is constant for each aquifer in the system unless
        # brine and CO2 have different relative permeabilities, residual
        # saturations or viscosities
        # For simple calculations, the mobility ratio Lambda is calculated at
        # the endpoint relative permeability values for the two phases.

        Sres = self.parameters.aquBrineResSaturation

        # Maximum value for brine relative permeability
        brinePermeability = np.ones(self.nSL)

        # Maximum value for CO2 relative permeability
        CO2Permeability = np.zeros(self.nSL)
        for ind in range(self.nSL):
            CO2Permeability[ind] = 1 - Sres[ind]

        brineViscosity = self.brineViscosityAquifer
        CO2Viscosity = self.CO2ViscosityAquifer 

        self.Lambda = np.ones(self.nSL)*CO2Permeability*brineViscosity/(
            brinePermeability*CO2Viscosity)

    def get_interface(self):
        """ Calculate height of the interface between brine and CO2. """
        Sres = self.parameters.aquBrineResSaturation

        # Set only for aquifers (not including the reservoir layer)
        for j in range(1, self.nSL):

            if self.CO2Mass[j-1] > 0:
                x = 2 * np.pi * (self.thicknessH[j]) * (1 - Sres[j]) * 1.0 * (
                    self.parameters.wellRadius ** 2) / self.CO2Mass[j - 1]
                
                # ADDED THIS
                warnings.simplefilter("error", RuntimeWarning)
                
                try:
                    if x < 2/self.Lambda[j]:
                        tempInterface = (1 - Sres[j]) * self.thicknessH[j]
                    
                    elif x >= 2*self.Lambda[j]:
                        tempInterface = 0.
                    else:
                        tempInterface = (1 / (self.Lambda[j] - 1)) * (
                            np.sqrt(2 * self.Lambda[j] / x) - 1) * self.thicknessH[j]
                except:
                    tempInterface = 0.

                self.interface[j] = min(
                    [tempInterface, (1 - Sres[j]) * self.thicknessH[j], 
                     self.interface[j - 1]])

    def CC_shale(self):

        PermWellShale = self.parameters.shalePermeability
        ShaleThickness = self.parameters.shaleThickness

        CO2MobilityShale = self.CO2MobilityShale
        brineMobilityShale = self.brineMobilityShale

        CC = (PermWellShale / ShaleThickness) * (CO2MobilityShale + brineMobilityShale)

        return CC

    def GG_shale(self):

        PermWellShale = self.parameters.shalePermeability
        gg = self.g

        self.CO2DensityShale = np.zeros(self.nSL)
        self.brineDensityShale = np.zeros(self.nSL)
        
        for ii in range(self.nSL):
            
            if ii == self.nSL-1:
                self.CO2DensityShale[ii] = self.CO2DensityAquifer[ii]
                self.brineDensityShale[ii] = self.brineDensityAquifer[ii]
            else:
                self.CO2DensityShale[ii] = (
                    self.CO2DensityAquifer[ii] + self.CO2DensityAquifer[ii + 1]
                    ) * 0.5
                self.brineDensityShale[ii] = (
                    self.brineDensityAquifer[ii] + self.brineDensityAquifer[ii + 1]
                    ) * 0.5
        
        CO2MobilityShale = self.CO2MobilityShale
        brineMobilityShale = self.brineMobilityShale

        GG = PermWellShale * gg * (
            (CO2MobilityShale * self.CO2DensityShale) 
            + (brineMobilityShale * self.brineDensityShale))

        return GG

    def FF_aquifer(self):

        SCO2 = self.interface*(1/self.thicknessH)
        gg = self.g

        CO2Density = self.parameters.CO2Density
        brineDensity = self.parameters.brineDensity

        AquiferThickness = self.thicknessH

        FF = (AquiferThickness / 2) * gg * (
            SCO2 * CO2Density+(1 - SCO2) * brineDensity)

        return FF

    def WV_aquifer(self):

        SCO2 = self.interface * (1 / self.thicknessH)

        CO2MobilityAq = self.CO2MobilityAq
        brineMobilityAq = self.brineMobilityAq

        EffMobility = SCO2 * CO2MobilityAq + (1 - SCO2) * brineMobilityAq

        cf_ave = self.parameters.compressibility

        tt = self.parameters.timePoint * units.day()
        EffTime = 0.92 * tt

        rw = self.parameters.wellRadius
        
        # Reservoir perm not defined separately
        PermAquiferHor = np.append(self.parameters.aquiferPermeability[0], 
                                   self.parameters.aquiferPermeability)
        
        denominator = (4 * EffMobility * PermAquiferHor * EffTime)
        uLeak = np.zeros((len(denominator)))
        for i in range(len(denominator)):
            if denominator[i] == 0:
                uLeak[i] = 0
            else:
                uLeak[i] = ((rw ** 2) * cf_ave) / denominator[i]

        GLeak = well_function(uLeak)

        AquiferThickness = self.thicknessH
        
        denominator = (4 * np.pi * EffMobility * AquiferThickness * PermAquiferHor)
        WV = np.zeros((len(denominator)))
        for i in range(len(denominator)):
            if denominator[i] == 0:
                WV[i] = 0
            else:
                WV[i] = GLeak[i] / denominator[i]

        return WV

def read_data(filename):
    """ Routine used for reading pressure and saturation data files."""
    # Check whether the file with given name exists
    if os.path.isfile(filename):
        data = np.genfromtxt(filename)
        return data
    return None


def well_function(x):
    """
    Return the approximation of the well function with even number
    of terms in expansion. Expansions with an odd number of terms
    diverge to plus infinity without crossing zero.
    """
    W = np.zeros(len(x))
    for i, u in list(enumerate(x)):
        if u <= 1.0:
            demoninators = [
                (2 * scm.factorial(2)), (3 * scm.factorial(3)), 
                (4 * scm.factorial(4)), (5 * scm.factorial(5)), 
                (6 * scm.factorial(6)), (7 * scm.factorial(7)), 
                (8 * scm.factorial(8))]
            
            if 0 in demoninators:
                W[i] = np.inf
            else:
                # Use this statement to catch cases where np.log(u) encounters 
                # an issue (divide by zero in log) - this way, the warning is 
                # not printed.
                warnings.simplefilter("error", RuntimeWarning)
                try:
                    W[i] = (-0.577216 - np.log(u) + u 
                            - (u ** 2 / (2 * scm.factorial(2))) 
                            + (u ** 3 / (3 * scm.factorial(3))) 
                            - (u ** 4 /(4 * scm.factorial(4))) 
                            + (u ** 5 / (5 * scm.factorial(5))) 
                            - (u ** 6 / (6 * scm.factorial(6))) 
                            + (u ** 7 / (7 * scm.factorial(7))) 
                            - (u ** 8 / (8 * scm.factorial(8))))
                except:
                    # For cases where np.log(u) encounters an error (divide by zero in log)
                    W[i] = np.inf
            
        elif u <= 9:
            uu = np.linspace(1.0, 9.0, num=9)
            Wu = np.array([0.219, 0.049, 0.013, 0.0038, 0.0011,
                           0.00036, 0.00012, 0.000038, 0.000012])
            Wfun = interp1d(uu, Wu, kind='linear')
            W[i] = Wfun(u)
        else:
            W[i] = 0.000001
    return W

if __name__ == "__main__":
    xx = [5.0, 6, 4]
    w_val = well_function(xx)
    print(w_val)
