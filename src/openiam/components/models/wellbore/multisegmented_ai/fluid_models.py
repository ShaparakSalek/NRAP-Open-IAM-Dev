# -*- coding: utf-8 -*-
"""
FluidModels object used by the Multisegmented Wellbore AI component.
"""
import os
import numpy as np
import pandas as pd

class FluidModels():
    """    
    Fluid models to calculate properties of CO2 and brine as a function of 
    pressure, temperature and salinity (for brine only). 
    
    Models were developed using a random forest machine learning algorithm and 
    random sampled data based on American Society of Mechanical Engineers steam
    table formulations (Meyer et al. 1993) for brine and the equation of state 
    for CO2 by Span and Wagner (1996) and Fenghour et al.(1998)
    """
    def __init__(self, BasicPath):
        
        from pickle import load
        import joblib
        
        self.BasicPath0 = BasicPath
       
        self.BrineFluidModel = joblib.load(
            os.path.join(self.BasicPath0, 'BrineFluid', 
                         'MultipleFluidProps_Brine.joblib'))
        self.Brine_scalerX = load(
            open(os.path.join(self.BasicPath0, 'BrineFluid', 
                              'scaler_x1_scaler.pkl'), 'rb'))
        self.Brine_scalerY = load(
            open(os.path.join(self.BasicPath0, 'BrineFluid', 
                              'scaler_y1_scaler.pkl'), 'rb'))
        
        import CoolProp.CoolProp
        self.CO2FluidModelTesting = CoolProp.CoolProp

    def CO2Model(self,Temp_C,P_MPA):
           
        FluidProps = pd.DataFrame(columns=["Temp,C", "Press,MPa"])    
        FluidProps["Temp,C"]    = Temp_C
        FluidProps["Press,MPa"] = P_MPA
        
        yPred = np.ndarray((FluidProps.shape[0],2))
        for ii in range(len(FluidProps["Temp,C"])):
            viscosity_ii = self.CO2FluidModelTesting.PropsSI(
                "V", "T", FluidProps["Temp,C"][ii] + 273.15, "P", 
                FluidProps["Press,MPa"][ii] * 1e6, 'CarbonDioxide') # "Pa-s", viscosity
            density_ii = self.CO2FluidModelTesting.PropsSI(
                "D", "T", FluidProps["Temp,C"][ii]+273.15, "P", 
                FluidProps["Press,MPa"][ii]*1e6, 'CarbonDioxide') # "kg/m^3", density
            yPred[ii,:] = [density_ii,viscosity_ii]
           
        return yPred
       
    def BrineModel(self, Temp_C, P_MPA, Salinity):
        
        FluidProps = pd.DataFrame(columns = ["Temp,C", "Press,MPa", "Salinity"])    
        FluidProps["Temp,C"]    = Temp_C
        FluidProps["Press,MPa"] = P_MPA
        FluidProps["Salinity"] = Salinity
        
        x_data_bot    = self.Brine_scalerX.transform(FluidProps)
        y_reduced_bot = self.BrineFluidModel.predict(x_data_bot)
        yPred         = self.Brine_scalerY.inverse_transform(y_reduced_bot)
        yPred[:,1] *= 1e-3 # viscosity unit conversion
        
        return yPred
