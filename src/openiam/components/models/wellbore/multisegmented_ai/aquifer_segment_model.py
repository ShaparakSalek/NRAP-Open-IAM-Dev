# -*- coding: utf-8 -*-
"""
DLAquiferSegmentModels object used by the Multisegmented Wellbore AI component.
"""
import os


class DLAquiferSegmentModels():
    """        
    S.Baek & Maruti working on. Aquifer segment model
    """
    
    def __init__(self, componentPath):
        
        import time
        import logging
        import joblib
        
        TargetVars = ['BrineRegression', 
                      'CO2Regression', 
                      'SaturationRegression',]
       
        t1 = time.time()
        for TargetVar in TargetVars:
            t0 = time.time()
        
            if TargetVar == 'BrineRegression':
                self.model_rf_qB = joblib.load(open(os.path.join(
                    componentPath, TargetVar, 'Brine_RF_lb_3_depth_15_132.joblib'), 'rb'))
            elif TargetVar == 'CO2Regression':
                self.model_rf_qC = joblib.load(open(os.path.join(
                    componentPath, TargetVar, 'CO2_RF_lb_3_depth_15_132.joblib'), 'rb'))
            elif TargetVar == 'SaturationRegression':
                self.model_rf_sC = joblib.load(open(os.path.join(
                    componentPath, TargetVar, 'RFmodel_lb_3_depth_13_132.joblib'), 'rb'))
                
            # print('* loading dl models: {:.0f} secs'.format(time.time()-t0))    
        
        # print('** total loading dl models: {:.0f} secs'.format(time.time()-t1))  
