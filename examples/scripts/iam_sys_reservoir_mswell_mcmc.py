'''
This example performs MCMC on the analytical reservoir model to produce
posterior distributions for the reservoir permeability (logResPerm),
reservoir porosity (reservoirPorosity) and reservoir thickness
(reservoirThickness) - parameters constrained by measured values of pressure
and |CO2| saturation at leaking well. The posterior parameter values are
then run through the multisegmented wellbore model in order to propagate
the uncertainty to the brine and |CO2| leakage from the same well
for default multisegmented wellbore component.

Example of run:
$ python iam_sys_reservoir_mswell_mcmc.py
'''

import sys
import os
import matplotlib.pyplot as plt
import numpy as np

from openiam.components.iam_base_classes import SystemModel
from openiam.components.analytical_reservoir_component import AnalyticalReservoir
from openiam.components.multisegmented_wellbore_component import MultisegmentedWellbore


if __name__=='__main__':
    # For multiprocessing in Spyder
    __spec__ = None
    # Define keyword arguments of the system model
    num_years = 5
    time_array = 365.25*np.arange(0.0, num_years+1)
    sm_model_kwargs = {'time_array': time_array}  # time is given in days
    nwalkers = 10  # Number of MCMC chains
    nsamples_per_walker = 100  # Number of samples per chain (walker)
    burnin = 20  # Number of samples to throw away at the beginning of chains

    # Create system model
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    # Add reservoir component model
    ares = sm.add_component_model_object(AnalyticalReservoir(name='ares', parent=sm))

    # Add parameters of reservoir component model
    ares.add_par('logResPerm', min=-13.5, max=-12.5, value=-13.)
    ares.add_par('reservoirPorosity', min=0.1, max=0.3, value=0.25)
    ares.add_par('reservoirThickness', min=20., max=40., value=30)

    # Add observations of reservoir component model to be used by the next component
    ares.add_obs('pressure')  # by default observations at all time points are added
    ares.add_obs('CO2saturation')

    # Create synthetic observations by running the model forward with
    # specified parameter 'values' defined above in add_par calls.
    # Observation data will be read in from files in the future
    sm.forward()
    pres_list = []
    sat_list = []
    for ind in range(len(time_array)):
        pres_list.append(ares.obs['pressure'+'_'+str(ind)].sim)
        sat_list.append(ares.obs['CO2saturation'+'_'+str(ind)].sim)

    # index + value option should be used if one is interested in analysis of
    # particular data points. User can also specify what are the measured
    # values of observations to compare them versus simulated.
    ares.add_obs('pressure', index=list(range(num_years+1)), value=pres_list)
    ares.add_obs('CO2saturation', index=list(range(num_years+1)), value=sat_list)
    ares.add_obs_to_be_linked('pressure')
    ares.add_obs_to_be_linked('CO2saturation')

    # Run MCMC to produce parameter posterior distributions based on consistency
    # with ares model observations 'pressure' and 'CO2saturation'
    ss = sm.emcee(nwalkers=nwalkers, nsamples=nsamples_per_walker, burnin=burnin)

    # Create trace plots. This is done somewhat manually here, but could
    # be easily encapsulated in a SystemModel method
    # Perm plots
    f1, ax = plt.subplots(3, 1, sharex=True, figsize=(10, 10))
    ax[0].plot(ss.chain[:, burnin:, 0].T, color='k', alpha=0.4)
    ax[0].axhline(-13, color='gray', lw=2)
    ax[0].set_ylabel('logResPerm')
    ax[0].set_ylim(-14, -12)
    # Porosity plots
    ax[1].plot(ss.chain[:, burnin:, 1].T, color='k', alpha=0.4)
    ax[1].axhline(0.25, color='gray', lw=2)
    ax[1].set_ylabel('Porosity')
    ax[1].set_ylim(0.1, 0.3)
    # Thickness plots
    ax[2].plot(ss.chain[:, burnin:, 2].T, color='k', alpha=0.4)
    ax[2].axhline(0.3, color='gray', lw=2)
    ax[2].set_xlabel('Step number')
    ax[2].set_ylabel('reservoirThickness')
    ax[2].set_ylim(20., 40.)
    f1.savefig('mcmc_traces.png')

    # Create panels plot to explore the histograms and correlations of
    # posterior parameter distributions
    sc = sm.create_sampleset(ss.chain[:, burnin:, :].reshape((-1, len(sm.pars))))
    sc.panels(bins=50, figsize=(10, 10), fontsize=16, ms=1,
              filename='mcmc_posteriors_panels.png')

    # Add multisegmented wellbore component so that reservoir uncertainty
    # can be propagated through a "default" wellbore to estimate
    # the uncertainty in CO2 and brine leakage.
    ms = sm.add_component_model_object(MultisegmentedWellbore(name='ms', parent=sm))

    # Add parameters of wellbore component linked to the parameters
    # of the same name from reservoir model. In this example all wellbore
    # component parameters are default parameters, i.e., linked to the
    # corresponding parameters of the reservoir model
    ms.add_par_linked_to_par('numberOfShaleLayers', ares.default_pars['numberOfShaleLayers'])
    ms.add_par_linked_to_par('shale1Thickness', ares.default_pars['shaleThickness'])
    ms.add_par_linked_to_par('shale2Thickness', ares.default_pars['shaleThickness'])
    ms.add_par_linked_to_par('shale3Thickness', ares.default_pars['shaleThickness'])
    ms.add_par_linked_to_par('aquifer1Thickness', ares.default_pars['aquiferThickness'])
    ms.add_par_linked_to_par('aquifer2Thickness', ares.default_pars['aquiferThickness'])
    ms.add_par_linked_to_par('reservoirThickness', ares.pars['reservoirThickness'])
    ms.add_par_linked_to_par('datumPressure', ares.default_pars['datumPressure'])

    # Add keyword arguments linked to the output provided by reservoir model
    ms.add_kwarg_linked_to_obs('pressure', ares.linkobs['pressure'])
    ms.add_kwarg_linked_to_obs('CO2saturation', ares.linkobs['CO2saturation'])
    ms.add_obs('brine_aquifer1')
    ms.add_obs('CO2_aquifer1')

    # Run posterior parameter distributions with multisegmented wellbore model
    # added to the system model
    sc.run(cpus=4, verbose=False)

    # Run prior (uniform) distribution through coupled model to determine
    # prior uncertainty in leakage. The nsamples equals the number of
    # mcmc samples minus the burnin for each chain
    sprior = sm.lhs(siz=nwalkers*nsamples_per_walker - nwalkers*burnin)
    sprior.run(cpus=2, verbose=False)

    # Plot resulting CO2 and brine leakage histograms to explore which show
    # the resulting leakage uncertainty due to reservoir permeability,
    # porosity, and thickness uncertainty
    f2, ax = plt.subplots(1, 2, figsize=(12, 6))
    # '_4' is used to indicate the index of time point of interest
    ax[0].hist(sprior.recarray['ms.brine_aquifer1_4'], bins=30, label='Prior')
    ax[0].hist(sc.recarray['ms.brine_aquifer1_4'], bins=30, label='Posterior')
    ax[0].legend()
    ax[1].hist(sprior.recarray['ms.CO2_aquifer1_4'], bins=30, label='Prior')
    ax[1].hist(sc.recarray['ms.CO2_aquifer1_4'], bins=30, label='Posterior')
    ax[1].legend()

    ax[0].set_ylabel("Count")
    ax[0].set_xlabel(r"Brine leakage rate, [kg/s]")
    ax[1].set_xlabel(r"CO$_2$ leakage rate, [kg/s]")
    f2.show()
    f2.savefig('mcmc_leakage_panels.png')
