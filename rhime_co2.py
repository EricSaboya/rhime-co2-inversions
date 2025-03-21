# ---------------------------------------------------------------------------------------
# rhime.py
# Created 23 July 2024 
# Author: Eric Saboya, School of Geographical Sciences, University of Bristol 
# ---------------------------------------------------------------------------------------
# RHIME CO2 inverse models
# > rhime_inversions (mole fraction sectoral inversions)
# 
# 
# ---------------------------------------------------------------------------------------

import os
import sys
import pickle
import numpy as np
import xarray as xr
from pathlib import Path
from typing import Optional
from collections import namedtuple

sys.path.append("/user/work/wz22079/projects/rhime-co2-inversions/")
import utils
from model_error_methods import model_error_method_parser
import inversion_setup as setup
import calculate_basis_functions as cbf
from get_co2_data import get_mf_obs_sims

#from tracers_co2 import get_14c_obs_sims
from inversion_mcmc import inferpymc, inferpymc_postprocessouts


def rhime_inversions(obs_dict: dict = None,
                     flux_dict: dict = None,
                     bc_dict: dict = None,
                     fp_dict: dict = None,
                     basis_dict: dict = None,
                     mcmc_dict: dict = None,
                     use_bc: bool = True,
                     model_error_method: str = None, 
                     sigma_freq: str = None,
                     outputname: str = None,
                     outputpath: str = None,
                     country_file: str = None,
                    ):
    """
    -------------------------------------------------------
    Regional Hierarchical Inverse Modelling Environment 

    RHIME is a regional hierarchical Bayesian inverse 
    model designed for inferrring fluxes from atmospheric
    trace gas measurements and a priori flux data.
    -------------------------------------------------------
    Args:
        obs_dict (dict | defaults to None):
            Dictionary of parameters for retrieving 
            observations.
        flux_dict (dict | defaults to None):
            Dictionary of parameters for retrieving
            fluxes.
        bc_dict (dict | defaults to None):
            Dictionary of parameters for retrieving
            boundary conditions.
        fp_dict (dict | defaults to None):
            Dictionary of parameters for retreiving 
            footprints.   
        basis_dict (dict | defaults to None):
            Dictionary of basis function parameters
        mcmc_dict (dict | defaults to None):
            Dictionary of MCMC function parameters
        use_bc (bool | defaults to True):
            Option to use/solve for boundary conditions 
            values in the inversion.
        model_error_method (str | defaults to None):
            Specify model error calculation method.
            One of:
            > "residual"
            > "fixed"
            > "None" 
        sigma_freq (str | defaults to None):
            Period over which the model error is estimated.
        bc_freq (str | defaults to None):
            Period over which the baseline is estimated.
            > "monthly": estimates BCs over one calendar month 
            > "30D": estimates BCs over 30-days
            > None: one BC for entire inversion period 
        outputname (str | defaults to None):
            Unique identifier for output/run name
        outputpath (str | defaults to None):
            Path of where output should be saved
        country_file (str | defaults to None):
            Path to countryfile mask to be used 
    -------------------------------------------------------    
    """
    # Get CO2 obs and create forward simulations 
    # (w./ Hall term for each sector)
    (data_dict, 
     sites, 
     inlet, 
     fp_height, 
     instrument, 
     averaging_period) = get_mf_obs_sims(flux_dict, 
                                         fp_dict, 
                                         obs_dict, 
                                         bc_dict, 
                                         use_bc,
                                        )
   
    # Update site parameters to remove sites with 
    # no data available during period of interest 
    if sites is None:
        sites = obs_dict["site"]
        inlet = obs_dict["inlet"]
        fp_height = fp_dict["fp_height"]
        instrument = obs_dict["instrument"]
        averaging_period = obs_dict["averaging_period"]

    basis_dict["site"] = sites 
    basis_dict["source"] = flux_dict["source"]
    basis_dict["domain"] = fp_dict["domain"]
    basis_dict["start_date"] = obs_dict["start_date"]
                
    # Calculate basis functions for each flux sector
    (fp_data, 
     tempdir, 
     basis_dir, 
     bc_basis_dir) = cbf.basis_functions_wrapper(data_dict,
                                                 basis_dict,
                                                 use_bc=use_bc,
                                                 outputname=outputname,
                                                 outputpath=outputpath,
                                                )

    # Apply data filtering
    if obs_dict["filters"] is not None:
        fp_data = utils.filtering(fp_data, obs_dict["filters"])

    # Calculate model errors
    fp_data = model_error_method_parser(fp_data, 
                                        model_error_method,
                                       )

    # Remove any sites that return empty data array post-filtering
    s_dropped = []
    for site in sites:
        if fp_data[site].time.values.shape[0] == 0:
            s_dropped.append(site)
            del fp_data[site]
            
    if len(s_dropped) != 0:
        sites = [s for i, s in enumerate(sites) if s not in s_dropped]
        print(f"\nDropping {s_dropped} sites as no data passed the filtering.\n")

    # Append model domain region to site attributes
    for i, site in enumerate(sites):
        fp_data[site].attrs["Domain"] = fp_dict["domain"]
            
    # Mole fraction multi-sector inversions 
    error = np.zeros(0)    # Observational uncertainty 
    Hbc = np.zeros(0)      # Basis function dosage for model domain boundary [region]
    Hx = np.zeros(0)       # Basis function dosage [region[sector], t]
    Hxerr = np.zeros(0)    # Coefficient of variability of the dosage in each basis function [region[sector], t]
    Y = np.zeros(0)
    Ymodelerror = np.zeros(0) # Model error 
    siteindicator = np.zeros(0)

    for i, site in enumerate(sites):
        # Select variables to drops NaNs from 
        drop_vars = []
        for var in ["H", "Herr", "H_bc", "mf", "mf_variability", "mf_repeatability", "y_model_err"]:
            if var in fp_data[site].data_vars:
                drop_vars.append(var)
                        
        # PyMC does not like NaNs, so drop them for the variables used below
        fp_data[site] = fp_data[site].dropna("time", subset=drop_vars)
        
        # Propagate mole fraction observational uncertainties for repeatability and variability
        # these terms should be added in quadrature
        myerror = np.zeros(0)
        if "mf_repeatability" in fp_data[site]:
            if len(myerror)==0:
                myerror = np.concatenate((myerror, fp_data[site]["mf_repeatability"].values**2))                        
            elif len(myerror) == len(fp_data[site]["mf_repeatability"].values):
                myerror += fp_data[site]["mf_repeatability"].values**2
            else:
                raise(f"Error array length does not match length of mole fraction repeatability values for {site}.")
        if "mf_variability" in fp_data[site]:
            if len(myerror)==0:
                myerror = np.concatenate((myerror, fp_data[site]["mf_variability"].values**2))
            elif len(myerror) == len(fp_data[site]["mf_variability"].values):
                myerror += fp_data[site]["mf_variability"].values**2
            else:
                raise(f"Error array length does not match length of mole fraction variability values for {site}.")  
        error = np.concatenate((error, np.sqrt(myerror)))
        
        # Concatenate observational mole fractions for each site to Y 
        Y = np.concatenate((Y, fp_data[site]["mf"].values))
        Ymodelerror = np.concatenate((Ymodelerror, fp_data[site]["y_model_err"].values))
        siteindicator = np.concatenate((siteindicator, np.ones_like(fp_data[site].mf.values) * i))
                
        if i == 0:
            Ytime = fp_data[site]["time"].values
        else:
            Ytime = np.concatenate((Ytime, fp_data[site]["time"].values))
            
        if use_bc is True:
            bc_freq = bc_dict["bc_freq"]
            
            if bc_freq == "monthly":
                Hmbc = setup.monthly_bcs(obs_dict["start_date"], 
                                         obs_dict["end_date"], 
                                         site, 
                                         fp_data,
                                        )
            elif bc_freq is None:
                Hmbc = fp_data[site]["H_bc"].values
            else:
                Hmbc = setup.create_bc_sensitivity(obs_dict["start_date"], 
                                                   obs_dict["end_date"], 
                                                   site, 
                                                   fp_data,
                                                   bc_freq,
                                                  )
        elif use_bc is False:
            Hmbc = np.zeros(0)
            
        if i == 0:
            Hbc = np.copy(Hmbc)
            Hx = fp_data[site].H.values
            Hxerr = fp_data[site].Herr.values
        else:
            Hbc = np.hstack((Hbc, Hmbc))
            Hx = np.hstack((Hx, fp_data[site].H.values))
            Hxerr = np.hstack((Hxerr, fp_data[site].Herr.values))
    
    # Mask source regions in Hx
    basis_region_mask = np.zeros_like(fp_data[site]["region"].values)
    count = 0 
    for emi in flux_dict["source"]:
        count += 1
        for i in range(len(fp_data[site]["region"].values)):
            if emi in fp_data[site]["region"].values[i]:
                basis_region_mask[i] = count
    basis_region_mask = basis_region_mask.astype(int)

    sigma_freq_index = setup.sigma_freq_indicies(Ytime, sigma_freq)

    # Aligning BC units 
    # if obs_dict["species"].lower()=="co2":
    #     if Hbc.mean() < 1e-2:
    #         Hbc = Hbc * 1e6
    
        
    # Run PyMC MCMC inversion (no tracer)
    print("******************************************")
    print("Running MCMC inversion without tracer ... ")
    print("******************************************")
    mcmc_results = inferpymc(Hx=Hx,
                             Hxerr=Hxerr,
                             basis_region_mask=basis_region_mask,
                             Y=Y,
                             error=error, 
                             Ymodelerror=Ymodelerror,
                             siteindicator=siteindicator,
                             sigma_freq_index=sigma_freq_index,
                             Hbc=Hbc,
                             xprior=mcmc_dict["xprior"],
                             bcprior=mcmc_dict["bcprior"],
                             sigprior=mcmc_dict["sigprior"],
                             nit=int(mcmc_dict["nit"]),
                             burn=int(mcmc_dict["burn"]),
                             tune=int(mcmc_dict["tune"]),
                             nchain=int(mcmc_dict["nchain"]),
                             sigma_per_site=True,
                             offsetprior=mcmc_dict["offsetprior"],
                             add_offset=mcmc_dict["add_offset"],
                             verbose=False,
                             save_trace=False,
                             use_bc=use_bc,
                            )

    print("**************************************")
    print("Processing MCMC inversion outputs ... ")
    print("**************************************")
    inferpymc_postprocessouts(mcmc_results=mcmc_results,
                              use_bc=use_bc,
                              mcmc_dict=mcmc_dict,
                              Hx=Hx,
                              Y=Y,
                              error=error,
                              Ymodelerror=Ymodelerror, 
                              Ytime=Ytime,
                              siteindicator=siteindicator,
                              sigma_freq_index=sigma_freq_index,
                              domain=fp_dict["domain"],
                              species=obs_dict["species"],
                              sites=sites,
                              start_date=obs_dict["start_date"],
                              end_date=obs_dict["end_date"],
                              outputname=outputname,
                              outputpath=outputpath,
                              country_unit_prefix=None,
                              emissions_name=flux_dict["source"],
                              emissions_store=flux_dict["store"],
                              Hbc=Hbc,
                              fp_data=fp_data,
                              country_file=country_file,
                              rerun_file=None,
                             )

    print("==== INVERSION COMPLETE ====")

# ------------------------------------------------------------------------------# 


