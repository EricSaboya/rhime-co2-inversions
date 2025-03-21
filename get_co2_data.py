# ---------------------------------------------------------------------------------------
# get_co2_data.py
# Created: 15 May 2024
# Author: Eric Saboya, School of Geographical Sciences, University of Bristol
# ---------------------------------------------------------------------------------------
# Functions for retrieving mole fraction datasets and creating simulations 
# ---------------------------------------------------------------------------------------

import os
import sys
import pint
import glob
import numpy as np
import pandas as pd
import xarray as xr
import openghg
import datetime as dt

from openghg.types import SearchError
from openghg.retrieve import get_flux, get_bc, get_obs_surface, get_footprint

sys.path.append("/user/work/wz22079/projects/rhime-co2-inversions/")

from data_checks import create_unit_registry, check_obs_units


def get_fps_data(fp_dict: dict)->dict:
    """
    -------------------------------------------------------
    Retrieve footprints for each site
    -------------------------------------------------------
    Args:
        fp_dict (dict): 
            Dictionary of footprints data specifications
            Keys: 'site', 
                  'fp_height', 
                  'domain', 
                  'species',
                  'start_date', 
                  'end_date', 
                  'store'
    
    Returns:
        Dictionary of OpenGHG footprint objects for each 
        site
    -------------------------------------------------------
    """
    fps_out_dict = {}
    sites = fp_dict["site"]

    for i, site in enumerate(sites):
        try:
            get_fps_data = get_footprint(site=fp_dict["site"][i],                            
                                         height=fp_dict["fp_height"][i],
                                         domain=fp_dict["domain"],
                                         start_date=fp_dict["start_date"],
                                         end_date=fp_dict["end_date"],
                                         store=fp_dict["store"],
                                         species=fp_dict["species"],
                                        )
            fps_out_dict[site] = get_fps_data
            
        except SearchError:
            print(f"\nNo footprint data found for {site}.\n")
            continue  # skip this site
    return fps_out_dict

def get_flux_data(flux_dict: dict)->dict:
    """
    -------------------------------------------------------
    Retrieve flux fields and apply any scaling factors to 
    the flux values
    -------------------------------------------------------
    Args:
        flux_dict (dict)
            Dictionary of flux data specifications
            Keys: 'species', 
                  'source', 
                  'domain', 
                  'flux_sf',
                  'start_date', 
                  'end_date', 
                  'store'
    
    Returns:
        Dictionary of OpenGHG flux data objects for each
        source 
    -------------------------------------------------------
    """
    flux_data_dict = {}
    for i, source in enumerate(flux_dict["source"]):
        get_flux_data = get_flux(species=flux_dict["species"],
                                 domain=flux_dict["domain"],
                                 source=source,
                                 start_date=flux_dict["start_date"],
                                 end_date=flux_dict["end_date"],
                                 store=flux_dict["store"],
                                )
        
        # Apply any scale factors to flux fields 
        if "flux_sf" in list(flux_dict.keys()): 
            if flux_dict["flux_sf"] is not None:
                if type(flux_dict["flux_sf"]) is dict:
                    print(f"Warning! Applying a multiplicative scale factor of {flux_dict["flux_sf"][source]} to {source}.")
                    get_flux_data.data.flux.values *= flux_dict["flux_sf"][source]
                elif type(flux_dict["flux_sf"]) in [int, float]:
                    print(f"Warning! Applying a multiplicative scale factor of {flux_dict["flux_sf"]} to {source}.")
                    get_flux_data.data.flux.values *= flux_dict["flux_sf"]
                elif type(flux_dict["flux_sf"]) is list:
                    print(f"Warning! Applying a multiplicative scale factor of {flux_dict["flux_sf"][i]} to {source}.")
                    get_flux_data.data.flux.values *= flux_dict["flux_sf"][i]
                else:
                    raise KeyError("Use either a dict, float or list for flux scale factors.")
        flux_data_dict[source] = get_flux_data        
    return flux_data_dict

def get_mf_bc_data(bc_dict: dict):
    """
    -------------------------------------------------------
    Retrieve boundary conditions data and apply any 
    scaling factors to these data
    -------------------------------------------------------
    Args:
        bc_dict (dict):
            Dictionary of BC data specifications
            Keys: 'species', 
                  'domain', 
                  'bc_input', 
                  'bc_freq', 
                  'bc_sf', 
                  'store'
                  'start_date', 
                  'end_date', 
    
    Returns:
        OpenGHG BC object for time of interest
    -------------------------------------------------------
    """
    get_bc_data = get_bc(species=bc_dict["species"],
                         domain=bc_dict["domain"],
                         bc_input=bc_dict["bc_input"],
                         start_date=bc_dict["start_date"],
                         end_date=bc_dict["end_date"],
                         store=bc_dict["store"],
                        )
    
    # Apply any scale factors to BCs
    if "bc_sf" in list(bc_dict.keys()):
        if bc_dict["bc_sf"] is not None:
            if type(bc_dict["bc_sf"]) is dict:
                for key in bc_dict["bc_sf"].keys():
                    print(f"Warning! Applying a multiplicative scale factor of {bc_dict["bc_sf"][key]} to {key}.")
                    get_bc_data.data[key].values *= bc_dict["bc_sf"][key]
            elif type(bc_dict["bc_sf"]) in [int, float]:
                for key in get_bc_data.data.keys():
                    print(f"Warning! Applying a multiplicative scale factor of {bc_dict["bc_sf"]} to {key}.")
                    get_bc_data.data[key].values *= bc_dict["bc_sf"]
            else:
                raise KeyError("Use a dict or float value for the BC scale factors.")
    return get_bc_data

def get_mf_obs(obs_dict: dict,
               dataset_source: str,
              )->dict:
    """
    -------------------------------------------------------
    Retrieves surface observation mole fractions and
    calculates the associated observational uncertainty
    -------------------------------------------------------
    Args:
       obs_dict (dict):
            Dictionary of observations data specifications
            Keys: 'site', 
                  'inlet', 
                  'species', 
                  'instrument',
                  'start_date', 
                  'end_date', 
                  'calibration_scale',
                  'data_level', 
                  'store', 
                  'averaging_period', 
                  'filters',
                  
       dataset_source (str):
           Specific key word used for PARIS 
           pseudo-observations

    Returns:
        Dictionary of OpenGHG surface observation 
        objects for each measurement station of interest
    -------------------------------------------------------
    """
    mf_obs_out_dict = {}
    sites = obs_dict["site"]
    
    for i, site in enumerate(sites):
        try:
            # Retrieves AVERAGED mole fractions data
            site_data = get_obs_surface(site=obs_dict["site"][i],
                                        species=obs_dict["species"],
                                        inlet=obs_dict["inlet"][i],
                                        start_date=obs_dict["start_date"],
                                        end_date=obs_dict["end_date"],
                                        icos_data_level=obs_dict["data_level"][i],
                                        average=obs_dict["averaging_period"][i],
                                        instrument=obs_dict["instrument"][i],
                                        calibration_scale=obs_dict["calibration_scale"],
                                        store=obs_dict["store"],
                                        keep_missing=True,
                                        dataset_source=dataset_source,
                                       )
            
            # Calculate observational uncertainty over averaging period 
            site_data_no_ave = get_obs_surface(site=obs_dict["site"][i],
                                               species=obs_dict["species"],
                                               inlet=obs_dict["inlet"][i],
                                               start_date=obs_dict["start_date"],
                                               end_date=obs_dict["end_date"],
                                               icos_data_level=obs_dict["data_level"][i],
                                               instrument=obs_dict["instrument"][i],
                                               calibration_scale=obs_dict["calibration_scale"],
                                               store=obs_dict["store"],
                                               keep_missing=True,
                                               dataset_source=dataset_source,
                                              )

            # Observation uncertainty is defined as the sum in quadrature of:
            # 1. The root-squared sum of the 1 min mf variabilities reported in the 
            #    uploaded data in the averaging period and divided by the number of 
            #    observations in the averaging period. This represents the "error" of 
            #    the sampled mean.
            #
            # 2. The variance of the mole fraction concentrations reported in the 
            #    uploaded data over the averaging period. This accounts for the 
            #    variability in mole fractions over the averaging period.

            
            # TO CHECK: variability for optical instruments, repeatability for ECD?
            if "mf_variability" in list(site_data_no_ave[site].keys()) and "mf_repeatability" in list(site_data_no_ave[site].keys()):
                # Sampled mean error of mf_variability (not averaged)
                ds_resampled_rep = np.sqrt((site_data_no_ave[site]["mf_variability"]**2).resample(time=obs_dict["averaging_period"][i]).sum())/site_data_no_ave[site]["mf_variability"].resample(time=obs_dict["averaging_period"][i]).count()

                # Propagate sampled mean error term of mf_variability with mf_repeatability
                site_data[site]["mf_repeatability"] = np.sqrt(site_data[site]["mf_repeatability"] **2 + ds_resampled_rep **2)
                site_data[site]["mf_variability"] = site_data[site]["mf_variability"]**2 
               
            elif "mf_variability" in list(site_data_no_ave[site].keys()) and "mf_repeatability" not in list(site_data_no_ave[site].keys()):
                # In this instance, there is no repeatability. We use the sampled mean error of the variability in lieu
                ds_resampled_rep = np.sqrt((site_data_no_ave[site]["mf_variability"]**2).resample(time=obs_dict["averaging_period"][i]).sum())/site_data_no_ave[site]["mf_variability"].resample(time=obs_dict["averaging_period"][i]).count()
                
                site_data[site]["mf_variability"] = site_data[site]["mf_variability"]**2
                site_data[site]["mf_repeatability"] = ds_resampled_rep

            elif "mf_variability" not in list(site_data_no_ave[site].keys()) and "mf_repeatability" in list(site_data_no_ave[site].keys()):
                site_data[site]["mf_variability"] = site_data[site]["mf_variability"]**2

            elif "mf_variability" not in list(site_data_no_ave[site].keys()) and "mf_repeatability" not in list(site_data_no_ave[site].keys()):
                print(f"{site} mole fractions have no associated uncertainties. Using 0.1% of mf value as the repeatability")
                pseudo_err = site_data[site]["mf"] * 0.001 
                # Calculates the sampled mean error from pseudo uncertainties 
                ds_resampled_rep = np.sqrt((pseudo_err**2).resample(time=obs_dict["averaging_period"][i]).sum())/site_data_no_ave[site]["mf"].resample(time=obs_dict["averaging_period"][i]).count()

                site_data[site]["mf_repeatability"] = ds_resampled_rep
                site_data[site]["mf_variability"] = site_data[site]["mf_variability"]**2 # Variance of mf over averaging period


            # Check units and convert to mol/mol if in ppm 
            site_data = check_obs_units(site_data)

            # Add data to output dictionary
            mf_obs_out_dict[site] = site_data

        except SearchError:
            print(f"\nNo obs data found for {site} \n")
            continue  # skip this site
        except AttributeError:
            print(f"\nNo data found for {site} between {obs_dict["start_date"]} and {obs_dict["end_date"]}.\n")
            continue  # skip this site
        else:
            if site_data is None:
                print(f"\nNo data found for {site} between {obs_dict["start_date"]} and {obs_dict["end_date"]}.\n")
                continue  # skip this site
    return mf_obs_out_dict

def get_mf_obs_sims(flux_dict: dict, 
                    fp_dict: dict, 
                    obs_dict: dict, 
                    bc_dict: dict, 
                    use_bc: bool,
                   ):
    """
    -------------------------------------------------------
    Function to produce forward simulations of 
    CO2 mole fractions. 
    -------------------------------------------------------
    Args:
        flux_dict (dict): 
            Dictionary of flux data specifications
            Keys: 'species', 
                  'source', 
                  'domain', 
                  'start_date', 
                  'end_date',
                  'store'
                  
        fp_dict (dict): 
            Dictionary of footprints data specifications
            Keys: 'site', 
                  'fp_height', 
                  'domain', 
                  'species',
                  'start_date', 
                  'end_date', 
                  'store'
                  
        obs_dict (dict): 
            Dictionary of observations data specifications
            Keys: 'site', 
                  'inlet', 
                  'species', 
                  'instrument',
                  'start_date', 
                  'end_date', 
                  'calibration_scale',
                  'data_level', 
                  'store', 
                  'averaging_period', 
                  'filters'
                  
        bc_dict (dict): 
            Dictionary of Boundary Conditions data specifications
            Keys: 'species', 
                  'domain', 
                  'bc_input', 
                  'bc_freq',
                  'start_date', 
                  'end_date', 
                  'store',
                  
        use_bc (bool): 
            True: Retrieves BCs as specified in bc_dict.
            False: Asssumes Obs are assumed to be baseline-subtracted 
                   And BC data are NOT used


    Returns:
        data_dict: 
            Dictionary containing modelled simulations, fluxes, BCs, 
        sites:
            List of returned sitenames that have available data over 
            period of interest
        inlet:
            List of inlet heights for sites that have available data  
            over period of interest 
        fp_height:
            List of footprints heights for sites that have available
            data over period of interest
        instrument:
            List of instrument names for sites that have available
            data over period of interest
        averaging_period:
            List of averaging periods for sites that have available
            data over period of interest
    -------------------------------------------------------
    """ 
    data_dict = {}
    sites = fp_dict["site"]
    
    # Get CO2 flux fields
    data_dict[".flux"] =  get_flux_data(flux_dict)
    
    # Get CO2 boundary conditions        
    if use_bc is True:
        data_dict[".bc"] = get_mf_bc_data(bc_dict)
    else:
        data_dict[".bc"] = None
    
    # Get dataset source entry (if exists, used for verification games)
    if "dataset_source" in obs_dict.keys():
        dataset_source = obs_dict["dataset_source"]
    else:
        dataset_source = None
        
    # Get mole fraction observations 
    obs_mf_dict = get_mf_obs(obs_dict, 
                             dataset_source=dataset_source,
                            )
    
    # Get footprint datasets
    fp_data_dict = get_fps_data(fp_dict)
   
    site_indices_to_keep = []
    for i, site in enumerate(sites):
        if (site in obs_mf_dict.keys()) and (site in fp_data_dict.keys()):
            # ____ Create CO2 mole fraction forward simulations ____
            try:

                site_ind = np.where(np.array(obs_dict["site"])==site)[0][0]
                
                model_scenario = openghg.analyse.ModelScenario(site=site,
                                                               species=obs_dict["species"],
                                                               inlet=obs_dict["inlet"][site_ind],
                                                               start_date=obs_dict["start_date"],
                                                               end_date=obs_dict["end_date"],
                                                               obs=obs_mf_dict[site],
                                                               footprint=fp_data_dict[site],
                                                               flux=data_dict[".flux"],
                                                               bc=data_dict[".bc"],
                                                              )

                if len(flux_dict["source"]) == 1:
                    scenario_combined = model_scenario.footprints_data_merge()
                    if str(obs_dict["species"]).lower() == "co2":
                        Hall_t, Hall_yxt = model_scenario._calc_modelled_obs_HiTRes(sources=flux_dict["source"],
                                                                                    averaging=obs_dict["averaging_period"][site_ind],
                                                                                    output_fpXflux=True,
                                                                                   )
                        scenario_combined["Hall"] = Hall_yxt
                    else:
                        Hall_t, Hall_yxt = model_scenario.calc_modelled_obs(sources=flux_dict["source"],
                                                                            averaging=obs_dict["averaging_period"][site_ind],
                                                                            output_fpXflux=True,
                                                                           )
                        scenario_combined["Hall"] = Hall_yxt


                elif len(flux_dict["source"]) > 1:
                    model_scenario_dict = {}
                    for source in flux_dict["source"]:
                        scenario_sector = model_scenario.footprints_data_merge(sources=source, 
                                                                               recalculate=True,
                                                                              )
                        if fp_dict["species"].lower() == "co2":
                            Hall_t, Hall_yxt = model_scenario._calc_modelled_obs_HiTRes(sources=source,
                                                                                        averaging=obs_dict["averaging_period"][site_ind],
                                                                                        output_fpXflux=True,
                                                                                       )
                        
                            model_scenario_dict[f"mf_mod_high_res_{source}"] = scenario_sector["mf_mod_high_res"]
                            model_scenario_dict[f"Hall_{source}"] = Hall_yxt
                    
                        elif fp_dict["species"].lower() != "co2":
                            Hall_t, Hall_yxt = model_scenario.calc_modelled_obs(sources=source,
                                                                                averaging=obs_dict["averaging_period"][site_ind],
                                                                                output_fpXflux=True,
                                                                               )
                        
                            model_scenario_dict[f"mf_mod_{source}"] = scenario_sector["mf_mod"]
                            model_scenario_dict[f"Hall_{source}"] = Hall_yxt

                    scenario_combined = model_scenario.footprints_data_merge(recalculate=True)
                    for k, v in model_scenario_dict.items():
                        scenario_combined[k] = v
            
                data_dict[site] = scenario_combined
                # data_dict[site].bc_mod.values *= 1e-3 # convert from ppb (default in openghg) to ppm
                data_dict[site].bc_mod.values *= 1e-9 # convert from ppb (default in openghg) to mol/mol
                site_indices_to_keep.append(site_ind)

            except SearchError:
                print(f"\nError in reading in BC or flux data for {site}. \n")

    
    if len(site_indices_to_keep) == 0:
        raise SearchError("No site data found. Exiting process.")

    # If data was not extracted correctly for any sites, 
    # drop these from the rest of the inversion
    if len(site_indices_to_keep) < len(sites):
        sites = [obs_dict["site"][s] for s in site_indices_to_keep]
        inlet = [obs_dict["inlet"][s] for s in site_indices_to_keep]
        fp_height = [fp_dict["fp_height"][s] for s in site_indices_to_keep]
        instrument = [obs_dict["instrument"][s] for s in site_indices_to_keep]
        averaging_period = [obs_dict["averaging_period"][s] for s in site_indices_to_keep]
    
        return data_dict, sites, inlet, fp_height, instrument, averaging_period

    else:
        return data_dict, None, None, None, None, None
