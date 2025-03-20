# *****************************************************************************
# data_preprocessing.py
# Created: 15 May 2024
# Author: Eric Saboya, School of Geographical Sciences, University of Bristol
# *****************************************************************************
# Functions for retrieving CO2 datasets needed for performing inversions 
# *****************************************************************************

import os
import glob
import numpy as np
import pandas as pd
import xarray as xr
import openghg
import datetime as dt

from openghg.types import SearchError

def get_obs_site_info():
    """
    Returns dictionary of sites in EUROPE model
    domain with their inlet heights as defined
    in the ICOS/DECC/ ObsPack
    """
    site_dict = {"MHD": "24m",
                 "BSD": "248m",
                 "HFD": "100m",
                 "RGL": "90m",
                 "TAC": "185m",
                 "WAO": "10m",
                 "BIR": "75m",
                 "CBW": "207m",
                 "CMN": "8m",
                 "GAT": "341m",
                 "HEI": "30m",
                 "HEL": "110m",
                 "HPB": "131m",
                 "HTM": "150m",
                 "HUN": "115m",
                 "JFJ": "13.9m",
                 "KIT": "200m",
                 "KRE": "250m",
                 "LIN": "98m",
                 "LUT": "60m",
                 "NOR": "100m",
                 "OPE": "120m",
                 "OXK": "163m",
                 "PAL": "12m",
                 "SAC": "100m",
                 "SSL": "35m",
                 "STE": "252m",
                 "TOH": "147m",
                 "TRN": "180m",
                 "UTO": "57m",
                 "WES": "14m"}   
    return site_dict

def get_fp_site_info():
    """
    Returns dictionary of sites in EUROPE model
    domain with their NAME inlet heights 
    """
    site_dict = {"MHD": "10m",
                 "BSD": "248m",
                 "HFD": "100m",
                 "RGL": "90m",
                 "TAC": "185m",
                 "WAO": "20m",
                 "BIR": "75m",
                 "CBW": "200m",
                 "CMN": "500m",
                 "GAT": "341m",
                 "HEI": "30m",
                 "HEL": "110m",
                 "HPB": "130m",
                 "HTM": "150m",
                 "HUN": "115m",
                 "JFJ": "1000m",
                 "KIT": "200m",
                 "KRE": "250m",
                 "LIN": "98m",
                 "LUT": "60m",
                 "NOR": "100m",
                 "OPE": "120m",
                 "OXK": "163m",
                 "PAL": "10m",
                 "SAC": "100m",
                 "SSL": "10m",
                 "STE": "252m",
                 "TOH": "147m",
                 "TRN": "180m",
                 "UTO": "57m",
                 "WES": "14m"}
    return site_dict

def get_hourly_times(start_time=None, end_time=None):
    """
    Creates array of hourly timestamps values
    between start_time and end_time. Defaults to 2021 values
    
    Args:
        start_time (str): First timestamp with format (DD-MM-YYYY)
        end_time (str): Last timestamp with format (DD-MM-YYYY)

    Returns:
        array of hourly timestamps between start_time and end_time
    """
    import datetime as dt
    
    if start_time == None and end_time == None:
        start_date_string = "01-01-2021 00:00:00.000"
        end_date_string = "01-01-2022 00:00:00.000"
    else:
        start_date_string = start_time + " 00:00:00.000"
        end_date_string = end_time + " 00:00:00.000"
    
    start_datetime = dt.datetime.strptime(start_date_string, '%d-%m-%Y %H:%M:%S.%f')
    end_datetime = dt.datetime.strptime(end_date_string, '%d-%m-%Y %H:%M:%S.%f')
    
    # Create a range of dates
    index = pd.date_range(start=start_datetime, 
                          end=end_datetime,
                          freq="1H")
    
    return index[0:-1]

def align_obs(site, start_date, end_date, obs_store):
    """
    Retrieve CO2 measurements data and align data
    to a temporally uniform array of hourly 
    values. Keeps missing data points. 

    Args:
        site (str): Site name 
        start_date (str): First timestamp with format (DD-MM-YYYY)
        end_date (str): Last timestamp with format (DD-MM-YYYY)
        obs_store (str): Object store containing CO2 obs data
    
    """
    site_dict = get_obs_site_info()
    index = get_hourly_times(start_date, end_date)
    
    ds_obs = openghg.retrieve.get_obs_surface(store=obs_store, 
                                              species="co2", 
                                              site=site, 
                                              inlet=site_dict[site],
                                              start_date=start_date, 
                                              end_date=end_date,
                                              average="1H",
                                              keep_missing=True,
                                             )

    co2_obs_array = np.zeros(len(index)) + np.nan
    co2_obs_n_array = np.zeros(len(index)) + np.nan
    co2_obs_var_array = np.zeros(len(index)) + np.nan

    co2_obs = ds_obs[site.lower()]["mf"].values
    co2_obs_n = ds_obs[site.lower()]["mf_number_of_observations"].values
    co2_obs_var = ds_obs[site.lower()]["mf_variability"].values
    co2_t = ds_obs[site.lower()]["time"].values
    
    for t in range(len(co2_t)):
        dt = np.array(index) - co2_t[t]
        dt = dt.astype(float)
        ind = np.intersect1d(np.where(dt>-1000), np.where(dt<1000))[0]
        co2_obs_array[ind] = co2_obs[t]
        co2_obs_n_array[ind] = co2_obs_n[t]
        co2_obs_var_array[ind] = co2_obs_var[t]
        
    return co2_obs_array, co2_obs_n_array, co2_obs_var_array 

def align_background(site, start_date, end_date, bg_dir=None):
    """
    Get pre-calculated baseline CO2 mole fractions 
    and align with hourly values over start_date and end_date.
    Args:
        site (str): Site name 
        start_date (str): First timestamp with format (DD-MM-YYYY)
        end_date (str): Last timestamp with format (DD-MM-YYYY)
        bg_dir (str): Directory containing baseline mole fractions files 

    """
    index = get_hourly_times(start_date, end_date)
    co2_bg_array = np.zeros(len(index)) + np.nan

    if bg_dir is None:
        bg_dir="/user/work/wz22079/projects/CO2/paris_verification_games_mk2/baselines/"
    bg_search = glob.glob(os.path.join(bg_dir, f"*{site}*v23*"))
    
    if len(bg_search) == 1:
        bg = xr.open_dataset(bg_search[0])
        bg_t = bg["time"].values
        bg_c = bg["mf_baseline"].values
        
        for t in range(len(bg_t)):
            dt = np.array(index) - bg_t[t]
            dt = dt.astype(float)
            ind = np.intersect1d(np.where(dt>-1000), np.where(dt<1000))[0]
            co2_bg_array[ind] = bg_c[t]
    else:
        raise("Multiple background files found. Improve search specification.")
          
    return index, co2_bg_array



def align_forward_sims(flux_dict, fp_dict, obs_dict, bc_dict, use_bc):
    """
    Function to produce forward simulations of 
    CO2 mole fraction baseline perturbations. 
    
    Args:
        flux_dict (dict): Dictionary of flux data specifications 
        fp_dict (dict): Dictionary of footprints data specifications
        obs_dict (dict): Dictionary of observations data specifications
        bc_dict (dict): Dictionary of Boundary Conditions data specifications
    """
    # *************************************************************** # 
    #   
    # ____ Get CO2 Flux Data ____
    # Get fluxes from flux_dict
    #
    flux_data_dict = {}
    for source in flux_dict["source"]:
        get_flux_data = openghg.retrieve.get_flux(species=flux_dict["species"],
                                                  domain=flux_dict["domain"],
                                                  source=source,
                                                  start_date=flux_dict["start_date"],
                                                  end_date=flux_dict["end_date"],
                                                  store=flux_dict["store"],
                                                 )
        flux_data_dict[source] = get_flux_data

     
    data_dict = {}
    sites = fp_dict["site"]
    site_indices_to_keep = []
    
    for i, site in enumerate(sites):
    
    # ___ Get CO2 footprints data for each site ____ 
        try:
            get_fps_data = openghg.retrieve.get_footprint(site=fp_dict["site"][i],                            
                                                          height=fp_dict["fp_height"][i],
                                                          domain=fp_dict["domain"],
                                                          start_date=fp_dict["start_date"],
                                                          end_date=fp_dict["end_date"],
                                                          store=fp_dict["store"],
                                                          species=fp_dict["species"],
                                                         )
        except SearchError:
            print(f"\nNo footprint data found for {site}.\n",)
            continue  # skip this site
            
    # ___ Get CO2 boundary conditions ____        
        if use_bc is True:
            get_bc_data = openghg.retrieve.get_bc(species=bc_dict["species"],
                                                  domain=bc_dict["domain"],
                                                  bc_input=bc_dict["bc_input"],
                                                  start_date=bc_dict["start_date"],
                                                  end_date=bc_dict["end_date"],
                                                  store=bc_dict["store"],
                                                 )
        else:
            get_bc_data = None    
    
    # ____ Get CO2 observations for each site ____    
        try:
            site_data = openghg.retrieve.get_obs_surface(site=obs_dict["site"][i],
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
                                                        )
            
            # Calculate observational uncertainty over averaging period 
            site_data_no_ave = openghg.retrieve.get_obs_surface(site=obs_dict["site"][i],
                                                                species=obs_dict["species"],
                                                                inlet=obs_dict["inlet"][i],
                                                                start_date=obs_dict["start_date"],
                                                                end_date=obs_dict["end_date"],
                                                                icos_data_level=obs_dict["data_level"][i],
                                                                instrument=obs_dict["instrument"][i],
                                                                calibration_scale=obs_dict["calibration_scale"],
                                                                store=obs_dict["store"],
                                                                keep_missing=True,
                                                                )

            # Observation uncertainty is defined as the sum in quadrature of:
            # 1. The root-squared sum of the 1 min mf variabilities reported in the uploaded data in the averaging period
            #    and divided by the number of observations in the averaging period. This represents the "error" of the 
            #    sampled mean.
            # 2. The variance of the mole fraction concentrations reported in the uploaded data over the averaging period.
            #    This accounts for the variability in mole fractions over the averaging period
            # 
            
            # Sum in quadrature of "mf_variability" values over averaging period
            ds_resampled_uncert = np.sqrt((site_data_no_ave[site]["mf_variability"]**2).resample(time=obs_dict["averaging_period"][i]).sum()) 
            / site_data_no_ave[site]["mf_variability"].resample(time=obs_dict["averaging_period"][i]).count()

            # Average variance of mole fractions across averaging period
            ds_resampled_var = site_data_no_ave[site]["mf"].resample(time=obs_dict["averaging_period"][i]).std(skipna=False, keep_attrs=True)**2

            site_data[site]["mf_variability"] = np.sqrt(ds_resampled_uncert**2 + ds_resampled_var**2)
    
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

        try:
        # ____ Create CO2 mole fraction forward simualtions ____
            model_scenario = openghg.analyse.ModelScenario(site=site,
                                                           species=obs_dict["species"],
                                                           inlet=obs_dict["inlet"][i],
                                                           start_date=obs_dict["start_date"],
                                                           end_date=obs_dict["end_date"],
                                                           obs=site_data,
                                                           footprint=get_fps_data,
                                                           flux=flux_data_dict,
                                                           bc=get_bc_data,
                                                          )

            if len(flux_dict["source"]) == 1:
                scenario_combined = model_scenario.footprints_data_merge()
                
            elif len(flux_dict["source"]) > 1:
                model_scenario_dict = {}
                for source in flux_dict["source"]:
                    scenario_sector = model_scenario.footprints_data_merge(sources=source, 
                                                                           recalculate=True)
                    if fp_dict["species"].lower() == "co2":
                        model_scenario_dict["mf_mod_high_res_" + source] = scenario_sector["mf_mod_high_res"]
                    
                    elif fp_dict["species"].lower() != "co2":
                        model_scenario_dict["mf_mod_" + source] = scenario_sector["mf_mod"]

                scenario_combined = model_scenario.footprints_data_merge(recalculate=True)
                for k, v in model_scenario_dict.items():
                    scenario_combined[k] = v
                    
            data_dict[site] = scenario_combined
            data_dict[site].bc_mod.values *= 1e-3
            site_indices_to_keep.append(i)

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






def get_model_params():
    obs_dict = {"species": "co2",
                "site": ["CBW", "MHD"],
                "inlet": ["200m", "24m"],
                "averaging_period": ["1H", "1H"],
                "instrument": ["", ""],
                "data_level": ["2", "2"],
                "store": "obs_nir_2024_01_25_store_zarr",
                "calibration_scale": None,
                "start_date": "2021-01-01",
                "end_date": "2021-01-10",
               }
    
    flux_dict = {"species": "co2", 
                 "domain": "EUROPE",
                 "source": ["paris-aten-bio", "paris-aten-fossil"], 
                 "start_date": "2021-01-01",
                 "end_date": "2021-01-10",
                 "store": "co2_store_zarr",
                }
    
    fp_dict = {"species": "co2",
               "domain": "EUROPE",
               "site" : ["CBW", "MHD"],
               "fp_height": ["200m", "10m"],
               "start_date": "2021-01-01",
               "end_date": "2021-01-10",
               "store": "shared_store_zarr",            
              }

    bc_dict = {"species": "co2",
               "domain": "EUROPE",
               "bc_input" : "camsv23_3h",
               "start_date": "2021-01-01",
               "end_date": "2021-01-10",
               "store": "co2_store_zarr",            
              }

    return obs_dict, flux_dict, fp_dict, bc_dict
    


    