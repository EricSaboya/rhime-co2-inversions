# ---------------------------------------------------------------------------------------
# Filename: example_run_rhime_202101.py
# Created 6 January 2025
# Author: Eric Saboya, School of Geographical Sciences, University of Bristol 
# Project: Testing RHIME
# ---------------------------------------------------------------------------------------
# About:
#   Example script for running multi-sector RHIME CO2 inversions 
#   Reads inversion inputs from inputs.py and saves inputs as a pickle dictionary 
# ---------------------------------------------------------------------------------------
# Summary:
#   Example file for running RHIME inversions
# ---------------------------------------------------------------------------------------

import sys 
import pickle
import numpy as np
sys.path.append("/user/work/wz22079/projects/rhime-co2-inversions/")
import rhime_co2

def inputs():
    """
    Function defining input dictionaries used for the inversions 
    Edit as needed 
    """
    # Dictionary containing observations data specifications
    obs_dict = {"species": "co2",
                "site": ["HFD", "MHD", "RGL", "TAC"],
                "inlet": ["100m", "24m", "90m", "185m"],
                "averaging_period": ["4H", "4H", "4H", "4H"],
                "instrument": ["picarro", "multiple", "picarro", "picarro"],
                "data_level": [None, "2", None, None],
                "store": "obs_nir_2024_01_25_store_zarr",
                "calibration_scale": None,
                "start_date": "2014-01-01",
                "end_date": "2014-02-01",
                "filters": ["daytime"],
               }
    
    # Dictionary containing CO2 fluxes data specifications 
    flux_dict = {"species": "co2", 
                 "domain": "EUROPE",
                 "source": ["edgar-fossil-hrly-flat", "vprm-gee", "vprm-resp"], 
                 "start_date": "2014-01-01",
                 "end_date": "2014-02-01",
                 "store": "uk_co2_zarr_store",
                 "flux_sf": {"edgar-fossil-hrly-flat": 1.0, 
                             "vprm-gee": 1.0,
                             "vprm-resp": 1.0,
                            }, 
                 "sector_dict": {"fossil": "edgar-fossil-hrly-flat",
                                 "gee": "vprm-gee",
                                 "resp": "vprm-resp",
                                },
                }

    # Footprints dictionary CO2 data specifications
    fp_dict = {"species": "co2",
               "domain": "EUROPE",
               "site" : ["HFD", "MHD", "RGL", "TAC"], 
               "fp_height": ["100m", "10m", "90m", "185m"],
               "start_date": "2014-01-01",
               "end_date": "2014-02-01",
               "store": "uk_co2_footprints_202406",            
              }
    
    # Boundary conditions dictionary specifications
    bc_dict = {"species": "co2",
               "domain": "EUROPE",
               "bc_input" : "camsv22-co2",
               "bc_freq": "monthly",
               "start_date": "2014-01-01",
               "end_date": "2014-02-01",
               "store": "uk_co2_zarr_store",
               "bc_sf": None,
              }

    # Basis functions dictionary 
    basis_dict = {"fp_basis_case": None,
                  "basis_directory": None, 
                  "fp_basis_algorithm": "weighted",
                  "nbasis": [50, 50, 50],
                  "bc_basis_case": "NESW",
                  "bc_basis_directory": "/group/chemistry/acrg/LPDM/bc_basis_functions/",
                 }

    # MCMC dict
    mcmc_inputs_dict = {"xprior": {"edgar-fossil-hrly-flat": {"pdf": "truncatednormal", "mu": 1.2, "sigma": 1.2, "lower":0.0},
                                   "vprm-gee": {"pdf": "truncatednormal", "mu": 1.0, "sigma": 2.0, "lower": 0.0},
                                   "vprm-resp": {"pdf": "truncatednormal", "mu": 1.0, "sigma": 2.0, "lower":0.0},
                                  },
                        "bcprior": {"pdf": "truncatednormal", "lower": 0.0, "mu":1.0, "sigma": 0.05},
                        "sigprior": {"pdf": "uniform", "lower": 0.1, "upper": 3.0}, 
                        "add_offset": False, 
                        "offsetprior": None, 
                        "nit": 5500,
                        "burn": 1000, 
                        "tune": 2000,
                        "nchain": 2,
                        "sigma_per_site": True
                       }
                        
    return obs_dict, flux_dict, fp_dict, bc_dict, basis_dict, mcmc_inputs_dict    


def main():    
    obs_dict, flux_dict, fp_dict, bc_dict, basis_dict, mcmc_dict = inputs()        
    use_bc = True
    model_error_method = "residual"
    sigma_freq = None
        
    outputname = f"UK_TFEB2025_"
    outputpath = "/group/chemistry/acrg/ES/co2/"
    country_file = "/group/chemistry/acrg/LPDM/countries/country_EUROPE_EEZ_PARIS_gapfilled.nc"
    
    rhime_co2.rhime_inversions(obs_dict=obs_dict,
                               flux_dict=flux_dict,
                               bc_dict=bc_dict,
                               fp_dict=fp_dict,
                               basis_dict=basis_dict,
                               mcmc_dict=mcmc_dict,
                               use_bc=use_bc,
                               model_error_method=model_error_method,
                               sigma_freq=sigma_freq,
                               outputname=outputname,
                               outputpath=outputpath,
                               country_file=country_file,
                              )
                     
                                
    input_dict = {"obs_inputs": obs_dict,
                  "flux_inputs": flux_dict,
                  "footprint_inputs": fp_dict,
                  "boundary_condition_inputs": bc_dict,
                  "basis_function_inputs": basis_dict,
                  "mcmc_inputs": mcmc_dict,
                 }

    species = obs_dict["species"]
    domain = flux_dict["domain"]
    start_date = obs_dict["start_date"]
    print("Saving inputs ... ")
    with open(f"{outputpath}/{species}_{domain}_{outputname}_{start_date}_INPUTS.pkl", "wb") as f:
        pickle.dump(input_dict, f)


if __name__ == "__main__":
    main()
