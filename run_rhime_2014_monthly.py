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
sys.path.append("/home/users/alice.ramsden/rhime-co2-inversions/")
import rhime_co2
import argparse

def inputs(start_date,end_date):
    """
    Function defining input dictionaries used for the inversions 
    Edit as needed 
    """
    
    # Dictionary containing observations data specifications
    obs_dict = {"species": "co2",
                "site": ["HFD", "MHD", "RGL", "TAC"],
                "inlet": ["100m", "24m", "90m", "185m"],
                "averaging_period": ["4H", "4H", "4H", "4H"],
                "instrument": ["picarro", None, "picarro", "picarro"],
                "data_level": [None, None, None, None],
                "store": "user",
                "calibration_scale": None,
                "start_date": start_date,#"2014-04-01",
                "end_date": end_date,#"2014-05-01",
                "filters": ["daytime"],
               }
    '''
    # Dictionary containing CO2 fluxes data specifications 
    flux_dict = {"species": "co2", 
                 "domain": "EUROPE",
                 "source": ["ff-edgar-mth","gpp-jules-mth", "rtot-jules-mth"], 
                "start_date": start_date,#"2014-04-01",
                "end_date": end_date,#"2014-05-01",
                 "store": "user",
                 "flux_sf": {"ff-edgar-mth": 1.0, 
                             "gpp-jules-mth": 1.0,
                             "rtot-jules-mth": 1.0,
                            }, 
                 "sector_dict": {"fossil": "ff-edgar-mth",
                                 "gee": "gpp-jules-mth",
                                 "rtot": "rtot-jules-mth",
                                },
                }
    '''
    # Dictionary containing CO2 fluxes data specifications 
    flux_dict = {"species": "co2", 
                 "domain": "EUROPE",
                 "source": ["ff-edgar-mth","gpp-jules-mth","rtot-jules-mth"], 
                "start_date": start_date,#"2014-04-01",
                "end_date": end_date,#"2014-05-01",
                 "store": "user",
                 "flux_sf": {"ff-edgar-mth": 1.0,
                             "gpp-jules-mth": 1.0,
                             "rtot-jules-mth": 1.0
                            }, 
                 
                 "sector_dict": {"fossil": "ff-edgar-mth",
                                 "gpp": "gpp-jules-mth",
                                 "rtot": "rtot-jules-mth"
                                },
                }

    # Footprints dictionary CO2 data specifications
    fp_dict = {"species": "co2",
               "domain": "EUROPE",
               "site" : ["HFD", "MHD", "RGL", "TAC"], 
               "fp_height": ["100m", "10m", "90m", "185m"],
                "start_date": start_date,#"2014-04-01",
                "end_date": end_date,#"2014-05-01",
               "store": "user",            
              }
    
    # Boundary conditions dictionary specifications
    bc_dict = {"species": "co2",
               "domain": "EUROPE",
               "bc_input" : "camsv19",
               "bc_freq": "monthly",
                "start_date": start_date,#"2014-04-01",
                "end_date": end_date,#"2014-05-01",
               "store": "user",
               "bc_sf": None,
              }

    # Basis functions dictionary 
    basis_dict = {"fp_basis_case": None,
                  "basis_directory": None, 
                  "fp_basis_algorithm": "weighted",
                  "nbasis": [50, 50],# 50],
                  "bc_basis_case": "NESW",
                  "bc_basis_directory": "/data/users/alice.ramsden/LPDM_co2/bc_basis_functions",
                 }

    # MCMC dict
    mcmc_inputs_dict = {"xprior": {"ff-edgar-mth": {"pdf": "truncatednormal", "mu": 1.2, "sigma": 1.2, "lower":0.0},
                                   "rtot-jules-mth": {"pdf": "truncatednormal", "mu": 1, "sigma": 2.0, "lower":0.0},
                                   "gpp-jules-mth": {"pdf": "truncatednormal", "mu": 1, "sigma": 2.0, "lower":0.0}},
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
    
    '''
    mcmc_inputs_dict = {"xprior": {"ff-edgar-mth": {"pdf": "truncatednormal", "mu": 1.2, "sigma": 1.2, "lower":0.0},
                                   "gpp-jules-mth": {"pdf": "truncatednormal", "mu": 1.0, "sigma": 2.0, "lower": 0.0},
                                   "rtot-jules-mth": {"pdf": "truncatednormal", "mu": 1.0, "sigma": 2.0, "lower":0.0},
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
    
    '''
                        
    return obs_dict, flux_dict, fp_dict, bc_dict, basis_dict, mcmc_inputs_dict    


def main(start_date,end_date):    
    obs_dict, flux_dict, fp_dict, bc_dict, basis_dict, mcmc_dict = inputs(start_date,end_date)        
    use_bc = True
    model_error_method = "residual"
    sigma_freq = None
        
    outputname = f"co2_4site_3sector_monthly"
    outputpath = "/home/users/alice.ramsden/data/co2_inversions/"
    country_file = "/home/users/alice.ramsden/data/LPDM_co2/countries/country_EUROPE_EEZ_PARIS_gapfilled.nc"
    
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
                               save_merged_data=True,
                               read_merged_data=False,
                               merged_data_name='/data/scratch/alice.ramsden/merged_data/co2_4site_3sector_monthly'
                              )

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Running MCMC inversion for CO2')
    parser.add_argument("start_date", help="Start date with format YYYY-MM-DD",nargs="?")                  
    parser.add_argument("end_date", help="End date with format YYYY-MM-DD",nargs="?")
    
    args = parser.parse_args()
    
    command_line_args = {}
    if args.start_date:
        command_line_args["start_date"] = args.start_date
    if args.end_date:
        command_line_args["end_date"] = args.end_date
    
    main(**command_line_args)
