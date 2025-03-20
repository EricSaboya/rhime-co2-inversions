import pint
import numpy as np

def create_unit_registry():
    """
    -------------------------------------------------------
    Create customised pint UnitRegistry for atmospheric 
    trace gas species
    -------------------------------------------------------
    """
    ureg = pint.UnitRegistry()
    ureg.define("ppb = 1e-9 mol/mol")
    ureg.define("ppt = 1e-12 mol/mol")
    ureg.define("m2 = m*m")
    
    return ureg

def check_obs_units(obs_data):
    """
    -------------------------------------------------------
    Apply unit corrections to atmospheric mole fraction
    observations retrieved from an OpenGHG objectstore
    -------------------------------------------------------
    Args:
        obs_data (OpenGHG.object)
            Retrieved obs data from OpenGHG store

    Returns:
        obs_data (OpenGHG.object)
            Obs data but presented in mol/mol
    -------------------------------------------------------
    """
    ureg = create_unit_registry()
    openghg_pint_unit_mapping = {"1e-6": "ppm",
                                 "1e-9": "ppb",
                                 "1e-12": "ppt",
                                }
    
    for key in obs_data.data.keys():
        if "units" in list(obs_data.data[key].attrs.keys()):
            i_unit = obs_data.data[key].attrs["units"]
            
            if i_unit in openghg_pint_unit_mapping.keys():
                conversion = ureg[openghg_pint_unit_mapping[i_unit]].to_reduced_units()
                
            elif i_unit == "mol/mol":
                conversion  = ureg[i_unit].to_reduced_units()

            else:
                raise KeyError("Unit not recognised.")

            obs_data.data[key].values *= conversion
            obs_data.data[key].attrs["units"] = "mol/mol"

    return obs_data

def check_footprint_units(fp_data):
    """
    -------------------------------------------------------
    Apply unit corrections to LPDM footprints
    retrieved from an OpenGHG objectstore
    -------------------------------------------------------
    Args:
        fp_data (OpenGHG.object)
            Retrieved footprints data from OpenGHG store

    Returns:
        fp_data (OpenGHG.object)
            Footprint data but presented in m2 * s / mol
    -------------------------------------------------------    
    """
    ureg = create_unit_registry()

    for key in fp_data.data.keys():
        if "units" in list(fp_data.data[key].attrs.keys()):
            if "mol" in fp_data.data[key].attrs["units"]:
                i_unit = fp_data.data[key].attrs["units"]
                conversion = ureg[i_unit].to_reduced_units()

                fp_data.data[key] *= conversion
                fp_data.data[key].attrs["units"] = conversion
                
    return fp_data




    



# def check_units(data_dict: dict) -> dict:
#     """
#     -------------------------------------------------------
#     Function to check data units in obs, forward sims, 
#     and BCs are consistent. 

#     Defaults to converting 'units' to ppm
#     -------------------------------------------------------
#     Args:
#         data_dict (dict):
#             Dictionary of xr.DataArray outputs from
#             get_co2_data.py

    
#     Returns:
#         data_dict (dict)
#             Dictionary where units are checked, amended,
#             and returned
#     -------------------------------------------------------
#     """
#     # Get list of sites
#     sites = []
#     for key in data_dict.keys():
#         if "." not in key:
#             sites.append(key)

#     # Get list of flux source keys
#     sources = list(data_dict[".flux"].keys())

#     # Get list of variables to check
#     myvars = {"mf": [220, 1000],
#               "mf_variability": [0.01, 30],
#               "mf_repeatability": [0.01, 30],
#               "mf_mod_high_res": [0.01, 30],
#               "bc_mod": [220, 1000],
#              }
    
#     # Keep data in PPM
#     for key in myvars.keys():
#         for site in sites:
#             check_var = np.nanmean(data_dict[site][key])

#             if (check_var>=myvars[key][0]) and (check_var<myvars[key][1]):
#                 print(f"CO2 {key} for {site} (very likely) in ppm")
                
#             else: 
#                 if key != "mf_mod_high_res":
#                     if (check_var*1e6> myvars[key][0]) and (check_var*1e6 <myvars[key][1]):
#                         print(f"Multiplying {key} variable by E6 ... ")
#                         data_dict[site][key].values *= 1e6
                        
#                     elif check_var * 1e6 > myvars[key][1]:
#                         raise ValueError(f"Check input data for {key}")
                        
#                 else:
#                     if (check_var*1e6> myvars[key][0]) and (check_var*1e6 <myvars[key][1]):
#                         print(f"Multiplying {key} variable by E6 ... ")
#                         data_dict[site][key].values *= 1e6
#                         for source in sources:
#                             print(f"Multiplying {key}, {source} variable by 10^6 ... ")
#                             data_dict[site][f"mf_mod_high_res_{source}"].values * 1e6
                
#                     elif check_var * 1e6 > myvars[key][1]:
#                          raise ValueError(f"Check input data for {key}")
                        
#     return data_dict