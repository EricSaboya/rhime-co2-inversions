# ---------------------------------------------------------------------------------------
# model_error_methods.py
# Created 25 Sept. 2024 
# Author: Eric Saboya, School of Geographical Sciences, University of Bristol 
# ---------------------------------------------------------------------------------------
# Methods for calculating model errors 
# > "residual"
# > "percentile"
# ---------------------------------------------------------------------------------------

import numpy as np 
import xarray as xr

def model_error_method_parser(data_dict: dict, 
                              method: str,
                             ):
    """
    Parser function for calculating model errors

    Args:
    data_dict (dict):
        Output dictionary from basis fucntions 
    method (str):
        Model error method of interest 
    """
    sites = []
    for key in data_dict.keys():
        if "." not in key:
            sites.append(key)      

    if method == None:
        data_dict = no_model_error(sites, data_dict)        
    elif method == "fixed":
        data_dict = fixed_method(sites, data_dict)
    elif method == "residual":
        data_dict = residual_method(sites, data_dict)
    return data_dict


def no_model_error(sites: list,
                   data_dict: dict,
                  ):
    """
    Model error of zero applied to data 
    """
    for site in sites:
        y_epsilon_m = data_dict[site].mf * 0.0
        data_dict[site]["y_model_err"] = y_epsilon_m
    return data_dict

def fixed_method(sites: list,
                 data_dict: dict,
                ):
    """
    Model error based on the mean obs-sim 
    difference 
    """
    for site in sites:
        y_epsilon_m = np.abs(np.nanmean(fp_data[site].mf - fp_data[site].mf_mod_high_res - fp_data[site].bc_mod))
        data_dict[site]["y_model_err"] = (data_dict[site].mf * 0.0) + y_epsilon_m
    return data_dict

def residual_method(sites: list, 
                    data_dict: dict,
                   ):
    """
    This method is explained in "Modeling of Atmospheric Chemistry" by Brasseur
    and Jacobs in Box 11.2 on p.499-500, following "Comparative inverse analysis of satellitle (MOPITT)
    and aircraft (TRACE-P) observations to estimate Asian sources of carbon monoxide", by Heald, Jacob,
    Jones, et.al. (Journal of Geophysical Research, vol. 109, 2004).

    Roughly, we assume that the observations y are equal to the modelled observations y_mod (mf_mod + bc_mod),
    plus a bias term b, and instrument, representation, and model error:

    y = y_mod + b + err_I + err_R + err_M

    Assuming the errors are mean zero, we have

    (y - y_mod) - mean(y - y_mod) = err_I + err_R + err_M  (*)

    where the mean is taken over all observations.

    Calculating the RMS of the LHS of (*) gives us an estimate for

    sqrt(sigma_I^2 + sigma_R^2 +  sigma_M^2),

    where sigma_I is the standard deviation of err_I, and so on.

    Thus a rough estimate for sigma_M is the RMS of the LHS of (*), possibly with the RMS of
    the instrument/observation and averaging errors removed (this isn't implemented here).

    Note: in the "non-robust" case, we are computing the standard deviation of y - y_mod. The mean on the LHS
    of equation (*) could be taken over a subset of the observation, in which case the value calculated is not
    a standard deviation. We wrote the derivation this way to match Brasseur and Jacobs.

    Args:
        ds_dict: dictionary of combined scenario datasets, keyed by site codes.
        robust: if True, use the "median absolute deviation" (https://en.wikipedia.org/wiki/Median_absolute_deviation)
            instead of the standard deviation. MAD is a measure of spread, similar to standard deviation, but
            is more robust to outliers.
        by_site: if True, return array with one mininum error value per site

    Returns:
        np.ndarray: estimated value(s) for model error.

    """
    for site in sites:
        y = data_dict[site].mf
        y_mod = data_dict[site].bc_mod + data_dict[site].mf_mod_high_res
        # Observational error
        y_epsilon_o = (y-y_mod) - np.nanmean(y-y_mod)
        # Instrumental error
        y_epsilon_i = np.sqrt(data_dict[site].mf_variability**2 + data_dict[site].mf_repeatability**2)
        # Model error 
        y_epsilon_m = np.sqrt(np.abs(y_epsilon_o**2 - y_epsilon_i**2))
        data_dict[site]["y_model_err"] = y_epsilon_m
    return data_dict
