# *****************************************************************************
# tracers_co2.py
# Created: 11 Feb. 2025
# Author: Eric Saboya, School of Geographical Sciences, University of Bristol
# *****************************************************************************
# Data functions for CO2 tracers (14C, APO) 
# *****************************************************************************
import os
import sys
import numpy as np
import xarray as xr

from openghg.retrieve import get_flux

sys.path.append("/user/work/wz22079/projects/rhime-co2-inversions/")
from get_co2_data import get_flux_data

###############################################################################
# RADIOCARBON FUNCTIONS
###############################################################################

def d14c_nuc()->float:
    """
    -------------------------------------------------------
    Returns nuclear power plant 14CO2 signature, as
    defined by Basu et al. (2016), in per mil
    -------------------------------------------------------
    """
    d13c = -8.0               # Assumed atmospheric value
    n = (975/(d13c+1000))**2  # Mass-dependent fractionation correction
    r_std = 1.176e-12         # Modern 14C standard value
    return 1000 * (n/r_std)

def d14c_fossil()->float:
    """
    -------------------------------------------------------
    Return fossil fuel combustion 'big delta' signature 
    in per mil
    -------------------------------------------------------
    """
    return -1000

def d14c_bio_to_atm(start_date: str,
                    end_date: str,
                    additive_bias: float = 0,
                   ):
    """
    -------------------------------------------------------
    Returns flux field of 'big delta' signature of 14CO2 
    fluxes from terrestrial biosphere carbon pool to 
    atmosphere.
    i.e. Respiration D14CO2 'big delta' fluxes, in per mil
    -------------------------------------------------------
    Args:
        start_date (str):
            Date for the start of period of interest 
            e.g. '2022-01-01'

        end_date (str):
            Date for the start of period of interest 
            e.g. '2022-02-01'
            
        additive_bias (float, defaults to 0):
            Additive bias to apply to the resp 14C fluxes
            across the period of interest

    Returns:
        OpenGHG flux field aligned with observations
    -------------------------------------------------------
    """
    d14c_resp = get_flux(species="14co2",
                         domain="EUROPE",
                         source="c14-resp-lpjguess-hrly-flat",
                         start_date=start_date,
                         end_date=end_date,
                        )

    if additive_bias is not None:
        d14c_resp.data.flux += additive_bias
    
    return d14c_resp
                               
#def d14c_ocean_to_atm(start_date: str,
#                      end_date: str,
#                      average: str,
#                     ):
#    """
#    Return timeseries of 'big delta' signature of 14C fluxes 
#    from ocean carbon pool to atmosphere.
#    Args:
#        start_date (str):
#            Date for the start of period of interest 
#            e.g. '2022-01-01'
#
#        end_date (str):
#            Date for the start of period of interest 
#            e.g. '2022-02-01'
#
#        average (str):
#            Data averaging 
#
#            
#    """
#    #### CHECK ALL OF THESE ONCE DATA UPLOADED TO AN OBJECT STORE !!!!
#    d14c_ocean = get_obs_surface(species="d14co2",
#                                 site="TBD",
#                                 store="TBD",
#                                 start_date=start_date,
#                                 end_date=end_date,
#                                 average=average,
#                                 keep_missing=True,
#                                )
#    return d14c_ocean


def get_14c_bg(start_date: str,
               end_date: str,
               method: str = "MHD",
               additive_bias: float = 0,
              ):
    """
    -------------------------------------------------------
    Create background D14CO2 ('big delta') values for time 
    period of interest using selected method.
    -------------------------------------------------------
    Args:
        start_date (str):
            Date for the start of period of interest 
            e.g. '2022-01-01'

        end_date (str):
            Date for the start of period of interest 
            e.g. '2022-02-01'
            
        method (str, Defaults to MHD):
            Method to calculate the background D14CO2 
            values

        bg_bias (float, Defaults to 0):
            Additive bias to apply to the background
            value across the period of interest

    Returns:
        Mean monthly background value
    -------------------------------------------------------    
    """
    # Convert start and end dates to datetime format
    t_start = dt.datetime.strptime(start_date, "%Y-%m-%d")
    t_end = dt.datetime.strptime(end_date, "%Y-%m-%d")

    if method == "MHD":
        # Get processed clean/marine air MHD observations
        c14_mhd_path = "/user/work/wz22079/projects/CO2/radiocarbon/data/icos/"
        c14_mhd_fname = "uheiiup_l2_2025_1_mhd_24m_int_14day_clean.c14.nc"
        c14_mhd_marine = xr.open_dataset(os.path.join(c14_mhd_path, c14_mhd_fname))

        # Calculate mean value across the period of interest
        #   If no data found in that period, keep extending by +/-
        #   7 days until a value is found
        mnthly_bg = c14_mhd_marine.sel(time=slice(t_start, t_end))
        c14_bg_mean = np.nanmean(mnthly_bg.D14CO2.values)
        
        if np.isnan(c14_bg_mean):
            while np.isnan(c14_bg_mean):
                t_start = t_start - dt.timedelta(days=7)
                t_end = t_end + dt.timedelta(days=7)
                mnthly_bg = c14_mhd_marine.sel(time=slice(t_start, t_end))
                c14_bg_mean = np.nanmean(mnthly_bg.D14CO2.values)

        if additive_bias is not None:
            return c14_bg_mean += additive_bias

    else:
        raise ValueError("Other options for D14CO2 ('big delta') background values are not yet available")
        


def get_14c_obs_sims(data_dict: dict,
                     c14_dict: dict,
                     nuclear_correction: bool = False,
                     emisource_to_sector: dict,
                     start_date: str,
                     end_date: str,
                    ):
    """
    -------------------------------------------------------
    Creates 14CO2 isofluxes for different source sectors 
    -------------------------------------------------------
    Args:
        data_dict (dict):
            Dictionary of xr.DataArray outputs from
            get_co2_data.py
            
        c14_dict (dict):
            Dictionary of c14 input values
            
        nuclear_correction (bool, defaults to False):
            Option to include a correction for nuclear 
            power plant emissions 
            
        emisource_to_sector (dict):
            Dictionary aligning flux sources with 
            sector names e.g. {'fossil': 'edgar-fossil'}

        start_date (str):
            Date for the start of period of interest 
            e.g. '2022-01-01'

        end_date (str):
            Date for the start of period of interest 
            e.g. '2022-02-01'


    Returns:
        data_dict (dict):
            Dictionary of xr.DataArray outputs from
            get_co2_data.py including 14C values

        resp_isoflux (OpenGHG.dataobject):
            Data object of respiration fluxes
            from LPJ-GUESS
    -------------------------------------------------------
    """
    # Check whether all essential sectors are included in 
    # emisource_to_sector input dictionary
    essential_sectors = ["fossil", "resp", "gpp"]
    input_sectors = list(emisource_to_sector.keys())
    
    for s in essential_sectors:
        if s not in list(emisource_to_sector.keys()):
            raise KeyError(f"Missing {s} sector from input values!")

    # Get site names
    sites = [key for key in list(data_dict.keys()) if key[0] != "."]
    
    # Get respiration 14CO2 signature values
    try:
        # Apply additive bias if included in the input dict.
        resp_bias = c14_dict["c14_resp_sig"]["additive_bias"]
        c14_resp = d14c_bio_to_atm(start_date, 
                                   end_date,
                                   resp_bias,
                                  )
    except:
        c14_resp = d14c_bio_to_atm(start_date, 
                                   end_date,
                                  )

    for site in sites:
        for i, key in enumerate(emisource_to_sector.keys()):
            if key == "fossil":    
                # Calculate fossil fuel isoflux
                ff_sector = emisource_to_sector["fossil"]
                ff_isoflux = d14c_fossil() * data_dict[site][f"mf_mod_high_res_{ff_sector}"]

                if i==0:
                    isoflux_s = ff_isoflux
                else:
                    isoflux_s += ff_isoflux

            elif key == "resp":
                # Calculate respiration isoflux
                resp_sector = emisource_to_sector["resp"]
                resp_isoflux = np.nansum(np.nansum(c14_resp.data.flux * data_dict[site][f"Hall_{resp_sector}"], axis=1), axis=1)

                if i==0:
                    isoflux_s = resp_isoflux
                else:
                    isoflux_s += resp_isoflux
 
            elif key == "gpp":
                # Disequilibrium term
                gpp_sector = emisource_to_sector["gpp"]
                ymod_total = data_dict[site]["mf_mod_high_res"] + data_dict[site]["bc_mod"]
                disequil_isoflux = ymod_total - data_dict[site][f"mf_mod_high_res_{gpp_sector}"]
            
            elif key == "ocean":
                # Ocean to atmosphere isoflux
                print("Warning! Not currently able to apply ocean D14C correction term!")
                # ocean_sector = emisource_to_sector["ocean"]
                # if i==0:
                    # isoflux_s = resp_isoflux
                # else:
                    # isoflux_s+= resp_isoflux

                
        # Get background D14CO2 values and create BG isoflux
        try:
            bg_bias = c14_dict["c14_bg"]["additive_bias"]
            c14_bg_mean = get_14c_bg(start_date,
                                     end_date,
                                     method="MHD",
                                     bg_bias=bg_bias,
                                    )
        except:
            c14_bg_mean = get_14c_bg(start_date,
                                     end_date,
                                     method="MHD",
                                    )

        bg_isoflux = c14_bg_mean * data_dict["bc_mod"].values
        

        # Apply nuclear power plant correction
        if nuclear_correction is False:
            d14co2_atm_sim = (isoflux_s + bg_isoflux)/disequil_isoflux

        else:
            print("Not yet an option")
            d14co2_atm_sim = (isoflux_s + bg_isoflux)/disequil_isoflux
        
        data_dict[site]["ymod_c14"] = d14co2_atm_sim + (data_dict[site]["bc_mod"].values * 0)
        data_dict[site]["ymodBC_c14"] = bg_isoflux + (data_dict[site]["bc_mod"].values * 0)

    return data_dict, resp_isoflux
