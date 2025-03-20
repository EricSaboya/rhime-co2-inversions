# *****************************************************************************
# sensitivity.py
# Created: 28 Feb 2025
# Author: Eric Saboya, School of Geographical Sciences, University of Bristol
# *****************************************************************************
# fp_sensitivity and bc_sensitivity functions
# *****************************************************************************


import sys
import glob
import json
import numpy as np
import pandas as pd
import xarray as xr
import dask.array as da
from pathlib import Path
from types import SimpleNamespace
from typing import Optional, Union

import convert
import calculate_basis_functions as cbf
from importlib import reload

sys.path.append("/user/work/wz22079/projects/CO2/inversions_mk2/")

from utils import combine_datasets, load_json, synonyms

def fp_sensitivity(data_dict: dict, 
                   domain: str, 
                   basis_case: str, 
                   basis_directory=None, 
                   verbose=True
                  )-> dict:
    """
    -------------------------------------------------------
    The fp_sensitivity function adds a sensitivity matrix, 
    H, to each site xarray dataframe in data_dict.
    
    Basis function data is in an array: lat, lon, no. regions.
    In each 'region'element of array there is a lat-lon grid 
    with 1 in region and 0 outside region.

    Region numbering must start from 1
    -------------------------------------------------------
    Args:
        data_dict (dict):
            Output from get_mf_obs_sims() function. 
            Dictionary of datasets.
        
          domain (str):
            Model domain name (str)
        
        basis_case (str, defaults to None):
            Basis case to read in. Examples of basis cases are 
            "NESW","stratgrad".
            String if only one basis case is required. 
            Dict if there are multiple
            sources that require separate basis cases. 
            n which case, keys in dict should
            reflect keys in emissions_name dict used in flux_dict.
      
        basis_directory (str, defaults to None):
            basis_directory can be specified if files are not 
            in the default directory. Must point to a directory 
            which contains subfolders organized by domain.

    Returns:
        dict (xarray.Dataset):
            Same format as data_dict with sensitivity matrix 
            and basis function grid added.
    -------------------------------------------------------
    """
    # List of sites
    sites = [key for key in list(data_dict.keys()) if key[0] != "."]

    # List of flux sectors
    flux_sources = list(data_dict[".flux"].keys())

    # Reads in fp basis function: array w/ dim[sector, lat, lon, time]
    basis_func = cbf.basis(domain=domain,
                           basis_case=basis_case,
                           basis_directory=basis_directory
                          )
        
    if "sector" not in basis_func.coords:
        print(("No sector info in basis function file, so assuming single basis grid "+
               "applies to all sectors"))
        
        for i, source in enumerate(flux_sources):
            if i == 0:
                basis_func_new = np.expand_dims(basis_func["basis"].astype(float), axis=0)
            else:
                basis_func_new = np.concatenate((basis_func_new,
                                                 np.expand_dims(basis_func["basis"].astype(float),axis=0)),
                                                 axis=0)

        basis_func["basis"] = xr.DataArray(data=basis_func_new,
                                           dims=["sector", "lat", "lon", "time"],
                                           coords={"sector": flux_sources,
                                                   "lat": basis_func.lat.values,
                                                   "lon": basis_func.lon.values,
                                                   "time": basis_func.time.values})

    for site in sites:
        for si, source in enumerate(flux_sources):
            if source in basis_func["sector"].values:
                source_ind = np.where(basis_func["sector"].values == source)[0]
                basis_func_source = basis_func["basis"][source_ind][0]
            else:
                print(f"Using %s as the basis case for {source}" %basis_func["sector"].values[0])
                basis_func_source = basis_func["basis"][0]

            
            if "fp_HiTRes" in list(data_dict[site].keys()): 
                site_bf = xr.Dataset({"fp_HiTRes": data_dict[site]["fp_HiTRes"],
                                      "fp": data_dict[site]["fp"]})
            else:
                site_bf = xr.Dataset({"fp": data_dict[site]["fp"]})
                

            if len(flux_sources) == 1:
                H_all_si = data_dict[site]["Hall"]
            elif len(flux_sources) > 1:
                H_all_si = data_dict[site][f"Hall_{source}"]

            H_all_v = H_all_si.values.reshape((len(site_bf.lat)*len(site_bf.lon), len(site_bf.time)))

            if "region" in list(basis_func.dims.keys()):
                if "time" in basis_func.basis.dims:
                    basis_func = basis_func.isel(time=0)

                site_bf = xr.merge([site_bf, basis_func_source])

                H = np.zeros((len(site_bf.region), len(site_bf.time)))
                Herr = np.zeros((len(site_bf.region), len(site_bf.time)))
                
                base_v = site_bf.basis.values.reshape((len(site_bf.lat) * len(site_bf.lon), len(site_bf.region)))

                for i in range(len(site_bf.region)):
                    H[i, :] = np.nansum(H_all_v * base_v[:, i, np.newaxis], axis=0)

                    s_ln = np.nanstd(np.log(np.abs(H_all_v * base_v[:, i, np.newaxis])), axis=0)
                    coeff_var_ln = np.sqrt(np.exp(s_ln**2)-1)
                    # coeff_var = np.nanstd(H_all_v * base_v[:, i, np.newaxis], axis=0)/np.nanmean(H_all_v * base_v[:, i, np.newaxis], axis=0)
                    
                    Herr[i, :] = np.abs(np.nan_to_num(coeff_var_ln))
                
                if source == "all":
                    if (sys.version_info < (3,0)):
                        region_name = site_bf.region
                    else:
                        region_name = site_bf.region.decode("ascii")
                else:
                    if (sys.version_info < (3,0)):
                        region_name = [source + "-" + reg for reg in site_bf.region.values]
                    else:
                        region_name = [source + "-" + reg.decode("ascii") for reg in site_bf.region.values]

                sens_coords = [("region", region_name), ("time", data_dict[site].coords["time"])]
                sensitivity = xr.DataArray(H, coords=sens_coords)
                sensitivity_err = xr.DataArray(Herr, coords=sens_coords)

            else:
                print("Warning: Using basis functions without a region dimension may be deprecated shortly.")

                site_bf = combine_datasets(site_bf, basis_func_source, method="nearest")

                H = np.zeros((int(np.max(site_bf.basis)), len(site_bf.time)))
                Herr = np.zeros((int(np.max(site_bf.basis)), len(site_bf.time)))

                basis_scale = xr.Dataset({"basis_scale": (["lat", "lon", "time"], np.zeros(np.shape(site_bf.basis)))},
                                         coords = site_bf.coords)
                site_bf = site_bf.merge(basis_scale)

                base_v = np.ravel(site_bf.basis.values[:, : ,0])
                for i in range(int(np.max(site_bf.basis))):
                    wh_ri = np.where(base_v == i + 1)
                    H[i, :] = np.nansum(H_all_v[wh_ri[0], :], axis=0)

                    s_ln = np.nanstd(np.log(np.abs(H_all_v[wh_ri[0], :])), axis=0)
                    coeff_var_ln = np.sqrt(np.exp(s_ln**2)-1)

                    # coeff_var = np.nanstd(H_all_v[wh_ri[0], :], axis=0)/np.nanmean(H_all_v[wh_ri[0], :], axis=0)
                    Herr[i, :] = np.abs(np.nan_to_num(coeff_var_ln))

                if source == "all":
                    region_name = list(range(1, np.max(site_bf.basis.values) + 1))
                else:
                    region_name = [source + "-" + str(reg) for reg in range(1, int(np.max(site_bf.basis.values) + 1))]

                sens_coords = {"region": (["region"], region_name),
                               "time": (["time"], data_dict[site].coords["time"].data),
                              }
                sens_dims = ["region", "time"]
                sensitivity = xr.DataArray(H, coords=sens_coords, dims=sens_dims)
                sensitivity_err = xr.DataArray(Herr, coords=sens_coords, dims=sens_dims)                


            if si == 0:
                concat_sensitivity = sensitivity
                concat_sensitivity_err = sensitivity_err
            else:
                concat_sensitivity = xr.concat((concat_sensitivity, sensitivity), dim="region")
                concat_sensitivity_err = xr.concat((concat_sensitivity_err, sensitivity_err), dim="region")

            sub_basis_cases = 0
            if source in basis_func["sector"].values:
                source_ind = np.where(basis_func["sector"].values == source)[0]
                basis_case_key = basis_func["sector"][source_ind]
                    
            elif "all" in basis_case.keys():
                source_ind = 0
                basis_case_key = "all"

        data_dict[site]["H"] = concat_sensitivity
        data_dict[site]["Herr"] = concat_sensitivity_err
        data_dict[".basis"] = basis_func["basis"]

    return data_dict


def bc_sensitivity(data_dict, 
                   domain, 
                   basis_case, 
                   bc_basis_directory=None,
                  ):
    """
    -------------------------------------------------------
    The bc_sensitivity adds H_bc to the sensitivity matrix,
    to each site xarray dataframe in fp_and_data.
    -------------------------------------------------------
    Args:
        data_dict (dict):
            Output from get_mf_obs_sims() function. 
            Dictionary of datasets.
        
        domain (str):
            Model domain name (str)
     
        basis_case (str):
           Basis case to read in. 
           Examples of basis cases are "NESW","stratgrad".
       
        bc_basis_directory (str):
           bc_basis_directory can be specified if files are 
           not in the default directory. Must point to a 
           directory which contains subfolders organized
           by domain. (optional)

    Returns:
      dict (xarray.Dataset):
        Same format as data_dict with HBc 
        sensitivity matrix added.
    -------------------------------------------------------
    """ 
    sites = [key for key in list(data_dict.keys()) if key[0] != "."]

    basis_func = cbf.basis_boundary_conditions(domain=domain, 
                                               basis_case=basis_case, 
                                               bc_basis_directory=bc_basis_directory
                                              )
    # sort basis_func into time order
    ind = basis_func.time.argsort()
    timenew = basis_func.time[ind]
    basis_func = basis_func.reindex({"time": timenew})

    species_info = load_json(filename="species_info.json")

    species = data_dict[sites[0]].attrs["species"]
    species = synonyms(species, species_info)

    for site in sites:
        # ES commented out line below as .bc not attribute. Also assume openghg adds all relevant particle data to file.
        #        if fp_and_data[site].bc.chunks is not None:
        for particles in [
            "particle_locations_n",
            "particle_locations_e",
            "particle_locations_s",
            "particle_locations_w",
        ]:
            data_dict[site][particles] = data_dict[site][particles].compute()

        # compute any chemical loss to the BCs, use lifetime or else set loss to 1 (no loss)
        if "lifetime" in species_info[species].keys():
            lifetime = species_info[species]["lifetime"]
            lifetime_hrs_list_or_float = convert.convert_to_hours(lifetime)

            # calculate the lifetime_hrs associated with each time point in fp_and_data
            # this is because lifetime can be a list of monthly values

            time_month = data_dict[site].time.dt.month
            if type(lifetime_hrs_list_or_float) is list:
                lifetime_hrs = [lifetime_hrs_list_or_float[item - 1] for item in time_month.values]
            else:
                lifetime_hrs = lifetime_hrs_list_or_float

            loss_n = np.exp(-1 * data_dict[site].mean_age_particles_n/lifetime_hrs).rename("loss_n")
            loss_e = np.exp(-1 * data_dict[site].mean_age_particles_e/lifetime_hrs).rename("loss_e")
            loss_s = np.exp(-1 * data_dict[site].mean_age_particles_s/lifetime_hrs).rename("loss_s")
            loss_w = np.exp(-1 * data_dict[site].mean_age_particles_w/lifetime_hrs).rename("loss_w")
            
        else:
            loss_n = data_dict[site].particle_locations_n.copy()
            loss_e = data_dict[site].particle_locations_e.copy()
            loss_s = data_dict[site].particle_locations_s.copy()
            loss_w = data_dict[site].particle_locations_w.copy()
            loss_n[:] = 1
            loss_e[:] = 1
            loss_s[:] = 1
            loss_w[:] = 1

        DS_particle_loc = xr.Dataset(
            {
                "particle_locations_n": data_dict[site]["particle_locations_n"],
                "particle_locations_e": data_dict[site]["particle_locations_e"],
                "particle_locations_s": data_dict[site]["particle_locations_s"],
                "particle_locations_w": data_dict[site]["particle_locations_w"],
                "loss_n": loss_n,
                "loss_e": loss_e,
                "loss_s": loss_s,
                "loss_w": loss_w,
            }
        )
        #                                 "bc":fp_and_data[site]["bc"]})

        DS_temp = combine_datasets(DS_particle_loc, data_dict[".bc"].data, method="ffill")

        DS = combine_datasets(DS_temp, basis_func, method="ffill")

        DS = DS.transpose("height", "lat", "lon", "region", "time")

        part_loc = np.hstack(
            [
                DS.particle_locations_n,
                DS.particle_locations_e,
                DS.particle_locations_s,
                DS.particle_locations_w,
            ]
        )

        loss = np.hstack([DS.loss_n, DS.loss_e, DS.loss_s, DS.loss_w])

        vmr_ed = np.hstack([DS.vmr_n, DS.vmr_e, DS.vmr_s, DS.vmr_w])

        bf = np.hstack([DS.bc_basis_n, DS.bc_basis_e, DS.bc_basis_s, DS.bc_basis_w])

        H_bc = np.zeros((len(DS.coords["region"]), len(DS["particle_locations_n"]["time"])))

        for i in range(len(DS.coords["region"])):
            reg = bf[:, :, i, :]
            H_bc[i, :] = np.nansum((part_loc * loss * vmr_ed * reg), axis=(0, 1))

        sensitivity = xr.Dataset(
            {"H_bc": (["region_bc", "time"], H_bc)},
            coords={"region_bc": (DS.coords["region"].values), "time": (DS.coords["time"])},
        )

        data_dict[site] = data_dict[site].merge(sensitivity)

    return data_dict