# utils.py 

import os
import re
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

from openghg.analyse import ModelScenario
from openghg.analyse import combine_datasets as openghg_combine_datasets

import convert
import calculate_basis_functions as cbf
from importlib import reload

def open_ds(path, chunks=None, combine=None):
    """
    Function efficiently opens xarray datasets.
    -----------------------------------
    Args:
      path (str):
      chunks (dict, optional):
        size of chunks for each dimension
        e.g. {'lat': 50, 'lon': 50}
        opens dataset with dask, such that it is opened 'lazily'
        and all of the data is not loaded into memory
        defaults to None - dataset is opened with out dask
      combine (str, optional):
        Way in which the data should be combined (if using chunks), either:
        'by_coords': order the datasets before concatenating (default)
        'nested': concatenate datasets in the order supplied

    Returns:
      ds (xarray)
    -----------------------------------
    """
    if chunks is not None:
        combine = "by_coords" if combine is None else combine
        ds = xr.open_mfdataset(path, chunks=chunks, combine=combine)
    else:
        # use a context manager, to ensure the file gets closed after use
        with xr.open_dataset(path) as ds:
            ds.load()

    return ds


def read_netcdfs(files, dim="time", chunks=None, verbose=True):
    """
    The read_netcdfs function uses xarray to open sequential netCDF files and
    and concatenates them along the specified dimension.
    Note: this function makes sure that file is closed after open_dataset call.
    -----------------------------------
    Args:
      files (list):
        List of netCDF filenames.
      dim (str, optional):
        Dimension of netCDF to use for concatenating the files.
        Default = "time".
      chunks (dict):
        size of chunks for each dimension
        e.g. {'lat': 50, 'lon': 50}
        opens dataset with dask, such that it is opened 'lazily'
        and all of the data is not loaded into memory
        defaults to None - dataset is opened with out dask

    Returns:
      xarray.Dataset:
        All files open as one concatenated xarray.Dataset object
    -----------------------------------
    """
    if verbose:
        print("Reading and concatenating files ...")
        for fname in files:
            print(fname)

    datasets = [open_ds(p, chunks=chunks) for p in sorted(files)]

    # reindex all of the lat-lon values to a common one to prevent floating point error differences
    with xr.open_dataset(files[0]) as temp:
        fields_ds = temp.load()
    fp_lat = fields_ds["lat"].values
    fp_lon = fields_ds["lon"].values

    datasets = [
        ds.reindex(indexers={"lat": fp_lat, "lon": fp_lon}, method="nearest", tolerance=1e-5)
        for ds in datasets
    ]

    combined = xr.concat(datasets, dim)

    return combined

def indexesMatch(dsa, dsb):
    """
    Check if two datasets need to be reindexed_like for combine_datasets
    -----------------------------------
    Args:
      dsa (xarray.Dataset) :
        First dataset to check
      dsb (xarray.Dataset) :
        Second dataset to check

    Returns:
      boolean:
        True if indexes match, False if datasets must be reindexed
    -----------------------------------
    """

    commonIndicies = [key for key in dsa.indexes.keys() if key in dsb.indexes.keys()]

    # test if each comon index is the same
    for index in commonIndicies:
        # first check lengths are the same to avoid error in second check
        if not len(dsa.indexes[index]) == len(dsb.indexes[index]):
            return False

        # check number of values that are not close (testing for equality with floating point)
        if index == "time":
            # for time iverride the default to have ~ second precision
            rtol = 1e-10
        else:
            rtol = 1e-5

        num_not_close = np.sum(
            ~np.isclose(
                dsa.indexes[index].values.astype(float),
                dsb.indexes[index].values.astype(float),
                rtol=rtol,
            )
        )
        if num_not_close > 0:
            return False

    return True

def combine_datasets(
    dataset_a: xr.Dataset,
    dataset_b: xr.Dataset,
    method: str | None = "nearest",
    tolerance: float | None = None,
) -> xr.Dataset:
    """Merges two datasets and re-indexes to the first dataset.

    If "fp" variable is found within the combined dataset,
    the "time" values where the "lat", "lon" dimensions didn't match are removed.

    NOTE: this is temporary solution while waiting for `.load()` to be added to openghg version of combine_datasets

    Args:
        dataset_a: First dataset to merge
        dataset_b: Second dataset to merge
        method: One of None, nearest, ffill, bfill.
                See xarray.DataArray.reindex_like for list of options and meaning.
                Defaults to ffill (forward fill)
        tolerance: Maximum allowed tolerance between matches.

    Returns:
        xarray.Dataset: Combined dataset indexed to dataset_a
    """
    return openghg_combine_datasets(dataset_a, dataset_b.load(), method=method, tolerance=tolerance)

# def combine_datasets(dsa, dsb, method="ffill", tolerance=None):
#     """
#     The combine_datasets function merges two datasets and re-indexes
#     to the FIRST dataset. If "fp" variable is found within the combined
#     dataset, the "time" values where the "lat","lon"dimensions didn't
#     match are removed.

#     Example:
#         ds = combine_datasets(dsa, dsb)
#     -----------------------------------
#     Args:
#       dsa (xarray.Dataset):
#         First dataset to merge
#       dsb (xarray.Dataset):
#         Second dataset to merge
#       method (str, optional):
#         One of {None, ‘nearest’, ‘pad’/’ffill’, ‘backfill’/’bfill’}
#         See xarray.DataArray.reindex_like for list of options and meaning.
#         Default = "ffill" (forward fill)
#       tolerance (int/float??):
#         Maximum allowed tolerance between matches.

#     Returns:
#       xarray.Dataset:
#         Combined dataset indexed to dsa
#     -----------------------------------
#     """
#     # merge the two datasets within a tolerance and remove times that are NaN (i.e. when FPs don't exist)

#     if not indexesMatch(dsa, dsb):
#         dsb_temp = dsb.reindex_like(dsa, method, tolerance=tolerance)
#     else:
#         dsb_temp = dsb

#     ds_temp = dsa.merge(dsb_temp)
#     if "fp" in list(ds_temp.keys()):
#         flag = np.where(np.isfinite(ds_temp.fp.mean(dim=["lat", "lon"]).values))
#         ds_temp = ds_temp[dict(time=flag[0])]
#     return ds_temp


def synonyms(search_string, info, alternative_label="alt"):
    """
     Check to see if there are other names that we should be using for
     a particular input. E.g. If CFC-11 or CFC11 was input,
     go on to use cfc-11, as this is used in species_info.json
     -----------------------------------
     Args:
       search_string (str):
         Input string that you're trying to match
       info (dict):
         Dictionary whose keys are the "default" values, and an
         variable that contains other possible names

    Returns:
         corrected string
     -----------------------------------
    """
    keys = list(info.keys())

    # First test whether site matches keys (case insensitive)
    out_strings = [k for k in keys if k.upper() == search_string.upper()]

    # If not found, search synonyms
    if len(out_strings) == 0:
        for k in keys:
            matched_strings = [s for s in info[k][alternative_label] if s.upper() == search_string.upper()]
            if len(matched_strings) != 0:
                out_strings = [k]
                break

    if len(out_strings) == 1:
        out_string = out_strings[0]
    else:
        out_string = None

    return out_string

def load_json(filename):
    """Load a JSON file from the internal data directory.

    Args:
        filename (str): Filename
    Returns:
        dict
    """
    # data_folder = Path(__file__).parent.joinpath("data")
    data_folder = Path("/user/home/wz22079/my_openghg/supplementary_data/openghg_defs/data")
    filepath = data_folder.joinpath(filename)
    return json.loads(filepath.read_text())

def filtering(datasets_in, filters, keep_missing=False):
    """
    Applies time filtering to entire dataset.
    Filters supplied in a list and then applied in order.
    For example if you wanted a daily, daytime average, you could do this:

        datasets_dictionary = filtering(datasets_dictionary,
                                    ["daytime", "daily_median"])

    The order of the filters reflects the order they are applied, so for
    instance when applying the "daily_median" filter if you only wanted
    to look at daytime values the filters list should be
    ["daytime","daily_median"]
    -----------------------------------
    Args:
      datasets_in (dict):
        Output from ModelScenario.footprints_merge(). Dictionary of datasets.
      filters (list):
        Filters to apply to the datasets.
          All options are:
            "daytime"           : selects data between 1100 and 1500 local solar time
            "daytime9to5"       : selects data between 0900 and 1700 local solar time
            "nighttime"         : Only b/w 23:00 - 03:00 inclusive
            "noon"              : Only 12:00 fp and obs used
            "daily_median"      : calculates the daily median
            "pblh"              : Only keeps times when pblh is > 50m away from the obs height
            "local_influence"   : Only keep times when localness is low
            "six_hr_mean"       :
            "local_lapse"       :
      keep_missing (bool) : Whether to reindex to retain missing data.

    Returns:
       Same format as datasets_in : Datasets with filters applied.
    -----------------------------------
    """
    if type(filters) is not list:
        filters = [filters]

    datasets = datasets_in.copy()

    def local_solar_time(dataset):
        """
        Returns hour of day as a function of local solar time
        relative to the Greenwich Meridian.
        """
        sitelon = dataset.release_lon.values[0]
        # convert lon to [-180,180], so time offset is negative west of 0 degrees
        if sitelon > 180:
            sitelon = sitelon - 360.0
        dataset["time"] = dataset.time + pd.Timedelta(minutes=float(24 * 60 * sitelon / 360.0))
        hours = dataset.time.to_pandas().index.hour
        return hours

    def local_ratio(dataset):
        """
        Calculates the local ratio in the surrounding grid cells
        """
        dlon = dataset.lon[1].values - dataset.lon[0].values
        dlat = dataset.lat[1].values - dataset.lat[0].values
        local_sum = np.zeros((len(dataset.mf)))

        for ti in range(len(dataset.mf)):
            release_lon = dataset.release_lon[ti].values
            release_lat = dataset.release_lat[ti].values
            wh_rlon = np.where(abs(dataset.lon.values - release_lon) < dlon / 2.0)
            wh_rlat = np.where(abs(dataset.lat.values - release_lat) < dlat / 2.0)
            if np.any(wh_rlon[0]) and np.any(wh_rlat[0]):
                local_sum[ti] = np.sum(
                    dataset.fp[
                        wh_rlat[0][0] - 2 : wh_rlat[0][0] + 3, wh_rlon[0][0] - 2 : wh_rlon[0][0] + 3, ti
                    ].values
                ) / np.sum(dataset.fp[:, :, ti].values)
            else:
                local_sum[ti] = 0.0

        return local_sum

    # Filter functions
    def daily_median(dataset, keep_missing=False):
        """Calculate daily median"""
        if keep_missing:
            return dataset.resample(indexer={"time": "1D"}).median()
        else:
            return dataset.resample(indexer={"time": "1D"}).median().dropna(dim="time")

    def six_hr_mean(dataset, keep_missing=False):
        """Calculate six-hour median"""
        if keep_missing:
            return dataset.resample(indexer={"time": "6H"}).mean()
        else:
            return dataset.resample(indexer={"time": "6H"}).mean().dropna(dim="time")

    def daytime(dataset, site, keep_missing=False):
        """Subset during daytime hours (11:00-15:00)"""
        hours = local_solar_time(dataset)
        ti = [i for i, h in enumerate(hours) if h >= 11 and h <= 15]

        if keep_missing:
            dataset_temp = dataset[dict(time=ti)]
            dataset_out = dataset_temp.reindex_like(dataset)
            return dataset_out
        else:
            return dataset[dict(time=ti)]

    def daytime9to5(dataset, site, keep_missing=False):
        """Subset during daytime hours (9:00-17:00)"""
        hours = local_solar_time(dataset)
        ti = [i for i, h in enumerate(hours) if h >= 9 and h <= 17]

        if keep_missing:
            dataset_temp = dataset[dict(time=ti)]
            dataset_out = dataset_temp.reindex_like(dataset)
            return dataset_out
        else:
            return dataset[dict(time=ti)]

    def nighttime(dataset, site, keep_missing=False):
        """Subset during nighttime hours (23:00 - 03:00)"""
        hours = local_solar_time(dataset)
        ti = [i for i, h in enumerate(hours) if h >= 23 or h <= 3]

        if keep_missing:
            dataset_temp = dataset[dict(time=ti)]
            dataset_out = dataset_temp.reindex_like(dataset)
            return dataset_out
        else:
            return dataset[dict(time=ti)]

    def noon(dataset, site, keep_missing=False):
        """Select only 12pm data"""
        hours = local_solar_time(dataset)
        ti = [i for i, h in enumerate(hours) if h == 12]

        if keep_missing:
            dataset_temp = dataset[dict(time=ti)]
            dataset_out = dataset_temp.reindex_like(dataset)
            return dataset_out
        else:
            return dataset[dict(time=ti)]

    def local_influence(dataset, site, keep_missing=False):
        """
        Subset for times when local influence is below threshold.
        Local influence expressed as a fraction of the sum of entire footprint domain.
        """
        if not dataset.filter_by_attrs(standard_name="local_ratio"):
            lr = local_ratio(dataset)
        else:
            lr = dataset.local_ratio

        pc = 0.1
        ti = [i for i, local_ratio in enumerate(lr) if local_ratio <= pc]
        if keep_missing is True:
            mf_data_array = dataset.mf
            dataset_temp = dataset.drop("mf")

            dataarray_temp = mf_data_array[dict(time=ti)]

            mf_ds = xr.Dataset(
                {"mf": (["time"], dataarray_temp)}, coords={"time": (dataarray_temp.coords["time"])}
            )

            dataset_out = combine_datasets(dataset_temp, mf_ds, method=None)
            return dataset_out
        else:
            return dataset[dict(time=ti)]

    def pblh(dataset, keep_missing=False):
        """
        Subset for times when observations are taken at a height more than
        50m away from (above or below) the PBLH.
        """

        ti = [
            i for i, pblh in enumerate(dataset.PBLH) if np.abs(float(dataset.inlet_height_magl) - pblh) > 50.0
        ]

        if len(ti) != 0:
            if keep_missing is True:
                mf_data_array = dataset.mf
                dataset_temp = dataset.drop("mf")

                dataarray_temp = mf_data_array[dict(time=ti)]

                mf_ds = xr.Dataset(
                    {"mf": (["time"], dataarray_temp)}, coords={"time": (dataarray_temp.coords["time"])}
                )

                dataset_out = combine_datasets(dataset_temp, mf_ds, method=None)
                return dataset_out
            else:
                return dataset[dict(time=ti)]

        else:
            print("PBLH filtering removed all datapoints so this filter is not applied to this site.")

    filtering_functions = {"daily_median": daily_median,        
                           "daytime": daytime,
                           "daytime9to5": daytime9to5,
                           "nighttime": nighttime,
                           "noon": noon,
                           "local_influence": local_influence,
                           "six_hr_mean": six_hr_mean,
                           "pblh": pblh,
                          }

    # Get list of sites
    sites = [key for key in list(datasets.keys()) if key[0] != "."]

    # Apply filtering
    for site in sites:
        for filt in filters:
            n_nofilter = datasets[site].time.values.shape[0]
            if filt in ["daily_median", "six_hr_mean", "pblh"]:
                datasets[site] = filtering_functions[filt](datasets[site], keep_missing=keep_missing)
            else:
                datasets[site] = filtering_functions[filt](datasets[site], site, keep_missing=keep_missing)
            n_filter = datasets[site].time.values.shape[0]
            n_dropped = n_nofilter - n_filter
            perc_dropped = np.round(n_dropped / n_nofilter * 100, 2)
            print(f"{filt} filter removed {n_dropped} ({perc_dropped} %) obs at site {site}")

    return datasets


def timeseries_HiTRes(flux_dict,
                      fp_HiTRes_ds=None,
                      fp_file=None,
                      output_TS=True,
                      output_fpXflux=True,
                      output_type="Dataset",
                      output_file=None,
                      verbose=False,
                      chunks=None,
                      time_resolution="1H",
                     ):
    """
    The timeseries_HiTRes function computes flux * HiTRes footprints.

    HiTRes footprints record the footprint at each 2 hour period back
    in time for the first 24 hours. Need a high time resolution flux
    to multiply the first 24 hours back of footprints. Need a residual
    flux to multiply the residual integrated footprint for the remainder
    of the 30 day period.
    -----------------------------------
    Args:
      fp_HiTRes_ds (xarray.Dataset)
        Dataset of high time resolution footprints. HiTRes footprints
        record the footprint at each timestep back in time for a given
        amount of time (e.g. hourly time steps back in time for the first
        24 hours).
      domain (str)
        Domain name. The footprint files should be sub-categorised by the domain.
      flux_dict (dict)
        This should be a dictionary of the form output in the format
        flux_dict: {'high_freq': flux_dataset, 'low_freq': flux_dataset}.
        This is because this function needs two time resolutions of fluxes as
        explained in the header.

        If there are multiple sectors, the format should be:
        flux_dict: {sector1 : {'high_freq' : flux_dataset, 'low_freq'  : flux_dataset},
                    sector2 : {'high_freq' : flux_dataset, 'low_freq'  : flux_dataset}}
      output_TS (bool)
        Output the timeseries. Default is True.
      output_fpXflux (bool)
        Output the sensitivity map. Default is True.
      verbose (bool)
        Show progress bar throughout loop
      chunks (dict)
       Size of chunks for each dimension
       e.g. {'lat': 50, 'lon': 50}
       opens dataset with dask, such that it is opened 'lazily'
       and all of the data is not loaded into memory
       defaults to None - dataset is opened with out dask

    Returns:
      xarray.Dataset or dict
        Same format as flux_dict['high_freq']:
        If flux_dict['high_freq'] is an xarray.Dataset then an xarray.Dataset is returned
        If flux_dict['high_freq'] is a dict of xarray.Datasets then a dict of xarray.Datasets
        is returned (an xarray.Dataset for each sector)

        If output_TS is True:
          Outputs the timeseries
        If output_fpXflux is True:
          Outputs the sensitivity map
    -----------------------------------
    """
    if verbose:
        print(f"\nCalculating timeseries with {time_resolution} resolution, this might take a few minutes")
        
    # Retrieve HiTRes footprint
    if fp_HiTRes_ds is None and fp_file is None:
        raise("Must provide either a footprint Dataset or footprint filename")
        
    elif fp_HiTRes_ds is None:
        fp_HiTRes_ds = read_netcdfs(fp_file, 
                                    chunks=chunks,
                                   )
        fp_HiTRes = fp_HiTRes_ds.fp_HiTRes
        
    else:
        if isinstance(fp_HiTRes_ds, xr.DataArray):
            fp_HiTRes = fp_HiTRes_ds
        elif fp_HiTRes_ds.chunks is None and chunks is not None:
            fp_HiTRes = fp_HiTRes_ds.fp_HiTRes.chunk(chunks)
        else:
            fp_HiTRes = fp_HiTRes_ds.fp_HiTRes

    # Resample fp to match the required time resolution
    fp_HiTRes = fp_HiTRes.resample(time=time_resolution).ffill()

    # Get H_back timestep and max number of hours back
    H_back_hour_diff = int(fp_HiTRes["H_back"].diff(dim="H_back").values.mean())
    max_H_back = int(fp_HiTRes["H_back"].values[-2])

    # Create time array to loop through, with the required resolution
    time_array = fp_HiTRes.time.values
    # Extract, as a dask array, to make the loop quicker
    fp_HiTRes = da.array(fp_HiTRes)

    # This is the number of hours over which the H_back will be resampled to match time_resolution
    H_resample = (
        int(time_resolution[0])
        if H_back_hour_diff == 1
        else 1
        if H_back_hour_diff == int(time_resolution[0])
        else None
    )
    if H_resample is None:
        print("Cannot resample H_back")
        return None
    # reverse the H_back coordinate to be chronological, and resample to match time_resolution
    fp_HiTRes = fp_HiTRes[:, :, :, ::-H_resample]

    # convert fluxes to dictionaries with sectors as keys
    flux = (
        {"total": flux_dict}
        if any([ff in list(flux_dict.keys()) for ff in ["high_freq", "low_freq"]])
        else flux_dict
    )
    flux = {
        sector: {freq: None if flux_freq is None else flux_freq for freq, flux_freq in flux_sector.items()}
        for sector, flux_sector in flux.items()
    }

    # extract the required time data
    flux = {
        sector: {
            freq: flux_freq.sel(
                time=slice(fp_HiTRes_ds.time[0] - np.timedelta64(max_H_back, "h"), fp_HiTRes_ds.time[-1])
            ).flux
            if flux_freq is not None
            else None
            for freq, flux_freq in flux_sector.items()
        }
        for sector, flux_sector in flux.items()
    }

    for sector, flux_sector in flux.items():
        if "high_freq" in flux_sector.keys() and flux_sector["high_freq"] is not None:
            # reindex the high frequency data to match the fp
            time_flux = np.arange(
                fp_HiTRes_ds.time[0].values - np.timedelta64(max_H_back, "h"),
                fp_HiTRes_ds.time[-1].values + np.timedelta64(time_resolution[0], time_resolution[1].lower()),
                time_resolution[0],
                dtype=f"datetime64[{time_resolution[1].lower()}]",
            )
            flux_sector["high_freq"] = flux_sector["high_freq"].reindex(time=time_flux, method="ffill")
        else:
            print(
                f"\nWarning: no high frequency flux data for {sector}, estimating a timeseries using the low frequency data"
            )
            flux_sector["high_freq"] = None

        if "low_freq" not in flux_sector.keys() or flux_sector["low_freq"] is None:
            print(f"\nWarning: no low frequency flux data for {sector}, resampling from high frequency data")
            flux_sector["low_freq"] = flux_sector["high_freq"].resample(time="MS").mean()

    # convert to array to use in numba loop
    flux = {
        sector: {
            freq: None
            if flux_freq is None
            else flux_freq.values
            if flux_freq.chunks is None
            else da.array(flux_freq)
            for freq, flux_freq in flux_sector.items()
        }
        for sector, flux_sector in flux.items()
    }

    # Set up a numpy array to calculate the product of the footprint (H matrix) with the fluxes
    if output_fpXflux:
        fpXflux = {
            sector: da.zeros((len(fp_HiTRes_ds.lat), len(fp_HiTRes_ds.lon), len(time_array)))
            for sector in flux.keys()
        }

    elif output_TS:
        timeseries = {sector: da.zeros(len(time_array)) for sector in flux.keys()}

    # month and year of the start of the data - used to index the low res data
    start = {dd: getattr(np.datetime64(time_array[0], "h").astype(object), dd) for dd in ["month", "year"]}

    # put the time array into tqdm if we want a progress bar to show throughout the loop
    iters = tqdm(time_array) if verbose else time_array
    ### iterate through the time coord to get the total mf at each time step using the H back coord
    # at each release time we disaggregate the particles backwards over the previous 24hrs
    for tt, time in enumerate(iters):
        # get 4 dimensional chunk of high time res footprint for this timestep
        # units : mol/mol/mol/m2/s
        fp_time = fp_HiTRes[:, :, tt, :]

        # get the correct index for the low res data
        # estimated using the difference between the current and start month and year
        current = {dd: getattr(np.datetime64(time, "h").astype(object), dd) for dd in ["month", "year"]}
        tt_low = current["month"] - start["month"] + 12 * (current["year"] - start["year"])

        # select the high res emissions for the corresponding 24 hours
        # if there aren't any high frequency data it will select from the low frequency data
        # this is so that we can compare emissions data with different resolutions e.g. ocean species
        emissions = {
            sector: flux_sector["high_freq"][:, :, tt : tt + fp_time.shape[2] - 1]
            if flux_sector["high_freq"] is not None
            else flux_sector["low_freq"][:, :, tt_low]
            for sector, flux_sector in flux.items()
        }
        # add an axis if the emissions is array is 2D so that it can be multiplied by the fp
        emissions = {
            sector: em_sec[:, :, np.newaxis] if len(em_sec.shape) == 2 else em_sec
            for sector, em_sec in emissions.items()
        }
        # select average monthly emissions for the start of the month
        emissions_end = {
            sector: flux_sector["low_freq"][:, :, tt_low] for sector, flux_sector in flux.items()
        }

        # Multiply the HiTRes footprint with the HiTRes emissions to give mf
        # we take all but the slice for H_back==24 as these are the hourly disaggregated fps
        # flux units : mol/m2/s;       fp units : mol/mol/mol/m2/s
        # --> mol/mol/mol/m2/s * mol/m2/s === mol / mol
        fpXflux_time = {sector: em_sec * fp_time[:, :, 1:] for sector, em_sec in emissions.items()}
        # multiply the monthly flux by the residual fp, at H_back==24
        fpXflux_end = {sector: em_end * fp_time[:, :, 0] for sector, em_end in emissions_end.items()}
        # append the residual emissions
        fpXflux_time = {
            sector: np.dstack((fp_fl, fpXflux_end[sector])) for sector, fp_fl in fpXflux_time.items()
        }

        for sector, fp_fl in fpXflux_time.items():
            if output_fpXflux:
                # Sum over time (H back) to give the total mf at this timestep
                fpXflux[sector][:, :, tt] = np.nansum(fp_fl, axis=2)

            elif output_TS:
                # work out timeseries by summing over lat, lon, & time (24 hrs)
                timeseries[sector][tt] = np.nansum(fp_fl)

    if output_fpXflux and output_TS:
        # if not already done then calculate the timeseries
        timeseries = {sector: fp_fl.sum(axis=(0, 1)) for sector, fp_fl in fpXflux.items()}

    if output_fpXflux:
        fpXflux = (
            {sec: (["lat", "lon", "time"], ff_sec) for sec, ff_sec in fpXflux.items()}
            if output_type.lower() == "dataset"
            else fpXflux
        )

        fpXflux = (
            xr.Dataset(
                fpXflux,
                coords={"lat": fp_HiTRes_ds.lat.values, "lon": fp_HiTRes_ds.lon.values, "time": time_array},
            )
            if output_type.lower() == "dataset"
            else {
                sector: xr.DataArray(
                    data=fpXflux_sector,
                    dims=["lat", "lon", "time"],
                    coords={
                        "lat": fp_HiTRes_ds.lat.values,
                        "lon": fp_HiTRes_ds.lon.values,
                        "time": time_array,
                    },
                )
                for sector, fpXflux_sector in fpXflux.items()
            }
            if output_type.lower() == "dataarray"
            else fpXflux
        )

        if output_type.lower() == "dataset":
            fpXflux = fpXflux if fpXflux.chunks is None else fpXflux.compute()
        else:
            fpXflux = {sec: ff if ff.chunks is None else ff.compute() for sec, ff in fpXflux.items()}

        if output_type.lower() == "dataarray" and list(flux.keys()) == ["total"]:
            fpXflux = fpXflux["total"]

    if output_fpXflux and not output_TS:
        return fpXflux

    else:
        # for each sector create a tuple of ['time'] (coord name) and the timeseries
        # if the output required is a dataset
        timeseries = (
            {sec: (["time"], ts_sec) for sec, ts_sec in timeseries.items()}
            if output_type.lower() == "dataset"
            else timeseries
        )
        timeseries = (
            xr.Dataset(timeseries, coords={"time": time_array})
            if output_type.lower() == "dataset"
            else {
                sector: xr.DataArray(data=ts_sector, dims=["time"], coords={"time": time_array})
                for sector, ts_sector in timeseries.items()
            }
            if output_type.lower() == "dataarray"
            else timeseries
        )

        if output_type.lower() == "dataset":
            timeseries = timeseries if timeseries.chunks is None else timeseries.compute()
        else:
            timeseries = {sec: tt if tt.chunks is None else tt.compute() for sec, tt in timeseries.items()}

        if output_type.lower() == "dataarray" and list(flux.keys()) == ["total"]:
            timeseries = timeseries["total"]

        if output_file is not None and output_type.lower() == "dataset":
            print(f"Saving to {output_file}")
            timeseries.to_netcdf(output_file)
        elif output_file is not None:
            print(f"output type must be dataset to save to file")

        if output_fpXflux:
            return timeseries, fpXflux
        elif output_TS:
            return timeseries



# def fp_sensitivity(data_dict: dict, 
#                    domain: str, 
#                    basis_case: str, 
#                    basis_directory=None, 
#                    verbose=True
#                   ):
#     """
#     The fp_sensitivity function adds a sensitivity matrix, H, to each
#     site xarray dataframe in fp_and_data.
#     Basis function data in an array: lat, lon, no. regions.
#     In each 'region'element of array there is a lat-lon grid with 1 in
#     region and 0 outside region.

#     Region numbering must start from 1
#     -----------------------------------
#     Args:
#       data_dict (dict):
#         Output from footprints_data_merge() function. Dictionary of datasets.
#       domain (str):
#         Domain name. The footprint files should be sub-categorised by the domain.
#       basis_case:
#         Basis case to read in. Examples of basis cases are "NESW","stratgrad".
#         String if only one basis case is required. Dict if there are multiple
#         sources that require separate basis cases. In which case, keys in dict should
#         reflect keys in emissions_name dict used in fp_data_merge.
#       basis_directory (str):
#         basis_directory can be specified if files are not in the default
#         directory. Must point to a directory which contains subfolders organized
#         by domain. (optional)

#     Returns:
#         dict (xarray.Dataset):
#           Same format as fp_and_data with sensitivity matrix and basis function grid added.
#     -----------------------------------
#     """
#     import calculate_basis_functions as cbf 
    
#     # List of sites
#     sites = [key for key in list(data_dict.keys()) if key[0] != "."]

#     # List of flux sectors
#     flux_sources = list(data_dict[".flux"].keys())

#     # Reads in fp basis function: array w/ dim[sector, lat, lon, time]
#     basis_func = cbf.basis(domain=domain,
#                            basis_case=basis_case,
#                            basis_directory=basis_directory
#                           )
        
#     if "sector" not in basis_func.coords:
#         print(("No sector info in basis function file, so assuming single basis grid "+
#                "applies to all sectors"))
        
#         for i, source in enumerate(flux_sources):
#             if i == 0:
#                 basis_func_new = np.expand_dims(basis_func["basis"].astype(float), axis=0)
#             else:
#                 basis_func_new = np.concatenate((basis_func_new,
#                                                  np.expand_dims(basis_func["basis"].astype(float),axis=0)),
#                                                 axis=0)

#         basis_func["basis"] = xr.DataArray(data=basis_func_new,
#                                            dims=["sector", "lat", "lon", "time"],
#                                            coords={"sector": flux_sources,
#                                                    "lat": basis_func.lat.values,
#                                                    "lon": basis_func.lon.values,
#                                                    "time": basis_func.time.values})

#     for site in sites:
#         for si, source in enumerate(flux_sources):
#             if source in basis_func["sector"].values:
#                 source_ind = np.where(basis_func["sector"].values == source)[0]
#                 basis_func_source = basis_func["basis"][source_ind][0]
#             else:
#                 print(f"Using %s as the basis case for {source}" %basis_func["sector"].values[0])
#                 basis_func_source = basis_func["basis"][0]

            
#             if "fp_HiTRes" in list(data_dict[site].keys()): 
#                 site_bf = xr.Dataset({"fp_HiTRes": data_dict[site]["fp_HiTRes"],
#                                       "fp": data_dict[site]["fp"]})
#             else:
#                 site_bf = xr.Dataset({"fp": data_dict[site]["fp"]})
                

#             if len(flux_sources) == 1:
#                 H_all_si = data_dict[site]["Hall"]
#             elif len(flux_sources) > 1:
#                 H_all_si = data_dict[site][f"Hall_{source}"]

#             H_all_v = H_all_si.values.reshape((len(site_bf.lat) * len(site_bf.lon), len(site_bf.time)))

#             if "region" in list(basis_func.dims.keys()):
#                 if "time" in basis_func.basis.dims:
#                     basis_func = basis_func.isel(time=0)

#                 site_bf = xr.merge([site_bf, basis_func_source])

#                 H = np.zeros((len(site_bf.region), len(site_bf.time)))
#                 base_v = site_bf.basis.values.reshape((len(site_bf.lat) * len(site_bf.lon), len(site_bf.region)))

#                 for i in range(len(site_bf.region)):
#                     H[i, :] = np.nansum(H_all_v * base_v[:, i, np.newaxis], axis=0)
                    
#                 if source == "all":
#                     if (sys.version_info < (3, 0)):
#                         region_name = site_bf.region
#                     else:
#                         region_name = site_bf.region.decode("ascii")
#                 else:
#                     if (sys.version_info < (3, 0)):
#                         region_name = [source + "-" + reg for reg in site_bf.region.values]
#                     else:
#                         region_name = [source + "-" + reg.decode("ascii") for reg in site_bf.region.values]

#                 sensitivity = xr.DataArray(H,
#                                            coords=[("region", region_name),
#                                                    ("time", data_dict[site].coords["time"])])

#             else:
#                 print("Warning: Using basis functions without a region dimension may be deprecated shortly.")

#                 site_bf = combine_datasets(site_bf, basis_func_source, method="nearest")

#                 H = np.zeros((int(np.max(site_bf.basis)), len(site_bf.time)))
#                 Herr = np.zeros((int(np.max(site_bf.basis)), len(site_bf.time)))

#                 basis_scale = xr.Dataset({"basis_scale": (["lat", "lon", "time"], np.zeros(np.shape(site_bf.basis)))},
#                                          coords = site_bf.coords)
#                 site_bf = site_bf.merge(basis_scale)

#                 base_v = np.ravel(site_bf.basis.values[:, : ,0])
#                 for i in range(int(np.max(site_bf.basis))):
#                     wh_ri = np.where(base_v == i + 1)
#                     H[i, :] = np.nansum(H_all_v[wh_ri[0], :], axis=0)
                    
#                     coeff_var = np.nanstd(H_all_v[wh_ri[0], :], axis=0)/np.nanmean(H_all_v[wh_ri[0], :], axis=0)
#                     Herr[i, :] = coeff_var

#                 if source == "all":
#                     region_name = list(range(1, np.max(site_bf.basis.values) + 1))
#                 else:
#                     region_name = [source + "-" + str(reg) for reg in range(1, int(np.max(site_bf.basis.values) + 1))]

#                 sens_coords = {"region": (["region"], region_name),
#                                "time": (["time"], data_dict[site].coords["time"].data),
#                               }
                
#                 # sensitivity = xr.DataArray(H, coords=[("region", region_name), ("time", data_dict[site].coords["time"].data)])
#                 # sensitivity_err = xr.DataArray(Herr, coords=[("region", region_name), ("time", data_dict[site].coords["time"].data)])                
#                 sensitivity = xr.DataArray(H, coords=sens_coords)
#                 sensitivity_err = xr.DataArray(Herr, coords=sens_coords)                


#             if si == 0:
#                 concat_sensitivity = sensitivity
#                 concat_sensitivity_err = sensitivity_err
#             else:
#                 concat_sensitivity = xr.concat((concat_sensitivity, sensitivity), dim="region")
#                 concat_sensitivity_err = xr.concat((concat_sensitivity_err, sensitivity_err), dim="region")

#             sub_basis_cases = 0
#             if source in basis_func["sector"].values:
#                 source_ind = np.where(basis_func["sector"].values == source)[0]
#                 basis_case_key = basis_func["sector"][source_ind]
                    
#             elif "all" in basis_case.keys():
#                 source_ind = 0
#                 basis_case_key = "all"

#         data_dict[site]["H"] = np.abs(np.nan_to_num(concat_sensitivity))
#         data_dict[site]["Herr"] = concat_sensitivity_err
#         data_dict[".basis"] = basis_func["basis"]

#     return data_dict


# def bc_sensitivity(data_dict, domain, basis_case, bc_basis_directory=None):
#     """
#     The bc_sensitivity adds H_bc to the sensitivity matrix,
#     to each site xarray dataframe in fp_and_data.
#     -----------------------------------
#     Args:
#       data_dict (dict):
#         Output from ModelScenario.footprints_data_merge() function. Dictionary of datasets.
#      domain (str):
#        Domain name. The footprint files should be sub-categorised by the domain.
#      basis_case (str):
#        Basis case to read in. Examples of basis cases are "NESW","stratgrad".
#      bc_basis_directory (str):
#        bc_basis_directory can be specified if files are not in the default
#        directory. Must point to a directory which contains subfolders organized
#        by domain. (optional)

#     Returns:
#       dict (xarray.Dataset):
#         Same format as fp_and_data with sensitivity matrix added.
#     -----------------------------------
#     """ 
#     sites = [key for key in list(data_dict.keys()) if key[0] != "."]

#     basis_func = cbf.basis_boundary_conditions(domain=domain, 
#                                                basis_case=basis_case, 
#                                                bc_basis_directory=bc_basis_directory
#                                               )
#     # sort basis_func into time order
#     ind = basis_func.time.argsort()
#     timenew = basis_func.time[ind]
#     basis_func = basis_func.reindex({"time": timenew})

#     species_info = load_json(filename="species_info.json")

#     species = data_dict[sites[0]].attrs["species"]
#     species = synonyms(species, species_info)

#     for site in sites:
#         # ES commented out line below as .bc not attribute. Also assume openghg adds all relevant particle data to file.
#         #        if fp_and_data[site].bc.chunks is not None:
#         for particles in [
#             "particle_locations_n",
#             "particle_locations_e",
#             "particle_locations_s",
#             "particle_locations_w",
#         ]:
#             data_dict[site][particles] = data_dict[site][particles].compute()

#         # compute any chemical loss to the BCs, use lifetime or else set loss to 1 (no loss)
#         if "lifetime" in species_info[species].keys():
#             lifetime = species_info[species]["lifetime"]
#             lifetime_hrs_list_or_float = convert.convert_to_hours(lifetime)

#             # calculate the lifetime_hrs associated with each time point in fp_and_data
#             # this is because lifetime can be a list of monthly values

#             time_month = data_dict[site].time.dt.month
#             if type(lifetime_hrs_list_or_float) is list:
#                 lifetime_hrs = [lifetime_hrs_list_or_float[item - 1] for item in time_month.values]
#             else:
#                 lifetime_hrs = lifetime_hrs_list_or_float

#             loss_n = np.exp(-1 * data_dict[site].mean_age_particles_n / lifetime_hrs).rename("loss_n")
#             loss_e = np.exp(-1 * data_dict[site].mean_age_particles_e / lifetime_hrs).rename("loss_e")
#             loss_s = np.exp(-1 * data_dict[site].mean_age_particles_s / lifetime_hrs).rename("loss_s")
#             loss_w = np.exp(-1 * data_dict[site].mean_age_particles_w / lifetime_hrs).rename("loss_w")
#         else:
#             loss_n = data_dict[site].particle_locations_n.copy()
#             loss_e = data_dict[site].particle_locations_e.copy()
#             loss_s = data_dict[site].particle_locations_s.copy()
#             loss_w = data_dict[site].particle_locations_w.copy()
#             loss_n[:] = 1
#             loss_e[:] = 1
#             loss_s[:] = 1
#             loss_w[:] = 1

#         DS_particle_loc = xr.Dataset(
#             {
#                 "particle_locations_n": data_dict[site]["particle_locations_n"],
#                 "particle_locations_e": data_dict[site]["particle_locations_e"],
#                 "particle_locations_s": data_dict[site]["particle_locations_s"],
#                 "particle_locations_w": data_dict[site]["particle_locations_w"],
#                 "loss_n": loss_n,
#                 "loss_e": loss_e,
#                 "loss_s": loss_s,
#                 "loss_w": loss_w,
#             }
#         )
#         #                                 "bc":fp_and_data[site]["bc"]})

#         DS_temp = combine_datasets(DS_particle_loc, data_dict[".bc"].data, method="ffill")

#         DS = combine_datasets(DS_temp, basis_func, method="ffill")

#         DS = DS.transpose("height", "lat", "lon", "region", "time")

#         part_loc = np.hstack(
#             [
#                 DS.particle_locations_n,
#                 DS.particle_locations_e,
#                 DS.particle_locations_s,
#                 DS.particle_locations_w,
#             ]
#         )

#         loss = np.hstack([DS.loss_n, DS.loss_e, DS.loss_s, DS.loss_w])

#         vmr_ed = np.hstack([DS.vmr_n, DS.vmr_e, DS.vmr_s, DS.vmr_w])

#         bf = np.hstack([DS.bc_basis_n, DS.bc_basis_e, DS.bc_basis_s, DS.bc_basis_w])

#         H_bc = np.zeros((len(DS.coords["region"]), len(DS["particle_locations_n"]["time"])))

#         for i in range(len(DS.coords["region"])):
#             reg = bf[:, :, i, :]
#             H_bc[i, :] = np.nansum((part_loc * loss * vmr_ed * reg), axis=(0, 1))

#         sensitivity = xr.Dataset(
#             {"H_bc": (["region_bc", "time"], H_bc)},
#             coords={"region_bc": (DS.coords["region"].values), "time": (DS.coords["time"])},
#         )

#         data_dict[site] = data_dict[site].merge(sensitivity)

#     return data_dict

def areagrid(lat, lon):
    """
    Calculates grid of areas (m2) given arrays of latitudes and longitudes
    -------------------------------------
    Args:
      lat (array):
        1D array of latitudes
      lon (array):
        1D array of longitudes

    Returns:
      area (array):
        2D array of areas of of size lat x lon
    -------------------------------------
    Example:
      import utils.areagrid
      lat=np.arange(50., 60., 1.)
      lon=np.arange(0., 10., 1.)
      area=utils.areagrid(lat, lon)
    """

    re = 6367500.0  # radius of Earth in m

    dlon = abs(np.mean(lon[1:] - lon[0:-1])) * np.pi / 180.0
    dlat = abs(np.mean(lat[1:] - lat[0:-1])) * np.pi / 180.0
    theta = np.pi * (90.0 - lat) / 180.0

    area = np.zeros((len(lat), len(lon)))

    for latI in range(len(lat)):
        if theta[latI] == 0.0 or np.isclose(theta[latI], np.pi):
            area[latI, :] = (re**2) * abs(np.cos(dlat / 2.0) - np.cos(0.0)) * dlon
        else:
            lat1 = theta[latI] - dlat / 2.0
            lat2 = theta[latI] + dlat / 2.0
            area[latI, :] = (re**2) * (np.cos(lat1) - np.cos(lat2)) * dlon

    return area

def get_country(domain, country_file=None):
    if country_file is None:
        if not os.path.exists(os.path.join(openghginv_path, "countries/")):
            os.makedirs(os.path.join(openghginv_path, "countries/"))
            raise FileNotFoundError(
                "Country definition file not found." f" Please add to {openghginv_path}/countries/"
            )
        else:
            country_directory = os.path.join(openghginv_path, "countries/")

        filenames = glob.glob(os.path.join(country_directory, f"country_{domain}.nc"))
        filename = filenames[0]
    else:
        filename = country_file

    with xr.open_dataset(filename) as f:
        lon = f.variables["lon"][:].values
        lat = f.variables["lat"][:].values

        # Get country indices and names
        if "country" in f.variables:
            country = f.variables["country"][:, :]
        elif "region" in f.variables:
            country = f.variables["region"][:, :]
        else:
            raise ValueError(f"Variables 'country' or 'region' not found in country file {filename}.")

        #         if (ukmo is True) or (uk_split is True):
        #             name_temp = f.variables['name'][:]
        #             f.close()
        #             name=np.asarray(name_temp)

        #         else:
        name = f.variables["name"].values.astype(str)

    result = dict(
        lon=lon,
        lat=lat,
        lonmax=np.max(lon),
        lonmin=np.min(lon),
        latmax=np.max(lat),
        latmin=np.min(lat),
        country=np.asarray(country),
        name=name,
    )
    return SimpleNamespace(**result)