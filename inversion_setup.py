# inversion_setup.py
# Created: 17 June 2024
# Author: Atmospheric Chemistry Research Group, University of Bristol
# *****************************************************************************
# Inversion setup functions for MCMC 
# *****************************************************************************

import numpy as np
import pandas as pd 
import xarray as xr
from openghg.retrieve import get_obs_surface 

def offset_matrix(siteindicator):
    """
    Set up a matrix that can be used to add an offset to each site.
    This will anchor to the first site (i.e. first site has no offset)
    -----------------------------------
    Args:
      siteindicator (array):
        Array of values used for indicating the indices associated
        with each site used in the inversion
    -----------------------------------
    """
    b = np.zeros((int(len(siteindicator)), int(max(siteindicator)) + 1))
    for i in range(int(max(siteindicator) + 1)):
        b[siteindicator == i, i] = 1.0

    return b


def sigma_freq_indicies(ytime, sigma_freq):
    """
    Create an index that splits times
    into given periods
    -----------------------------------
    Args:
      ytime (array of datetime64):
        concatanted array of time values for observations
      sigma_freq (str):
        either "monthly", a pandas format string ("30D"), or None
        this is the period of time to divide the time array into

    Returns:
      output (array):
        index array that defines periods against time
    -----------------------------------
    """
    ydt = pd.to_datetime(ytime)
    output = np.zeros(shape=len(ytime)).astype(int)
    if sigma_freq is None:
        # output already all 0's as expected for this setting
        pass
    elif sigma_freq.lower() == "monthly":
        months = ydt.month
        years = ydt.year
        months_u = np.unique(months)
        years_u = np.unique(years)

        # incrementally set sigma indicies for each month in each year
        count = 0
        for y in years_u:
            for m in months_u:
                indicies = (years == y) & (months == m)
                output[indicies] = count
                count += 1
    else:
        # divide the time between t0 and ti by sigma_freq, then floor
        # to calculate number of integer intervals the calculation is
        # performed in seconds as division by pd time_delta is not allowed
        time_delta = pd.to_timedelta(sigma_freq)
        fractional_freq_time = (ydt - np.amin(ydt)).total_seconds() / time_delta.total_seconds()
        output[:] = np.floor(fractional_freq_time.values).astype(int)

    return output

def opends(fn):
    """
    Open a netcdf dataset with xarray
    -----------------------------------
    Args:
      fn (str):
        Netcdf file to be opened

    Returns:
        xarray.Dataset:
            netcdf file as  dataset
    -----------------------------------
    """
    with xr.open_dataset(fn) as load:
        ds = load.load()

        return ds

def monthly_bcs(start_date, end_date, site, fp_data):
    """
    Creates a sensitivity matrix (H-matrix) for the boundary
    conditions, which will map monthly boundary condition
    scalings to the observations. This is for a single site.
    -----------------------------------
    Args:
      start_date (str):
        Start time of inversion "YYYY-mm-dd"
      end_date (str):
        End time of inversion "YYYY-mm-dd"
      site (str):
        Site that you're creating it for
      fp_data (dict):
        Output from utils..bc_sensitivity

    Returns:
      hmbc (array):
        Sensitivity matrix by month for observations
    -----------------------------------
    """
    allmonth = pd.date_range(start_date, end_date, freq="MS")[:-1]
    nmonth = len(allmonth)
    curtime = pd.to_datetime(fp_data[site].time.values).to_period("M")
    pmonth = pd.to_datetime(fp_data[site].resample(time="MS").mean().time.values)
    hmbc = np.zeros((4 * nmonth, len(fp_data[site].time.values)))
    count = 0
    for cord in range(4):
        for m in range(0, nmonth):
            if allmonth[m] not in pmonth:
                count += 1
                continue
            mnth = allmonth[m].month
            yr = allmonth[m].year
            mnthloc = np.where(np.logical_and(curtime.month == mnth, curtime.year == yr))[0]
            hmbc[count, mnthloc] = fp_data[site].H_bc.values[cord, mnthloc]
            count += 1

    return hmbc

def create_bc_sensitivity(start_date, end_date, site, fp_data, freq):
    """
    Creates a sensitivity matrix (H-matrix) for the boundary
    conditions, which will map boundary condition scalings to
    the observations. This is for a single site. The frequency
    that the boundary condition sensitivity is specified over
    must be given in days. Currently only works for a
    boundary condition from each cardinal direction.
    -----------------------------------
    Args:
      start_date (str):
        Start time of inversion "YYYY-mm-dd"
      end_date (str):
        End time of inversion "YYYY-mm-dd"
      site (str):
        Site that you're creating it for
      fp_data (dict):
        Output from ModelScenario()
      freq (str):
        Length-scale over which boundary condition sensitivities are
        specified over. Specified as in pandas, e.g. "30D".

    Returns:
      hmbc (array):
        Sensitivity matrix by for observations to boundary conditions
    -----------------------------------
    """
    dys = int("".join([s for s in freq if s.isdigit()]))
    alldates = pd.date_range(
        pd.to_datetime(start_date), pd.to_datetime(end_date) + pd.DateOffset(days=dys), freq=freq
    )
    ndates = np.sum(alldates < pd.to_datetime(end_date))
    curdates = fp_data[site].time.values
    hmbc = np.zeros((4 * ndates, len(fp_data[site].time.values)))
    count = 0
    for cord in range(4):
        for m in range(0, ndates):
            dateloc = np.where(
                np.logical_and(
                    curdates >= alldates[m].to_datetime64(), curdates < alldates[m + 1].to_datetime64()
                )
            )[0]
            if len(dateloc) == 0:
                count += 1
                continue
            hmbc[count, dateloc] = fp_data[site].H_bc.values[cord, dateloc]
            count += 1

    return hmbc