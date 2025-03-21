# inversion_mcmc.py
# Created: 17 June 2024
# Author: Atmospheric Chemistry Research Group, University of Bristol
# *****************************************************************************
# About
#   Functions for performing RHIME inversion.
# *****************************************************************************
import os
import re
import sys
import getpass
import pymc as pm
import arviz as az
import numpy as np
import pandas as pd
import xarray as xr
import pytensor.tensor as pt

from scipy import stats
from pathlib import Path
from typing import Optional, Union

import convert
import utils
from inversion_setup import offset_matrix

# from openghg_inversions.hbmcmc.hbmcmc_output import define_output_filename
# from openghg_inversions.config.version import code_version

def get_function_dict():
    """
    Returns dictionary of PyMC distributions
    """
    functiondict = {"uniform": pm.Uniform,
                    "flat": pm.Flat,
                    "halfflat": pm.HalfFlat,
                    "normal": pm.Normal,
                    "truncatednormal": pm.TruncatedNormal,
                    "halfnormal": pm.HalfNormal,
                    "skewnormal": pm.SkewNormal,
                    "beta": pm.Beta,
                    "kumaraswamy": pm.Kumaraswamy,
                    "exponential": pm.Exponential,
                    "laplace": pm.Laplace,
                    "studentt": pm.StudentT,
                    "halfstudentt": pm.HalfStudentT,
                    "cauchy": pm.Cauchy,
                    "halfcauchy": pm.HalfCauchy,
                    "gamma": pm.Gamma,
                    "inversegamma": pm.InverseGamma,
                    "weibull": pm.Weibull,
                    "lognormal": pm.Lognormal,
                    "chisquared": pm.ChiSquared,
                    "wald": pm.Wald,
                    "pareto": pm.Pareto,
                    "exgaussian": pm.ExGaussian,
                    "vonmises": pm.VonMises,
                    "triangular": pm.Triangular,
                    "gumbel": pm.Gumbel,
                    "rice": pm.Rice,
                    "logistic": pm.Logistic,
                    "logitnormal": pm.LogitNormal,
                    "interpolated": pm.Interpolated,
                   }
    return functiondict


def parseprior(name, 
               prior_params, 
               shape=(),
               sigma_sf=0,
              ):
    """
    Parses all continuous distributions for PyMC 3.8:
    https://docs.pymc.io/api/distributions/continuous.html
    This format requires updating when the PyMC distributions update,
    but is safest for code execution
    -----------------------------------
    Args:
      name (str):
        name of variable in the pymc model
      prior_params (dict):
        dict of parameters for the distribution,
        including 'pdf' for the distribution to use
      shape (array):
        shape of distribution to be created.
        Default shape = () is the same as used by PyMC3
      sigma_sf (float):
          Scale factor for sigma value of the pdf
    -----------------------------------
    """
    functiondict = get_function_dict()

    pdf = prior_params["pdf"]
    # Get a dictionary of the pdf arguments
    params = {x: prior_params[x] for x in prior_params if x != "pdf"}
    if "sigma" in params.keys():
        params["sigma"] = np.sqrt(params["sigma"]**2 + sigma_sf**2)
    return functiondict[pdf.lower()](name, shape=shape, **params)


def inferpymc(Hx,
              Hxerr,
              basis_region_mask,
              Y,
              error,
              Ymodelerror,
              siteindicator,
              sigma_freq_index,
              Hbc,
              xprior,
              bcprior,
              sigprior,
              nit=1e4,
              burn=1500,
              tune=3000,
              nchain=2,
              sigma_per_site=True,
              offsetprior=None,
              add_offset=False,
              verbose=False,
              save_trace=False,
              use_bc=True,
             ):
    """
    Uses PyMC for Bayesian inference of: 
    - Scaling of flux field
    - Scaling of boundary conditions (optional)
    - Model error scaling value
    This uses a Normal likelihood but the (hyper)prior PDFs can selected by user.
    -----------------------------------
    Args:
      Hx (array):
        Transpose of the sensitivity matrix to map emissions to measurement.
        This is the same as what is given from fp_data[site].H.values, where
        fp_data is the output from e.g. footprint_data_merge, but where it
        has been stacked for all sites.
      Hbc (array):
        Same as above but for boundary conditions
      basis_region_mask (array):
        Mask to indicate which Hx elements are for each flux sector
      Y (array):
        Measurement vector containing all measurements
      error (arrray):
        Measurement error vector, containg a value for each element of Y.
      siteindicator (array):
        Array of indexing integers that relate each measurement to a site
      sigma_freq_index (array):
        Array of integer indexes that converts time into periods
      xprior (dict):
        Dictionary containing information about the prior PDF for emissions.
        The entry "pdf" is the name of the analytical PDF used, see
        https://docs.pymc.io/api/distributions/continuous.html for PDFs
        built into pymc3, although they may have to be coded into the script.
        The other entries in the dictionary should correspond to the shape
        parameters describing that PDF as the online documentation,
        e.g. N(1,1**2) would be: xprior={pdf:"normal", "mu":1, "sd":1}.
        Note that the standard deviation should be used rather than the
        precision. Currently all variables are considered iid.
      bcprior (dict):
        Same as above but for boundary conditions.
      sigprior (dict):
        Same as above but for model error.
      offsetprior (dict):
        Same as above but for bias offset. Only used is addoffset=True.
      sigma_per_site (bool):
        Whether a model sigma value will be calculated for each site independantly (True) or all sites together (False).
        Default: True
      add_offset (bool):
        Add an offset (intercept) to all sites but the first in the site list. Default False.
      min_model_error (float):
        Minimum model error to impose on species baseline 
      verbose:
        When True, prints progress bar

    Returns:
      outs (array):
        MCMC chain for emissions scaling factors for each basis function.
      bcouts (array):
        MCMC chain for boundary condition scaling factors.
      sigouts (array):
        MCMC chain for model error.
      Ytrace (array):
        MCMC chain for modelled obs.
      YBCtrace (array):
        MCMC chain for modelled boundary condition.
      convergence (str):
        Passed/Failed convergence test as to whether mutliple chains
        have a Gelman-Rubin diagnostic value <1.05
      step1 (str):
        Type of MCMC sampler for emissions and boundary condition updates.
        Currently it's hardwired to NUTS (probably wouldn't change this
        unless you're doing something obscure).
      step2 (str):
        Type of MCMC sampler for model error updates.
        Currently it's hardwired to a slice sampler. This parameter is low
        dimensional and quite simple with a slice sampler, although could
        easily be changed.

    TO DO:
       - Allow non-iid variables
    -----------------------------------
    """
    burn = int(burn)     # No. of iterations to discard in MCMC
    hx = Hx.T            # Transpose of matrix of dosages with dimemsions of [t, region[sector]]
    hxerr = Hxerr.T      # Transpose of matrix of coefficient of variability of dosages with dimensions of [t, region[sector]]
    nx = hx.shape[1]     # No. basis functions for all sectors stacked 
    ny = len(Y)          # No. time steps
    nit = int(nit)       # Total no. of iterations in MCMC 

    if use_bc is True:
        hbc = Hbc.T
        nbc = hbc.shape[1]
        
    nflux = len(list(xprior.keys()))   # No. flux sectors

    # Convert siteindicator into a site indexer
    if sigma_per_site:
        _sites = siteindicator.astype(int)
        nsites = np.amax(_sites) + 1
    else:
        _sites = np.zeros_like(siteindicator).astype(int)
        nsites = 1
    nsigmas = np.amax(sigma_freq_index) + 1

    if add_offset:
        B = offset_matrix(siteindicator)

    with pm.Model() as model:
        # ---- Flux hyperparameter ---- #
        all_distributions_x = []
        nsec_inds = 0
        hxerr_tmean = np.nanmean(hxerr, axis=0)
        for i, key in enumerate(xprior.keys()):
            sector_len = len(np.squeeze(np.where(basis_region_mask == i + 1)))
            sector_mask = np.squeeze(np.where(basis_region_mask == i + 1))
            print(f"Sector {key} is using {sector_len} basis functions")

            # group_distributions = [parseprior(f"{key}_{j}", xprior[key], sigma_sf=hxerr_tmean[int(j+nsec_inds)]) for j in range(sector_len)]
            group_distributions = [parseprior(f"{key}_{j}", xprior[key]) for j in range(sector_len)]
            all_distributions_x.extend(group_distributions)
            nsec_inds += sector_len
    

        hx_dot_x = pt.dot(hx, all_distributions_x)
        
        # ---- Model error hyperparameter ---- # 
        sig = parseprior("sig", sigprior, shape=(nsites, nsigmas))

        # ---- Boundary Conditions ---- # 
        if use_bc is True:
            xbc = [parseprior(f"xbc_{j}", bcprior) for j in range(nbc)]
            all_distributions_x.extend(xbc)

            if add_offset:
                offset = parseprior("offset", offsetprior, shape=nsites-1)
                offset_vec = pt.concatenate((np.array([0]), offset), axis=0)
                mu = hx_dot_x + pt.dot(hbc, xbc) + pt.dot(B, offset_vec) 
            else:
                mu = hx_dot_x + pt.dot(hbc, xbc)

        else:
            if add_offset:
                offset = parseprior("offset", offsetprior, shape=nsites-1)
                offset_vec = pt.concatenate((np.array([0]), offset), axis=0)
                mu = hx_dot_x + pt.dot(B, offset_vec) 
            else:
                mu = hx_dot_x 

        model_error = Ymodelerror * sig[_sites, sigma_freq_index]    
        epsilon = pt.sqrt(error**2 + model_error**2 + 0.25)
        y = pm.Normal("y", mu=mu, sigma=np.sqrt(epsilon), observed=Y, shape=ny)

        step1 = pm.NUTS(vars=all_distributions_x)
        step2 = pm.Slice(vars=[sig])
        
        trace = pm.sample(nit, 
                          tune=int(tune), 
                          chains=nchain,
                          step=[step1, step2], 
                          progressbar=False, 
                          cores=nchain) #step=pm.Metropolis())#  #target_accept=0.8,

    xouts = {}
    bcouts = {}

    # ---- Calculate Gelman-Rubin Statistic for each posterior ----- 
    GelmanRubin_dict = {}
    for key in trace.posterior.keys():
        posterior_chain_means = []
        posterior_chain_var = []
        for i in range(nchain):
            posterior_chain_means.append(np.nanmean(trace.posterior[key].values[i, burn:nit]))
            posterior_chain_var.append(np.nanstd(trace.posterior[key].values[i, burn:nit])**2)
            
        # Calculate Gelman-Rubin Convergence Diagnostic 
        B = np.nanstd(posterior_chain_means)**2
        W = np.nanmean(posterior_chain_var)
        L = nit-burn
        GR_stat = (((L-1)/L)*W + (1/L)*B)/W

        tolerance = 0.05
        if np.abs((1-GR_stat)<tolerance):
            GelmanRubin_dict[key] = {"convergence": "passed", 
                                     "GR_value": GR_stat,
                                    }

            # For MCMC posteriors that converge, append 0th chain 
            # values to outs dictionary 
            if key.split("_")[0] in xprior.keys():
                xouts[key] = trace.posterior[key][0, burn:nit]
            elif "xbc" in key:
                bcouts[key] = trace.posterior[key][0, burn:nit]
            elif "sig" in key:
                sigouts = trace.posterior["sig"][0, burn:nit]
            elif "offset" in key:
                offsetouts = trace.posterior["offset"][0, burn:nit]
            
        else:
            GelmanRubin_dict[key] = {"convergence": "failed", 
                                     "GR_value": GR_stat,
                                    }
            print(f"{key} did not achieve a Gelman-Rubin statistic within a 5% tolerance")

            # For MCMC posteriors that don't converge append values of 0.0 to
            # dictionary 
            if key.split("_")[0] in xprior.keys():
                xouts[key] = trace.posterior[key][0, burn:nit] * 0.0
            elif "xbc" in key:
                bcouts[key] = trace.posterior[key][0, burn:nit] * 0.0
            elif "sig" in key:
                sigouts = trace.posterior["sig"][0, burn:nit] * 0.0
            elif "offset" in key:
                offsetouts = trace.posterior["offset"][0, burn:nit] * 0.0


    # Include option for calculating Monte Carlo standard Errors?


    # X TERMS 
    hx_dot_x_posterior = []
    xouts_posterior = {}
    
    # hx_dot_x_posterior has shape [nsector, nbasis, time, nit-burn]
    for i, key in enumerate(xprior.keys()):
        sector_mask = np.squeeze(np.where(basis_region_mask == i+1))
        hx_dot_x_posterior_sector = []
        x_posterior_sector = []
        # Combine each basis function per flux sector with posterior scale factors
        for j in sector_mask:
            k = np.mod(j, len(sector_mask))
            xtrace_posterior_sector_key = f"{key}_{k}"
            hx_sector = np.reshape(hx[:, j], (hx[:, j].shape[0], 1))
            x_posterior = np.reshape(xouts[xtrace_posterior_sector_key].values, (xouts[xtrace_posterior_sector_key].values.shape[0], 1))
            x_posterior_sector.append(xouts[xtrace_posterior_sector_key].values)
            
            hx_dot_x_posterior_sector.append([np.dot(hx_sector, x_posterior.T)])
            
        # Stack all basis functions per flux sector
        # hx_dot_x_posterior should therefore just be [sector1, sector2, ..., sectorN]
        # where ith sector has dimensions [time, basis function]
        hx_dot_x_posterior.append(np.vstack(hx_dot_x_posterior_sector))
        xouts_posterior[key] = np.vstack(x_posterior_sector)

    # CALCULATE PERTURBED TRACE 
    for i in range(len(xprior.keys())):
        if i == 0:
            YPERTtrace = np.sum(hx_dot_x_posterior[0], axis=0)
        elif i!=0:
            YPERTtrace += np.sum(hx_dot_x_posterior[i], axis=0) 
        
    # OFFSETS 
    if add_offset:
        offset_trace = np.hstack([np.zeros((int(nit-burn), 1)), offsetouts])
        OFFSETtrace = np.dot(B, offset_trace.T)

    # BOUNDARY CONDITIONS
    if use_bc: 
        bcouts_array = np.zeros((int(nit-burn), nbc))
        bcouts_array[:, 0] = bcouts["xbc_0"].values
        bcouts_array[:, 1] = bcouts["xbc_1"].values
        bcouts_array[:, 2] = bcouts["xbc_2"].values
        bcouts_array[:, 3] = bcouts["xbc_3"].values
        
        YBCtrace = np.dot(Hbc.T, bcouts_array.T)
        
    # Sum over basis functions, then sectors for regional fluxes 
    if add_offset:
        Ytrace = np.squeeze(YPERTtrace + OFFSETtrace)
    else:
        Ytrace = np.squeeze(YPERTtrace)
        
    
    # Collect outputs into dictionary 
    result = {"xouts": xouts_posterior,
              "sigouts": sigouts,
              "Ytrace": Ytrace,
              "convergence": GelmanRubin_dict, 
              "step1": step1,
              "step2": step2,
             }
    if add_offset:
        result["offset_outs"] = offsetouts
        result["OFFSETtrace"] = OFFSETtrace
        
    if use_bc:
        result["bcouts"] = bcouts_array
        result["YBCtrace"] = YBCtrace
        
    return result


def inferpymc_postprocessouts(mcmc_results,
                              use_bc,
                              mcmc_dict,
                              Hx,
                              Y,
                              error,
                              Ymodelerror,
                              Ytime,
                              siteindicator,
                              sigma_freq_index,
                              domain,
                              species,
                              sites,
                              start_date,
                              end_date,
                              outputname,
                              outputpath,
                              country_unit_prefix,
                              emissions_name,
                              emissions_store,
                              Hbc: Optional[np.ndarray] = None,
                              fp_data=None,
                              country_file=None,
                              rerun_file=None,
                             ):
    """
    Takes the output from inferpymc function, along with some other input
    information, and places it all in a netcdf output. This function also
    calculates the mean posterior emissions for the countries in the
    inversion domain and saves it to netcdf.
    Note that the uncertainties are defined by the highest posterior
    density (HPD) region and NOT percentiles (as the tdMCMC code).
    The HPD region is defined, for probability content (1-a), as:
        1) P(x \in R | y) = (1-a)
        2) for x1 \in R and x2 \notin R, P(x1|y)>=P(x2|y)
    -------------------------------
    Args:
      xouts (array):
        MCMC chain for emissions scaling factors for each basis function.
      bcouts (array):
        MCMC chain for boundary condition scaling factors.
      sigouts (array):
        MCMC chain for model error.
      convergence (str):
        Passed/Failed convergence test as to whether mutliple chains
        have a Gelman-Rubin diagnostic value <1.05
      Hx (array):
        Transpose of the sensitivity matrix to map emissions to measurement.
        This is the same as what is given from fp_data[site].H.values, where
        fp_data is the output from e.g. footprint_data_merge, but where it
        has been stacked for all sites.
      Hbc (array):
        Same as above but for boundary conditions
      Y (array):
        Measurement vector containing all measurements
      error (arrray):
        Measurement error vector, containg a value for each element of Y.
      Ytrace (array):
        Trace of modelled y values calculated from mcmc outputs and H matrices
      YBCtrace (array):
        Trace of modelled boundary condition values calculated from mcmc outputs and Hbc matrices
      step1 (str):
        Type of MCMC sampler for emissions and boundary condition updates.
      step2 (str):
        Type of MCMC sampler for model error updates.
      xprior (dict):
        Dictionary containing information about the prior PDF for emissions.
        The entry "pdf" is the name of the analytical PDF used, see
        https://docs.pymc.io/api/distributions/continuous.html for PDFs
        built into pymc3, although they may have to be coded into the script.
        The other entries in the dictionary should correspond to the shape
        parameters describing that PDF as the online documentation,
        e.g. N(1,1**2) would be: xprior={pdf:"normal", "mu":1, "sd":1}.
        Note that the standard deviation should be used rather than the
        precision. Currently all variables are considered iid.
      bcprior (dict):
        Same as above but for boundary conditions.
      sigprior (dict):
        Same as above but for model error.
      offsetprior (dict):
        Same as above but for bias offset. Only used is addoffset=True.
      Ytime (pandas datetime array):
        Time stamp of measurements as used by the inversion.
      siteindicator (array):
        Numerical indicator of which site the measurements belong to,
        same length at Y.
      sigma_freq_index (array):
        Array of integer indexes that converts time into periods
      domain (str):
        Inversion spatial domain.
      species (str):
        Species of interest
      sites (list):
        List of sites in inversion
      start_date (str):
        Start time of inversion "YYYY-mm-dd"
      end_date (str):
        End time of inversion "YYYY-mm-dd"
      outputname (str):
        Unique identifier for output/run name.
      outputpath (str):
        Path to where output should be saved.
      country_unit_prefix ('str', optional)
        A prefix for scaling the country emissions. Current options are:
        'T' will scale to Tg, 'G' to Gg, 'M' to Mg, 'P' to Pg.
        To add additional options add to acrg_convert.prefix
        Default is none and no scaling will be applied (output in g).
      burn (int):
        Number of iterations burned in MCMC
      tune (int):
        Number of iterations used to tune step size
      nchain (int):
        Number of independent chains run
      sigma_per_site (bool):
        Whether a model sigma value was be calculated for each site independantly (True)
        or all sites together (False).
      fp_data (dict, optional):
        Output from footprints_data_merge + sensitivies
      emissions_name (list, optional):
        Update: Now a list with "source" values as used when adding emissions data to
        the OpenGHG object store.
      basis_directory (str, optional):
        Directory containing basis function file
      country_file (str, optional):
        Path of country definition file
      add_offset (bool):
        Add an offset (intercept) to all sites but the first in the site list. Default False.
      rerun_file (xarray dataset, optional):
        An xarray dataset containing the ncdf output from a previous run of the MCMC code.

    Returns:
        netdf file containing results from inversion
    -------------------------------
    TO DO:
        - Look at compressability options for netcdf output
        - I'm sure the number of inputs can be cut down or found elsewhere.
        - Currently it can only work out the country total emissions if
          the a priori emissions are constant over the inversion period
          or else monthly (and inversion is for less than one calendar year).
    """

    print("Post-processing PyMC output")

    # Get parameters for output file
    nit = mcmc_dict["nit"] - mcmc_dict["burn"]  # No. of MCMC iterations used for inferring posterior
    nx = Hx.shape[0] # No. basis functions for all sectors stacked 
    ny = len(Y) # No. of time steps 
    nui = np.arange(2)
    steps = np.arange(nit)
    nparam = np.arange(nx)
    nmeasure = np.arange(ny)

    # ---- OFFSET PARAMETERS ---- #
    if mcmc_dict["add_offset"] == True:
        noff = mcmc_results["offset_outs"].shape[0]
        nOFF = np.arange(noff)
        OFFSETtrace = mcmc_results["OFFSETtrace"]

        YmodmuOFF = np.mean(OFFSETtrace, axis=1)            # mean
        YmodmedOFF = np.median(OFFSETtrace, axis=1)         # median
        YmodmodeOFF = np.zeros(shape=OFFSETtrace.shape[0])  # mode

        for i in range(0, OFFSETtrace.shape[0]):
            # if sufficient no. of iterations use a KDE to calculate mode
            # else, mean value used in lieu
            if np.nanmax(OFFSETtrace[i, :]) > np.nanmin(OFFSETtrace[i, :]):
                xes_off = np.linspace(np.nanmin(OFFSETtrace[i, :]), np.nanmax(OFFSETtrace[i, :]), 200)
                kde = stats.gaussian_kde(OFFSETtrace[i, :]).evaluate(xes_off)
                YmodmodeOFF[i] = xes_off[kde.argmax()]
            else:
                YmodmodeOFF[i] = np.mean(OFFSETtrace[i, :])

        Ymod95OFF = az.hdi(OFFSETtrace.T, 0.95)
        Ymod68OFF = az.hdi(OFFSETtrace.T, 0.68)

    # ---- Y-BC HYPERPARAMETER ---- #
    if use_bc:
        nbc = Hbc.shape[0]
        nBC = np.arange(nbc)
        YBCtrace = mcmc_results["YBCtrace"]
        bcouts = mcmc_results["bcouts"]

        YmodmuBC = np.mean(YBCtrace, axis=1)                # mean
        YmodmedBC = np.median(YBCtrace, axis=1)             # median
        YmodmodeBC = np.zeros(shape=YBCtrace.shape[0])      # mode

        for i in range(0, YBCtrace.shape[0]):
            # if sufficient no. of iterations use a KDE to calculate mode
            # else, mean value used in lieu
            if np.nanmax(YBCtrace[i, :]) > np.nanmin(YBCtrace[i, :]):
                xes_bc = np.linspace(np.nanmin(YBCtrace[i, :]), np.nanmax(YBCtrace[i, :]), 200)
                kde = stats.gaussian_kde(YBCtrace[i, :]).evaluate(xes_bc)
                YmodmodeBC[i] = xes_bc[kde.argmax()]
            else:
                YmodmodeBC[i] = np.mean(YBCtrace[i, :])

        Ymod95BC = az.hdi(YBCtrace.T, 0.95)
        Ymod68BC = az.hdi(YBCtrace.T, 0.68)
        YaprioriBC = np.sum(Hbc, axis=0)

    # ---- Y-VALUES HYPERPARAMETER (XOUTS * H) ---- #
    Ytrace = mcmc_results["Ytrace"]
    Ymodmu = np.mean(Ytrace, axis=1)
    Ymodmed = np.median(Ytrace, axis=1)
    Ymodmode = np.zeros(shape=Ytrace.shape[0])

    for i in range(0, Ytrace.shape[0]):
        # if sufficient no. of iterations use a KDE to calculate mode
        # else, mean value used in lieu
        if np.nanmax(Ytrace[i, :]) > np.nanmin(Ytrace[i, :]):
            xes = np.arange(np.nanmin(Ytrace[i, :]), np.nanmax(Ytrace[i, :]), 0.05)
            kde = stats.gaussian_kde(Ytrace[i, :]).evaluate(xes)
            Ymodmode[i] = xes[kde.argmax()]
        else:
            Ymodmode[i] = np.mean(Ytrace[i, :])

    Ymod95 = az.hdi(Ytrace.T, 0.95)
    Ymod68 = az.hdi(Ytrace.T, 0.68)

    if use_bc:
        Yapriori = np.sum(Hx.T, axis=1) + np.sum(Hbc.T, axis=1)
    else:
        Yapriori = np.sum(Hx.T, axis=1)

    sitenum = np.arange(len(sites))

    if fp_data is None and rerun_file is not None:
        lon = rerun_file.lon.values
        lat = rerun_file.lat.values
        site_lat = rerun_file.sitelats.values
        site_lon = rerun_file.sitelons.values
        bfds = rerun_file.basisfunctions
    else:
        lon = fp_data[sites[0]].lon.values
        lat = fp_data[sites[0]].lat.values
        site_lat = np.zeros(len(sites))
        site_lon = np.zeros(len(sites))
        for si, site in enumerate(sites):
            site_lat[si] = fp_data[site].release_lat.values[0]
            site_lon[si] = fp_data[site].release_lon.values[0]
        bfds = fp_data[".basis"]

    # ---- Calculate mean and mode posterior scale map and flux field ---- # 
    # NB. Basis field [sector, lat, lon, time]
    scalemap_mu_dict = {}   # Mean scale map dictionary
    scalemap_mode_dict = {} # Mode scale map dictionary 

    xouts = mcmc_results["xouts"]
    xoutsave = np.zeros(shape=(len(nparam), nit))
    nbasis_sectors = []
    
    for i, flux_sector in enumerate(xouts.keys()):
        nbasis_sec = xouts[flux_sector].shape[0]
        
        xoutsave[int(np.sum(nbasis_sectors)): int(np.sum(nbasis_sectors)+nbasis_sec), :] = xouts[flux_sector]
        nbasis_sectors.append(nbasis_sec)
        
        scalemap_mu = np.zeros_like(bfds[i].values)    # Mean scale map
        scalemap_mode = np.zeros_like(bfds[i].values)  # Mode scale map
        
        for npm in range(1+int(bfds[i].values.max())):
            scalemap_mu[np.squeeze(bfds[i].values)==npm] = np.mean(xouts[flux_sector][npm-1, :])

            if np.nanmax(xouts[flux_sector][npm-1, :]) > np.nanmin(xouts[flux_sector][npm-1, :]):
                xes = np.arange(np.nanmin(xouts[flux_sector][npm-1, :]), np.nanmax(xouts[flux_sector][npm-1, :]), 0.01)
                kde = stats.gaussian_kde(xouts[flux_sector][npm-1, :]).evaluate(xes)
                scalemap_mode[np.squeeze(bfds[i].values) == npm] = xes[kde.argmax()]
            else:
                scalemap_mode[np.squeeze(bfds[i].values) == npm] = np.mean(xouts[npm-1, :])

        scalemap_mu_dict[flux_sector] = scalemap_mu[:, :, 0]
        scalemap_mode_dict[flux_sector] = scalemap_mode[:, :, 0]

    xprior = mcmc_dict["xprior"]

    scalemap_mu_flux = np.zeros_like(np.squeeze(bfds.values))
    scalemap_mode_flux = np.zeros_like(np.squeeze(bfds.values))
    for i, key in enumerate(xprior.keys()):
        scalemap_mu_flux[i] = scalemap_mu_dict[key]
        scalemap_mode_flux[i] = scalemap_mode_dict[key]
    
    
    # Get Fluxes 
    if rerun_file is not None:
        flux_array_all = np.expand_dims(rerun_file.fluxapriori.values, 2)
    else:
        if emissions_name is None:
            raise ValueError("Emissions name not provided.")
        else:
            apriori_flux = np.zeros_like(np.squeeze(bfds.values))
            for i, key in enumerate(xprior.keys()):
                t_flux = fp_data[".flux"][key].data["time"]
                if len(t_flux) > 1:
                    # Use time-averaged flux field 
                    apriori_flux[i,:,:] = xr.DataArray.mean(fp_data[".flux"][key].data.flux, dim="time")
                else:
                    apriori_flux[i,:,:] = fp_data[".flux"][key].data.flux.values[:,:,:]            

    
    # if flux_array_all.shape[2] == 1:
    #     print("\nAssuming flux prior is annual and extracting first index of flux array.")
    #     apriori_flux = flux_array_all[:, :, 0]
    # else:
    #     print("\nAssuming flux prior is monthly.")
    #     print(f"Extracting weighted average flux prior from {start_date} to {end_date}")
    #     allmonths = pd.date_range(start_date, end_date, freq="1h").values[:-1]
    #     allmonths -= 1  # to align with zero indexed array

    #     apriori_flux = np.zeros_like(flux_array_all[:, :, 0])

    #     # calculate the weighted average flux across the whole inversion period
    #     for m in np.unique(allmonths):
    #         apriori_flux += flux_array_all[:, :, m] * np.sum(allmonths == m) / len(allmonths)

    aposteriori_flux_mode = scalemap_mode_flux * apriori_flux
    aposteriori_flux_mean = scalemap_mu_flux * apriori_flux

    # Basis functions to save
    bfarray = bfds.values - 1

    # Calculate country totals
    area = utils.areagrid(lat, lon)
    if not rerun_file:
        c_object = utils.get_country(domain, country_file=country_file)
        cntryds = xr.Dataset(
            {"country": (["lat", "lon"], c_object.country), "name": (["ncountries"], c_object.name)},
            coords={"lat": (c_object.lat), "lon": (c_object.lon)},
        )
        cntrynames = cntryds.name.values
        cntrygrid = cntryds.country.values
    else:
        cntrynames = rerun_file.countrynames.values
        cntrygrid = rerun_file.countrydefinition.values

    nsector = len(xprior.keys())
    cntrymean = np.zeros((nsector, len(cntrynames)))
    cntrymedian = np.zeros((nsector, len(cntrynames)))
    cntrymode = np.zeros((nsector, len(cntrynames)))
    cntry68 = np.zeros((nsector, len(cntrynames), len(nui)))
    cntry95 = np.zeros((nsector, len(cntrynames), len(nui)))
    cntrysd = np.zeros((nsector, len(cntrynames)))
    cntryprior = np.zeros((nsector, len(cntrynames)))
        
    molarmass = convert.molar_mass(species)

    unit_factor = convert.prefix(country_unit_prefix)
    if country_unit_prefix is None:
        country_unit_prefix = ""
    country_units = country_unit_prefix + "g"
    if rerun_file is not None:
        obs_units = rerun_file.Yobs.attrs["units"].split(" ")[0]
    else:
        obs_units = str(fp_data[site].mf.attrs['units'])


    fluxsector = []
    for i, key in enumerate(xprior.keys()):
        fluxsector.append(key)
        for ci, cntry in enumerate(cntrynames):
            cntrytottrace = np.zeros(len(steps))
            cntrytotprior = 0
            
            for bf in range(int(np.max(bfarray[i]))):   # bfarray runs from 0 to n-1
                bothinds = np.logical_and(cntrygrid == ci, bfarray[i,:,:,0] == bf)
                cntrytottrace += (
                    np.sum(area[bothinds].ravel() * apriori_flux[i, bothinds].ravel() * 3600 * 24 * 365 * molarmass)
                    * xouts[key][bf, :]
                    / unit_factor
                )
                cntrytotprior += (
                    np.sum(area[bothinds].ravel() * apriori_flux[i, bothinds].ravel() * 3600 * 24 * 365 * molarmass)
                    / unit_factor
                )
            cntrymean[i, ci] = np.mean(cntrytottrace)
            cntrymedian[i, ci] = np.median(cntrytottrace)

            if np.nanmax(cntrytottrace) > np.nanmin(cntrytottrace):
                xes = np.linspace(np.nanmin(cntrytottrace), np.nanmax(cntrytottrace), 200)
                kde = stats.gaussian_kde(cntrytottrace).evaluate(xes)
                cntrymode[i, ci] = xes[kde.argmax()]
            else:
                cntrymode[i, ci] = np.mean(cntrytottrace)

            cntrysd[i, ci] = np.std(cntrytottrace)
            cntry68[i, ci, :] = az.hdi(cntrytottrace, 0.68)
            cntry95[i, ci, :] = az.hdi(cntrytottrace, 0.95)
            cntryprior[i, ci] = cntrytotprior

    # Make convergence results suitable for saving 
    conv_xpdf_gr = []       # Gelman-Rubin value for x
    conv_xpdf_result = []   # Gelman-Rubin result for x 

    conv_sig_gr = []        # Gelman-Rubin value for model error 
    conv_sig_result = []    # Gelman-Rubin result for model error

    conv_bc_gr = []         # Gelman-Rubin value for BCs 
    conv_bc_result = []     # Gelman-Rubin result for BCs

    for i, key in enumerate(mcmc_results['convergence'].keys()):
        if "sig" not in key or "xbc" not in key:
            conv_xpdf_result.append(mcmc_results['convergence'][key]['convergence'])
            conv_xpdf_gr.append(mcmc_results['convergence'][key]['GR_value'])
        elif "sig" in key:
            conv_sig_result.append(mcmc_results['convergence'][key]['convergence'])
            conv_sig_gr.append(mcmc_results['convergence'][key]['GR_value'])
        elif "xbc" in key:
            conv_bc_result.append(mcmc_results['convergence'][key]['convergence'])
            conv_bc_gr.append(mcmc_results['convergence'][key]['GR_value'])

       
    # Make output netcdf file
    data_vars = {
        "Yobs": (["nmeasure"], Y),
        "Yerror": (["nmeasure"], error),
        "Ytime": (["nmeasure"], Ytime),
        "Ymodelerror_prior": (["nmeasure"], Ymodelerror),
        
        "Yapriori": (["nmeasure"], Yapriori),
        "Ymodmean": (["nmeasure"], Ymodmu),
        "Ymodmedian": (["nmeasure"], Ymodmed),
        "Ymodmode": (["nmeasure"], Ymodmode),
        "Ymod95": (["nmeasure", "nUI"], Ymod95),
        "Ymod68": (["nmeasure", "nUI"], Ymod68),
        
        "xtrace": (["nparam", "steps"], xoutsave),
        
        "sigtrace": (["steps", "nsigma_site", "nsigma_time"], mcmc_results["sigouts"].values),
        "siteindicator": (["nmeasure"], siteindicator),
        "sigmafreqindex": (["nmeasure"], sigma_freq_index),
        "sitenames": (["nsite"], sites),
        "sitelons": (["nsite"], site_lon),
        "sitelats": (["nsite"], site_lat),
        "fluxapriori": (["fluxsector", "lat", "lon"], np.squeeze(apriori_flux)),
        "fluxaposteriori_mean": (["fluxsector", "lat", "lon"], np.squeeze(aposteriori_flux_mean)),
        "fluxaposteriori_mode": (["fluxsector", "lat", "lon"], np.squeeze(aposteriori_flux_mode)),
        "scalingmean": (["fluxsector", "lat", "lon"], np.squeeze(scalemap_mu_flux)),
        "scalingmode": (["fluxsector", "lat", "lon"], np.squeeze(scalemap_mode_flux)),
        "basisfunctions": (["fluxsector", "lat", "lon"], np.squeeze(bfarray)),
        "countrymean": (["fluxsector", "countrynames"], cntrymean),
        "countrymedian": (["fluxsector", "countrynames"], cntrymedian),
        "countrymode": (["fluxsector", "countrynames"], cntrymode),
        "countrysd": (["fluxsector", "countrynames"], cntrysd),
        "country68": (["fluxsector", "countrynames", "nUI"], cntry68),
        "country95": (["fluxsector", "countrynames", "nUI"], cntry95),
        "countryapriori": (["fluxsector", "countrynames"], cntryprior),
        "countrydefinition": (["lat", "lon"], cntrygrid),
        "xsensitivity": (["nmeasure", "nparam"], Hx.T),
    }

    coords = {       
        "stepnum": (["steps"], steps),
        "paramnum": (["nlatent"], nparam),
        "measurenum": (["nmeasure"], nmeasure),
        "UInum": (["nUI"], nui),
        "nsites": (["nsite"], sitenum),
        "nsigma_time": (["nsigma_time"], np.unique(sigma_freq_index)),
        "nsigma_site": (["nsigma_site"], np.arange(mcmc_results["sigouts"].shape[1]).astype(int)),
        "lat": (["lat"], lat),
        "lon": (["lon"], lon),
        "fluxsector": (["fluxsector"], fluxsector),
        "nbasis": (["nbasis"], np.array(nbasis_sectors)),
        "countrynames": (["countrynames"], cntrynames),
    }

    if use_bc:
        data_vars.update(
            {
                "YaprioriBC": (["nmeasure"], YaprioriBC),
                "YmodmeanBC": (["nmeasure"], YmodmuBC),
                "YmodmedianBC": (["nmeasure"], YmodmedBC),
                "YmodmodeBC": (["nmeasure"], YmodmodeBC),
                "Ymod95BC": (["nmeasure", "nUI"], Ymod95BC),
                "Ymod68BC": (["nmeasure", "nUI"], Ymod68BC),
                "bctrace": (["steps", "nBC"], bcouts),
                "bcsensitivity": (["nmeasure", "nBC"], Hbc.T),
            }
        )
        coords["numBC"] = (["nBC"], nBC)

    if mcmc_dict["add_offset"] == True:
        data_vars.update(
            {"Yoffmean": (["nmeasure"], YmodmuOFF),
             "Yoffmedian": (["nmeasure"], YmodmedOFF),
             "Yoffmode": (["nmeasure"], YmodmodeOFF),
             "Yoff68": (["nmeasure", "nUI"], Ymod68OFF),
             "Yoff95": (["nmeasure", "nUI"], Ymod95OFF),
            }
        )


    outds = xr.Dataset(data_vars, coords=coords)

    outds.fluxaposteriori_mode.attrs["units"] = "mol/m2/s"
    outds.fluxaposteriori_mean.attrs["units"] = "mol/m2/s"
    outds.fluxapriori.attrs["units"] = "mol/m2/s"
    
    outds.Yobs.attrs["units"] = obs_units + " " + "mol/mol"
    outds.Yapriori.attrs["units"] = obs_units + " " + "mol/mol"
    outds.Ymodmean.attrs["units"] = obs_units + " " + "mol/mol"
    outds.Ymodmedian.attrs["units"] = obs_units + " " + "mol/mol"
    outds.Ymodmode.attrs["units"] = obs_units + " " + "mol/mol"
    outds.Ymod95.attrs["units"] = obs_units + " " + "mol/mol"
    outds.Ymod68.attrs["units"] = obs_units + " " + "mol/mol"
    outds.Yerror.attrs["units"] = obs_units + " " + "mol/mol"
    
    outds.countrymean.attrs["units"] = country_units
    outds.countrymedian.attrs["units"] = country_units
    outds.countrymode.attrs["units"] = country_units
    outds.country68.attrs["units"] = country_units
    outds.country95.attrs["units"] = country_units
    outds.countrysd.attrs["units"] = country_units
    outds.countryapriori.attrs["units"] = country_units
    outds.xsensitivity.attrs["units"] = obs_units + " " + "mol/mol"
    outds.sigtrace.attrs["units"] = obs_units + " " + "mol/mol"

    outds.Yobs.attrs["longname"] = "observations"
    outds.Yerror.attrs["longname"] = "measurement error"
    outds.Ytime.attrs["longname"] = "time of measurements"
    outds.Yapriori.attrs["longname"] = "a priori simulated measurements"
    outds.Ymodmean.attrs["longname"] = "mean of posterior simulated measurements"
    outds.Ymodmedian.attrs["longname"] = "median of posterior simulated measurements"
    outds.Ymodmode.attrs["longname"] = "mode of posterior simulated measurements"
    outds.Ymod68.attrs["longname"] = " 0.68 Bayesian credible interval of posterior simulated measurements"
    outds.Ymod95.attrs["longname"] = " 0.95 Bayesian credible interval of posterior simulated measurements"
    outds.xtrace.attrs["longname"] = "trace of unitless scaling factors for emissions parameters"
    outds.sigtrace.attrs["longname"] = "trace of model error parameters"
    outds.siteindicator.attrs["longname"] = "index of site of measurement corresponding to sitenames"
    outds.sigmafreqindex.attrs["longname"] = "perdiod over which the model error is estimated"
    outds.sitenames.attrs["longname"] = "site names"
    outds.sitelons.attrs["longname"] = "site longitudes corresponding to site names"
    outds.sitelats.attrs["longname"] = "site latitudes corresponding to site names"
    outds.fluxapriori.attrs["longname"] = "mean a priori flux over period"
    outds.fluxaposteriori_mode.attrs["longname"] = "mode posterior flux over period"
    outds.fluxaposteriori_mean.attrs["longname"] = "mean posterior flux over period"
    outds.scalingmean.attrs["longname"] = "mean scaling factor field over period"
    outds.scalingmode.attrs["longname"] = "mode scaling factor field over period"
    outds.basisfunctions.attrs["longname"] = "basis function field"
    outds.countrymean.attrs["longname"] = "mean of ocean and country totals"
    outds.countrymedian.attrs["longname"] = "median of ocean and country totals"
    outds.countrymode.attrs["longname"] = "mode of ocean and country totals"
    outds.country68.attrs["longname"] = "0.68 Bayesian credible interval of ocean and country totals"
    outds.country95.attrs["longname"] = "0.95 Bayesian credible interval of ocean and country totals"
    outds.countrysd.attrs["longname"] = "standard deviation of ocean and country totals"
    outds.countryapriori.attrs["longname"] = "prior mean of ocean and country totals"
    outds.countrydefinition.attrs["longname"] = "grid definition of countries"
    outds.xsensitivity.attrs["longname"] = "emissions sensitivity timeseries"

    if use_bc:
        outds.YmodmeanBC.attrs["units"] = obs_units + " " + "mol/mol"
        outds.YmodmedianBC.attrs["units"] = obs_units + " " + "mol/mol"
        outds.YmodmodeBC.attrs["units"] = obs_units + " " + "mol/mol"
        outds.Ymod95BC.attrs["units"] = obs_units + " " + "mol/mol"
        outds.Ymod68BC.attrs["units"] = obs_units + " " + "mol/mol"
        outds.YaprioriBC.attrs["units"] = obs_units + " " + "mol/mol"
        outds.bcsensitivity.attrs["units"] = obs_units + " " + "mol/mol"

        outds.YaprioriBC.attrs["longname"] = "a priori simulated boundary conditions"
        outds.YmodmeanBC.attrs["longname"] = "mean of posterior simulated boundary conditions"
        outds.YmodmedianBC.attrs["longname"] = "median of posterior simulated boundary conditions"
        outds.YmodmodeBC.attrs["longname"] = "mode of posterior simulated boundary conditions"
        outds.Ymod68BC.attrs[
            "longname"
        ] = " 0.68 Bayesian credible interval of posterior simulated boundary conditions"
        outds.Ymod95BC.attrs[
            "longname"
        ] = " 0.95 Bayesian credible interval of posterior simulated boundary conditions"
        outds.bctrace.attrs[
            "longname"
        ] = "trace of unitless scaling factors for boundary condition parameters"
        outds.bcsensitivity.attrs["longname"] = "boundary conditions sensitivity timeseries"

    if mcmc_dict["add_offset"] == True:
        outds.Yoffmean.attrs["units"] = obs_units + " " + "mol/mol"
        outds.Yoffmedian.attrs["units"] = obs_units + " " + "mol/mol"
        outds.Yoffmode.attrs["units"] = obs_units + " " + "mol/mol"
        outds.Yoff95.attrs["units"] = obs_units + " " + "mol/mol"
        outds.Yoff68.attrs["units"] = obs_units + " " + "mol/mol"

        outds.Yoffmean.attrs["longname"] = "mean of posterior simulated offset between measurements"
        outds.Yoffmedian.attrs["longname"] = "median of posterior simulated offset between measurements"
        outds.Yoffmode.attrs["longname"] = "mode of posterior simulated offset between measurements"
        outds.Yoff68.attrs["longname"] = " 0.68 Bayesian credible interval of posterior simulated offset between measurements"
        outds.Yoff95.attrs["longname"] = " 0.95 Bayesian credible interval of posterior simulated offset between measurements"


    outds.attrs["Start date"] = start_date
    outds.attrs["End date"] = end_date
    outds.attrs["Latent sampler"] = str(mcmc_results['step1'])[19:-25]
    outds.attrs["Hyper sampler"] = str(mcmc_results['step2'])[19:-25]
    outds.attrs["Burn in"] = str(int(mcmc_dict["burn"]))
    outds.attrs["Tuning steps"] = str(int(mcmc_dict["tune"]))
    outds.attrs["Number of chains"] = str(int(mcmc_dict["nchain"]))
    outds.attrs["Error for each site"] = str(mcmc_dict["sigma_per_site"])
    outds.attrs["Emissions Prior"] = str(mcmc_dict["xprior"])
    outds.attrs["Model Error Prior"] = str(mcmc_dict["sigprior"]) 
    
    # outds.attrs["Emissions Prior"] = "".join(["{0},{1},".format(k, v) for k, v in xprior.items()])[:-1]
    # outds.attrs["Model error Prior"] = "".join(["{0},{1},".format(k, v) for k, v in sigprior.items()])[:-1]
    if use_bc:
        # outds.attrs["BCs Prior"] = "".join(["{0},{1},".format(k, v) for k, v in bcprior.items()])[:-1]
        outds.attrs["Boundary Conditions prior"] = str(mcmc_dict["bcprior"])
    
    if mcmc_dict["add_offset"] == True:
        # outds.attrs["Offset Prior"] = "".join(["{0},{1},".format(k, v) for k, v in offsetprior.items()])[:-1]
        outds.attrs["Offset Prior"] = str(mcmc_dict["offsetprior"])
    
    outds.attrs["Creator"] = getpass.getuser()
    outds.attrs["Date created"] = str(pd.Timestamp("today"))

    outds.attrs["xPDF GelmanRubin"] = conv_xpdf_gr
    outds.attrs["bcPDF GelmanRubin"] = conv_bc_gr
    outds.attrs["sigPDF GelmanRubin"] = conv_sig_gr
    outds.attrs["xPDF result"] = conv_xpdf_result
    outds.attrs["bcPDF result"] = conv_bc_result
    outds.attrs["sigPDF result"] = conv_sig_result



    # outds.attrs["Repository version"] = code_version()

    # variables with variable length data types shouldn't be compressed
    # e.g. object ("O") or unicode ("U") type
    do_not_compress = []
    dtype_pat = re.compile(r"[<>=]?[UO]")  # regex for Unicode and Object dtypes
    for dv in outds.data_vars:
        if dtype_pat.match(outds[dv].data.dtype.str):
            do_not_compress.append(dv)

    # setting compression levels for data vars in outds
    comp = dict(zlib=True, complevel=5)
    encoding = {var: comp for var in outds.data_vars if var not in do_not_compress}

    # output_filename = define_output_filename(outputpath, species, domain, outputname, start_date, ext=".nc")

    output_filename = os.path.join(outputpath, f"{species}_{domain}_{outputname}_{start_date}.nc")
    
    Path(outputpath).mkdir(parents=True, exist_ok=True)
    outds.to_netcdf(output_filename, encoding=encoding, mode="w")
