# ---------------------------------------------------------------------------------------
# postprocessing_output.py
# Created 3 Dec 2024 
# Author: Eric Saboya, School of Geographical Sciences, University of Bristol 
# ---------------------------------------------------------------------------------------
# Postprocessing individual output files from RHIME-CO2
# Functions to plot:
# > mf timeseries of posterior, prior and observations for each site
# > emissions values of prior and posterior for each sector (time series and flux map)
# ---------------------------------------------------------------------------------------

import os 
import pickle
import numpy as np
import pandas as pd
import xarray as xr
import datetime as dt
 
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
import matplotlib.dates as mdates
import matplotlib.text as text

from matplotlib import ticker
from matplotlib import gridspec
from matplotlib import rcParams
from matplotlib.lines import Line2D

from matplotlib.colors import ListedColormap, LinearSegmentedColormap

def update_rcParams(key, val):
    if key in rcParams:
        rcParams[key] = val

update_rcParams('font.size', 11)
update_rcParams('font.family', 'serif')
update_rcParams('xtick.major.size', 8)
update_rcParams('xtick.labelsize', 'large')
update_rcParams('xtick.direction', "in")
update_rcParams('xtick.minor.visible', False) #
update_rcParams('xtick.top', False) #
update_rcParams('ytick.major.size', 8)
update_rcParams('ytick.labelsize', 'large')
update_rcParams('ytick.direction', "in")
update_rcParams('ytick.minor.visible', False) #
update_rcParams('ytick.right', False) # 
update_rcParams('xtick.minor.size', 4)
update_rcParams('ytick.minor.size', 4)
update_rcParams('xtick.major.pad', 10)
update_rcParams('ytick.major.pad', 10)
update_rcParams('legend.numpoints', 1)
#update_rcParams('use.tex',True)

import matplotlib as mpl 
mpl.rcParams['axes.linewidth'] = 2

# --------------------------------------------------------------------------------------- # 

def processing_mf_outputs(ds)->dict:
    """
    -------------------------------------------------------
    Process mole fraction outputs into a dictionary 
    with site names as keys
    -------------------------------------------------------
    Args:
        ds (DataArray):
            RHIME output netcdf file
            
    Returns:
        Dictionary containing a priori, a posteriori, observed mfs
        and their uncertainties.
    -------------------------------------------------------
    """
    # Mole fraction comparisons
    sitenames = ds["sitenames"].values
    site_dict = {}

    for i, sname in enumerate(sitenames):
        mydict = {}
        site_ind = np.where(ds["siteindicator"].values==i)

        mydict["t"] = ds["Ytime"][site_ind].values
        mydict["yobs"] = ds["Yobs"][site_ind].values
        mydict["yobs_err"] = ds["Yerror"][site_ind].values
        mydict["y_apriori"] = ds["Yapriori"][site_ind].values
        mydict["y_apriori_bc"] = ds["YaprioriBC"][site_ind].values
    
        mydict["y_apost"] = ds["Ymodmean"][site_ind].values
        mydict["y_apost_68"] = ds["Ymod68"][site_ind].values    
    
        mydict["y_apost_bc"] = ds["YmodmeanBC"][site_ind].values    
        mydict["y_apost_bc_68"] = ds["Ymod68BC"][site_ind].values    

        site_dict[sname] = mydict

    return site_dict


def plot_mf_timeseries(site_dict: dict,
                       period: list,
                       show_fig: bool,
                       save_fig: bool,
                       save_dir: str,
                       model_run_name: str,
                      )->None:
    """
    Plot mole fractions time series for 
    inversion period
    -------------------------------------------------------
    Args:
        site_dict (dict):
            Dictionary output from `processing_mf_outputs`
        period (list):
            List containing start_date and end_date
        show_fig (bool):
            Show the plot
        save_fig (bool)
            Option to save the figure
        save_dir (str):
            Directory where figures are saved
        model_run_name (str):
            Name that indicates which model setup was run
    -------------------------------------------------------
    """

    start_date = dt.datetime.strptime(period[0], "%Y-%m-%d")
    end_date = dt.datetime.strptime(period[1], "%Y-%m-%d")
        
    
    for site in site_dict.keys():

        fsave_name = f"RHIME_CO2_mf-timeseries_{model_run_name}_{start_date}_{site}.jpg"
        
        fig, ax = plt.subplots(figsize=(12,5))
        
        # Plot CO2 mf observations 
        ax.errorbar(site_dict[site]["t"], 
                    site_dict[site]["yobs"], 
                    yerr=site_dict[site]["yobs_err"], 
                    fmt="o-", 
                    color="k", 
                    ecolor="k", 
                    label="Observed", 
                    alpha=0.7,)
        
        # Plot CO2 a priori simulated
        ax.plot(site_dict[site]["t"], 
                site_dict[site]["y_apriori"], 
                "^:", 
                color="b",
                label="A priori simulated")

        # Plot CO2 a posteriori values 
        ax.plot(site_dict[site]["t"], 
                site_dict[site]["y_apost"], 
                "s-.", 
                color="r",
                label="A posteriori")
        
        ax.fill_between(site_dict[site]["t"], 
                        site_dict[site]["y_apost_68"][:,0], 
                        site_dict[site]["y_apost_68"][:,1], 
                        color="r", 
                        alpha=0.35)
    
        ax.set_title(site)
        ax.legend(loc=2)
        ax.set_ylabel(r"CO$_2$ mole fraction (mol mol$^{-1}$)")

        ax.set_xlim((start_date, end_date))
        fig.tight_layout()
        if save_fig:
            plt.savefig(os.path.join(save_dir, fsave_name), dpi=300)
        
        else:
            plt.show()

        if show_fig is False:
            plt.clf()
            plt.close()

def plot_mf_wobg_timeseries(site_dict: dict,
                            period: list,
                            save_dir: str,
                            model_run_name: str,
                           ):
    """
    Plot mole fractions (background subtracted)
    for each station
    Args:
        site_dict (dict):
            Dictionary output from `processing_mf_outputs`
        period (list):
            List containing datetime values for 
            start_date and end_date
        save_dir (str):
            Directory where figures are saved
        model_run_name (str):
            Name that indicates which model setup was run

    """
    start_date, end_date = period[0], period[1]

    for site in site_dict.keys():
        fsave_name = f"RHIME_CO2_mf-wobg-timeseries_{model_run_name}_{start_date}.jpg"
        
        fig, ax = plt.subplots(figsize=(12,5))

        # Plot observations - a priori baselines 
        ax.errorbar(site_dict[site]["t"], 
                    site_dict[site]["yobs"]-site_dict[site]["y_apriori_bc"], 
                    yerr=site_dict[site]["yobs_err"], 
                    fmt="^:", 
                    color="k", 
                    ecolor="k", 
                    label="Observed", 
                    alpha=0.7)
        
        # Plot a priori baseline deviations
        ax.plot(site_dict[site]["t"], 
                site_dict[site]["y_apriori"]-site_dict[site]["y_apriori_bc"], 
                "s--", 
                color="b",
                label="apriori")

        # Plot a posteriori baseline deviations
        ax.plot(site_dict[site]["t"], 
                site_dict[site]["y_apost"]-site_dict[site]["y_apost_bc"], 
                "o-", 
                color="r",
                label="aposteriori")
    
        ax.set_title(site)
        ax.legend(loc=2)
        ax.set_ylabel(r"Baseline-perturbed CO$_2$ mole fraction (ppm)")

        ax.set_xlim((start_date, end_date))
        fig.tight_layout()
        plt.savefig(os.path.join(save_dir, fsave_name), dpi=300)
        plt.clf()
        plt.close()

def plot_mf_bg_timeseries(site_dict: dict,
                          period: list,
                          save_dir: str,
                          model_run_name: str,
                         ):
    """
    Plot modelled background mole fractions
    for each station
    Args:
        site_dict (dict):
            Dictionary output from `processing_mf_outputs`
        period (list):
            List containing datetime values for 
            start_date and end_date
        save_dir (str):
            Directory where figures are saved
        model_run_name (str):
            Name that indicates which model setup was run 

    """
    start_date, end_date = period[0], period[1]

    for site in site_dict.keys():
        fsave_name = f"RHIME_CO2_mf-bg-timeseries_{model_run_name}_{start_date}.jpg"
        
        fig, ax = plt.subplots(figsize=(12,5))
        
        # Plot a priori baseline 
        ax.plot(site_dict[site]["t"], 
                site_dict[site]["y_apriori_bc"], 
                "s--", 
                color="b",
                label="apriori")

        # Plot a posteriori baseline deviations
        ax.plot(site_dict[site]["t"], 
                site_dict[site]["y_apost"]-site_dict[site]["y_apost_bc"], 
                "o-", 
                color="r",
                label="aposteriori")
    
        ax.set_title(site)
        ax.legend(loc=2)
        ax.set_ylabel(r"Baseline CO$_2$ mole fraction (ppm)")

        ax.set_xlim((start_date, end_date))
        fig.tight_layout()
        plt.savefig(os.path.join(save_dir, fsave_name), dpi=300)
        plt.clf()
        plt.close()


def plot_uk_flux_maps(ds
                      save_dir,
                      model_run_name,
                     ):
    """
    Plot fluxes of UK totals for each sector
    Args:
        ds (DataArray):
            RHIME output file
        save_dir (str):
            Directory where figures are saved
        model_run_name (str):
            Name that indicates which model setup was run
        
    """

    start_date = ds.attrs['Start date']
    fsave_name = f"RHIME_CO2_uk-fluxes_{model_run_name}_{start_date}.jpg"
    
    import cartopy
    import cartopy.crs as ccrs
    import matplotlib.ticker as mticker
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

    # economist_cols = ["#EBEDFA","#D6DBF5","#475ED1","#2E45B8","#1F2E7A","#141F52"]
    economist_cols = ['#2E45B8','#475ED1','#EBEDFA','#F6423C','#E3120B']
    economist_cmap = LinearSegmentedColormap.from_list("economist_cmap", economist_cols)

    zissou_cols = ["#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00"]
    zissou_cmap = LinearSegmentedColormap.from_list("zissou_cmap", zissou_cols)

    lat, lon = ds['lat'].values, ds['lon'].values
    xs, ys = np.meshgrid(lon, lat)
    diff = ds['fluxaposteriori_mean'].values-ds['fluxapriori'].values
    nsectors = len(ds['fluxsector'])

    # Make the plot
    fig, ax = plt.subplots(nsectors, 3, 
                           figsize=(12, 10),
                           subplot_kw={'projection': ccrs.PlateCarree()})

    
    for i, fname in enumerate(ds['fluxsector'].values):
        
        # Column 1: a priori fluxes
        k0 = ax[i,0].pcolormesh(xs, ys, ds['fluxapriori'][i,:,:], cmap=zissou_cmap)
        plt.colorbar(k0, ax=ax[i,0], orientation='vertical', label="mol/m2/s")
        ax[i,0].text(-0.07, 0.55, f'{fname}', 
                     va='bottom', 
                     ha='center', 
                     rotation='vertical', 
                     rotation_mode='anchor', 
                     transform=ax[i,0].transAxes)
        
        # Column 2: a posteriori fluxes
        k1 = ax[i,1].pcolormesh(xs, ys, ds['fluxaposteriori_mean'][i,:,:], cmap=zissou_cmap)
        plt.colorbar(k1, ax=ax[i,1], orientation='vertical', label="mol/m2/s")

        # Column 3: a posteriori - a priori fluxes 
        k2 = ax[i,2].pcolormesh(xs, ys, diff[i,:,:], cmap=economist_cmap)
        plt.colorbar(k2, ax=ax[i,2], orientation='vertical', label="mol/m2/s")


    ax[0,0].set_title("a priori")
    ax[0,1].set_title("a posteriori")
    ax[0,2].set_title("a posteriori - a priori")

    for i in range(nsectors):
        for j in range(3):
            ax[i,j].coastlines(lw=1.0)
            ax[i,j].set_xlim((-10.5, 2))
            ax[i,j].set_ylim((49.5, 61))
    
    fig.tight_layout()
    plt.savefig(os.path.join(save_dir, fsave_name), dpi=300)
    plt.clf()
    plt.close()


def plot_uk_flux_scalings_map(ds):
    """
    """
    start_date = ds.attrs['Start date']
    fsave_name = f"RHIME_CO2_uk-fluxscalings_{model_run_name}_{start_date}.jpg"
    
    import cartopy
    import cartopy.crs as ccrs
    import matplotlib.ticker as mticker
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

    # economist_cols = ["#EBEDFA","#D6DBF5","#475ED1","#2E45B8","#1F2E7A","#141F52"]
    economist_cols = ['#2E45B8','#475ED1','#EBEDFA','#F6423C','#E3120B']
    economist_cmap = LinearSegmentedColormap.from_list("economist_cmap", economist_cols)

    zissou_cols = ["#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00"]
    zissou_cmap = LinearSegmentedColormap.from_list("zissou_cmap", zissou_cols)

    lat, lon = ds['lat'].values, ds['lon'].values
    xs, ys = np.meshgrid(lon, lat)
    nsectors = len(ds['fluxsector'])

    fname = ds['fluxsector'].values

    # Should include a fourth subplot with total scaling 
    if nsectors+1 <= 3:
        nrows=1
        ncols=3
    elif nsectors+1 == 4:
        nrows=2
        ncols=2
    elif nsectors+1 == 5:
        nrows=2
        ncols=3
    elif nsectors+1 ==6:
        nrows=2
        ncols=3
    else:
        nrows=3
        ncols=3
    
    # Make the plot
    fig, ax = plt.subplots(nrows, ncols, 
                           figsize=(13, 4),
                           subplot_kw={'projection': ccrs.PlateCarree()})

    count = 0
    for i in range(nrows):
        for j in range(ncols):
            
            k0 = ax[i,j].pcolormesh(xs, ys, ds['scalingmean'][count,:,:], cmap=economist_cmap, vmin=0.75, vmax=1.25)
            plt.colorbar(k0, ax=ax[i,j], orientation='vertical', label='Mean Scaling')        
            ax[i,j].set_title(f"{fname[count]}")
            count += 1
    
    for i in range(nrows):
        for j in range(ncols):
            ax[i,j].coastlines(lw=1.0)
            ax[i,j].set_xlim((-10.5, 2))
            ax[i,j].set_ylim((49.5, 61))

fig.tight_layout()
plt.show()



def plot_uk_emissions_timeseries():
    """
    """
    fig, ax = plt.subplots(figsize=(12,4))
    ax.errorbar(t_s1, fossil_mean[0], yerr=fossil_err[0], fmt='o', color='k', ecolor='k', ms=7, label="Fossil a posteriori", capsize=5)
    ax.errorbar(t_s1+dt.timedelta(hours=10), fossil_pri[0], fmt='D', color='k', ecolor='k', ms=7, label="Fossil a priori")

    ax.errorbar(t_s2, gee_mean[0], yerr=gee_err[0], fmt='o', color='b', ecolor='b', ms=5, label="GEE a posteriori", capsize=5)
    ax.errorbar(t_s2+dt.timedelta(hours=10), gee_pri[0], fmt='D', color='b', ecolor='b', ms=7, label="GEE a priori")


    ax.errorbar(t_s3, resp_mean[0], yerr=resp_err[0], fmt='o', color='g', ecolor='g', ms=7, label="Respiration a posteriori", capsize=5)
    ax.errorbar(t_s3+dt.timedelta(hours=10), resp_pri[0], fmt='D', color='g', ecolor='g', ms=7, label="Respiration a priori")

    ax.set_xlim((dt.datetime(2014,1,1), dt.datetime(2014,2,1)))
    ax.legend(loc=3)
    ax.set_ylabel(r"UK CO$_2$ emissions (Tg yr$^{-1}$)")
    fig.tight_layout()
    plt.show()





def postprocessing(fname: str,
                   model_run_name: str, 
                   save_dir: str,
                   mf_timeseries: bool,
                   mf_wobg_timeseries: bool,
                   mf_bg_timeseries: bool,
                   uk_flux_maps: bool,
                   uk_flux_scalings_map: bool,
                   
                  ):
    """
    """
    # Read output file and processing dictionary 
    ds = xr.open_dataset(fname)
    site_dict = processing_mf_outputs(ds)

    start_date = ds.attrs['Start date']
    end_date = ds.attrs['End date']
    
    if mf_timeseries is True:
        plot_mf_timeseries(site_dict, period, save_dir, model_run_name)

    if mf_wobg_timeseries is True:
        