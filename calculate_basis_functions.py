# ---------------------------------------------------------------------------------------
# calculate_basis_functions.py
# Created: 15 May 2024
# Author: Eric Saboya, School of Geographical Sciences, University of Bristol
# ---------------------------------------------------------------------------------------
# Functions for calculating basis functions for CO2 data sets
# ---------------------------------------------------------------------------------------

import os
import glob
import uuid
import getpass
import scipy.optimize
import numpy as np
import xarray as xr
import pandas as pd

from pathlib import Path
from typing import Optional
from collections import namedtuple
from importlib import reload


# *****************************************************************************
# BASIS FUNCTION ALGORITHMS 
# *****************************************************************************

class quadTreeNode:
    def __init__(self, xStart, xEnd, yStart, yEnd):
        self.xStart = xStart
        self.xEnd = xEnd
        self.yStart = yStart
        self.yEnd = yEnd

        self.child1 = None  # top left
        self.child2 = None  # top right
        self.child3 = None  # bottom left
        self.child4 = None  # bottom right

    def isLeaf(self):
        if self.child1 or self.child2 or self.child3 or self.child4:
            return False
        else:
            return True

    def createChildren(self, grid, limit):
        value = np.sum(grid[self.xStart : self.xEnd, self.yStart : self.yEnd])  # .values

        # stop subdividing if finest resolution or bucket level reached
        if value < limit or (self.xEnd - self.xStart < 2) or (self.yEnd - self.yStart < 2):
            return

        dx = self.xEnd - self.xStart
        dy = self.yEnd - self.yStart

        # create 4 children for subdivison
        self.child1 = quadTreeNode(self.xStart, self.xStart + dx // 2, self.yStart, self.yStart + dy // 2)
        self.child2 = quadTreeNode(
            self.xStart + dx // 2, self.xStart + dx, self.yStart, self.yStart + dy // 2
        )
        self.child3 = quadTreeNode(
            self.xStart, self.xStart + dx // 2, self.yStart + dy // 2, self.yStart + dy
        )
        self.child4 = quadTreeNode(
            self.xStart + dx // 2, self.xStart + dx, self.yStart + dy // 2, self.yStart + dy
        )

        # apply recursion on all child nodes
        self.child1.createChildren(grid, limit)
        self.child2.createChildren(grid, limit)
        self.child3.createChildren(grid, limit)
        self.child4.createChildren(grid, limit)

    def appendLeaves(self, leafList):
        # recursively append all leaves/end nodes to leafList
        if self.isLeaf():
            leafList.append(self)
        else:
            self.child1.appendLeaves(leafList)
            self.child2.appendLeaves(leafList)
            self.child3.appendLeaves(leafList)
            self.child4.appendLeaves(leafList)

def quadTreeGrid(grid, limit):
    """
    -------------------------------------------------------
    Apply quadtree division algorithm
    -------------------------------------------------------
    Args:
      grid (array):
        2d numpy array to apply quadtree division to
      limit (float):
        Use value as bucket level for defining maximum subdivision

    Returns:
      outputGrid (array):
        2d numpy grid, same shape as grid, with values correpsonding to
        each  box from boxList
      boxList: (list of lists)
        Each sublist describes the corners of a quadtree leaf
    -------------------------------------------------------
    """
    # start with a single node the size of the entire input grid:
    parentNode = quadTreeNode(0, grid.shape[0], 0, grid.shape[1])
    parentNode.createChildren(grid, limit)

    leafList = []
    boxList = []
    parentNode.appendLeaves(leafList)

    outputGrid = np.zeros_like(grid)

    for i, leaf in enumerate(leafList):
        outputGrid[leaf.xStart : leaf.xEnd, leaf.yStart : leaf.yEnd] = i
        boxList.append([leaf.xStart, leaf.xEnd, leaf.yStart, leaf.yEnd])

    return outputGrid, boxList

def quadtreebasisfunction(
    emissions_name,
    fp_all,
    sites,
    start_date,
    domain,
    species,
    outputname,
    outputdir=None,
    nbasis=50,
):
    """
    -------------------------------------------------------
    Creates a basis function with nbasis grid cells using a quadtree algorithm.
    The domain is split with smaller grid cells for regions which contribute
    more to the a priori (above basline) mole fraction. This is based on the
    average footprint over the inversion period and the a priori emissions field.
    Output is a netcdf file saved to /Temp/<domain> in the current directory
    if no outputdir is specified or to outputdir if specified.
    The number of basis functions is optimised using dual annealing. Probably
    not the best or fastest method as there should only be one minima, but doesn't
    require the Jacobian or Hessian for optimisation.
    -------------------------------------------------------
    Args:
      emissions_name (list):
        List of "source" key words as used for retrieving specific emissions
        from the object store.
      fp_all (dict):
        Output from footprints_data_merge() function. Dictionary of datasets.
      sites (list):
        List of site names (This could probably be found elsewhere)
      start_date (str):
        String of start date of inversion
      domain (str):
        The inversion domain
      species (str):
        Atmospheric trace gas species of interest (e.g. 'co2')
      outputname (str):
        Identifier or run name
      outputdir (str, optional):
        Path to output directory where the basis function file will be saved.
        Basis function will automatically be saved in outputdir/DOMAIN
        Default of None makes a temp directory.
      nbasis (int/list):
        Desired number of basis function regions.
        If int: same nbasis value used for all flux sectors 
        If list: specifiy nbasis value per flux sector 

    Returns:
        If outputdir is None, then returns a Temp directory. The new basis function is saved in this Temp directory.
        If outputdir is not None, then does not return anything but saves the basis function in outputdir.
    -------------------------------------------------------
    """
    if emissions_name == None:
        raise ValueError("emissions_name needs to be specified")
    
    # No. of flux sectors
    nsectors = len(emissions_name)
    
    # If one nbasis value provided, we assume that is 
    # the No. of basis functions for each flux sector 
    if isinstance(nbasis, int):
        nbasis = [nbasis] * nsectors
        
    basis_per_sector = {}
    
    # Calculate mean combined footprint of all sites being used 
    meanfp = np.zeros((fp_all[sites[0]].fp.shape[0], fp_all[sites[0]].fp.shape[1]))
    div = 0
    for site in sites:
        meanfp += np.sum(fp_all[site].fp.values, axis=2)
        div += fp_all[site].fp.shape[2]
    meanfp /= div

    print("Using absolute values of fluxes to calculate basis functions")
    for i in range(0, nsectors):
        flux_i = fp_all[".flux"][emissions_name[i]].data.flux.values
        absflux = np.absolute(flux_i)
        meanflux = np.squeeze(absflux)
        
        if meanflux.shape != meanfp.shape:
            meanflux = np.mean(meanflux, axis=2)
        
        fps = meanfp * meanflux
        print(f"Calculating basis functions for {emissions_name[i]} using {nbasis[i]} basis functions ...")

        def qtoptim(x):
            basisQuad, boxes = quadTreeGrid(fps, x)
            return (nbasis[i] - np.max(basisQuad) - 1) ** 2

        cost = 1e6
        pwr = 0
        while cost > 3.:
            optim = scipy.optimize.dual_annealing(qtoptim, np.expand_dims([0, 100 / 10**pwr], axis=0))
            cost = np.sqrt(optim.fun)
            pwr += 1
            if pwr > 10:
                raise Exception("Quadtree did not converge after max iterations.")
        basisQuad, boxes = quadTreeGrid(fps, optim.x[0])
        
        basis_per_sector[emissions_name[i]] = np.expand_dims(basisQuad + 1, axis=2) 


    lon = fp_all[sites[0]].lon.values
    lat = fp_all[sites[0]].lat.values     
    time = [pd.to_datetime(start_date)]

    base = []
    for key in basis_per_sector.keys():
        base.append(basis_per_sector[key])
    base = np.array(base)

    newds = xr.Dataset({"basis" : (["sector", "lat", "lon", "time"], base)}, 
                        coords={"time": (["time"], time), 
                                "lat": (["lat"], lat), 
                                "lon": (["lon"], lon),
                                "sector": (["sector"], emissions_name)
                       })  
 
    newds.lat.attrs["long_name"] = "latitude" 
    newds.lon.attrs["long_name"] = "longitude" 
    newds.lat.attrs["units"] = "degrees_north"
    newds.lon.attrs["units"] = "degrees_east"    
    newds.attrs["creator"] = getpass.getuser()
    newds.attrs["date created"] = str(pd.Timestamp.today())
    
    if outputdir is None:
        cwd = os.getcwd()
        tempdir = os.path.join(cwd, f"Temp_{str(uuid.uuid4())}")
        os.mkdir(tempdir)    
        os.mkdir(os.path.join(tempdir,f"{domain}/"))
        newds.to_netcdf(os.path.join(tempdir, domain, f"quadtree_{species}-{outputname}_{domain}_{start_date.split('-')[0]}.nc"), mode="w")
        return tempdir
    else:
        basisoutpath = os.path.join(outputdir, domain)
        if not os.path.exists(basisoutpath):
            os.makedirs(basisoutpath)
        newds.to_netcdf(os.path.join(basisoutpath, f"quadtree_{species}-{outputname}_{domain}_{start_date.split('-')[0]}.nc"), mode="w")

        return outputdir


# BUCKET BASIS FUNCTIONS
def load_landsea_indices():
    """
    -------------------------------------------------------
    Load UKMO array with indices that separate
    land and sea regions in EUROPE domain
    -------------------------------------------------------
    land = 1
    sea = 0
    -------------------------------------------------------
    """
    landsea_indices = xr.open_dataset("/user/work/wz22079/country_masks/country-EUROPE-UKMO-landsea-2023.nc")
    return landsea_indices["country"].values

def bucket_value_split(grid, bucket, offset_x=0, offset_y=0):
    """
    Algorithm that will split the input grid (e.g. fp * flux)
    such that the sum of each basis function region will
    equal the bucket value or by a single array element.

    The number of regions will be determined by the bucket value
    i.e. smaller bucket value ==> more regions
         larger bucket value ==> fewer regions
    ------------------------------------
    Args:
        grid: np.array
            2D grid of footprints * flux, or whatever
            grid you want to split. Could be: population
            data, spatial distribution of bakeries, you chose!

        bucket: float
            Maximum value for each basis function region

    Returns:
        array of tuples that define the indices for each basis function region
        [(ymin0, ymax0, xmin0, xmax0), ..., (yminN, ymaxN, xminN, xmaxN)]
    """

    if np.sum(grid) <= bucket or grid.shape == (1, 1):
        return [(offset_y, offset_y + grid.shape[0], offset_x, offset_x + grid.shape[1])]

    else:
        if grid.shape[0] >= grid.shape[1]:
            half_y = grid.shape[0] // 2
            return bucket_value_split(grid[0:half_y, :], bucket, offset_x, offset_y) + bucket_value_split(
                grid[half_y:, :], bucket, offset_x, offset_y + half_y
            )

        elif grid.shape[0] < grid.shape[1]:
            half_x = grid.shape[1] // 2
            return bucket_value_split(grid[:, 0:half_x], bucket, offset_x, offset_y) + bucket_value_split(
                grid[:, half_x:], bucket, offset_x + half_x, offset_y
            )

# Optimize bucket value to number of desired regions
def get_nregions(bucket, grid):
    """Returns no. of basis functions for bucket value"""
    return np.max(bucket_split_landsea_basis(grid, bucket))

def optimize_nregions(bucket, grid, nregion, tol, max_iter=2000):
    """
    Optimize bucket value to obtain nregion basis functions within +/- tol.
    """
    step_factor = 0.005  # Initial step factor
    iteration = 0

    while iteration < max_iter:
        current_nregions = get_nregions(bucket, grid)
        error = current_nregions - nregion

        if abs(error) <= tol:
            return bucket

        # Dynamically adjust step size based on the magnitude of the error
        dynamic_step = step_factor * (1 + abs(error) / max(nregion, 1))

        if error < 0:
            bucket *= (1 - dynamic_step)
        else:
            bucket *= (1 + dynamic_step)

        iteration += 1

    # Return the last computed bucket if convergence wasn't achieved
    return bucket

# def optimize_nregions(bucket, grid, nregion, tol):
#     """
#     Optimize bucket value to obtain nregion basis functions
#     within +/- tol.
#     """
#     # print(bucket, get_nregions(bucket, grid))
#     if get_nregions(bucket, grid) <= nregion + tol and get_nregions(bucket, grid) >= nregion - tol:
#         return bucket

#     if get_nregions(bucket, grid) < nregion + tol:
#         bucket = bucket * 0.995
#         return optimize_nregions(bucket, grid, nregion, tol)

#     elif get_nregions(bucket, grid) > nregion - tol:
#         bucket = bucket * 1.005
#         return optimize_nregions(bucket, grid, nregion, tol)

def bucket_split_landsea_basis(grid, bucket, offset_x=0, offset_y=0):
    """
    Same as bucket_split_basis but includes
    land-sea split. i.e. basis functions cannot overlap sea and land
     ------------------------------------
    Args:
        grid: np.array
            2D grid of footprints * flux, or whatever
            grid you want to split. This could be population
            data, distribution of bakeries. You choose!

        bucket: float
            Maximum value for each basis function region

    Returns:
        2D array with basis function values

    """
    landsea_indices = load_landsea_indices()
    myregions = bucket_value_split(grid, bucket)

    mybasis_function = np.zeros(shape=grid.shape)

    for i in range(len(myregions)):
        ymin, ymax = myregions[i][0]+offset_y, myregions[i][1]+offset_y
        xmin, xmax = myregions[i][2]+offset_x, myregions[i][3]+offset_x

        inds_y0, inds_x0 = np.where(landsea_indices[ymin:ymax, xmin:xmax] == 0)
        inds_y1, inds_x1 = np.where(landsea_indices[ymin:ymax, xmin:xmax] == 1)

        count = np.max(mybasis_function)

        if len(inds_y0) != 0:
            count += 1
            for i in range(len(inds_y0)):
                mybasis_function[inds_y0[i]+ymin-offset_y, inds_x0[i]+xmin-offset_x] = count

        if len(inds_y1) != 0:
            count += 1
            for i in range(len(inds_y1)):
                mybasis_function[inds_y1[i]+ymin-offset_y, inds_x1[i]+xmin-offset_x] = count

    return mybasis_function

def nregion_landsea_basis(grid, bucket=1, nregion=100, tol=1, offset_x=0, offset_y=0):
    """
    Obtain basis function with nregions (for land-sea split)
    ------------------------------------
    Args:
        grid: np.array
            2D grid of footprints * flux, or whatever
            grid you want to split. Could be: population
            data, spatial distribution of bakeries, you chose!

        bucket: float
            Initial bucket value for each basis function region.
            Defaults to 1

        nregion: int
            Number of desired basis function regions
            Defaults to 100

        tol: int
            Tolerance to find number of basis function regions.
            i.e. optimizes nregions to +/- tol
            Defaults to 1

    Returns:
        basis_function np.array
        2D basis function array

    """
    bucket_opt = optimize_nregions(bucket, grid, nregion, tol)
    basis_function = bucket_split_landsea_basis(grid, bucket_opt, offset_x, offset_y)
    return basis_function


def bucketbasisfunction(emissions_name: (str, list),
                        data_dict: dict,
                        sites: (str, list),
                        start_date: str,
                        domain: str,
                        species: str,
                        outputname: str,
                        outputdir: str,
                        nbasis: (int, list),
                       ):
    """
    Basis functions calculated using a weighted region approach
    where each basis function / scaling region contains approximately
    the same value
    -----------------------------------
    Args:
      emissions_name (str/list):
        List of keyword "source" args used for retrieving 
        emissions files from the object store.
      data_all (dict):
        data_dict dictionary object as produced from get_co2_data
      sites (str/list):
        List of measurements sites being used
      start_date (str):
        Start date of period of inference
      domain (str):
        Name of model domain
      species (str):
        Name of atmospheric trace gas species
      outputname (str):
        Name of inversion run
      outputdir (str):
        Directory where inversion run outputs are saved
      nbasis (int/list):
        Desired number of basis function regions.
        If int: same nbasis value used for all flux sectors 
        If list: specifiy nbasis value per flux sector 

    """
    if emissions_name is None:
        raise ValueError(" 'emissions_name' needs to be specified \n")
    
    # No. of flux sectors
    nsectors = len(emissions_name)
    
    # If one nbasis value provided, we assume that is 
    # the no. of basis functions for each flux sector 
    if isinstance(nbasis, int):
        nbasis = [nbasis] * nsectors
        
    basis_per_sector = {}
    
    # Calculate mean combined footprint of all sites being used 
    meanfp = np.zeros((data_dict[sites[0]].fp.shape[0], data_dict[sites[0]].fp.shape[1]))
    div = 0
    for site in sites:
        # meanfp += np.sum(data_dict[site].fp.values, axis=2)
        # div += data_dict[site].fp.shape[2]
        meanfp += data_dict[site].fp.mean(dim="time").values
    meanfp /= len(sites)

    # Calculate basis function per flux sector 
    print("Using absolute values of fluxes to calculate basis functions")
    for i in range(0, nsectors):
        flux_i = data_dict[".flux"][emissions_name[i]].data.flux
        absflux = np.absolute(flux_i)
                
        if absflux.shape != meanfp.shape:
            meanflux = xr.DataArray.mean(absflux, dim="time")
            print(meanflux.shape)
            
            if meanflux.shape == meanfp.shape:
                fps = meanfp * meanflux.values
            else:
                raise ValueError("Footprint and Flux dimensions do not match.")

        # Check whether a significant proportion of the fps domain has zero values
        # If so, apply basis functions to region where fluxes exist
        fps_nonzero_inds = np.where(fps!=0)
        prop = len(fps_nonzero_inds[0])/len(fps.ravel())
        print(f"Proportion of non-zero flux*fp grid cells in domain {prop}")
        
        if prop > 0.55:
            print(f"Calculating basis functions for {emissions_name[i]} using {nbasis[i]} basis functions ...")
            print("Calculating basis functions over entire model domain")
            
            # Use median grid value as starting point for buckets value 
            starting_bucket_value = np.nanmedian(fps)
            bucket_basis_i = nregion_landsea_basis(fps, starting_bucket_value, nbasis[i])
            basis_per_sector[emissions_name[i]] = np.expand_dims(bucket_basis_i, axis=2)
            
        else:
            print(f"Calculating basis functions for {emissions_name[i]} using {nbasis[i]} basis functions ...")
            print("Calculating basis functions over an inner domain only.")
            
            # Find sub-domain where fluxes exist
            i_min, i_max = np.nanmin(fps_nonzero_inds[0]), np.nanmax(fps_nonzero_inds[0])
            j_min, j_max = np.nanmin(fps_nonzero_inds[1]), np.nanmax(fps_nonzero_inds[1])

            n, m = fps.shape[0], fps.shape[1]
            
            # Inner region where values exist
            fps_inner = fps[i_min:i_max+1, j_min:j_max+1]
            starting_bucket_value = np.nanmedian(fps_inner)
            
            # Use median grid value as starting point for buckets value 
            bucket_basis_i = nregion_landsea_basis(fps_inner, 
                                                   starting_bucket_value, 
                                                   nbasis[i], 
                                                   1, 
                                                   j_min, 
                                                   i_min,
                                                  )

            new_basis_grid = np.zeros(fps.shape)
            
            # region 1
            for k in range(0, i_min):
                for j in range(0, j_min):
                    new_basis_grid[k,j] = 1 + bucket_basis_i.max()
            # region 2
            for k in range(i_min, i_max):
                for j in range(0, j_min):
                    new_basis_grid[k,j] = 2 + bucket_basis_i.max()
            # region 3
            for k in range(i_max, n):
                for j in range(0, j_min):
                    new_basis_grid[k,j] = 3 + bucket_basis_i.max()
            # region 4
            for k in range(0, i_min):
                for j in range(j_min, j_max):
                    new_basis_grid[k,j] = 4 + bucket_basis_i.max()
            # region 5
            for k in range(i_max, n):
                for j in range(j_min, j_max):
                    new_basis_grid[k,j] = 5 + bucket_basis_i.max()
            # region 6
            for k in range(0, i_min):
                for j in range(j_max, m):
                    new_basis_grid[k,j] = 6 + bucket_basis_i.max()
            # region 7
            for k in range(i_min, i_max):
                for j in range(j_max, m):
                    new_basis_grid[k,j] = 7 + bucket_basis_i.max()
            # region 8
            for k in range(i_max, n):
                for j in range(j_max, m):
                    new_basis_grid[k,j] = 8 + bucket_basis_i.max()
            # Inner region
            for k in range(i_min, i_max):
                for j in range(j_min, j_max):
                    new_basis_grid[k,j] = bucket_basis_i[k-i_min, j-j_min]

            # Add sector basis function to dict 
            basis_per_sector[emissions_name[i]] = np.expand_dims(new_basis_grid, axis=2)
        
    lon = data_dict[sites[0]].lon.values
    lat = data_dict[sites[0]].lat.values
    time = [pd.to_datetime(start_date)]
    
    base = []
    for key in basis_per_sector.keys():
        base.append(basis_per_sector[key])
    base = np.array(base)

    # Create xarray dataset with basis function per sector
    newds = xr.Dataset({"basis": (["sector", "lat", "lon", "time"], base)},
                        coords = {"time": (["time"], time),
                                  "lat": (["lat"], lat),
                                  "lon": (["lon"], lon),
                                  "sector": (["sector"], emissions_name),
                                 })

    newds.lat.attrs["long_name"] = "latitude"
    newds.lon.attrs["long_name"] = "longitude"
    newds.lat.attrs["units"] = "degrees_north"
    newds.lon.attrs["units"] = "degrees_east"
    newds.attrs["creator"] = getpass.getuser()
    newds.attrs["date created"] = str(pd.Timestamp.today())

    if outputdir is None:
        # cwd = os.getcwd()
        cwd = "/user/home/wz22079/my_openghg/openghg_inversions/scratch/"
        tempdir = os.path.join(cwd, f"Temp_{str(uuid.uuid4())}")
        os.mkdir(tempdir)
        os.mkdir(os.path.join(tempdir, f"{domain}/"))
        newds.to_netcdf(os.path.join(tempdir, domain, f"weighted_co2-{outputname}_{domain}_{start_date.split('-')[0]}{start_date.split('-')[1]}.nc"), mode="w")
        return tempdir

    else:
        basisoutpath = os.path.join(outputdir, domain)
        if not os.path.exists(basisoutpath):
            os.makedirs(basisoutpath)
        newds.to_netcdf(os.path.join(basisoutpath, f"weighted_co2-{outputname}_{domain}_{start_date.split('-')[0]}{start_date.split('-')[1]}.nc"), mode='w')
        return outputdir


# *****************************************************************************
# BASIS FUNCTION UTILITIES 
# *****************************************************************************

def basis(domain: str, 
          basis_case: str, 
          basis_directory: Optional[str] = None,
         ):
    """
    The basis function reads in the all matching files for the
    basis case and domain as an xarray Dataset.

    Expect filenames of the form:
        [basis_directory]/domain/"basis_case"_"domain"*.nc
        e.g. [/data/shared/LPDM/basis_functions]/EUROPE/sub_transd_EUROPE_2014.nc

    TODO: More info on options for basis functions.
    -----------------------------------
    Args:
      domain (str):
        Domain name. The basis files should be sub-categorised by the domain.
      basis_case (str):
        Basis case to read in.
        Examples of basis cases are "voroni","sub-transd","sub-country_mask",
        "INTEM".
      basis_directory (str, optional):
        basis_directory can be specified if files are not in the default
        directory. Must point to a directory which contains subfolders
        organized by domain.

    Returns:
      xarray.Dataset:
        combined dataset of matching basis functions
    -----------------------------------
    """
    from utils import read_netcdfs
    
    openghginv_path = "/user/home/wz22079/my_openghg/openghg_inversions/scratch/"
    
    if basis_directory is None:
        if not os.path.exists(os.path.join(openghginv_path, "basis_functions/")):
            os.makedirs(os.path.join(openghginv_path, "basis_functions/"))
        basis_directory = os.path.join(openghginv_path, "basis_functions/")

    file_path = os.path.join(basis_directory, domain, f"{basis_case}*.nc")
    files = sorted(glob.glob(file_path))

    if len(files) == 0:
        raise IOError(
            f"\nError: Can't find basis function files for domain '{domain}' \
                          and basis_case '{basis_case}' "
        )

    basis_ds = read_netcdfs(files)

    return basis_ds


def basis_boundary_conditions(domain: str, 
                              basis_case: str, 
                              bc_basis_directory: Optional[str] = None,
                             ):
    """
    The basis_boundary_conditions function reads in all matching files
    for the boundary conditions basis case and domain as an xarray Dataset.

    Expect filesnames of the form:
        [bc_basis_directory]/domain/"basis_case"_"domain"*.nc
        e.g. [/data/shared/LPDM/bc_basis_directory]/EUROPE/NESW_EUROPE_2013.nc

    TODO: More info on options for basis functions.
    -----------------------------------
    Args:
      domain (str):
        Domain name. The basis files should be sub-categorised by the domain.
      basis_case (str):
        Basis case to read in. Examples of basis cases are "NESW","stratgrad".
      bc_basis_directory (str, optional):
        bc_basis_directory can be specified if files are not in the default directory.
        Must point to a directory which contains subfolders organized by domain.

    Returns:
      xarray.Datset:
        Combined dataset of matching basis functions
    -----------------------------------
    """
    from utils import read_netcdfs
    
    openghginv_path = "/user/home/wz22079/my_openghg/openghg_inversions/scratch/"

    if bc_basis_directory is None:
        if not os.path.exists(os.path.join(openghginv_path, "bc_basis_functions/")):
            os.makedirs(os.path.join(openghginv_path, "bc_basis_functions/"))
        bc_basis_directory = os.path.join(openghginv_path, "bc_basis_functions/")

    file_path = os.path.join(bc_basis_directory, domain, f"{basis_case}_{domain}*.nc")

    files = sorted(glob.glob(file_path))
    file_no_acc = [ff for ff in files if not os.access(ff, os.R_OK)]
    files = [ff for ff in files if os.access(ff, os.R_OK)]

    if len(file_no_acc) > 0:
        print(
            "Warning: unable to read all boundary conditions basis function files which match this criteria:"
        )
        [print(ff) for ff in file_no_acc]

    if len(files) == 0:
        raise IOError(
            "\nError: Can't find boundary condition basis function files for domain '{0}' "
            "and basis_case '{1}' ".format(domain, basis_case)
        )

    basis_ds = read_netcdfs(files)

    return basis_ds


BasisFunction = namedtuple("BasisFunction", ["description", "algorithm"])
basis_functions = {"quadtree": BasisFunction("quadtree algorithm", quadtreebasisfunction),
                   "weighted": BasisFunction("weighted by data algorithm", bucketbasisfunction),
                  }



def basis_functions_wrapper(data_dict: dict, 
                            basis_dict: dict, 
                            use_bc: bool,
                            outputname: Optional[str] = None,
                            outputpath: Optional[str] = None,
                           ):
    """
    Wrapper function for selecting basis function
    algorithm.
    
    Args:
        data_dict (dict): 
            Dictionary of observations and forward simulations
            datasets created from get_co2_data.py
        basis_dict (dict): 
            Dictionary of basis functions specifications
        use_bc (bool):
            True: Retrieves BCs as specified in bc_dict.
            False: Asssumes Obs are assumed to be baseline-subtracted 
                   And BC data are NOT used
        outputname (str):
            Name of output
        outputpath (str):
            Path to save output                  
    """
    import utils
    from sensitivity import fp_sensitivity, bc_sensitivity
        
    # Extract inputs from basis_dict 
    fp_basis_case=basis_dict["fp_basis_case"]
    basis_directory=basis_dict["basis_directory"]
    basis_algorithm=basis_dict["fp_basis_algorithm"]
    nbasis=basis_dict["nbasis"]
    bc_basis_case=basis_dict["bc_basis_case"]
    bc_basis_directory=basis_dict["bc_basis_directory"]
    domain=basis_dict["domain"]
    source=basis_dict["source"]
    site=basis_dict["site"]
    start_date=basis_dict["start_date"]
    
    if fp_basis_case is not None:
        if basis_algorithm:
            print(f"Basis algorithm {basis_algorithm} and basis case {fp_basis_case} supplied; using {fp_basis_case}.")
        
        basis_data_array = basis(domain=domain, 
                                 basis_case=fp_basis_case, 
                                 basis_directory=basis_directory,
                                ).basis

        tempdir = None
        
    elif basis_algorithm is None:
        raise ValueError("One of `fp_basis_case` or `basis_algorithm` must be specified.")

    else:
        try:
            basis_function = basis_functions[basis_algorithm]
            
        except KeyError as e:
            raise ValueError(
                "Basis algorithm not recognised. Please use either 'quadtree' or 'weighted', or input a basis function file"
            ) from e
        print(f"Using {basis_function.description} to derive basis functions.")

        
        tempdir = basis_function.algorithm(source,
                                           data_dict, 
                                           site,
                                           start_date,
                                           domain,
                                           "CO2",
                                           outputname,
                                           outputpath,
                                           nbasis,
                                          )

        fp_basis_case = "weighted_co2-" + outputname
        basis_directory = tempdir
                                                                                               
    fp_data = fp_sensitivity(data_dict, 
                             domain,
                             fp_basis_case,
                             basis_directory,
                            )

    if use_bc is True:
        fp_data = bc_sensitivity(data_dict,
                                 domain=domain,
                                 basis_case=bc_basis_case,
                                 bc_basis_directory=bc_basis_directory,
                                )

    return fp_data, tempdir, basis_directory, bc_basis_directory




