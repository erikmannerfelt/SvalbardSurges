import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import socket
import os
import warnings

# The PROJ installation points to the wrong directory for the proj.db file, which needs to be fixed on this computer
if socket.gethostname() == "DESKTOP-09DFBN6":
    os.environ["PROJ_DATA"] =  "C:\\Users\\eliss\\anaconda3\\envs\\SvalbardSurges\\Lib\\site-packages\\pyproj\\proj_dir\\share\\proj"

import xarray as xr
import variete
import variete.vrt.vrt
import rasterio as rio
import geoutils as gu
import geopandas as gpd


# Catch a deprecation warning that arises from skgstat when importing xdem
with warnings.catch_warnings():
    import numba
    warnings.simplefilter("ignore", numba.NumbaDeprecationWarning)
    import xdem


def main():

    # load IS2 data
    data = xr.open_dataset("nordenskiold_land-is2.nc")

    # set bounds for subsetting data
    bounds = {
        "bottom": 8617861,
        "right": 552643,
        "top": 8637186,
        "left": 540154
    }

    # subset IS2 data by bounds
    is2_subset = subset_is2(data, bounds, "scheelebreen")

    # plot the cropped IS2 data
    #plot_pts(is2_subset, "h_te_best_fit")

    # load DEM cropped to bounds
    dem_subset = load_dem(bounds, "scheelebreen")

    # load shapefile
    scheelebreen_shp = load_shp("Scheelebreen")

    # clip DEM to glacier area outlines
    cropped_dem = mask_dem(dem_subset, scheelebreen_shp)

    # get elevation difference between IS2 and reference DEM
    is2_dh = IS2_DEM_difference(cropped_dem, is2_subset, "scheelebreen")

    # plot elevation differences !!! functions need improvement (choose which variable to plot)
    #is2_dh = xr.open_dataset("cache/scheelebreen-is2-dh.nc")
    #plot_pts(is2_dh, "dh")

    hypsometric_binning("cache/scheelebreen-is2-dh.nc", "cache/S0_DTM5_2011_25163_33_scheelebreen_cropped.vrt")

def subset_is2(data, bounds, label):
    # subset IS2 data by bounds, create label for them

    cache_path = Path(f"cache/{label}-is2.nc")

    # if subset already exists open dataset
    if cache_path.is_file():
        return xr.open_dataset(cache_path)

    # subset data based on bounds
    subset = data.where(
        (data.easting > bounds["left"]) & (data.easting < bounds["right"]) & (data.northing > bounds["bottom"]) & (
                    data.northing < bounds["top"]), drop=True)

    #
    cache_path.parent.mkdir(exist_ok = True)
    subset.to_netcdf(cache_path)

    return subset

def plot_pts(data, var):
    # creates scatter plot of data

    plt.scatter(data.easting, data.northing, c=data[var])
    plt.gca().set_aspect('equal')
    plt.show()

def load_shp(glacier_name):
    # loads shp based on glacier name

    # paths
    file_name = Path("GAO_SfM_1936_1938_v3.shp")
    dir_name = Path("C:/Users/eliss/SvalbardSurges/GAO")
    file_path = dir_name/file_name

    # load shapefile converted to EPSG:32633
    gao = gpd.read_file(file_path).to_crs(32633)

    # filter by glacier name
    gao_glacier = gao.query(f"NAME=='{glacier_name}'")

    # if there's more glaciers with the same name select the largest one (??)

    return gao_glacier

def load_dem(bounds, label):
    # loads DEM

    # paths
    file_name = Path("S0_DTM5_2011_25163_33.tif")
    dir_name = Path("C:/Users/eliss/SvalbardSurges/DEMs")
    file_path = dir_name/("NP_"+file_name.stem)/file_name
    vrt_warped_filepath = Path(f"cache/{file_name.stem}_{label}_warped.vrt")
    vrt_cropped_filepath = Path(f"cache/{file_name.stem}_{label}_cropped.vrt")

    # convert bounds (dict) to bounding box (list)
    bbox = rio.coords.BoundingBox(**bounds)

    # warp vrt (virtual raster), dst coord system EPSG:32633 (WGS-84)
    variete.vrt.vrt.vrt_warp(vrt_warped_filepath, file_path, dst_crs=32633)

    # crop warped vrt to bbox
    variete.vrt.vrt.build_vrt(vrt_cropped_filepath, vrt_warped_filepath, output_bounds=bbox)

    # create DEM object (??)
    dem = xdem.DEM(vrt_cropped_filepath, load_data=False)

    # plot DEM
    #dem.show()
    #plt.show()

    return dem

def mask_dem(dem, shp):
    # masks DEM data by the glacier area outlines (.shp)
    # run before sampling raster

    # rasterize the shapefile to fit the DEM
    shp_rasterized = gu.Vector(shp).create_mask(dem)
    shp_rasterized.show(cmap="Purples")
    #plt.show()

    # extract values inside the glacier area outlines
    dem.load()
    dem.set_mask(~shp_rasterized)


    #dem.show(cmap="Purples")
    #plt.show()

    return dem

def IS2_DEM_difference(dem, is2, label):
    # get elevation difference between IS2 data and a reference DEM

    # create empty list for delta_h (elevation differences)
    delta_h = []

    # for each IS2 measurement point
    for i in range(len(is2.easting)):
        # value of easting, northing and height from IS2 data
        e = is2.easting.values[i]
        n = is2.northing.values[i]
        elev_is2 = is2.h_te_best_fit.values[i]

        # value of elevation (DEM) for current pt in the
        elev_dem = dem.value_at_coords(e, n)

        # subtract the values
        elev_difference = elev_dem - elev_is2

        # if the elevation difference has a nodata value (outside of GAO boundaries), append nodata
        # else append elevation difference as float
        # maybe there's a smoother way of doing this - for example remove all the pts with nodata value from dataset (now there is gonna be a lot of nan values...)
        if elev_difference == dem.nodata:
            delta_h.append(np.nan)
        else:
            delta_h.append(float(elev_difference))

    # create new variable and assign elevation differences to existing IS2 dataset
    is2 = is2.assign(dh = delta_h)

    # values are assigned as coordinates, fixes them to be assigned as values
    is2_dh = is2.reset_index("dh").reset_coords("dh")

    # save as netcdf file
    cache_path = Path(f"cache/{label}-is2-dh.nc")
    is2_dh.to_netcdf(cache_path)

    return is2_dh

def hypsometric_binning(is2_path, ref_dem_path):
    # not working
    # open datasets with xarray
    is2 = xr.open_dataset(is2_path) # ICESat-2 data
    ref_dem = xr.open_dataset(ref_dem_path) # cropped DEM


    # write crs
    is2.rio.write_crs("epsg:32633", inplace=True)
    ref_dem.rio.write_crs("epsg:32633", inplace=True)

    # rename easting, northing to x, y (only necessary in IS2 data)
    is2 = is2.rename({'easting': 'x','northing': 'y'})

    # create a subset only containing variables x, y, dh
    #rast.reset_index("dh").reset_coords("dh") # resets dh values from coords to variable
    is2 = is2[['x', 'y', 'dh']]

    # convert tu ndarray
    is2 = is2.to_array()
    is2 = is2.to_numpy()
    ref_dem = ref_dem.to_array()
    ref_dem = ref_dem.to_numpy()

    #do hypsometric binning
    xdem.volume.hypsometric_binning(is2, ref_dem)

    return





if __name__ == "__main__":
    main()
