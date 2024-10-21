import pandas as pd

from . import gww_requests
from shapely.geometry import shape
import copy
import hydromt
from hydromt._utils.log import _setuplog as setuplog
import geopandas as gpd
from pyflwdir import FlwdirRaster
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import differential_evolution
import xarray as xr
from matplotlib.colors import LightSource


# default logger
logger = setuplog("./raster.log", log_level=10)

def power_law(h, a, b, h0):
    return a*(h-h0)**b

def inv_power_law(A, a, b, h0):
    return (A/a)**(1/b) + h0

def int_power_law(h, a, b, h0):
    return a/(b+1)*(h - h0)**(b+1)
# helper functions

def inv_int_power_law(V, a, b, h0):
    return (V/(a/(b+1)))**(1/(b+1)) + h0

def cell_surface(xi, yi, R=6371e3):
    """

    Parameters
    ----------
    xi
    yi
    R

    Returns
    -------

    """
    """
    compute surface from x and y coordinates in latlon
    """
    resx = xi[0][1] - xi[0][0]
    #     print(resx)
    resy = yi[1][0] - yi[0][0]
    lat1 = (yi - 0.5 * resy) / 180 * np.pi
    lat2 = (yi + 0.5 * resy) / 180 * np.pi
    lon1 = (xi - 0.5 * resx) / 180 * np.pi
    lon2 = (xi + 0.5 * resx) / 180 * np.pi

    dx = np.arccos(np.sin(lat1) * np.sin(lat1) + np.cos(lat1) * np.cos(lat1) * np.cos(lon2 - lon1)) * R
    dy = np.abs(2 * np.pi * R / 360 * resy)
    return dx * dy

def check_vars(ds, vars):
    for var in vars:
        if var not in ds.data_vars:
            raise ValueError(f'variable "{var}" is expected in Dataset but not found, available variables are {ds.data_vars}')

def wrap_rmse(pars_iter, *args):
    return rmse(pars_iter, *args)

def rmse(pars, func, X, Y):
    """mean of sum of squares between evaluation of function with provided parameters and X input, and Y as dependent variable.

    Parameters
    ----------
    pars : list or tuple
        parameter passed as *args to func
    func: function def
        receiving X and *pars as input and returning predicted Y as result
    X: tuple with lists or array-likes
        indepent variable(s).
    Y: list or array-like
        dependent variable, predicted by func

    Returns
    -------

    Y_pred : list or array-like
        predicted Y from X and pars by func
    """
    Y_pred = func(X, *pars)
    err = np.sqrt(np.sum((Y_pred - Y) ** 2))
    return err




class Bathymetry:
    """

    Hypsometry class that can
    - mask reservoir area
    - get bathymetry from geo functions (e.g. messager et al., 2016)
    - create hypso relationships
    - apply relationships on passed surface area time series

    """

    def __init__(
            self,
            polygon,
            crs=4326,
            dam_loc=None,
            downstream_reservoir=0.,
            above_reservoir=0.,
    ):
        """
        Bathymetry data with functionalities to find reservoir surface areas, mask these, and fill these with
        geomorphological logic using stream directions as input.

        Parameters
        ----------
        polygon : Polygon
            approximate area of the reservoir, only used to find a reasonable bounding box for topography data
        crs : int, optional
            coordinate reference system EPSG code belonging to anticipated retrieved topography dataset (default 4326)
        dam_loc : tuple (float, float), optional
            exact x, y coordinate of dam location (default None). If provided, algorithms will place the dam at set
            location. If not provided, a search will be performed to find the most likely dam location
        downstream_reservoir : float, optional
            distance downstream of dam location, to find representative in-stream elevation downstream of dam and
            mask existing elevation values, and fill in new elevation values (default: 0.). Recommended value 500 m.
        above_reservoir : float, optional
            amount of meters above the masked and interpolated area, used to extend water level - area relationship.
            This may be needed to ensure the hyposmetry curve extends far enough above the water surface during taking
            af the static elevation data. Particularly important when water levels were very low during the static
            terrain data acquisition (default: 0.)
        """
        self.mask = None
        self.ds_topography = None
        self.flwdir = None
        self.downstream_reservoir = downstream_reservoir
        self.above_reservoir = above_reservoir
        self.path = None
        self.dam_loc = dam_loc
        # response = gww_requests.get_reservoir(self.reservoir_id)
        # polygon = shape(response.json()["geometry"])
        # set polygon as a mask
        self.polygon = polygon
        self.crs = crs
        self.logger = logger


    @property
    def x(self):
        """
        Returns
        -------
        array-like
            x-axis of topography dataset
        """
        return self.ds_topography.x if self.ds_topography is not None else None

    @property
    def y(self):
        """
        Returns
        -------
        array-like
            y-axis of topography dataset
        """
        return self.ds_topography.y if self.ds_topography is not None else None

    @property
    def cell_area(self):
        """
        Returns
        -------
        array-like
            2d array with surface area per grid cell (in m2) of topography dataset
        """
        if self.x is not None:
            xi, yi = np.meshgrid(self.x, self.y)
            # get surface area per cell
            return cell_surface(xi, yi)
        else:
            return None


    @property
    def mask_stream(self):
        """
        Returns
        -------
        DataArray
            mask on cells with path on main stream
        """
        return self.ds_topography["path_boolean"]

    @property
    def gdf_polygon(self):
        """
        A Geopandas data frame representation of the provided polygon. Useful for plotting or saving to a GIS compatible
        vector file.
        
        Returns
        -------
        gdf : geopandas.GeoDataFrame
            contains one row with the provided Polygon as geodataframe representation. Can be used for plotting or easy
            export to geojson or GIS compatible files.
        """
        return gpd.GeoDataFrame(index=[1], crs=self.crs, geometry=[self.polygon])


    @property
    def downstream_basin(self):
        """
        Returns
        -------
        np.array
            boolean values indicating all cells that lie within or downstream of reservoir polygon
        """
        return self.flwdir.accuflux(self.mask.values) > 0

    @property
    def upstream_area(self):
        """

        Returns
        -------
        np.array
            boolean array showing all areas upstream of the estimated dam wall location
        """
        if self.dam_loc is None:
            dam_loc = self.get_downstream_limit()
            x, y = dam_loc.x.values[0], dam_loc.y.values[0]
        else:
            x, y = self.dam_loc
        return self.flwdir.basins(xy=(x, y))


    def read_topography(self, fn, **kwargs):
        """
        Read and set topography in attribute "ds_topography". The file must be xarray compatible (e.g. zarr or netcdf)
        and must at minimum contain the following data variables:
        "elevtn": 2D-elevation (float) [m]
        "flwdir": pyflwdir compatible flow directions (integer type) [-]
        "lndslp": slope of the land surface following stream direction (float) [-]
        "strord": stream order (int) [-]
        "uparea": upstream area (float) [m2]

        Parameters
        ----------
        fn : str,
            Path to file (e.g. netcdf or zarr) containing required data variables
        **kwargs : dict, optional
            additional keyword arguments to pass as flags (True/False) to ``self.set_topography``. these may be

            - ``add_flow_direction`` (default: True)
            - ``add_derivatives`` (default: True)
            - ``add_mask`` (default: True)

            It is recommended not to set these flags to False, so normally this is not required.

        Returns
        -------
        None

        """
        ds = xr.open_dataset(fn)
        self.set_topography(ds, **kwargs)


    def get_topography_from_hydromt(
            self,
            data_catalog,
            add_flow_direction=True,
            add_derivatives=True,
            add_mask=True,
            path_or_key="merit_hydro",
            **kwargs
    ):
        """
        Retrieve topography (including flow directions if indicated) from hydromt Data Catalogue. User must provide
        key name of dataset to retrieve. Topography is added as xarray.Dataset to object under self.ds_topography

        Parameters
        ----------
        add_flow_direction: bool, optional
            if True (default) then also flow directions are added to dataset
        add_derivatives: bool, optional
            if True (defaukt) also on-the-fly add the main stream (boolean field), slope, distance to most downstream
            point to topography dataset.
        add_mask: bool, optional
            add a mask raster (.mask) as burned layer of polygon of reservoir
        path_or_key: str, optional
            name (as in hydromt data catalogue) of dataset to retrieve as topography
        **kwargs: dict, optional
            set of key word arguments to pass to hydromt.data_catalogue.get_rasterdataset (e.g. buffer=30 for some
            buffering around reservoir is recommended)

        Returns
        -------


        """
        ds = data_catalog.get_rasterdataset(data_like=path_or_key, geom=self.gdf_polygon, **kwargs)
        ds = ds.compute()
        self.set_topography(
            ds,
            add_flow_direction=add_flow_direction,
            add_derivatives=add_derivatives,
            add_mask=add_mask
        )


    def set_topography(self, ds, add_flow_direction=True, add_derivatives=True, add_mask=True):
        """

        Parameters
        ----------
        ds
        add_flow_direction : bool, optional
            if True (default), convert data variable in dataset with name "flwdir" into a pyflwdir object, which is
            required to enable following stream directions. Needed for filling operations hence should usually
            be left as ``True``.
        add_derivatives : bool, optional
            if True (default) add several derivatives of elevation to dataset including the main flow path over the
            reservoir and the distance of each pixel to the outlet location within the topography dataset
        add_mask : bool, optional
            if True (default) add a masked topography, where cells that are on the reservoir surface are removed.
            this masked topography is filled with other methods to establish bathymetry in these area, so normally
            this is needed and ``add_mask`` should therefore be True.
        Returns
        -------
        None

        """
        if ds.raster.crs is None and self.crs != None:
            ds.raster.set_crs(self.crs)
        check_vars(ds, ["elevtn", "lndslp", "flwdir", "uparea", "strord"])
        self.ds_topography = ds
        if add_flow_direction:
            self.set_flow_direction()
        if add_mask:
            mask = self.ds_topography["elevtn"] * 0.
            mask = mask.raster.rasterize(self.gdf_polygon)
            self.mask = mask
        if add_derivatives:
            self.set_path()
            self.set_distance()
            self.set_mask_topography()

    def get_mask_topography(self):
        """
        makes a masked representation of topography as bathymetry layer by excluding areas upstream of dam wall and
        within a certain elevation limit. Can then be used to gradually fill up elevation through functionality in this class.

        Returns
        -------
        xr.DataArray
            2d array containing elevation with areas masked that are likely the water surface of the reservoir.
        """
        if self.dam_loc is None:
            # try and find downstream location
            xy_down = self.get_downstream_limit()
        else:
            xy_down = self.ds_topography.sel(
                x=[self.dam_loc[0]],
                y=[self.dam_loc[1]],
                method="nearest"
            )
        xy_up = self.get_upstream_limit()
        elevation_max = xy_up.elevtn.values[0]
        # mask out intermediate values
        elevtn_max = xy_up.elevtn.values[0]
        # mask out intermediate values
        bas = self.flwdir.basins(xy=(xy_down.x.values[0], xy_down.y.values[0]))
        # now get elevation and remove area within reservoir
        elevtn = copy.deepcopy(self.ds_topography["elevtn"].values)
        # filter out intermediate areas from elevation
        nodata_filter = np.all([bas, elevtn < elevtn_max], axis=0)
        elevtn[nodata_filter] = np.nan
        # add masked elev as "bathymetry" layer
        return elevtn

    def set_mask_topography(self):
        """
        Masks out areas in the "elevtn" data variable of the topography dataset that are likely on the reservoir
        surface and sets the masked elevation values on a new data variable "bathymetry".

        Returns
        -------
        None, sets additional variable on ``self.ds_topography["bathymetry"]``.
        """
        elevtn = self.get_mask_topography()
        self.ds_topography["bathymetry"] = (("y", "x"), elevtn)

    def set_flow_direction(self, name="flwdir"):
        """
        Convert flow directions from ordinary map to pyflwdir compatible flow direction map

        Parameters
        ----------
        name: str, optional
            name of DataArray in .ds_topography to convert (default: "flwdir")

        Returns
        -------
        None, sets pyflwdir object on ``self.flwdir``
        """
        assert(self.ds_topography is not None), "Topography not yet downloaded, first run `get_topography`"
        da = self.ds_topography[name]
        # try to set the spatial reference
        da.raster.set_crs(self.ds_topography.spatial_ref.crs_wkt)
        flwdir = hydromt.gis.flw.flwdir_from_da(da.astype(np.uint8), ftype="d8")
        self.flwdir = flwdir

        # self.flwdir = hydromt.flw.flwdir_from_da(self.ds_topography[name].astype(np.uint8), ftype="d8")

    def get_path(self, xy):
        """
        Returns the main flow path (i.e. with largest upstream area) from a provided xy location as boolean array

        Parameters
        ----------
        xy : tuple (float, float)
            Location for which to provide main flow path upstream

        Returns
        -------
        array-like (boolean)
            2d array with main flow path upstream of xy
        """
        # cache largest upstream per pixel
        largest_up_areacell = self.flwdir.main_upstream(uparea=self.ds_topography.uparea.values)
        path_up_from_start = FlwdirRaster.path(
            self.flwdir,
            xy=xy,
            unit='m',
            direction='up'
        )
        path_idx = list(path_up_from_start)[0]
        path = np.array(
            np.zeros(
                (len(self.ds_topography.y),
                 len(self.ds_topography.x))
            ),
            dtype="bool"
        )
        np.put(path, path_idx, True)
        return path

    def set_path(self, name="path_boolean"):
        """
        Sets main path across reservoir from most downstream point in the entire area of topography dataset
        on a new data variable "path_boolean"

        Returns
        -------
        None, upstream largest path set on ``self.ds_topography["path_boolean"]``
        """

        assert(self.ds_topography is not None), "Topography not yet downloaded, first run `get_topography`"
        uparea_masked = self.ds_topography.uparea.where(
            self.flwdir.accuflux(self.mask).astype(bool)
        )
        most_downstream = uparea_masked.where(
            uparea_masked == uparea_masked.max(),
            drop=True
        )

        # most_downstream = self.ds_topography.uparea.where(
        #     self.ds_topography.uparea == self.ds_topography.uparea.max(),
        #     drop=True
        # )
        # find downstream end
        xy = (most_downstream.x.values[0], most_downstream.y.values[0])
        path = self.get_path(xy)
        self.ds_topography[name] = (("y", "x"), path)


    def set_distance(self, name="distance"):
        """
        For each cell in the entire array get the distance to the start of path coordinate in meters
        Distance is added as new data variable in ``self.ds_topography``.

        Parameters
        ----------
        name: str, optional
            name of layer to add distance (default: "distance")

        Returns
        -------
        None, upstream largest path set on ``self.ds_topography["distance"]``

        """
        assert(self.flwdir is not None), "Flow directions not yet derived, first run `get_flow_direction`"
        self.distance = FlwdirRaster.stream_distance(self.flwdir, unit="m")
        self.ds_topography[name] = (('y', 'x') , self.distance)

    def skeletonize_bathymetry(self, iter=5, verbose=True):
        """
        fills in missing bathymetry in self.ds_topography["bathymetry"] starting from highest
        Strahler streamorder found within missing part of grid iterating several times to lower
        streamorder (i.e. smaller tributaries). Bathymetry is filled in following Messager et al. 2016, using

        dh = x**b * tan(s)

        where x is distance from shore of the lake surface, b is an exponent and s is the local slope at the found
        shoreline (i.e. upstream of missing area).

        After this operation ``self.ds_topography["bathymetry"]`` will contain filled in bathymetry values over
        the main stream, and tributaries of the main stream, tributaries of the tributaries, and so on...
        up until strahler order equal to the order of the main stream minus the value set for ``iter``. Values in even
        smaller tributaries can be interpolated with a 2-dimensional interpolation strategy using ``self.interpolate``.

        Parameters
        ----------
        iter : int, optional
            amount of iterations to go through whilst searching for streams. E.g. if the highest
            streamorder found in missing areas is 10, and iter=7 (default), then the lowest streamorder for which bathymetry values
            will be skeletonized is 4.

        Returns
        -------
        None, skeletonized bathymetry values set in larger streamorders in bathymetry layer in
            ``self.ds_topography["bathymetry"]``
        """
        strords = 0
        # find highest stream order to start with
        cur_str_order = int(self.ds_topography.strord.max().values.max())
        # now wrap func and loop
        while (strords < iter and cur_str_order >= 1):
            # update the downstream bathymetry map (changes with every step)
            bathym_down = self.flwdir.downstream(self.ds_topography.bathymetry.values)
            self.ds_topography["bathymetry_down"] = (("y", "x"), bathym_down)
            # find all locations with current streamorder, where the downstream bathym
            # already exists but where the value itself is NaN
            xys_down = self.ds_topography.where(
                self.ds_topography.bathymetry
            ).where(
                self.ds_topography.strord == cur_str_order
            ).where(
                np.isfinite(
                    self.flwdir.downstream(self.ds_topography.bathymetry)
                )
            ).where(
                np.isnan(self.ds_topography.bathymetry)
            )
            # turn the found points into a stacked layer with dim points
            xys_down = xys_down.stack(points=("y", "x"))
            # only save the non NaN locations
            xys_down = xys_down.where(np.isfinite(xys_down), drop=True)
            if len(xys_down.points) > 0:
                if verbose:
                    print(f"found {len(xys_down.points)} connecting points on stream order {cur_str_order}")
                # found a stream order within unmapped reservoir
                strords += 1
                # loop through all found points and interpolate stream bathym
                for n in range(len(xys_down.points)):
                    xy_down = xys_down.isel(points=n)
                    # get the entire upstream path of position (up until the known bathym)
                    mask_stream = self.get_path((xy_down.x.values, xy_down.y.values))
                    xy_up = self.get_upstream_limit(mask_stream)
                    if len(xy_up.points) != 0:
                        self.set_stream_bathymetry(xy_up=xy_up, xy_down=xy_down, mask_stream=mask_stream)
                    else:
                        # when cell is already at the edge of the shore, then nothing will be found
                        if verbose:
                            print(f"No upstream point found for downstream point {n}")
            # move to one order smaller streams
            cur_str_order -= 1

    def plot_2d(self, ax=None, field="elevtn", mask_stream=True, **kwargs):
        """
        Make a simple spatial plot of a layer in ``self.ds_topography``.

        Parameters
        ----------
        ax : matplotlib axes object
            existing axes object to plot variable on. If not provided, an axes will be created.
        field: str, optional
            name of layer to plot in 2d
        mask_stream : bool, optional
            set whether stream pixels should be masked or not
        **kwargs : dict, optional
            keyword arguments to pass to plot.

        Returns
        -------
        ax : matplotlib axes
            axes containing plot.

        """
        if ax is None:
            ax = plt.axes()
        if mask_stream:
            self.ds_topography[field].where(self.mask_stream).plot(**kwargs)
        else:
            self.ds_topography[field].plot(**kwargs)
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlabel('lon [degrees east]')
        ax.set_ylabel('lat [degrees north]')
        return ax


    def plot_topo(self, ax=None, field="elevtn", mask_stream=True, ls_kwargs={"azdeg": 315, "altdeg": 45}, **kwargs):
        """
        Make a shaded relief plot of an elevation layer

        Parameters
        ----------
        ax : matplotlib axes object
            existing axes object to plot variable on. If not provided, an axes will be created.
        field: str, optional
            name of layer to plot in 2d
        mask_stream : bool, optional
            set whether stream pixels should be masked or not
        ls_kwargs : dict, optional
            dictionary of keyword - args to pass to ``matplotlib.colors.LightSource``, default
            ``{"azdeg": 315, "altdeg": 45}``.
        **kwargs : dict, optional
            keyword arguments to pass to ``LightSource.shade``.
        Returns
        -------
        ax : matplotlib axes
            axes containing plot.

        """
        ls = LightSource(**ls_kwargs)
        cmap = plt.cm.gist_earth
        if ax is None:
            ax = plt.axes()
        rgb = ls.shade(
            self.ds_topography[field].values,
            cmap=cmap,
            dx=np.sqrt(self.cell_area[0][0]),
            dy=np.sqrt(self.cell_area[0][0]),
            **kwargs
        )
        ax.imshow(rgb)
        ax.set_aspect('equal', adjustable='box')
        ax.set_title(f'{field}')
        return ax


    def plot_along_stream(self, ax=None, field="elevtn", **kwargs):
        """
        Plot values of topography layer along the main stream only

        Parameters
        ----------
        ax : matplotlib axes object
            existing axes object to plot variable on. If not provided, an axes will be created.
        field: str, optional
            name of layer to plot in 2d
        **kwargs : dict, optional
            keyword arguments to pass to ``matplotlib.pyplot.plot``.

        Returns
        -------
        ax : matplotlib axes
            axes containing plot.

        """
        # extract and mask values along stream
        da_along_stream = self.ds_topography[field].where(self.mask_stream)
        # squeeze y and x into one dimension
        da_along_stream = da_along_stream.stack(distance=["y", "x"])
        stream_distance = self.ds_topography.distance.values.flatten()
        # replace coordinates by distance along stream
        da_along_stream = da_along_stream.assign_coords({"distance": stream_distance})
        # filter out all nan
        da_along_stream = da_along_stream.dropna(dim="distance")
        # sort
        da_along_stream = da_along_stream.sortby("distance")
        # plot
        if ax is None:
            ax = plt.axes()
        da_along_stream.plot(ax=ax, **kwargs)
        ax.set_title(f"{field}")
        return ax

    def get_upstream_limit(self, mask_stream=None):
        """
        get upstream start location for interpolation with Messager et al. method.

        Parameters
        ----------
        mask_stream : array-like, boolean, optional
            boolean with stream locations to use to search for upstream limits. If not provided, the main stream
            will be used.
        Returns
        -------
        tuple (x, y) : float, float
            coordinate of the upstream select in-stream point for bathymetry interpolation

        """
        if mask_stream is None:
            # use the main stream to find the right spot, else, use defined stream mask
            mask_stream = self.mask_stream
        else:
            mask_stream = mask_stream.astype(bool)
        # first limit selection to points upstream of reservoir by masking on reservoir, stream and (not) downstream basin
        filtered = self.ds_topography.where(self.mask == False).where(mask_stream).where(
            self.downstream_basin == False)

        # find x meters above lowest and define point
        # filtered = filtered.where(
        #     self.flwdir.upstream_sum(self.ds_topography.elevtn.values) > float(filtered.elevtn.min()) + self.above_reservoir
        # )
        filtered = filtered.where(filtered.elevtn > filtered.elevtn.min() + self.above_reservoir)

        xy = filtered.where(filtered.distance == filtered.distance.min(), drop=True)
        return xy.stack(points=("y", "x"))


    def get_hypsometry(self, add_to_max_h=0., resolution=1.):
        """

        Parameters
        ----------
        add_to_max_h : float, optional
            when set to > 0. (default) extend the bathymetry that amount above the estimated dam elevation
        resolution : float, optional
            set the vertical resolution [m] to sample hypsometry elevation - area points. Default: 1.

        Returns
        -------
        Hypsometry
            Object containing area and level points and methods to fit relationships and estimate one variable from
            another using those relationships (see Hypsometry class for more information)

        """
        bathymetry = self.ds_topography["bathymetry"].where(self.upstream_area)
        h_min = self.get_downstream_limit().elevtn.values[0]
        h_max = self.get_upstream_limit().elevtn.values[0] + add_to_max_h
        water_level = np.arange(
            int(np.floor(h_min)),
            int(np.ceil(h_max)),
            resolution
        )

        # get area per water level
        surface_area = [float(((bathymetry <= h) * self.cell_area).sum(dim="x").sum(dim="y").values) for h in
                water_level]  # km2

        return Hypsometry(np.float64(water_level), np.float64(surface_area))


    def get_downstream_limit(self):
        """

        get downstream end location for interpolation and return topography data at that point.
        The dam location is searched as the largest slope within the stream mask. Then a x meters below (vertical) or downstream (horizontal) reservoir point is searched and returned

        Returns
        -------
        xy : xr.Dataset
            topography variables and coordinates for the downstream (single point) location

        """

        # find dam wall counting from within downstream area

        # first limit selection to points within and downstream of reservoir
        filtered = self.ds_topography.where(self.mask_stream).where(
            self.downstream_basin)
        # compute slope along the stream
        river_slope = self.flwdir.downstream(filtered.elevtn) - filtered.elevtn
        # then search for likely position of dam wall
        # dist_dam = filtered.where(np.abs(filtered.lndslp) == np.abs(filtered.lndslp).max(), drop=True).distance.values[0]
        dist_dam = filtered.where(np.abs(river_slope) == np.abs(river_slope).max(), drop=True).distance.values[0]

        # first get rid of anything above the dam minus some distance downstream
        filtered = filtered.where(filtered.distance <= dist_dam - self.downstream_reservoir) #1km downstream of dam wall

        # get the minimum distance (i.e. closest to outlet
        xy = filtered.where(filtered.distance == filtered.distance.max(), drop=True)
        return xy


    def set_stream_bathymetry(self, xy_up=None, xy_down=None, mask_stream=None):
        """
        Fill in bathymetry over one stream path, defined by an upstream and downstream point in the flow network.
        Infilling is done using the Messager et al. power law using the slope at the shore and the height difference
        between the upstream and downstream point as constraints to find the power as follows:

        dh = i * x**b

        Where:
        - dh is the difference in elevation between shore and the point of interest [m]
        - x is distance from shore [m]
        - b is a power that defines the decline in slope over distance [-]


        fitting the power

        Parameters
        ----------
        xy_up : xr.Dataset, optional
            topography data at upstream point in flow network. If not defined, will be derived from the masked
            topography following the main stream through the masked area.
        xy_down : xr.Dataset, optional
            topography data at downstream point in flow network. If not defined, will be derived from the masked
            topography following the main stream through the masked area.
        mask_stream :
            mask of cells over which interpolation should be performed. If not defined, the main stream
            will be used.

        Returns
        -------
        None, bathymetry values set on defined stream in bathymetry layer in ``self.ds_topography["bathymetry"]``

        """
        if mask_stream is None:
            # fill in the main stream
            mask_stream = self.mask_stream
        if xy_up is None:
            xy_up = self.get_upstream_limit()
        if xy_down is None:
            xy_down = self.get_downstream_limit()
        distance = self.ds_topography.where(
            self.ds_topography.distance < xy_up.distance.values
        ).where(
            self.ds_topography.distance >= xy_down.distance.values
        ).where(mask_stream).distance
        # add or update a layer with downstream bathymetry values
        # reduce by minimal distance found
        distance = distance - distance.min()
        # turn around so that distance at upstream end becomes zero and dam side maximum
        distance = distance.max() - distance

        # define parameter b
        log_slope = np.log(xy_up.lndslp.values)
        log_hdiff = np.log(xy_up.elevtn.values - xy_down.bathymetry_down.values)
        log_dist = np.log(xy_up.distance.values - xy_down.distance.values)

        b = (log_hdiff - log_slope) / log_dist

        # apply
        depth = xy_up.lndslp.values * distance ** b
        # depth.plot()
        # subtract elevation value
        elev_stream = xy_up.elevtn.values - depth

        bath = self.ds_topography["bathymetry"].values
        # bath.plot()

        bath[np.isfinite(elev_stream.values)] = elev_stream.values[np.isfinite(elev_stream.values)]
        self.ds_topography["bathymetry"][:] = bath

    def interpolate(self, method="cubic", **kwargs):
        """
        Linear interpolation of missing bathymetry records. Should be done only after performance of all steps for
        intermediate data such as ``self.skeletonize_bathymetry``.

        Parameters
        ----------
        method : str, optional
            interpolation method used to interpolate missing bathymetry data (after infilling stream elevations).
            passed to ``xr.Dataset.raster.interpolate_na``.
        **kwargs : dict, optional
            keyword arguments passed to ``xr.Dataset.raster.interpolate_na``.

        Returns
        -------
        None, missing bathymetry values in ``ds_topography["bathymetry"]`` filled with interpolated values

        """
        bath = self.ds_topography["bathymetry"]
        bath.raster.set_nodata(np.nan)
        # bath.raster.set_crs(self.ds_topography[])
        bath = bath.raster.interpolate_na(method=method, **kwargs)
        # remove the spatial_ref property (from rioxarray)
        if "spatial_ref" in bath.coords:
            bath = bath.drop("spatial_ref")
        self.ds_topography["bathymetry"] = bath


class Hypsometry():
    """
    class containing hypsometry points, derivation of relationships as properties, and derivation of volumes
    from surface area or water level.

    """
    def __init__(self, water_level, area):
        self.area = area
        self.water_level = water_level
        self.params = None


    @property
    def water_level(self):
        return self._water_level

    @water_level.setter
    def water_level(self, water_level):
        if water_level is None:
            self._water_level = None
        else:
            self._water_level = water_level


    @property
    def area(self):
        return self._area

    @area.setter
    def area(self, surface_area):
        if surface_area is None:
            self._area = None
        else:
            self._area = surface_area

    @property
    def volume(self):
        return self.volume_from_wl(self.water_level)


    @property
    def area_fit(self):
        return self.area_from_wl(self.water_level)


    @property
    def volume_fit(self):
        return self.volume


    @property
    def power_law_params(self):
        """
        parameters of fit through values in ``self.water_level`` and ``self.area`` following power law:


        A(h) = a * (h - h0)**b

        Where A [m2] is the surface area of a lake, h [m] is the water level of a lake, and a [-], h0 [m] and b[-] are
        fit parameters. First guess of the parameters is estimated by a log-transform linear fit, and improved using a
        curve fitting method.

        Returns
        -------
        dict, parameters of fit, with keynames ``a``, ``b`` and ``h0``.
        """

        def func_powerlaw(x, a, b, h0):
            return a * np.maximum(x - h0, 0) ** b

        if self.params is None:

            # def func_powerlaw(x, a, b):
            #     return a * x ** b
            # first guess
            X = self.water_level
            Y = self.area
            bounds = ([0.001, 10000], [0.5, 4], [-1000, self.water_level.min()])

            # x = np.log(self.water_level - self.water_level.min())
            # y = np.log(self.area)
            # y = y[np.isfinite(x)]
            # x = x[np.isfinite(x)]
            # slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
            # a = np.exp(intercept)
            # b = slope
            # h0 = self.water_level.min()
            # improvement with a curve_fit over power law a, b (h0 already known)
            # popt, pcov = curve_fit(
            #     func_powerlaw,
            #     self.water_level - self.water_level.min(),
            #     self.area,
            #     p0=[a, b]
            # )

            result = differential_evolution(
                wrap_rmse,
                args=(func_powerlaw, X, Y),
                bounds=bounds,
                workers=1,
                popsize=10,
                updating="deferred",
                seed=0,
            )
            # unravel parameters
            a, b, h0 = result.x
            self.params = {"a": a, "b": b, "h0": h0}
        return self.params
        #
        # return {
        #     "a": popt[0],
        #     "b": popt[1],
        #     "h0": h0
        # }

    @property
    def wl_area_table(self):
        """
        Create a pandas table representation of the water level and area pairs

        Returns
        -------
        pd.DataFrame
        """
        df = pd.DataFrame(
            {"water_level": self.water_level, "surface_area": self.area},
            index=np.arange(len(self.water_level)) + 1
        )
        df.index.name = "index"
        return df
    def area_from_wl(self, wl):
        """
        Returns computed surface area from water levels, using the power law fit parameters

        Parameters
        ----------
        wl : float or array-like with floats
            water level(s) to compute surface area(s) for.

        Returns
        -------
        area : float or array-like with floats
            surface area(s), same length as wl
        """
        return power_law(wl, **self.power_law_params)


    def wl_from_area(self, area):
        """
        Returns computed water level from surface area, using the power law fit parameters

        Parameters
        ----------
        area : float or array-like with floats
            surface area(s), to compute water level(s) for.

        Returns
        -------
        wl : float or array-like with floats
            water level(s), same length as area

        """
        return inv_power_law(area, **self.power_law_params)


    def volume_from_wl(self, wl):
        """
        Returns computed volume from water levels, using the power law fit parameters

        Parameters
        ----------
        wl : float or array-like with floats
            water level(s) to compute volume(s) for.

        Returns
        -------
        volume : float or array-like with floats
            volume(s), same length as wl

        """
        return int_power_law(
            wl,
            **self.power_law_params
        )

    def wl_from_volume(self, V):
        return inv_int_power_law(
            V,
            **self.power_law_params
        )
    def volume_from_area(self, area):
        """
        Returns computed volumes from area, using the power law fit parameters

        Parameters
        ----------
        area : float or array-like with floats
            area(s) to compute volumes(s) for.

        Returns
        -------
        volume : float or array-like with floats
            volume(s), same length as area

        """
        return int_power_law(
            self.wl_from_area(area),
            **self.power_law_params
        )

    def plot(self, x="water_level", y="area", ax=None, add_fit=True):
        """
        Plot the hypsometric points and fits for user provided variables.

        Parameters
        ----------
        x : str, optional
            variable to plot on x-axis, can be "water_level" (default), "area" or "volume"
        y : str, optional
            variable to plot on y-axis, can be "water_level" (default), "area" or "volume"
        ax : matplotlib.axes, optional
            axes object to plot in, if not provided, and axes object will be made.
        add_fit : bool, optional
            Also plot the fitted relationship (default: True)

        Returns
        -------
        ax : axes object

        """
        p = []
        if ax is None:
            ax = plt.axes()
        p1, = ax.plot(getattr(self, x), getattr(self, y), "x", label="samples", zorder=2)
        p.append(p1)
        if add_fit:
            if x == "water_level":
                # just use values as is
                xvals = self.water_level
            else:
                xvals = getattr(self, x + "_fit")
            if y == "water_level":
                # just use values as is
                yvals = self.water_level
            else:
                yvals = getattr(self, y + "_fit")
            p2, = ax.plot(
                xvals,
                yvals,
                "ro",
                zorder=1,
                label="fit a: {:1.1f}, b: {:1.2f}, h0: {:1.1f}".format(
                    self.power_law_params["a"],
                    self.power_law_params["b"],
                    self.power_law_params["h0"]
                )
            )
            p.append(p2)
        ax.grid()
        ax.set_xlabel(x.replace("_", " "))
        ax.set_ylabel(y.replace("_", " "))
        ax.legend(handles=p)
        return ax

    def to_file(self, fn, **kwargs):
        """
        Write water level and area to file

        Parameters
        ----------
        fn : str,
            Path to file

        """
        self.wl_area_table.to_csv(fn, **kwargs)
