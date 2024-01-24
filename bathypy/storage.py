from functools import partial
from typing import Callable, List

import ee
import numpy as np
from tqdm import tqdm
import pandas as pd
import geopandas as gpd
from eepackages import assets, utils


# functions related to storage detection and time series
def get_waterbody_images(
    waterbody, start, stop, scale, opt_missions=None, mostly_clean=True
):
    """
    Get mostly clean images on the waterbody

    :param waterbody: ee.Feature, waterbody polygon
    :param start: str, date stop '2018-01-01'
    :param stop: str, date stop '2020-01-01'
    :param scale: int, basic scale of computation
    :param opt_missions: Optional, list of satellite products to include ['L4', 'L5', 'L7', 'L8', 'S2']
    """
    geom = ee.Feature(waterbody).geometry()

    missions = ["L4", "L5", "L7", "L8", "S2"]

    if opt_missions:
        missions = opt_missions

    images = assets.getImages(
        geom,
        {
            "resample": True,
            "filter": ee.Filter.date(start, stop),
            "missions": missions,
            "scale": scale * 10,
        },
    )

    options = {"scale": scale * 5}

    g = geom.buffer(300, scale)

    if mostly_clean:
        return assets.getMostlyCleanImages(images, g, options).sort("system:time_start")
    else:
        return (
            assets.addQualityScore(images, g, options).filter(
                ee.Filter.gt("quality_score", 0)
            )
        ).sort("system:time_start")


def compute_surfacewaterarea_singleimage(image, waterbody, scale, water_occurrence):
    """
    Compute surface water area in an individual image

    :param image: ee.Image, source image
    :param waterbody: ee.Feature, waterbody polygon
    :param scale: int, basic scale of computation
    :param water_occurrence: ee.Image water occurrence layer
    """

    geom = ee.Feature(waterbody).geometry()
    fill_percentile = 75  # Low confidence in our prior
    ndwi_bands = ["green", "swir"]

    water_max_image = ee.Image().float().paint(waterbody.buffer(150), 1)
    max_area = waterbody.area(scale)

    t = image.get("system:time_start")

    image = image.updateMask(water_max_image).updateMask(
        image.select("swir").min(image.select("nir")).gt(0.001)
    )

    ndwi = image.normalizedDifference(ndwi_bands)

    # TODO: make otsu thresholding server-side compatible
    th = utils.computeThresholdUsingOtsu(ndwi, scale, geom, 0.5, 0.7, -0.2)
    # th = 0.2

    water = ndwi.gt(th)

    area = (
        ee.Image.pixelArea()
        .mask(water)
        .reduceRegion(**{"reducer": ee.Reducer.sum(), "geometry": geom, "scale": scale})
        .get("area")
    )

    # fill missing, estimate water probability as the lowest percentile.
    water_edge = ee.Algorithms.CannyEdgeDetector(ndwi, 0.5, 0.7)

    # image mask
    image_mask = ndwi.mask()
    image_mask = image_mask.focal_min(
        ee.Number(scale).multiply(1.5), "square", "meters"
    )

    # TODO: exclude non-water/missing boundsry
    water_edge = water_edge.updateMask(image_mask)

    # get water probability around edges
    # P(D|W) = P(D|W) * P(W) / P(D) ~=  P(D|W) * P(W)
    p = (
        water_occurrence.mask(water_edge)
        .reduceRegion(
            **{
                "reducer": ee.Reducer.percentile([fill_percentile]),
                "geometry": geom,
                "scale": scale,
            }
        )
        .values()
        .get(0)
    )

    # TODO: exclude edges belonging to cloud/water or cloud/land
    # TODO: use multiple percentiles (confidence margin)

    p = ee.Algorithms.If(ee.Algorithms.IsEqual(p, None), 101, p)

    water_fill = water_occurrence.gt(ee.Image.constant(p)).updateMask(
        water.unmask(0, False).Not()
    )

    # exclude false-positive, where we're sure in a non-water
    non_water = ndwi.lt(-0.15).unmask(0, False)
    water_fill = water_fill.updateMask(non_water.Not())

    fill = (
        ee.Image.pixelArea()
        .mask(water_fill)
        .reduceRegion(**{"reducer": ee.Reducer.sum(), "geometry": geom, "scale": scale})
        .get("area")
    )

    area_filled = ee.Number(area).add(fill)
    filled_fraction = ee.Number(fill).divide(area_filled)

    return (
        image.addBands(water_fill.rename("water_fill"))
        .addBands(water_edge.rename("water_edge"))
        .addBands(ndwi.rename("ndwi"))
        .addBands(water.rename("water"))
        .set(
            {
                "p": p,
                "area": area,
                "area_filled": area_filled,
                "filled_fraction": filled_fraction,
                "system:time_start": t,
                "ndwi_threshold": th,
                "waterOccurrenceExpected": water_occurrence.mask(water_edge),
            }
        )
    )


def compute_surfacewaterarea_JRC(waterbody, start, stop, scale):
    """
    Compute surface area from JRC water occurrence

    :param waterbody: ee.Feature, waterbody polygon
    :param start: str, date stop '2018-01-01'
    :param stop: str, date stop '2020-01-01'
    :param scale: int, basic scale of computation
    """
    geom = ee.Feature(waterbody).geometry()

    jrc_monthly = ee.ImageCollection("JRC/GSW1_0/MonthlyHistory")

    def compute_area(i):
        area = (
            ee.Image.pixelArea()
            .mask(i.eq(2))
            .reduceRegion(
                {"reducer": ee.Reducer.sum(), "geometry": geom, "scale": scale}
            )
        )
        return i.set({"area": area.get("area")})

    water = jrc_monthly.filterDate(start, stop).map(compute_area)

    return water


def compute_filled_waterarea(image, waterbody):
    """
    Compute the water mask from image and combines labels non-water:0, water_mask:1, filled_water_mask:2

    :param image: ee.Image, individual image
    :param waterbody: ee.Feature, waterbody polygon
    """
    water_occurrence = (
        ee.Image("JRC/GSW1_3/GlobalSurfaceWater")
        .select("occurrence")
        .unmask(0)
        .resample("bicubic")
        .divide(100)
    )

    water_occurrence = water_occurrence.mask(water_occurrence)

    scale = waterbody.geometry().area().sqrt().divide(200).max(10)

    image_area = compute_surfacewaterarea_singleimage(
        image, waterbody, scale, water_occurrence
    )

    # properties = ['MISSION', 'ndwi_threshold', 'quality_score', 'area_filled', 'filled_fraction', 'p', 'system:time_start', 'area']
    # properties_new = ["mission", "ndwi_threshold", "quality_score", "water_area_filled", "water_area_filled_fraction", "water_area_p", "water_area_time", "water_area_value"]

    # Combining filled and raw water masks as 0: no water, 1: water mask, 2: water_occurrence filled water mask
    combined_watermask = (
        image_area.select("water")
        .unmask(0, False)
        .add(
            image_area.select("water_fill")
            .unmask(0, False)
            .where(image_area.select("water_fill").eq(1), 2)
        )
    )

    return image_area.addBands(combined_watermask.rename("combined_water_mask"))


def time_distance(image, icesat_date):
    """
    Function that returns the distance between the icesat_date and the image date (in millis)

    :param image: ee.Image, individual image
    :param icesat_date: str: date e.g. '2018-01-01'
    """
    return image.set(
        "date_dist",
        ee.Number(image.get("system:time_start"))
        .subtract(ee.Date(icesat_date).millis())
        .abs(),
    )


def estimate_watermask_icesatrack(
    clean_images, icesat_date, waterbody, max_time_window
):
    """
    Function that returns a date-distance ordered list of the closest water masks associated to the icesatrack date

    :param image: ee.Image, individual image
    :param icesat_date: str: date e.g. '2018-01-01'
    :param waterbody: ee.Feature, waterbody polygon
    :param max_time_window: int, max distance to check for masks in days [centered window]

    """

    # filter images to be near icesat
    icesat_date = ee.Date(icesat_date)
    start = icesat_date.advance(-max_time_window, "day")
    stop = icesat_date.advance(max_time_window, "day")
    images = clean_images.filterDate(start, stop)

    water = images.map(lambda i: compute_filled_waterarea(i, waterbody))

    water = water.filter(
        ee.Filter.And(
            ee.Filter.neq("p", 101),
            ee.Filter.gt("ndwi_threshold", -0.15),
            ee.Filter.lt("ndwi_threshold", 0.5),
            ee.Filter.lt("filled_fraction", 0.6),
        )
    )

    # get time distance to mask
    water = water.map(lambda i: time_distance(i, icesat_date)).sort(
        "date_dist"
    )  # sorted from closest

    return water


def get_watermasks_icesat2(
    clean_images, icesatdataframe, waterbody, max_time_window=10, subsample_is2_frac=0.1
):
    """
    Function that loops in a icesatdataframe and per unique date extracts the matching water mask.

    TODO: OPTIMIZE LOOP

    :param clean_image: ee.ImageCollection, clean images from get_waterbody_images()
    :param icesatdataframe: dataframe of icesat2 tracks
    :param waterbody: ee.Feature, waterbody polygon
    :param max_time_window: int, max distance to check for masks in days [centered window]
    :param subsample_is2_frac: float, 0-1 fraction for subsampling icesat2 tracks (reduces comp effort)

    Two dataframes are returned:
    :output df_is2: icesat2 dataframe subsampled with values for: water_mask, water_occurrence and distance from water_mask border
    :output water_level: dataframe with unique dates, mask area (filled) and icesat2 median height for pixels labelled as water (and water_filled)
    """

    dates_ic2 = np.unique(icesatdataframe.date)

    samples_dates = []
    dates = []
    water_area_sampled = []

    for ic_date_i in tqdm(dates_ic2):
        try:
            water_area = estimate_watermask_icesatrack(
                clean_images, ic_date_i, waterbody, max_time_window
            )

            sampled_is2 = icesatdataframe[icesatdataframe.date == ic_date_i].sample(
                frac=subsample_is2_frac
            )

            # select closest water mask in time
            water_mask = (
                water_area.first().select("combined_water_mask").rename("water")
            )  ### NOTE: consider refining the date matching, now we only select closest, include perhaps also other filters e.g. quality

            # add distance to waterbody shoreline as additional variable
            distance_int = (
                water_mask.gt(0)
                .Not()
                .fastDistanceTransform(256)
                .reproject(ee.Projection("EPSG:3857").atScale(10))
                .sqrt()
                .rename("dist_m")
            )

            distance_ext = (
                water_mask.gt(0)
                .fastDistanceTransform(256)
                .reproject(ee.Projection("EPSG:3857").atScale(10))
                .sqrt()
                .multiply(-1)
                .rename("dist_m")
            )

            distance = distance_int.add(distance_ext)

            water_occurrence = (
                ee.Image("JRC/GSW1_3/GlobalSurfaceWater")
                .select("occurrence")
                .unmask(0)
                .resample("bicubic")
                .divide(100)
            )

            water_mask = water_mask.addBands(distance).addBands(water_occurrence)

            # get client side results by batch
            results_w = []
            results_d = []
            results_oc = []

            BATCH_SIZE = 40000

            for i in range(0, len(sampled_is2), BATCH_SIZE):
                coords = list(
                    zip(
                        sampled_is2.lon.values[i : (i + BATCH_SIZE)],
                        sampled_is2.lat.values[i : (i + BATCH_SIZE)],
                    )
                )
                features = ee.List(coords).map(
                    lambda o: ee.Feature(ee.Geometry.Point(coords=o))
                )
                features = ee.FeatureCollection(features)

                # NOTE: explore strategy to reduce the three properties together
                values_w = water_mask.unmask(-999, False).reduceRegions(
                    collection=features,
                    reducer=ee.Reducer.first().setOutputs(["water"]),
                    scale=10,
                    tileScale=1,
                )
                values_d = water_mask.unmask(-999, False).reduceRegions(
                    collection=features,
                    reducer=ee.Reducer.first().setOutputs(["dist_m"]),
                    scale=10,
                    tileScale=1,
                )
                values_oc = water_mask.unmask(-999, False).reduceRegions(
                    collection=features,
                    reducer=ee.Reducer.first().setOutputs(["occurrence"]),
                    scale=10,
                    tileScale=1,
                )

                results_w = results_w + values_w.aggregate_array("water").getInfo()
                results_d = results_d + values_d.aggregate_array("dist_m").getInfo()
                results_oc = (
                    results_oc + values_oc.aggregate_array("occurrence").getInfo()
                )

                sampled_is2_mask = sampled_is2.copy()
                sampled_is2_mask["water"] = results_w
                sampled_is2_mask["dist_m"] = results_d
                sampled_is2_mask["w_occurrence"] = results_oc
                samples_dates.append(sampled_is2_mask)
                dates.append(ic_date_i)
                water_area_sampled.append(water_mask.get("area_filled").getInfo())

        except:
            print("empty date")
            continue

    df_is2 = pd.concat(samples_dates)
    df_is2[df_is2.w_occurrence < 0] = np.nan

    water_level = pd.DataFrame(df_is2[df_is2.water > 0].groupby("date")["h"].median())
    areas = pd.DataFrame(water_area_sampled, index=dates, columns=["water_area"])
    water_level["water_area"] = areas.loc[water_level.index]["water_area"]

    return df_is2, water_level


def sample_from_DEM_datasets(region):
    """
    Obtain 5000 samples from GEDI, ALOS, NASADEM DEMs on the water occurrence region around a reservoir

    TODO: allow for a larger number of samples (batch, getInfo limits to 5000 features per call)
    TODO: use spatial correlation DEM - WATERoccurrence to bias the sampling scheme

    :param region: ee.Feature, waterbody polygon

    Three dataframes are returned:
    :output sample_gedi: geopandasDF sample_gedi samples (elevation, water_occurrence)
    :output sample_alos: geopandasDF sample_alos samples (elevation, water_occurrence)
    :output sample_nasadem: geopandasDF sample_nasadem samples (elevation, water_occurrence)
    """
    aoi_bounds = region.geometry().buffer(500)

    geoid = ee.Image("users/gena/egm96-15")
    alos = ee.ImageCollection("JAXA/ALOS/AW3D30/V3_2")
    srtm = ee.Image("NASA/NASADEM_HGT/001").select("elevation")
    glo30 = ee.ImageCollection("projects/sat-io/open-datasets/GLO-30")
    glo30elev = glo30.mosaic().setDefaultProjection(glo30.first().projection())

    gedi = ee.ImageCollection("LARSE/GEDI/GEDI02_A_002_MONTHLY").filterBounds(
        aoi_bounds
    )
    gedi = gedi.map(
        lambda x: x.float().set({"bandCount": x.bandNames().size()})
    ).filter(ee.Filter.eq("bandCount", 138))

    water_occurrence = (
        ee.Image("JRC/GSW1_3/GlobalSurfaceWater")
        .select("occurrence")
        .divide(100)
        .unmask(0)
        .resample("bilinear")
        .reduceNeighborhood(ee.Reducer.mean(), ee.Kernel.circle(20, "meters"))
    )
    water_occurrence = water_occurrence.updateMask(water_occurrence.unitScale(0, 0.04))
    waterOccurrenceMin = 0.01

    mask = water_occurrence.gt(waterOccurrenceMin)
    water_occurrence = water_occurrence.mask(mask)

    # stratified sampling based on water_occurrence probability (50000 samples, scale 20 meters)
    n = 5000
    scale = 20

    samplesGEDI = (
        water_occurrence.addBands(
            gedi.select(["elev_lowestmode", "delta_time"])
            .mean()
            .rename(["GEDI", "delta_time"])
            .subtract(
                ee.Image(geoid.resample("bicubic")).addBands(ee.Image.constant(0))
            )
        )
        .updateMask(water_occurrence.gt(0.01))
        .addBands(water_occurrence.multiply(100).rename("class").int())
        .stratifiedSample(
            numPoints=n,
            classBand="class",
            region=aoi_bounds,
            scale=scale,
            seed=42,
            dropNulls=True,
            tileScale=6,
            geometries=True,
        )
    )

    samplesALOS = (
        water_occurrence.addBands(
            alos.select("DSM")
            .map(
                lambda i: i.resample("bicubic")
                .rename("ALOS")
                .convolve(ee.Kernel.gaussian(60, 45, "meters"))
            )
            .mosaic()
            .updateMask(mask)
        )
        .addBands(water_occurrence.multiply(100).rename("class").int())
        .stratifiedSample(
            numPoints=n,
            classBand="class",
            region=aoi_bounds,
            scale=scale,
            seed=42,
            dropNulls=True,
            tileScale=6,
            geometries=True,
        )
    )

    samplesSRTM = (
        water_occurrence.addBands(
            srtm.resample("bicubic")
            .rename("NASADEM")
            .convolve(ee.Kernel.gaussian(60, 45, "meters"))
        )
        .addBands(water_occurrence.multiply(100).rename("class").int())
        .stratifiedSample(
            numPoints=n,
            classBand="class",
            region=aoi_bounds,
            scale=scale,
            seed=42,
            dropNulls=True,
            tileScale=6,
            geometries=True,
        )
    )

    samplesGLO = (
        water_occurrence.addBands(
            glo30elev.resample("bicubic")
            .rename("GLO")
            .convolve(ee.Kernel.gaussian(60, 45, "meters"))
        )
        .addBands(water_occurrence.multiply(100).rename("class").int())
        .stratifiedSample(
            numPoints=n,
            classBand="class",
            region=aoi_bounds,
            scale=scale,
            seed=42,
            dropNulls=True,
            tileScale=6,
            geometries=True,
        )
    )

    sample_gedi = samplesGEDI.randomColumn("r").sort("r").limit(n).getInfo()
    samples_alos = samplesALOS.randomColumn("r").sort("r").limit(n).getInfo()
    samples_srtm = samplesSRTM.randomColumn("r").sort("r").limit(n).getInfo()
    samples_glo = samplesGLO.randomColumn("r").sort("r").limit(n).getInfo()

    sample_gedi_df = pd.DataFrame(
        [
            (
                v["properties"]["GEDI"],
                v["properties"]["occurrence_mean"],
                v["properties"]["delta_time"],
            )
            for v in sample_gedi["features"]
        ],
        columns=["h", "water_occurrence", "delta_time"],
    )
    sample_alos_df = pd.DataFrame(
        [
            (v["properties"]["ALOS"], v["properties"]["occurrence_mean"])
            for v in samples_alos["features"]
        ],
        columns=["h", "water_occurrence"],
    )
    sample_nasadem_df = pd.DataFrame(
        [
            (v["properties"]["NASADEM"], v["properties"]["occurrence_mean"])
            for v in samples_srtm["features"]
        ],
        columns=["h", "water_occurrence"],
    )
    sample_glo_df = pd.DataFrame(
        [
            (v["properties"]["GLO"], v["properties"]["occurrence_mean"])
            for v in samples_glo["features"]
        ],
        columns=["h", "water_occurrence"],
    )

    coordsgedi = np.array(
        [v["geometry"]["coordinates"] for v in sample_gedi["features"]]
    )
    sample_gedi = gpd.GeoDataFrame(
        sample_gedi_df, geometry=gpd.points_from_xy(coordsgedi[:, 0], coordsgedi[:, 1])
    )

    coordsalos = np.array(
        [v["geometry"]["coordinates"] for v in samples_alos["features"]]
    )
    sample_alos = gpd.GeoDataFrame(
        sample_alos_df, geometry=gpd.points_from_xy(coordsalos[:, 0], coordsalos[:, 1])
    )

    coordsnasadem = np.array(
        [v["geometry"]["coordinates"] for v in samples_srtm["features"]]
    )
    sample_nasadem = gpd.GeoDataFrame(
        sample_nasadem_df,
        geometry=gpd.points_from_xy(coordsnasadem[:, 0], coordsnasadem[:, 1]),
    )

    coordsglo = np.array(
        [v["geometry"]["coordinates"] for v in samples_glo["features"]]
    )
    sample_glo = gpd.GeoDataFrame(
        sample_glo_df, geometry=gpd.points_from_xy(coordsglo[:, 0], coordsglo[:, 1])
    )

    return sample_gedi, sample_alos, sample_nasadem, sample_glo


def collect_reservoir_points(
    reservoir_id: str, dates: List[str]
) -> ee.FeatureCollection:
    """
    Collects a featurecollection of all icesatdata for the requested reservoir id and dates.

    input:
        reservoir_id (str): reservoir id
        dates (List): list of dates in YYYY-MM-dd format

    returns (ee.FeatureCollection) containing icesat2 data for reservoir and dates
    """
    pass


def get_watermasks_icesat2_remote(
    clean_images: ee.ImageCollection,
    icesat_fc: ee.FeatureCollection,
    waterbody: ee.Feature,
    max_time_window: int = 10,
) -> ee.FeatureCollection:
    """
    Function that loops in a icesatdataframe and per unique date extracts the matching water mask.

    :param clean_images (ee.ImageCollection): clean images from get_waterbody_images()
    :param icesat_fc (ee.FeatureCollection): icesat_data for date and reservoir
    :param waterbody: ee.Feature, waterbody polygon
    :param max_time_window: int, max distance to check for masks in days [centered window]
    :param subsample_is2_frac: float, 0-1 fraction for subsampling icesat2 tracks (reduces comp effort)

    Two dataframes are returned:
    :output df_is2: icesat2 dataframe subsampled with values for: water_mask, water_occurrence and distance from water_mask border
    :output water_level: dataframe with unique dates, mask area (filled) and icesat2 median height for pixels labelled as water (and water_filled)
    """

    # What about a pass around midnight?
    dates_fc: ee.FeatureCollection = icesat_fc.distinct("system:time_start")
    dates_fc = dates_fc.map(
        lambda f: ee.Feature(
            None, {"system:time_start": ee.Feature(f).get("system:time_start")}
        )
    )

    water_occurrence = (
        ee.Image("JRC/GSW1_3/GlobalSurfaceWater")
        .select("occurrence")
        .unmask(0)
        .resample("bicubic")
        .divide(100)
    )

    # Water occurrence S2 DW
    waterb_ee = ee.Geometry(waterbody.geometry()).buffer(300)

    start, stop = "2018-01-02", "2021-12-03"

    COL_FILTER = ee.Filter.And(ee.Filter.bounds(waterb_ee), ee.Filter.date(start, stop))

    # Dynamic world water occurence:
    dwCol_m_c = (
        ee.ImageCollection("GOOGLE/DYNAMICWORLD/V1")
        .filter(COL_FILTER)
        .select("water")
        .mean()
        .unmask(0)
        .rename("DW_wp")
    )

    def _get_icesat_statistics(date_f, clean_images, waterbody, max_time_window):
        date: ee.Number = ee.Number(ee.Feature(date_f).get("system:time_start"))
        water_area: ee.ImageCollection = estimate_watermask_icesatrack(
            clean_images, date, waterbody, max_time_window
        )
        is2_points: ee.FeatureCollection = icesat_fc.filter(
            ee.Filter.eq("system:time_start", date)
        )

        def get_and_sample_water_area(
            water_area: ee.ImageCollection, is_points: ee.FeatureCollection
        ):
            water_mask = (
                ee.Image(water_area.first())
                .select("combined_water_mask")
                .rename("water")
            )  ### NOTE: consider refining the date matching, now we only select closest, include perhaps also other filters e.g. quality

            # add distance to waterbody shoreline as additional variable

            distance_int = (
                water_mask.gt(0)
                .Not()
                .fastDistanceTransform(256)
                .sqrt()
                .reproject(ee.Projection("EPSG:3857").atScale(10))
                .multiply(10)
                .rename("dist_m")
            )
            distance_ext = (
                water_mask.gt(0)
                .fastDistanceTransform(256)
                .sqrt()
                .reproject(ee.Projection("EPSG:3857").atScale(10))
                .multiply(-10)
                .rename("dist_m")
            )

            distance = distance_int.add(distance_ext)

            # add distance and water occurence to mask
            water_mask = (
                water_mask.addBands(distance)
                .addBands(water_occurrence)
                .addBands(dwCol_m_c)
            )

            return (
                water_mask.unmask(-999, False)
                .reduceRegions(
                    collection=is_points,
                    reducer=ee.Reducer.first(),  # first reducer will take first image (only image) at point
                    scale=10,
                    tileScale=1,
                )
                .map(
                    lambda f: f.set("water_area_filled", water_mask.get("area_filled"))
                )
                .map(
                    lambda f: f.set(
                        "water_mask_time", water_mask.get("system:time_start")
                    )
                )
            )

        # select closest water mask in time
        # return empty FeatureCollection if no image was found
        output_fc: ee.FeatureCollection = ee.Algorithms.If(
            water_area.size().lt(1),
            ee.FeatureCollection([]),
            get_and_sample_water_area(water_area, is2_points),
        )

        return output_fc

    mappable: Callable = partial(
        _get_icesat_statistics,
        clean_images=clean_images,
        waterbody=waterbody,
        max_time_window=max_time_window,
    )
    return dates_fc.map(mappable).flatten()
