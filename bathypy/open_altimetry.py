from datetime import date
from typing import List, Tuple

from dateutil.parser import parse
import requests
from retry import retry


def get_overpasses(
    track_ids: List[int], start: date = date.min, stop: date = date.max
) -> List[Tuple[str, int]]:
    """Retrieves all overpasses for the given track ids

    :param track_ids: list of icesat2 track ids
    :type track_ids: list
    :return: A list of lists containing the icesat2 track ids and the
    corresponding dates
    :rtype: list
    """

    trackdates_url = (
        "https://openaltimetry.earthdatacloud.nasa.gov/data/icesat2/getTrackDate"
    )
    dates = api_get_call(trackdates_url)
    overpasses = []
    for track_id in track_ids:
        for track_date in dates["track_{}".format(track_id)].split(","):
            if start < parse(track_date).date() < stop:
                overpasses.append((track_date, track_id))
    return overpasses


def get_track_ids(coordinates):
    """Retrieves icesat2 track ids for a given bounding box"

    :param coordinates: tuple of min and max xy coordinates of a bounding box
    :type coordinates: tuple of floats
    :return: list of track ids
    :rtype: list
    """
    minx, miny, maxx, maxy = coordinates
    tracks_url = (
        "https://openaltimetry.earthdatacloud.nasa.gov/data/api/icesat2/getTracks"
        f"?minx={minx}"
        f"&miny={miny}"
        f"&maxx={maxx}"
        f"&maxy={maxy}"
        "&outputFormat=json"
    )
    data = api_get_call(tracks_url)
    return data["output"]["track"]


def get_icesat2_data_from_OA_api(
    icesat2_product: str,
    boundingbox: list,
    track_date: str,
    track_id: str,
    series_included=["Medium", "High"],
) -> dict:
    """Gets icesat2 data from the Open Altimetry API for a given bounding box, track
    date, and track id. The type of icesat2 data can be specified.

    :param icesat2_product: Which icesat2 product to get, i.e. 'atl03' or 'atl08'.
    :type icesat2_product: str
    :param boundingbox: Bounding box in the format of [min x , min y, max x, max y]
        in EPSG 4326 coordinates.
    :type boundingbox: list
    :param track_date: The date of the icesat2 track.
    :type track_date: str
    :param track_id: The id of the icesat2 track
    :type track_id: str
    :param series_included: The confidence level of the quality of the icesat2 data.
    :type series_included: list, optional
    :return: The icesat2 data as nested dict.
    :rtype: dict
    """
    series_included = "&".join(
        ["photonConfidence=" + series.lower() for series in series_included]
    )

    minx, miny, maxx, maxy = boundingbox
    OA_API_URL = (
        "https://openaltimetry.earthdatacloud.nasa.gov/data/api/icesat2/{}?"
        "&minx={}&miny={}&maxx={}&maxy={}&date={}&trackId={}&{}"
        "&beamName=gt3r&beamName=gt3l&beamName=gt2r&beamName=gt2l&beamName=gt1r&beamName=gt1l".format(
            icesat2_product,
            minx,
            miny,
            maxx,
            maxy,
            track_date,
            track_id,
            series_included,
        )
    )
    return api_get_call(OA_API_URL)


@retry(delay=5, tries=10)
def api_get_call(url):
    """Makes an api get request with retries.

    :param url: url string for making a get request
    :type url: str
    :return: request response as JSON
    :rtype: dict
    """
    with requests.get(url) as response:
        return response.json()
