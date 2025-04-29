from fastapi import APIRouter, HTTPException, Query, Depends
from typing import Optional

from ..schemas.geocode import QueryAddressName, QueryAddressGeoCoords, PROVIDER_AMAP, PROVIDER_GOOGLE
# from ...services.geocode_service import GeocodeService, get_geocode_service

router = APIRouter(prefix='/geocode', tags=['geocode'])

@router.post('/', response_model=QueryAddressGeoCoords)
def geocode_address(address: QueryAddressName):

    # querey coordinates from amap API
    if address.provider == PROVIDER_AMAP:
        # Mock query
        # lon, lat = getCoordFromName(address.name)
        # Transform coordinates from cgcs2000 to wgs84
        # lon, lat = transform(target: wgs84, source: cgcs2000, coords: (lon, lat))
        # return QueryAddressGeoCoords(coords=(lon, lat))
        return QueryAddressGeoCoords(coords=(0.0, 0.0))
    else:
        # lon, lat = getCoordFromName(address.name)
        return QueryAddressGeoCoords(coords=(0.0, 0.0))