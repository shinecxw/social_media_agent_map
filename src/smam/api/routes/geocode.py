from fastapi import APIRouter, HTTPException, Query, Depends
from typing import Optional

from ..models.geocode import GeocodeRequest, GeocodeResponse
from ...services.geocode_service import GeocodeService, get_geocode_service

router = APIRouter(prefix='/geocode', tags=['geocode'])


@router.get("/", response_model=GeocodeResponse)
async def geocode_address(
    address: str = Query(..., description="地址或地点名称"),
    city: Optional[str] = Query(None, description="城市名称"),
    country: Optional[str] = Query(None, description="国家名称"),
    geocode_service: GeocodeService = Depends(get_geocode_service)
):
    """
    根据地址进行地理编码，返回经纬度坐标
    """
    try:
        result = await geocode_service.geocode_address(address, city, country)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"地理编码失败: {str(e)}")


@router.post("/", response_model=GeocodeResponse)
async def geocode_address_post(
    request: GeocodeRequest,
    geocode_service: GeocodeService = Depends(get_geocode_service)
):
    """
    根据地址进行地理编码，返回经纬度坐标 (POST方法)
    """
    try:
        result = await geocode_service.geocode_address(
            request.address, request.city, request.country
        )
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"地理编码失败: {str(e)}")


@router.get("/reverse", response_model=GeocodeResponse)
async def reverse_geocode(
    latitude: float = Query(..., description="纬度"),
    longitude: float = Query(..., description="经度"),
    geocode_service: GeocodeService = Depends(get_geocode_service)
):
    """
    根据经纬度坐标进行反向地理编码，返回地址信息
    """
    try:
        result = await geocode_service.reverse_geocode(latitude, longitude)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"反向地理编码失败: {str(e)}")