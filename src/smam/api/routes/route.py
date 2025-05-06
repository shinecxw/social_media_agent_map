from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import Dict, Any
from smam.services.planner import get_amap_coordinates, gcj02_to_wgs84

router = APIRouter()

class PointRequest(BaseModel):
    name: str

@router.post("/point/")
async def get_point_coordinates(request: PointRequest):
    """获取地点的坐标信息并转换为WGS84坐标系"""
    try:
        # 获取高德地图坐标
        gcj_coords = get_amap_coordinates(request.name)
        if not gcj_coords:
            raise HTTPException(status_code=404, detail=f"找不到地点: {request.name}")
        
        # 转换为WGS84坐标
        wgs_lon, wgs_lat = gcj02_to_wgs84(*gcj_coords)
        
        return {
            "name": request.name,
            "gcj02": {"lon": gcj_coords[0], "lat": gcj_coords[1]},
            "wgs84": {"lon": wgs_lon, "lat": wgs_lat}
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))