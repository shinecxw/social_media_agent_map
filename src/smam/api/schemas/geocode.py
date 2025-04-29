from pydantic import BaseModel, Field, field_validator
from typing import Optional, List, Dict, Any

PROVIDER_AMAP = 0
PROVIDER_GOOGLE = 1

class QueryAddressName(BaseModel):
    """
    The queried address name
    """
    name: str
    provider: int = Field(default=PROVIDER_AMAP, description="Geocoding provider: 0 for Amap, 1 for Google")
    epsg: int = Field(default=4326, description="Coordinate system EPSG code")

    @field_validator('name')
    def validate_name(cls, value):
        if not value:
            raise ValueError('Name must not be empty')
        return value
    
    @field_validator('provider')
    def validate_provider(cls, value):
        if value not in [PROVIDER_AMAP, PROVIDER_GOOGLE]:
            raise ValueError('Provider must be either PROVIDER_AMAP or PROVIDER_GOOGLE')
        return value

class QueryAddressGeoCoords(BaseModel):
    """The response data of query address"""

    coords: tuple[float, float] # (longitude, latitude)

    @field_validator('coords')
    def validate_coords(cls, value):
        if len(value) != 2:
            raise ValueError('Coordinates must be a tuple of (longitude, latitude)')
        return value


# class GeocodeRequest(BaseModel):
#     """地理编码请求模型"""
#     address: str = Field(..., description="地址或地点名称")
#     city: Optional[str] = Field(None, description="城市名称")
#     country: Optional[str] = Field(None, description="国家名称")


# class GeocodeLocation(BaseModel):
#     """地理位置坐标模型"""
#     latitude: float = Field(..., description="纬度")
#     longitude: float = Field(..., description="经度")


# class GeocodeResponse(BaseModel):
#     """地理编码响应模型"""
#     query: str = Field(..., description="原始查询内容")
#     formatted_address: str = Field(..., description="格式化的地址")
#     location: GeocodeLocation = Field(..., description="位置坐标")
#     components: Dict[str, Any] = Field(default_factory=dict, description="地址组成部分")
#     confidence: float = Field(..., description="结果置信度", ge=0, le=1)
#     additional_info: Optional[Dict[str, Any]] = Field(None, description="附加信息")