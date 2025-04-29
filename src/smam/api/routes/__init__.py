from fastapi import APIRouter

from . import geocode

api_router = APIRouter()

api_router.include_router(geocode.router)