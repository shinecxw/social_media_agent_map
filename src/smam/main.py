from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from .api.routes import api_router

app = FastAPI(
    title="Social Media Agent Map",
    description="An agent enhanced map application used for navigation based on social media multi-modal data.",
    version="0.1.0",
)

# CORE setting
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(api_router)