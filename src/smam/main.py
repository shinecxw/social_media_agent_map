from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles
from smam.api.routes.route import router as route_router  # ✅ 路由必须注册
import os

app = FastAPI()

app.include_router(route_router, prefix="/api/route")  # ✅ 正确注册路由前缀

# 静态文件
static_path = os.path.join(os.path.dirname(__file__), "../static")
app.mount("/", StaticFiles(directory=static_path, html=True), name="static")
