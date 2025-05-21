import sys
import os
from pathlib import Path

# 修正 BASE_DIR 为项目根目录
BASE_DIR = Path(__file__).resolve().parent.parent  # 直接到 D:\git\smam\social_media_agent_map
sys.path.append(str(BASE_DIR))

from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles
from smam.api.routes.route import router as route_router

app = FastAPI()

app.include_router(route_router, prefix="/api/route")

# 修正静态文件路径
static_path = os.path.join(BASE_DIR, "static")
print(f"Computed static path: {static_path}")  # 调试输出
if not os.path.exists(static_path):
    print(f"Warning: Static directory '{static_path}' does not exist.")
else:
    app.mount("/", StaticFiles(directory=static_path, html=True), name="static")

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="127.0.0.1", port=8000)