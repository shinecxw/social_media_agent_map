# social_media_agent_map

An agent enhanced map application used for navigation based on social media multi-modal data.

## 项目简介

本项目是一个基于社交媒体多模态数据的智能地图导航应用。通过整合社交媒体数据和智能代理技术，为用户提供更丰富、更智能的导航体验。

### 环境要求

- Python 3.12
- 相关依赖包 (见 requirements.txt)

### 安装步骤


### 使用方法
# 激活虚拟环境（如 .venv）
. .venv/Scripts/activate

# 设置模块搜索路径
$env:PYTHONPATH="src"

# 启动服务
python -m uvicorn smam.main:app --reload
