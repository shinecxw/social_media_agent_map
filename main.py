import uvicorn
from src.smam.main import app

if __name__ == "__main__":
    uvicorn.run("src.smam.main:app", host="0.0.0.0", port=8000, reload=True)