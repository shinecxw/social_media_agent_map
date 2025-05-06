import os
import math
import requests
import networkx as nx
import geopandas as gpd
from shapely.geometry import Point
from shapely.geometry import LineString, MultiLineString

AMAP_KEY = "b3ad8416f6f2f251a3a53d3b10f93833"  # 替换为你自己的高德 KEY

# ========= 坐标系转换 ========= #
PI = math.pi
A = 6378245.0
EE = 0.00669342162296594323

def out_of_china(lng, lat):
    return not (73.66 < lng < 135.05 and 3.86 < lat < 53.55)

def transform_lat(x, y):
    ret = -100.0 + 2.0 * x + 3.0 * y + 0.2 * y*y + \
        0.1 * x * y + 0.2 * math.sqrt(abs(x))
    ret += (20.0 * math.sin(6.0 * x * PI) +
            20.0 * math.sin(2.0 * x * PI)) * 2.0 / 3.0
    ret += (20.0 * math.sin(y * PI) +
            40.0 * math.sin(y / 3.0 * PI)) * 2.0 / 3.0
    ret += (160.0 * math.sin(y / 12.0 * PI) +
            320 * math.sin(y * PI / 30.0)) * 2.0 / 3.0
    return ret

def transform_lng(x, y):
    ret = 300.0 + x + 2.0 * y + 0.1 * x * x + \
        0.1 * x * y + 0.1 * math.sqrt(abs(x))
    ret += (20.0 * math.sin(6.0 * x * PI) +
            20.0 * math.sin(2.0 * x * PI)) * 2.0 / 3.0
    ret += (20.0 * math.sin(x * PI) +
            40.0 * math.sin(x / 3.0 * PI)) * 2.0 / 3.0
    ret += (150.0 * math.sin(x / 12.0 * PI) +
            300.0 * math.sin(x / 30.0 * PI)) * 2.0 / 3.0
    return ret

def gcj02_to_wgs84(lng, lat):
    if out_of_china(lng, lat):
        return lng, lat
    dlat = transform_lat(lng - 105.0, lat - 35.0)
    dlng = transform_lng(lng - 105.0, lat - 35.0)
    radlat = lat / 180.0 * PI
    magic = math.sin(radlat)
    magic = 1 - EE * magic * magic
    sqrtmagic = math.sqrt(magic)
    dlat = (dlat * 180.0) / ((A * (1 - EE)) / (magic * sqrtmagic) * PI)
    dlng = (dlng * 180.0) / (A / sqrtmagic * math.cos(radlat) * PI)
    mglat = lat + dlat
    mglng = lng + dlng
    return lng * 2 - mglng, lat * 2 - mglat

# ========= 坐标查询 ========= #
def get_amap_coordinates(name):
    url = "https://restapi.amap.com/v3/place/text"
    params = {
        "key": AMAP_KEY,
        "keywords": name,
        "offset": 1,
        "page": 1
    }
    resp = requests.get(url, params=params).json()
    if resp.get("status") == "1" and resp.get("pois"):
        lng, lat = map(float, resp["pois"][0]["location"].split(","))
        print(f"📍 查询地址：{name} → 高德坐标：({lng}, {lat})")
        return lng, lat
    print(f"❌ 查询地址失败：{name}，响应：{resp}")
    return None

# ========= 图构建 ========= #
def build_graph(shapefile_path):
    gdf = gpd.read_file(shapefile_path)
    G = nx.Graph()
    for _, row in gdf.iterrows():
        geom = row.geometry
        if isinstance(geom, LineString):
            lines = [geom]
        elif isinstance(geom, MultiLineString):
            lines = list(geom.geoms)
        else:
            continue  # 忽略非线类型

        for line in lines:
            coords = list(line.coords)
            for i in range(len(coords) - 1):
                p1, p2 = coords[i], coords[i+1]
                dist = Point(p1).distance(Point(p2))
                G.add_edge(p1, p2, weight=dist)
    return G

def nearest_node(G, point):
    pt = Point(point)
    return min(G.nodes, key=lambda n: Point(n).distance(pt))

def path_to_geojson(coords):
    return {
        "type": "FeatureCollection",
        "features": [{
            "type": "Feature",
            "geometry": {
                "type": "LineString",
                "coordinates": coords
            },
            "properties": {}
        }]
    }

# ========= 单段路径推进 ========= #
def advance_segment(G, current_node, instructions):
    target_distance = 100  # 默认距离
    for step in instructions:
        if "distance" in step:
            try:
                target_distance = float(step["distance"].replace("m", "").strip())
            except:
                pass
        if "time" in step:
            try:
                minutes = int(step["time"].replace("分钟", "").strip())
                target_distance = minutes * 60  # 步行速度假设 1m/s
            except:
                pass

    total = 0
    path = [current_node]
    while total < target_distance:
        neighbors = list(G.neighbors(path[-1]))
        if not neighbors:
            break
        next_node = min(neighbors, key=lambda n: Point(n).distance(Point(path[-1])))
        dist = Point(path[-1]).distance(Point(next_node))
        if next_node in path or dist == 0:
            break
        path.append(next_node)
        total += dist
    return path

# ========= 核心路径模拟主函数 ========= #
def compute_full_path(data, shapefile_path):
    G = build_graph(shapefile_path)
    start_name = data["start"]["name"]
    start_coord_gcj = get_amap_coordinates(start_name)
    if not start_coord_gcj:
        raise Exception("无法获取起点坐标")

    start_wgs = gcj02_to_wgs84(*start_coord_gcj)
    current_node = nearest_node(G, start_wgs)

    full_coords = []
    instructions = []

    for seg_id, seg in data["segments"].items():
        path_instr = seg.get("path_instructions", [])
        segment_path = advance_segment(G, current_node, path_instr)
        full_coords.extend(segment_path if not full_coords else segment_path[1:])
        current_node = segment_path[-1]
        text = " → ".join([" ".join(v for v in step.values()) for step in path_instr])
        instructions.append({
            "segment": seg_id,
            "description": text,
            "point_count": len(segment_path)
        })

    return {
        "geojson": path_to_geojson(full_coords),
        "instructions": instructions
    }
