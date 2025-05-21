import math
from shapely.geometry import Point, LineString,MultiLineString
import requests
import networkx as nx
import geopandas as gpd
import time
AMAP_KEY = "465b45001f158d735671fff7d1ee10cd"  # 替换为你自己的高德 KEY

# ========= 坐标系转换 ========= #
PI = math.pi
A = 6378245.0
EE = 0.00669342162296594323
# ========= 方向角度映射 ========= #
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
            300 * math.sin(x / 30.0 * PI)) * 2.0 / 3.0
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
    return None  #118.794787,32.041737
# 方向映射表
direction_mapping = {
    "east": 90, "e": 90, "东": 90, "向东": 90,
    "northeast": 45, "ne": 45, "东北": 45,
    "north": 0, "n": 0, "北": 0, "向北": 0,
    "northwest": 315, "nw": 315, "西北": 315,
    "west": 270, "w": 270, "西": 270, "向西": 270,
    "southwest": 225, "sw": 225, "西南": 225,
    "south": 180, "s": 180, "南": 180, "向南": 180,
    "southeast": 135, "se": 135, "东南": 135
}
def nearest_node(G, point):
    pt = Point(point)
    return min(G.nodes, key=lambda n: Point(n).distance(pt))
def path_to_geojson(path):
    coordinates = [[p[0], p[1]] for p in path]  # 包括所有路径点
    return {
        "type": "FeatureCollection",
        "features": [{
            "type": "Feature",
            "geometry": {
                "type": "LineString",
                "coordinates": coordinates  # 确保返回完整路径的所有点
            },
            "properties": {}
        }]
    }
# 相对方向映射表
relative_direction_mapping = {
    "forward": "straight",
    "straight": "straight",
    "back": "back",
    "left_forward": -45,    # 左前
    "right_forward": 45,    # 右前
    "left_back": -135,      # 左后
    "right_back": 135,      # 右后
    "go_straight": "straight",
    "to_end": "end",
    "left": -90,           # 左转
    "right": 90            # 右转
}
# Haversine 公式计算两点间距离（单位：米）
def haversine(lon1, lat1, lon2, lat2):
    R = 6371000  # 地球半径（米）
    lon1, lat1, lon2, lat2 = map(math.radians, [lon1, lat1, lon2, lat2])
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = math.sin(dlat/2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon/2)**2
    c = 2 * math.asin(math.sqrt(a))
    return R * c


def parse_time_to_minutes(time_str):
    """将时间字符串（如 '5min' 或 '1.5hour'）转换为分钟数"""
    import re
    if not time_str or not isinstance(time_str, str):
        return None
    time_match = re.match(r'(\d+(?:\.\d+)?)(min|hour)?', time_str.lower())
    if not time_match:
        print(f"      警告: 无效时间格式 {time_str}")
        return None
    time_value = float(time_match.group(1))
    time_unit = time_match.group(2) or 'min'
    if time_unit == 'hour':
        time_value *= 60  # 转换为分钟
    return time_value


def process_relative_direction(current_direction, rel_direction):
    """处理相对方向"""
    if not current_direction or not rel_direction:
        return None
    rel_direction = rel_direction.lower().replace(" ", "_")
    if rel_direction in relative_direction_mapping:
        if isinstance(relative_direction_mapping[rel_direction], (int, float)):
            return (current_direction + relative_direction_mapping[rel_direction]) % 360
        elif relative_direction_mapping[rel_direction] == "straight":
            return current_direction
        elif relative_direction_mapping[rel_direction] == "back":
            return (current_direction + 180) % 360
        elif relative_direction_mapping[rel_direction] == "end":
            return "end"
    return None
def calculate_direction(point1, point2):
    """计算从 point1 到 point2 的方向角度（北向为 0 度，顺时针增加）"""
    if not (isinstance(point1, tuple) and isinstance(point2, tuple) and len(point1) >= 2 and len(point2) >= 2):
        print(f"      警告: 点格式不正确 - {point1}, {point2}")
        return 0
    
    dx = point2[0] - point1[0]
    dy = point2[1] - point1[1]
    
    # 防止极小值导致的计算误差
    if abs(dx) < 1e-10 and abs(dy) < 1e-10:
        print(f"      警告: 两点极度接近 - {point1}, {point2}")
        return 0
        
    angle = (90 - math.degrees(math.atan2(dy, dx))) % 360
    #print(f"      计算方向: 从 {point1} 到 {point2}，dx={dx:.6f}, dy={dy:.6f}，角度={angle:.2f}°")
    return angle
def find_nearest_road_point(G, point):
    """找到距离某点最近的路网节点"""
    min_dist = float('inf')
    nearest_point = None
    for node in G.nodes():
        dist = math.hypot(node[0] - point[0], node[1] - point[1])
        if dist < min_dist:
            min_dist = dist
            nearest_point = node
    return nearest_point


def parse_distance(step, transportation="walk"):
    """解析距离字符串，优先使用 distance 字段，若无则从 time 字段计算，若无则返回默认 100 米"""
    import re
    print(f"      调试: 进入 parse_distance, step = {step}, transportation = {transportation}")
 # 如果明确只有交通方式，没有距离/时间/方向/动作，则返回 None
    if transportation and not step.get("distance") and not step.get("direction") and not step.get("time") and not step.get("action"):
        return None
    
    # 首先尝试解析 distance
    distance = step.get("distance")
    print(f"      调试: distance = {distance}")
    if distance:
        match = re.match(r'(\d+(?:\.\d+)?)(m|km)?', distance.lower())
        if match:
            value = float(match.group(1))
            unit = match.group(2) or 'm'
            if unit == 'km':
                value *= 1000  # 转换为米
            print(f"      调试: 从 distance 解析到 {value} 米")
            return value
        print(f"      调试: distance 格式无效，返回默认 100 米")
        return 100.0

    # 如果没有 distance，尝试解析 time
    time = step.get("time")
    print(f"      调试: time = {time}")
    if time:
        time_value = parse_time_to_minutes(time)
        print(f"      调试: parse_time_to_minutes 返回 {time_value}")
        if time_value is not None:
            speed = get_speed_by_transportation(transportation)
            print(f"      调试: 速度 = {speed} 米/分钟")
            calculated_distance = time_value * speed
            print(f"      调试: 计算距离 = {calculated_distance:.2f} 米")
            return calculated_distance
        print(f"      调试: time 解析失败，返回默认 100 米")
        return 100.0

    # 如果没有 distance 和 time，返回默认 100 米
    print(f"      调试: 无 distance 和 time，返回默认 100 米")
    return 100.0
def interpolate_point(G, start, end, ratio):
    """在起点和终点之间按比例插值计算一个点，并添加到图中"""
    start_x, start_y = start
    end_x, end_y = end
    x = start_x + (end_x - start_x) * ratio
    y = start_y + (end_y - start_y) * ratio
    new_point = (x, y)
    
    original_weight = G[start][end]['weight'] if G.has_edge(start, end) else haversine(start[0], start[1], end[0], end[1])
    weight1 = original_weight * ratio
    weight2 = original_weight * (1 - ratio)
    
    G.add_edge(start, new_point, weight=weight1)
    G.add_edge(new_point, end, weight=weight2)
    
    if G.has_edge(start, end):
        G.remove_edge(start, end)
    
    return new_point
def get_direction_from_angle(angle):
    directions = ['north', 'northeast', 'east', 'southeast', 
                 'south', 'southwest', 'west', 'northwest']
    
    angle = angle % 360
    index = round(angle / 45) % 8
    return directions[index]
def filter_by_direction(current, neighbors, direction, last_direction=None, angle_tolerance=30):
    """按方向过滤邻居节点，优先选择同一道路"""
    if not direction or not neighbors:
        return neighbors

    target_angle = None
    # Handle the case when direction is already a numeric value (float/int)
    if isinstance(direction, (int, float)):
        target_angle = direction
    else:
        # Process string direction
        direction = direction.lower()
        if direction in direction_mapping:
            target_angle = direction_mapping[direction]
        elif direction in relative_direction_mapping:
            if isinstance(relative_direction_mapping[direction], (int, float)):
                target_angle = (last_direction + relative_direction_mapping[direction]) % 360 if last_direction is not None else None
            elif relative_direction_mapping[direction] == "straight":
                target_angle = last_direction
            elif relative_direction_mapping[direction] == "back":
                target_angle = (last_direction + 180) % 360 if last_direction is not None else None
    
    if target_angle is None:
        return neighbors

    def direction_score(nbr):
        dx = nbr[0] - current[0]
        dy = nbr[1] - current[1]
        angle = math.degrees(math.atan2(dy, dx)) % 360
        angle = (90 - angle) % 360
        angle_diff = min((angle - target_angle) % 360, (target_angle - angle) % 360)
        return angle_diff

    filtered = []
    for nbr in neighbors:
        angle_diff = direction_score(nbr)
        if angle_diff <= angle_tolerance:
            filtered.append((nbr, angle_diff))

    if not filtered:
        print(f"        无方向匹配，尝试放宽容差至 {angle_tolerance * 1.5}°")
        for nbr in neighbors:
            angle_diff = direction_score(nbr)
            if angle_diff <= angle_tolerance * 1.5:
                filtered.append((nbr, angle_diff))

    filtered.sort(key=lambda x: x[1])
    filtered_neighbors = [item[0] for item in filtered]
    print(f"        按方向 {direction} (角度={target_angle}°) 过滤后有 {len(filtered_neighbors)} 个邻居，偏差最小的方向: {filtered[0][1] if filtered else '无'}°")
    return filtered_neighbors
def get_speed_by_transportation(transport):
    """返回交通方式的速度（单位：米/分钟）"""
    speeds = {
        "walk": 84,    # 1.4 米/秒 * 60 = 84 米/分钟
        "bike": 240,   # 4 米/秒 * 60 = 240 米/分钟
        "car": 600     # 10 米/秒 * 60 = 600 米/分钟
    }
    speed = speeds.get(transport.lower(), 84)
    print(f"      调试: get_speed_by_transportation({transport}) 返回 {speed}")
    return speed


def find_nearest_exit(G, current_node, current_direction, max_distance=1000, roads_gdf=None, current_road_id=None, transport="walk"):
    """找到最近的交叉路口或出口"""
    path = [current_node]
    total_distance = 0
    visited = {current_node}

    print(f"      寻找最近出口，当前道路ID: {current_road_id}, 最大距离: {max_distance}米")

    while total_distance < max_distance:
        neighbors = list(G.neighbors(current_node))
        for neighbor in neighbors:
            neighbor_road_id = G.nodes[neighbor].get("road_id")
            if neighbor_road_id != current_road_id and neighbor not in visited:
                path.append(neighbor)
                total_distance += G[current_node][neighbor]['weight']
                new_direction = calculate_direction(current_node, neighbor)
                print(f"      找到出口: {neighbor}, 距离: {total_distance:.2f}米, 方向: {new_direction:.2f}°")
                return path, total_distance, new_direction

        if not neighbors:
            print(f"      无更多邻居，停止搜索")
            break

        next_node = min(neighbors, key=lambda n: G[current_node][n]['weight'])
        edge_weight = G[current_node][next_node]['weight']
        if total_distance + edge_weight > max_distance:
            break
        path.append(next_node)
        total_distance += edge_weight
        visited.add(next_node)
        current_node = next_node
        current_direction = calculate_direction(path[-2], path[-1])

    print(f"      达到最大距离 {max_distance}米，未找到出口")
    return [current_node], 0.0, current_direction
def calculate_direction(point1, point2):
    if not (isinstance(point1, tuple) and isinstance(point2, tuple)):
        return 0
    dx = point2[0] - point1[0]
    dy = point2[1] - point1[1]
    angle = (90 - math.degrees(math.atan2(dy, dx))) % 360
    #print(f"      计算方向: 从 {point1} 到 {point2}，dx={dx:.6f}, dy={dy:.6f}，角度={angle:.2f}°")
    return angle

def explore_all_paths_no_backtrack(G, start_node, target_distance, transport="walk"):
    all_paths = []
    visited_edges = set()  # 记录已访问的边 (u, v)
    stack = [(start_node, [start_node], 0.0, None)]  # (节点, 路径, 总距离, 方向)

    print(f"      探索所有满足距离 {target_distance}米的路径，起点: {start_node}")

    while stack:
        current_node, path, total_distance, prev_direction = stack.pop()

        if abs(total_distance - target_distance) < 1e-5:  # 近似匹配目标距离
            current_direction = calculate_direction(path[-2], path[-1]) if len(path) > 1 else prev_direction
            all_paths.append((path, total_distance, current_direction))
            print(f"      找到路径，终点: {path[-1]}，距离: {total_distance:.2f}米，方向: {current_direction:.2f}°")
            continue

        if total_distance > target_distance * 1.5:  # 防止过长路径
            continue

        neighbors = list(G.neighbors(current_node))
        for neighbor in neighbors:
            edge = tuple(sorted([current_node, neighbor]))  # 确保边唯一
            if edge not in visited_edges:
                edge_weight = G[current_node][neighbor]['weight']
                new_distance = total_distance + edge_weight
                new_path = path + [neighbor]
                new_direction = calculate_direction(path[-1], neighbor) if len(path) >= 1 else None

                if new_distance <= target_distance:
                    stack.append((neighbor, new_path, new_distance, new_direction))
                    visited_edges.add(edge)
                elif new_distance > target_distance:
                    ratio = (target_distance - total_distance) / edge_weight
                    interpolated_point = interpolate_point(G, current_node, neighbor, ratio)
                    new_path[-1:] = [interpolated_point]
                    new_direction = calculate_direction(new_path[-2], interpolated_point) if len(new_path) >= 2 else None
                    all_paths.append((new_path, target_distance, new_direction))
                    print(f"      截断路径，终点: {interpolated_point}，距离: {target_distance:.2f}米，方向: {new_direction:.2f}°")

        print(f"      当前节点 {current_node}，累计距离: {total_distance:.2f}米，路径长度: {len(path)}")

    print(f"      共生成 {len(all_paths)} 条路径")
    return all_paths if all_paths else [([start_node], 0.0, None)]



def find_intersection(G, current_node, current_road_id, max_distance=500):
    path = [current_node]
    total_distance = 0
    visited = {current_node}

    print(f"      开始寻找路口，当前道路ID: {current_road_id}, 最大距离: {max_distance}米")

    while total_distance < max_distance:
        neighbors = [
            n for n in G.neighbors(current_node)
            if n not in visited and G.nodes[n].get("road_id") == current_road_id
        ]

        # 检查当前节点是否为路口
        all_neighbors = list(G.neighbors(current_node))
        for neighbor in all_neighbors:
            neighbor_road_id = G.nodes[neighbor].get("road_id")
            if neighbor_road_id != current_road_id:
                print(f"      找到路口: {current_node}，连接道路 {current_road_id} 和 {neighbor_road_id}")
                return path, total_distance

        if not neighbors:
            print(f"      没有更多同一道路的邻居，停止搜索")
            return None, 0

        # 选择权重最小的邻居继续前进
        next_node = min(neighbors, key=lambda n: G[current_node][n]['weight'])
        edge_weight = G[current_node][next_node]['weight']
        total_distance += edge_weight
        path.append(next_node)
        visited.add(next_node)
        current_node = next_node

    print(f"      达到最大搜索距离 {max_distance}米，未找到路口")
    return None, 0


def find_path_end(G, current_node, direction, max_distance=5000, roads_gdf=None, current_road_id=None, transport="walk"):
    """
    寻找道路的最远端点，基于 shapefile 的几何信息生成完整路径，优先选择符合当前方向的端点。
    如果当前节点在道路端点，尝试移动到另一端，并通过路网图扩展到实际尽头。
    """
    path = [current_node]
    total_distance = 0
    current_direction = direction

    if current_road_id is None and current_node in G.nodes:
        current_road_id = G.nodes[current_node].get("road_id")
    if current_road_id is None:
        print(f"      错误: 当前节点 {current_node} 没有 road_id")
        return [current_node], 0, current_direction

    print(f"      寻找道路 {current_road_id} 的最远端点，当前方向: {current_direction:.2f}°，交通方式: {transport}")

    # 获取道路几何信息
    road_geom = None
    road_coords = []
    if roads_gdf is not None and current_road_id is not None:
        for idx, row in roads_gdf.iterrows():
            if int(row["OBJECTID"]) == current_road_id:
                road_geom = row.geometry
                if isinstance(road_geom, LineString):
                    road_coords = list(road_geom.coords)
                elif isinstance(road_geom, MultiLineString):
                    p = Point(current_node)
                    min_dist = float('inf')
                    for line in road_geom.geoms:
                        dist = line.distance(p)
                        if dist < min_dist:
                            min_dist = dist
                            road_coords = list(line.coords)
                    if not road_coords:
                        road_coords = list(road_geom.geoms[0].coords)
                print(f"      找到道路 {current_road_id}，坐标点数: {len(road_coords)}")
                print(f"      道路坐标: {road_coords}")
                break

    if not road_coords:
        print(f"      警告: 未找到道路 {current_road_id} 的几何信息，使用路网图扩展")
        return extend_path_to_end(G, current_node, current_direction, max_distance, current_road_id, transport)

    # 确定道路起点和终点
    road_start = tuple(road_coords[0])
    road_end = tuple(road_coords[-1])
    dist_to_start = haversine(current_node[0], current_node[1], road_start[0], road_start[1])
    dist_to_end = haversine(current_node[0], current_node[1], road_end[0], road_end[1])

    # 计算道路的整体方向（基于起点到终点）
    road_direction = calculate_direction(road_start, road_end) if dist_to_start > 1e-10 and dist_to_end > 1e-10 else current_direction
    print(f"      道路 {current_road_id} 整体方向: {road_direction:.2f}°")

    # 计算当前节点到起点的方向和偏差
    angle_to_start = calculate_direction(current_node, road_start) if dist_to_start > 1e-10 else road_direction
    angle_diff_start = min((angle_to_start - current_direction) % 360, (current_direction - angle_to_start) % 360)

    # 计算当前节点到终点的方向和偏差
    angle_to_end = calculate_direction(current_node, road_end) if dist_to_end > 1e-10 else road_direction
    angle_diff_end = min((angle_to_end - current_direction) % 360, (current_direction - angle_to_end) % 360)

    # 如果当前节点接近任一端点，优先选择另一端
    target_end = road_end
    if dist_to_end < 1.0:  # 接近终点，选择起点
        target_end = road_start
        print(f"      当前节点接近终点 {road_end}，选择起点 {road_start} 作为目标")
    elif dist_to_start < 1.0:  # 接近起点，选择终点
        target_end = road_end
        print(f"      当前节点接近起点 {road_start}，选择终点 {road_end} 作为目标")
    else:
        # 选择与当前方向偏差最小的端点
        if angle_diff_start < angle_diff_end or (angle_diff_start == angle_diff_end and dist_to_start < dist_to_end):
            target_end = road_start
        print(f"      选择目标端点: {target_end}, 方向偏差: {angle_diff_start if target_end == road_start else angle_diff_end:.2f}°")

    print(f"      道路起点: {road_start}, 终点: {road_end}")
    print(f"      距离起点: {dist_to_start:.2f}米，方向偏差: {angle_diff_start:.2f}°")
    print(f"      距离终点: {dist_to_end:.2f}米，方向偏差: {angle_diff_end:.2f}°")

    # 找到当前节点在道路坐标中的位置
    min_dist = float('inf')
    start_idx = 0
    for i, coord in enumerate(road_coords):
        dist = haversine(current_node[0], current_node[1], coord[0], coord[1])
        if dist < min_dist:
            min_dist = dist
            start_idx = i
    print(f"      当前节点最接近道路坐标索引: {start_idx}，坐标: {road_coords[start_idx]}")

    # 确定前进方向（基于目标端点）
    if target_end == road_start:
        direction_indices = range(start_idx, -1, -1)  # 向起点方向
    else:
        direction_indices = range(start_idx, len(road_coords))  # 向终点方向

    # 沿道路几何前进到目标端点
    current_node = tuple(current_node)  # 确保是 tuple
    for i in direction_indices:
        next_coord = tuple(road_coords[i])
        if next_coord == current_node:
            continue
        segment_dist = haversine(current_node[0], current_node[1], next_coord[0], next_coord[1])
        if total_distance + segment_dist > max_distance:
            print(f"      达到最大距离 {max_distance}米，停止搜索")
            break
        path.append(next_coord)
        total_distance += segment_dist
        if next_coord not in G.nodes:
            G.add_node(next_coord, road_id=current_road_id)
        G.add_edge(current_node, next_coord, weight=segment_dist)
        current_node = next_coord
        print(f"      移动到 {current_node}，累计距离: {total_distance:.2f}米")

    # 如果未到达道路尽头，使用路网图扩展
    if total_distance < max_distance:
        extended_path, ext_distance, new_direction = extend_path_to_end(
            G, current_node, current_direction, max_distance - total_distance, current_road_id, transport
        )
        if extended_path and len(extended_path) > 1:
            path.extend(extended_path[1:])
            total_distance += ext_distance
            current_node = extended_path[-1]
            current_direction = new_direction
            print(f"      从路网图扩展路径，新增距离: {ext_distance:.2f}米，新节点: {current_node}")

    # 更新方向
    if len(path) >= 2:
        current_direction = calculate_direction(path[-2], path[-1])
        print(f"      更新方向为 {current_direction:.2f}°")

    # 验证路径点格式
    for point in path:
        if not isinstance(point, tuple) or len(point) != 2:
            print(f"      错误: 路径点 {point} 不是有效的坐标元组")
            raise ValueError(f"无效路径点: {point}")

    print(f"      路径生成完成，路径长度: {len(path)}，总距离: {total_distance:.2f}米，最终方向: {current_direction:.2f}°")
    return path, total_distance, current_direction
def extend_path_to_end(G, current_node, direction, max_distance, current_road_id, transport):
    """从路网图扩展路径到道路尽头，不依赖方向限制"""
    path = [current_node]
    total_distance = 0
    current_direction = direction
    visited = {current_node}

    print(f"      开始从路网图扩展路径，当前方向: {current_direction:.2f}°，最大距离: {max_distance}米")

    while total_distance < max_distance:
        # 筛选同道路的邻居节点
        neighbors = [
            n for n in G.neighbors(current_node)
            if n not in visited and G.nodes[n].get("road_id") == current_road_id
        ]
        if not neighbors:
            print(f"      无更多同道路邻居，停止扩展")
            break

        # 选择任意一个邻居（无需方向筛选）
        next_node = neighbors[0]  # 可以优化为选择距离最近的邻居
        edge_weight = G[current_node][next_node]['weight']
        path.append(next_node)
        total_distance += edge_weight
        visited.add(next_node)
        current_node = next_node
        current_direction = calculate_direction(path[-2], path[-1])
        print(f"      扩展到 {current_node}，累计距离: {total_distance:.2f}米，方向: {current_direction:.2f}°")

    return path, total_distance, current_direction


def walk_along_road(G, current_node, target_distance, current_direction=None, roads_gdf=None, transport="walk"):
    """
    沿着初始方向直走，支持多段路，方向偏差不超过5°，直到满足目标距离。
    参数：
        G: 路网图 (networkx.Graph)
        current_node: 当前节点坐标 (tuple)
        target_distance: 目标距离 (float, 单位：米)
        current_direction: 当前方向 (float, 单位：度，可选)
        roads_gdf: 道路几何数据 (GeoDataFrame, 可选)
        transport: 交通方式 (str, 默认 "walk")
    返回：
        path: 生成的路径 (list of nodes)
        actual_distance: 实际移动距离 (float)
        final_direction: 最终方向 (float)
    """
    path = [current_node]
    total_distance = 0
    visited = {current_node}

    # 确定初始方向（如果未提供）
    if current_direction is None:
        neighbors = list(G.neighbors(current_node))
        if neighbors:
            closest_neighbor = min(neighbors, key=lambda n: G[current_node][n]['weight'])
            current_direction = calculate_direction(current_node, closest_neighbor)
            print(f"      初始化方向: 从最近节点确定为 {current_direction:.2f}°")
        else:
            print(f"      警告: 当前节点 {current_node} 没有邻居，无法确定初始方向")
            current_direction = 0.0

    print(f"      沿着方向 {current_direction:.2f}° 直走，目标距离: {target_distance}米，允许偏差: ±5°，交通方式: {transport}")

    while total_distance < target_distance:
        neighbors = [n for n in G.neighbors(current_node) if n not in visited]
        if not neighbors:
            print("      没有更多邻居，停止搜索")
            break

        # 过滤方向偏差不超过5°的邻居
        filtered_neighbors = []
        for neighbor in neighbors:
            angle = calculate_direction(current_node, neighbor)
            angle_diff = min((angle - current_direction) % 360, (current_direction - angle) % 360)
            if angle_diff <= 5.0:
                filtered_neighbors.append((neighbor, angle_diff))

        if not filtered_neighbors:
            print(f"      警告: 未找到方向偏差≤5°的邻居，当前方向 {current_direction:.2f}°，停止搜索")
            break

        # 选择方向偏差最小的邻居
        filtered_neighbors.sort(key=lambda x: x[1])
        next_node = filtered_neighbors[0][0]
        edge_weight = G[current_node][next_node]['weight']
        print(f"      边 {current_node} -> {next_node}: 权重 {edge_weight:.2f}米, 方向偏差 {filtered_neighbors[0][1]:.2f}°")

        # 检查是否需要截断
        if total_distance + edge_weight > target_distance:
            remaining_distance = target_distance - total_distance
            fraction = remaining_distance / edge_weight
            new_lon = current_node[0] + (next_node[0] - current_node[0]) * fraction
            new_lat = current_node[1] + (next_node[1] - current_node[1]) * fraction
            interpolated_node = (new_lon, new_lat)
            path.append(interpolated_node)
            total_distance = target_distance
            print(f"      截断边 {current_node} -> {next_node} 在 {interpolated_node}: {remaining_distance:.2f}米 (比例: {fraction:.2f})")
            break

        path.append(next_node)
        total_distance += edge_weight
        visited.add(next_node)
        current_node = next_node
        current_direction = calculate_direction(path[-2], current_node)
        print(f"      移动到 {current_node}，累计距离: {total_distance:.2f}米，剩余距离: {(target_distance - total_distance):.2f}米，更新方向: {current_direction:.2f}°")

    print(f"      直走路径生成，路径长度: {len(path)}, 实际距离: {total_distance:.2f}米，最终方向: {current_direction:.2f}°")
    return path, total_distance, current_direction

def calculate_path_distance(G, path):
    """计算路径的总距离"""
    total = 0.0
    for i in range(len(path) - 1):
        if G.has_edge(path[i], path[i+1]):
            total += G[path[i]][path[i+1]]['weight']
        elif G.has_edge(path[i+1], path[i]):
            total += G[path[i+1]][path[i]]['weight']
    return total

def walk_distance(G, current_node, target_distance, direction=None, current_direction=None, roads_gdf=None, transport="walk"):
    path = [current_node]
    total_distance = 0
    visited = set([current_node])
    # 添加一个变量记录回溯过的边
    backtracked_edges = set()
    strict_direction = direction is not None

    print(f"      按方向 {direction if direction is not None else current_direction} 移动，目标距离: {target_distance}米，交通方式: {transport}")

    # 计算 target_angle
    if isinstance(direction, str):
        target_angle = direction_mapping.get(direction.lower(), 0)
    elif isinstance(direction, (int, float)):
        target_angle = float(direction)
    else:
        target_angle = current_direction if current_direction is not None else 0
        print(f"      调试: direction 为 None，使用 current_direction={target_angle}°")

    while total_distance < target_distance:
        all_neighbors = [n for n in G.neighbors(current_node) if n not in visited]
        
        # 修改：过滤出回溯过的边
        filtered_neighbors = []
        for n in all_neighbors:
            edge = (current_node, n)
            if edge not in backtracked_edges:
                filtered_neighbors.append(n)
        
        all_neighbors = filtered_neighbors  # 使用过滤后的邻居列表
        
        if not all_neighbors:
            print(f"      没有未访问的邻居节点，累计距离: {total_distance:.2f}米")
            if len(path) > 1:
                print(f"      回溯到 {path[-2]}")
                # 添加回溯边到集合中
                backtracked_edges.add((path[-2], path[-1]))
                backtracked_edges.add((path[-1], path[-2]))  # 添加反向边也标记
                
                # 减去回溯边的距离
                last_edge_weight = G[path[-2]][path[-1]]['weight'] if path[-2] in G[path[-1]] else 0
                total_distance -= last_edge_weight
                current_node = path[-2]
                path.pop()
                visited.remove(path[-1])  # 修改：移除最后一个节点，而不是当前节点
                print(f"      回溯后累计距离: {total_distance:.2f}米")
                continue
            break

        neighbors = all_neighbors
        if direction and strict_direction:
            reference_angle = current_direction if current_direction is not None else (
                calculate_direction(path[-2], current_node) if len(path) > 1 else None
            )
            filtered_neighbors = filter_by_direction(
                current_node, all_neighbors, direction, reference_angle
            )
            neighbors = filtered_neighbors if filtered_neighbors else all_neighbors
            print(f"      方向 {direction} 过滤后邻居数: {len(neighbors)}")

        # 选择方向偏差最小的邻居
        best_neighbor = None
        min_angle_diff = float('inf')
        best_weight = float('inf')

        for neighbor in neighbors:
            edge_weight = G[current_node][neighbor]['weight']
            dx = neighbor[0] - current_node[0]
            dy = neighbor[1] - current_node[1]
            angle = (90 - math.degrees(math.atan2(dy, dx))) % 360
            angle_diff = min((angle - target_angle) % 360, (target_angle - angle) % 360)
            if angle_diff < min_angle_diff or (angle_diff == min_angle_diff and edge_weight < best_weight):
                best_neighbor = neighbor
                min_angle_diff = angle_diff
                best_weight = edge_weight

        if best_neighbor:
            remaining_distance = target_distance - total_distance
            print(f"      边 {current_node} -> {best_neighbor}: 权重 {best_weight:.2f}米, 方向偏差 {min_angle_diff:.2f}°, 当前累计距离: {total_distance:.2f}米, 剩余距离: {remaining_distance:.2f}米")
            
            if total_distance + best_weight <= target_distance:
                total_distance += best_weight
                path.append(best_neighbor)
                visited.add(best_neighbor)
                current_node = best_neighbor
            else:
                ratio = remaining_distance / best_weight if best_weight > 0 else 0
                interpolated_point = interpolate_point(G, current_node, best_neighbor, ratio)
                total_distance = target_distance
                path.append(interpolated_point)
                print(f"      截断边 {current_node} -> {best_neighbor} 在 {interpolated_point}: {remaining_distance:.2f}米 (比例: {ratio:.2f})")
                break
        else:
            print(f"      无有效邻居，终止路径，累计距离: {total_distance:.2f}米")
            break

    print(f"      完整路径长度: {len(path)}, 累计距离: {total_distance:.2f}米")
    return path, total_distance