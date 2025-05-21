from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import Dict, Any, List, Optional
import os
import math
from shapely.geometry import Point

from smam.services.planner import (
    get_amap_coordinates,
    gcj02_to_wgs84,
    build_graph,
    convert_and_print_coordinates,
    advance_segment,
    path_to_geojson,
    insert_point_on_edge,
    get_direction_from_angle
)
from smam.services.utils import haversine  # 仅导入 haversine

router = APIRouter()

SHAPEFILE_PATH = os.path.join(os.path.dirname(__file__), "D:/git/smam/social_media_agent_map/roads_nj/roads_nj.shp")

class PointRequest(BaseModel):
    name: str

class RouteRequest(BaseModel):
    start: Dict[str, str]
    end: Dict[str, str]
    segments: Dict[str, Any]

@router.post("/point/")
async def get_point_coordinates(request: PointRequest):
    """获取地点的坐标信息并转换为WGS84坐标系"""
    try:
        wgs_coords = convert_and_print_coordinates(request.name)
        if wgs_coords is None:
            raise HTTPException(status_code=404, detail=f"找不到地点: {request.name}")
        
        gcj_coords = get_amap_coordinates(request.name)
        return {
            "name": request.name,
            "gcj02": {"lon": gcj_coords[0], "lat": gcj_coords[1]},
            "wgs84": {"lon": wgs_coords[0], "lat": wgs_coords[1]}
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
@router.post("/simulate/")
async def simulate_route(request: RouteRequest):
    try:
        start_name = request.start.get("name")
        gcj_coords = get_amap_coordinates(start_name)
        if not gcj_coords:
            raise HTTPException(status_code=404, detail=f"无法获取起点 '{start_name}' 的坐标")

        start_wgs = gcj02_to_wgs84(*gcj_coords)
        G, roads_gdf = build_graph(SHAPEFILE_PATH)

        current_node = insert_point_on_edge(G, start_wgs, tolerance=9e-5)
        full_path = [current_node]
        instructions = []
        segment_paths = []
        current_direction = None
        possible_end_points = []
        circles = []
        circle_created = False  # 初始化圆形范围生成标志
        sorted_segments = sorted(request.segments.items(), key=lambda x: int(x[0].split('_')[-1]))

        # 获取最后一个segment的ID和编号
        last_segment_id = sorted_segments[-1][0] if sorted_segments else None
        last_seg_num = int(last_segment_id.split('_')[-1]) if last_segment_id else 0
        print(f"最后一个segment的ID是: {last_segment_id}, 编号是: {last_seg_num}")

        for seg_id, segment in sorted_segments:
            current_seg_num = int(seg_id.split('_')[-1])
            print(f"处理segment: {seg_id}, 编号: {current_seg_num}")
            
            path_instructions = segment.get("path_instructions", [])
            segment_result = advance_segment(G, current_node, path_instructions, current_direction, roads_gdf=roads_gdf)

            if not isinstance(segment_result, tuple) or len(segment_result) != 6:
                raise ValueError(f"无效的 advance_segment 返回值")

            segment_path, segment_distance, current_direction, all_paths, circle_needed, circle_radius = segment_result

            if circle_needed and circle_radius is not None:
                circle_geojson = create_circle_geojson(current_node, circle_radius)
                circles.append(circle_geojson)
                circle_created = True  # 设置标志
            elif segment_path and len(segment_path) > 1:  # 仅当路径实际移动时更新
                current_node = segment_path[-1]
                full_path.extend(segment_path[1:] if len(full_path) > 0 else segment_path)
                segment_paths.append({
                    "id": seg_id,
                    "path": segment_path,
                    "distance": segment_distance,
                    "point_count": len(segment_path),
                    "seg_num": current_seg_num
                })
                instr_text = format_instructions(path_instructions)
                instructions.append({
                    "segment": seg_id,
                    "description": instr_text,
                    "distance": segment_distance,
                    "point_count": len(segment_path)
                })
                
                # 记录该段中所有可能的终点
                for path_idx, (path, dist, dir) in enumerate(all_paths):
                    if len(path) >= 1 and not (path_idx == 0 and path == segment_path):
                        end_point = path[-1]
                        # 记录该点所属的段ID, 用于后续判断颜色
                        possible_end_points.append({
                            "coordinates": {"lon": end_point[0], "lat": end_point[1]},
                            "segment_id": seg_id,
                            "seg_num": current_seg_num,
                            "distance": dist,
                            "direction": dir
                        })

        if len(full_path) <= 1 and not circles:
            raise HTTPException(status_code=400, detail="无法生成有效路径")

        end_node = full_path[-1]

        # 定义颜色常量
        START_COLOR = "#0000FF"  # 起点 - 蓝色
        END_COLOR = "#FF0000"    # 最终终点 - 红色
        INTERMEDIATE_COLOR = "#00FF00"  # 中间点 - 绿色
        
        route_geojson = {
            "type": "FeatureCollection",
            "features": []
        }

        # 添加起点 - 蓝色
        route_geojson["features"].append({
            "type": "Feature",
            "geometry": {"type": "Point", "coordinates": [start_wgs[0], start_wgs[1]]},
            "properties": {"color": START_COLOR, "type": "startpoint", "name": start_name}
        })

        # 添加路径
        if len(full_path) > 1:
            geojson = path_to_geojson(full_path)
            geojson["features"][0]["properties"] = {"color": START_COLOR, "type": "path"}
            route_geojson["features"].append(geojson["features"][0])

        # 添加路径点 - 区分中间段终点(绿色)和最终终点(红色)
        added_points = set()  # 用于跟踪已添加的点坐标，防止重复
        
        # 处理所有段的终点
        for segment in segment_paths:
            seg_id = segment["id"]
            seg_num = segment["seg_num"]
            
            # 确保有路径点可以添加
            if len(segment["path"]) > 0:
                # 获取该段的最后一个点
                end_point = segment["path"][-1]
                point_key = f"{end_point[0]:.6f},{end_point[1]:.6f}"
                
                # 如果该点已经添加过，则跳过
                if point_key in added_points:
                    print(f"跳过重复点: {point_key}")
                    continue
                
                # 判断是否是最后一个segment
                is_last_segment = (seg_num == last_seg_num)
                point_color = END_COLOR if is_last_segment else INTERMEDIATE_COLOR
                point_type = "endpoint" if is_last_segment else "intermediate_point"
                
                print(f"添加{'最终' if is_last_segment else '中间'}点: {point_key}, 段{seg_id} (值:{seg_num}), 颜色: {point_color}")
                
                route_geojson["features"].append({
                    "type": "Feature",
                    "geometry": {"type": "Point", "coordinates": end_point},
                    "properties": {
                        "color": point_color, 
                        "type": point_type,
                        "segment_id": seg_id,
                        "visible": True        # 添加可见性属性，确保点可见
                    }
                })
                added_points.add(point_key)

        # 添加圆形范围
        for circle in circles:
            route_geojson["features"].append(circle)

        # 添加可能的终点路径点
        if possible_end_points:
            for end_point in possible_end_points:
                # 获取点坐标
                coords = [end_point["coordinates"]["lon"], end_point["coordinates"]["lat"]]
                point_key = f"{coords[0]:.6f},{coords[1]:.6f}"
                
                # 如果该点已经添加过，则跳过
                if point_key in added_points:
                    print(f"跳过可能终点重复点: {point_key}")
                    continue
                
                segment_id = end_point["segment_id"]
                seg_num = end_point["seg_num"]
                
                # 判断颜色：只有最后一个segment的终点才为红色，其他都为绿色
                is_last_segment = (seg_num == last_seg_num)
                point_color = END_COLOR if is_last_segment else INTERMEDIATE_COLOR
                
                print(f"添加可能{'最终' if is_last_segment else '中间'}终点: {point_key}, 段{segment_id} (值:{seg_num}), 颜色: {point_color}")
                
                route_geojson["features"].append({
                    "type": "Feature",
                    "geometry": {
                        "type": "Point",
                        "coordinates": coords
                    },
                    "properties": {
                        "color": point_color,
                        "type": "possible_endpoint",
                        "segment_id": segment_id,
                        "distance": end_point["distance"]
                    }
                })
                added_points.add(point_key)

        total_distance = calculate_path_distance(G, full_path)

        return {
            "start": {
                "name": start_name,
                "coordinates": {"lon": start_wgs[0], "lat": start_wgs[1]}
            },
            "end": {
                "name": request.end.get("name", "终点"),
                "coordinates": {"lon": end_node[0], "lat": end_node[1]}
            },
            "route": route_geojson,
            "distance": total_distance,
            "node_count": len(full_path),
            "segments": [
                {
                    "id": seg["id"],
                    "distance": seg["distance"],
                    "point_count": seg["point_count"],
                    "color": "#00CED1"
                } for i, seg in enumerate(segment_paths)
            ],
            "instructions": instructions,
            "possible_end_points": possible_end_points
        }

    except Exception as e:
        import traceback
        print(traceback.format_exc())
        raise HTTPException(status_code=500, detail=str(e))

def create_circle_geojson(center, radius, num_points=64):
    """
    创建一个圆形 Polygon feature，用于渲染地图上真实范围的圆。
    :param center: (lon, lat)
    :param radius: 单位：米
    :param num_points: 圆上点的数量，越多越圆
    """
    lon, lat = center
    earth_radius = 6371000  # 地球半径（米）

    coords = []
    for i in range(num_points + 1):
        angle = math.radians(float(i) / num_points * 360)
        dx = radius * math.cos(angle)
        dy = radius * math.sin(angle)

        dlat = dy / earth_radius
        dlon = dx / (earth_radius * math.cos(math.pi * lat / 180))

        new_lat = lat + dlat * 180 / math.pi
        new_lon = lon + dlon * 180 / math.pi
        coords.append([new_lon, new_lat])

    return {
        "type": "Feature",
        "geometry": {
            "type": "Polygon",
            "coordinates": [coords]
        },
        "properties": {
            "type": "circle",
            "color": "#00CED1",
            "opacity": 0.3
        }
    }

def calculate_path_distance(G, path):
    """计算路径的总距离"""
    total_distance = 0
    
    i = 0
    while i < len(path) - 1:
        p1, p2 = path[i], path[i + 1]
        
        if G.has_edge(p1, p2):
            weight = G[p1][p2]['weight']
            total_distance += weight
            print(f"边 {p1} -> {p2}: {weight:.2f}米")
        else:
            edge_data = G.get_edge_data(p1, p2)
            if edge_data and 'original_edge' in edge_data:
                orig_u, orig_v = edge_data['original_edge']
                total_length = haversine(orig_u[1], orig_u[0], orig_v[1], orig_v[0])
                segment_length = haversine(p1[1], p1[0], p2[1], p2[0])
                if total_length > 0:
                    ratio = segment_length / total_length
                    segment_weight = G[orig_u][orig_v]['weight'] * ratio
                    total_distance += segment_weight
                    print(f"截断边 {p1} -> {p2} 在边 {orig_u} -> {orig_v} 上: {segment_weight:.2f}米 (比例: {ratio:.2f})")
            else:
                dist = haversine(p1[1], p1[0], p2[1], p2[0])
                total_distance += dist
                print(f"未找到原始边，使用直线距离: {p1} -> {p2}: {dist:.2f}米")
        
        i += 1
    return total_distance

def format_instructions(instructions):
    """格式化路径指令为字符串"""
    if not instructions:
        return "No instructions provided"
    
    result = []
    for step in instructions:
        step_parts = []
        
        if "direction" in step:
            step_parts.append(f"Go {step['direction']}")
        
        if "transportation" in step:
            step_parts.append(f"by {step['transportation']}")
        
        if "distance" in step:
            step_parts.append(f"for {step['distance']}")
        
        if "time" in step:
            step_parts.append(f"({step['time']})")
        
        result.append(" ".join(step_parts))
    
    return " → ".join(result)