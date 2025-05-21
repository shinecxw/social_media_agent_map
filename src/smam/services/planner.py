import math
import requests
import networkx as nx
import geopandas as gpd
from shapely.geometry import LineString, MultiLineString, Point
from fastapi import HTTPException
import shapely.geometry as geom
import pandas as pd
from smam.services.utils import (
    process_relative_direction,
    calculate_direction,
    get_amap_coordinates,
    haversine,
    filter_by_direction,
    get_direction_from_angle,
    walk_distance,
    find_path_end,
    find_intersection,
    gcj02_to_wgs84,
    path_to_geojson,
    get_speed_by_transportation,
    walk_along_road,
    find_nearest_exit,
    direction_mapping,
    parse_time_to_minutes,
    parse_distance,
    explore_all_paths_no_backtrack,
    find_nearest_road_point,
    calculate_path_distance

)

def build_graph(shapefile_path):
    """构建路网图，使用几何长度计算权重，并记录节点所属道路"""
    gdf = gpd.read_file(shapefile_path)
    G = nx.Graph()

    for idx, row in gdf.iterrows():
        geom = row.geometry
        if "OBJECTID" not in row or not isinstance(row["OBJECTID"], (int, float)) or pd.isna(row["OBJECTID"]):
            raise ValueError(f"Shapefile 行 {idx} 缺少有效的数字 'OBJECTID' 字段，当前值: {row.get('OBJECTID')}")
        road_id = int(row["OBJECTID"])
        distance = row.get("distance", 0)
        if isinstance(geom, LineString):
            lines = [geom]
        elif isinstance(geom, MultiLineString):
            lines = list(geom.geoms)
        else:
            continue

        for line in lines:
            coords = list(line.coords)
            if len(coords) < 2:
                continue

            total_geom_length = sum(
                haversine(coords[i][0], coords[i][1], coords[i+1][0], coords[i+1][1])
                for i in range(len(coords) - 1)
            )

            for i in range(len(coords) - 1):
                p1, p2 = coords[i], coords[i + 1]
                lon1, lat1 = p1
                lon2, lat2 = p2
                segment_geom_length = haversine(lon1, lat1, lon2, lat2)
                segment_weight = segment_geom_length
                G.add_edge(p1, p2, weight=segment_weight, geometry=line, road_id=road_id)
                if p1 not in G.nodes:
                    G.add_node(p1, road_id=road_id)
                else:
                    G.nodes[p1]['road_id'] = road_id
                if p2 not in G.nodes:
                    G.add_node(p2, road_id=road_id)
                else:
                    G.nodes[p2]['road_id'] = road_id
                #print(f"子段 {p1} -> {p2}: 权重 {segment_weight:.2f}米，所属道路 ID: {road_id}")

    # 验证边是否有 road_id
    for u, v, data in G.edges(data=True):
        if 'road_id' not in data:
            print(f"警告: 边 {u} -> {v} 缺少 road_id，可能影响路径生成")

    print(f"路网节点数: {len(G.nodes)}，边数: {len(G.edges)}")
    components = list(nx.connected_components(G))
    print(f"路网连通分量数量: {len(components)}")
    if len(components) > 1:
        print("警告：路网不完全连通，可能导致路径生成失败")
        for i, comp in enumerate(components):
            print(f"连通分量 {i+1}：节点数 {len(comp)}")
    return G, gdf

def insert_point_on_edge(G, point, tolerance=1e-5):
    """在最近的边上插入点，并继承原始边的 road_id"""
    px, py = point
    min_dist = float('inf')
    best_edge = None
    best_proj = None

    point_geom = geom.Point(px, py)

    # 尝试找到最近的边
    for u, v, data in G.edges(data=True):
        line = geom.LineString([u, v])
        proj = line.interpolate(line.project(point_geom))
        dist = proj.distance(point_geom)
        if dist < min_dist:
            min_dist = dist
            best_edge = (u, v)
            best_proj = (proj.x, proj.y)

    if best_edge and best_proj and min_dist < tolerance:
        u, v = best_edge
        weight = G[u][v]['weight']
        road_id = G[u][v].get('road_id')
        if road_id is None:
            print(f"警告: 边 {u} -> {v} 没有 road_id，可能影响路径生成")

        G.remove_edge(u, v)

        total_length = math.hypot(v[0] - u[0], v[1] - u[1])
        d1 = math.hypot(best_proj[0] - u[0], best_proj[1] - u[1])
        d2 = math.hypot(v[0] - best_proj[0], v[1] - best_proj[1])

        w1 = weight * (d1 / total_length) if total_length > 0 else 0
        w2 = weight * (d2 / total_length) if total_length > 0 else 0

        G.add_edge(u, best_proj, weight=w1, road_id=road_id)
        G.add_edge(best_proj, v, weight=w2, road_id=road_id)
        G.add_node(best_proj, road_id=road_id)
        print(f"插入新节点 {best_proj}，继承 road_id: {road_id}")
        return best_proj

    # 如果未找到合适的边，选择最近的节点
    nearest_node = find_nearest_road_point(G, point)
    print(f"未找到合适的边插入点 {point}，选择最近节点 {nearest_node}")
    return nearest_node


def advance_segment(G, current_node, instructions, prev_direction=None, roads_gdf=None):
    if not instructions:
        print(f"      无指令，返回单点路径: {current_node}")
        return [current_node], 0.0, prev_direction, [], False, None

    path = [current_node]
    total_distance = 0.0
    current_direction = prev_direction
    current_road_id = G.nodes[current_node].get("road_id") if current_node in G.nodes else None
    all_paths = []
    circle_needed = False
    circle_radius = None

    if current_direction is None:
        neighbors = list(G.neighbors(current_node))
        if neighbors:
            closest_neighbor = min(neighbors, key=lambda n: G[current_node][n]['weight'])
            current_direction = calculate_direction(current_node, closest_neighbor)
            print(f"      初始化方向: 从最近邻居 {closest_neighbor} 确定为 {current_direction:.2f}°")
        else:
            print(f"      警告: 当前节点 {current_node} 无邻居，设置默认方向 0°")
            current_direction = 0.0

    print(f"      初始位置: {current_node}, 初始方向: {current_direction:.2f}°, 道路ID: {current_road_id}")

    for i, step in enumerate(instructions):
        print(f"    执行指令 {i+1}/{len(instructions)}: {step}")
        action = step.get("action")
        direction = step.get("direction")
        reference_direction = step.get("reference_direction")
        transportation = step.get("transport", "walk")
        target_time = step.get("time")
        target_distance = None  # 推迟 parse_distance 调用

        print(f"      目标距离: {target_distance}米, 交通方式: {transportation}, 指令方向: {direction or '未指定'}, 参照方向: {reference_direction or '未指定'}, 时间: {target_time or '未指定'}")
        if not isinstance(current_node, tuple) or len(current_node) != 2:
            print(f"      错误: 当前节点 {current_node} 不是有效坐标元组")
            raise ValueError(f"无效节点: {current_node}")

        # 第一层：动作 go straight
        if action == "go straight" and step.get("distance"):
            target_distance = parse_distance(step, transportation)
            segment_path, segment_distance, new_direction = walk_along_road(
                G, current_node, target_distance, current_direction, roads_gdf=roads_gdf, transport=transportation
            )
            circle_needed = False
            circle_radius = None
            if segment_path and len(segment_path) > 1:
                path.extend(segment_path[1:])
                total_distance = calculate_path_distance(G, path)
                current_node = segment_path[-1]
                current_road_id = G.nodes[current_node].get("road_id") if current_node in G.nodes else current_road_id
                current_direction = new_direction
                print(f"      go straight 路径生成，距离: {segment_distance:.2f}米, 新节点: {current_node}, 方向: {current_direction:.2f}°")
            continue

        # 第一层：动作 go to the end
        elif action == "go to the end":
            circle_needed = False
            circle_radius = None
            segment_path, segment_distance, new_direction = find_path_end(
                G, current_node, current_direction, max_distance=5000, roads_gdf=roads_gdf,
                current_road_id=current_road_id, transport=transportation
            )
            if segment_path and len(segment_path) > 1:
                path.extend(segment_path[1:])
                total_distance = calculate_path_distance(G, path)
                current_node = segment_path[-1]
                current_road_id = G.nodes[current_node].get("road_id") if current_node in G.nodes else current_road_id
                current_direction = new_direction
                print(f"      走到道路尽头，距离: {segment_distance:.2f}米, 新节点: {current_node}, 方向: {current_direction:.2f}°")
            else:
                print(f"      警告: 未生成有效路径，保持节点 {current_node}")
            continue

        # 第一层：动作 Turn（相对方向）
        elif action == "Turn" and isinstance(direction, str):
            new_direction = process_relative_direction(current_direction, direction)
            if new_direction is None:
                print(f"      警告: 无法处理相对方向 {direction}")
                continue
            new_node, _, turn_path, new_road_id = handle_turn(
                G, current_node, current_direction, direction, roads_gdf, transport=transportation
            )
            if turn_path and len(turn_path) > 1:
                path.extend(turn_path[1:])
                turn_distance = sum(G[path[i]][path[i+1]]['weight'] for i in range(len(path)-len(turn_path), len(path)-1))
                total_distance += turn_distance
                current_node = new_node
                current_road_id = new_road_id if new_road_id is not None else current_road_id
                current_direction = new_direction
                print(f"      交通方式 + 转弯 {direction}，距离: {turn_distance:.2f}米, 新节点: {current_node}, 方向: {current_direction:.2f}°, 新道路 ID: {current_road_id}")
            else:
                print(f"      交通方式 + 转弯 {direction}，更新节点: {new_node}, 方向: {new_direction:.2f}°, 新道路 ID: {new_road_id}")
                current_node = new_node
                current_road_id = new_road_id if new_road_id is not None else current_road_id
                current_direction = new_direction
            continue

        # 第一层：动作 OUT
        elif action == "OUT":
            circle_needed = False
            circle_radius = None
            segment_path, segment_distance, new_direction = find_nearest_exit(
                G, current_node, current_direction, max_distance=1000, roads_gdf=roads_gdf,
                current_road_id=current_road_id, transport=transportation
            )
            if segment_path and len(segment_path) > 1:
                path.extend(segment_path[1:])
                total_distance = calculate_path_distance(G, path)
                current_node = segment_path[-1]
                current_road_id = G.nodes[current_node].get("road_id") if current_node in G.nodes else current_road_id
                current_direction = new_direction
                print(f"      动作 OUT，找到最近出口，距离: {segment_distance:.2f}米, 新节点: {current_node}, 方向: {current_direction:.2f}°")
            else:
                print(f"      警告: 未找到有效出口，保持节点 {current_node}")
            continue

        # 第一层：只有交通方式（没有任何其他条件）
        elif transportation and not action and not direction and not reference_direction and not step.get("distance") and not target_time:
            print(f"      仅有交通方式: {transportation}，生成圆形范围")
            circle_radius = 2000  # 步行默认 2km
            if transportation == "bike":
                circle_radius = 5000  # 骑行 5km
            elif transportation == "car":
                circle_radius = 10000  # 驾车 10km
            circle_needed = True
            all_paths = []
            print(f"      生成圆形范围，半径: {circle_radius:.0f}米，中心点: {current_node}")
            continue

        # 第一层：在“有交通方式”下的第二层嵌套
        elif transportation:
            # 第二层：交通方式 + 绝对方向 + 时间（优先处理 direction 和 time 的组合）
            if direction and target_time and not step.get("distance") and not reference_direction and not action:
                abs_direction = None
                if isinstance(direction, str):
                    direction = direction.lower().replace("°", "")
                    abs_direction = direction_mapping.get(direction, None)
                    if abs_direction is None:
                        try:
                            abs_direction = float(direction)
                        except ValueError:
                            print(f"      警告: 方向 {direction} 格式错误，跳过")
                            continue
                elif isinstance(direction, (int, float)):
                    abs_direction = float(direction)
                else:
                    print(f"      警告: 无效方向格式 {direction}，跳过")
                    continue

                time_value = parse_time_to_minutes(target_time)
                if time_value is None:
                    print(f"      警告: 无法解析时间 {target_time}，跳过")
                    continue
                speed = get_speed_by_transportation(transportation)
                target_distance = time_value * speed
                print(f"      交通方式 + 绝对方向 + 时间: 方向 {direction}（{abs_direction}°），时间 {target_time}，转换为距离: {target_distance:.2f}米")

                segment_path, segment_distance = walk_distance(
                    G, current_node, target_distance, direction=abs_direction,
                    current_direction=current_direction, roads_gdf=roads_gdf, transport=transportation
                )
                if segment_path and len(segment_path) > 1:
                    path.extend(segment_path[1:])
                    total_distance = calculate_path_distance(G, path)
                    current_node = segment_path[-1]
                    current_road_id = G.nodes[current_node].get("road_id") if current_node in G.nodes else current_road_id
                    current_direction = calculate_direction(segment_path[-2], current_node) if len(segment_path) >= 2 else current_direction
                    print(f"      交通方式 + 绝对方向 + 时间: 按绝对方向 {abs_direction}° 移动，距离: {segment_distance:.2f}米, 新节点: {current_node}, 方向: {current_direction:.2f}°")
                else:
                    print(f"      警告: 无法生成路径，路径长度: {len(segment_path)}")
                continue

            # 第二层：交通方式 + 绝对方向 + 距离
            elif direction and step.get("distance") and not target_time and not reference_direction and not action:
                target_distance = parse_distance(step, transportation)
                abs_direction = None
                if isinstance(direction, str):
                    direction = direction.lower().replace("°", "")
                    abs_direction = direction_mapping.get(direction, None)
                    if abs_direction is None:
                        try:
                            abs_direction = float(direction)
                        except ValueError:
                            print(f"      警告: 方向 {direction} 格式错误，跳过")
                            continue
                elif isinstance(direction, (int, float)):
                    abs_direction = float(direction)
                else:
                    print(f"      警告: 无效方向格式 {direction}，跳过")
                    continue

                segment_path, segment_distance = walk_distance(
                    G, current_node, target_distance, direction=abs_direction,
                    current_direction=current_direction, roads_gdf=roads_gdf, transport=transportation
                )
                if segment_path and len(segment_path) > 1:
                    path.extend(segment_path[1:])
                    total_distance = calculate_path_distance(G, path)
                    current_node = segment_path[-1]
                    current_road_id = G.nodes[current_node].get("road_id") if current_node in G.nodes else current_road_id
                    current_direction = calculate_direction(segment_path[-2], current_node) if len(segment_path) >= 2 else current_direction
                    print(f"      交通方式 + 绝对方向 + 距离: 按绝对方向 {abs_direction}° 移动，距离: {segment_distance:.2f}米, 新节点: {current_node}, 方向: {current_direction:.2f}°")
                else:
                    print(f"      警告: 无法生成路径，路径长度: {len(segment_path)}")
                continue

            # 第二层：交通方式 + 绝对方向（无距离和时间时，使用默认距离）
            elif direction and not step.get("distance") and not target_time and not reference_direction and not action:
                abs_direction = None
                if isinstance(direction, str):
                    direction = direction.lower().replace("°", "")
                    abs_direction = direction_mapping.get(direction, None)
                    if abs_direction is None:
                        print(f"      警告: 方向 {direction} 未在映射表中，跳过")
                        continue
                elif isinstance(direction, (int, float)):
                    abs_direction = float(direction)
                else:
                    print(f"      警告: 无效方向格式 {direction}，跳过")
                    continue

                target_distance = 1000  # 默认距离 1000 米
                print(f"      交通方式 + 绝对方向: 方向 {direction}（{abs_direction}°），使用默认距离: {target_distance:.2f}米")

                segment_path, segment_distance = walk_distance(
                    G, current_node, target_distance, direction=abs_direction,
                    current_direction=current_direction, roads_gdf=roads_gdf, transport=transportation
                )
                if segment_path and len(segment_path) > 1:
                    path.extend(segment_path[1:])
                    total_distance = calculate_path_distance(G, path)
                    current_node = segment_path[-1]
                    current_road_id = G.nodes[current_node].get("road_id") if current_node in G.nodes else current_road_id
                    current_direction = calculate_direction(segment_path[-2], current_node) if len(segment_path) >= 2 else current_direction
                    print(f"      交通方式 + 绝对方向: 按绝对方向 {abs_direction}° 移动，距离: {segment_distance:.2f}米, 新节点: {current_node}, 方向: {current_direction:.2f}°")
                else:
                    print(f"      警告: 无法生成路径，路径长度: {len(segment_path)}")
                continue

            # 第二层：交通方式 + 距离（无 direction）
            elif step.get("distance") and not direction and not reference_direction and not action:
                target_distance = parse_distance(step, transportation)
                print(f"      处理距离指令: {target_distance}米，交通方式: {transportation}")
                if i == 0 or current_direction is None:
                    all_possible_paths = explore_all_paths_no_backtrack(G, current_node, target_distance, transportation)
                    if all_possible_paths:
                        all_paths.extend(all_possible_paths)
                        main_path, main_distance, main_direction = all_possible_paths[0]
                        path.extend(main_path[1:])
                        total_distance += main_distance
                        current_node = main_path[-1]
                        current_direction = main_direction
                        print(f"      只有距离 {target_distance}米，第一步或无方向，找到 {len(all_possible_paths)} 条路径，选用第一条，实际距离: {main_distance:.2f}米, 新节点: {current_node}, 方向: {current_direction:.2f}°")
                    else:
                        print(f"      警告: 未生成有效路径，保持节点 {current_node}")
                else:
                    segment_path, segment_distance = walk_distance(
                        G, current_node, target_distance, current_direction=current_direction,
                        roads_gdf=roads_gdf, transport=transportation
                    )
                    if segment_path and len(segment_path) > 1:
                        path.extend(segment_path[1:])
                        total_distance = calculate_path_distance(G, path)
                        current_node = segment_path[-1]
                        current_direction = calculate_direction(segment_path[-2], segment_path[-1]) if len(segment_path) >= 2 else current_direction
                        print(f"      只有距离 {target_distance}米，非第一步，沿方向 {current_direction:.2f}°，实际距离: {segment_distance:.2f}米, 新节点: {current_node}")
                continue

            # 第二层：交通方式 + 时间（无 direction 和 distance）
            elif target_time and not step.get("distance") and not direction and not reference_direction and not action:
                time_value = parse_time_to_minutes(target_time)
                if time_value is None:
                    print(f"      警告: 无法解析时间 {target_time}，跳过")
                    continue
                speed = get_speed_by_transportation(transportation)
                target_distance = time_value * speed
                print(f"      处理时间指令: {target_time}，交通方式: {transportation}，转换为距离: {target_distance:.2f}米")
                if i == 0 or current_direction is None:
                    all_possible_paths = explore_all_paths_no_backtrack(G, current_node, target_distance, transportation)
                    if all_possible_paths:
                        all_paths.extend(all_possible_paths)
                        main_path, main_distance, main_direction = all_possible_paths[0]
                        path.extend(main_path[1:])
                        total_distance += main_distance
                        current_node = main_path[-1]
                        current_direction = main_direction
                        print(f"      交通方式 + 时间: 第一步或无方向，转换为距离 {target_distance:.2f}米，找到 {len(all_possible_paths)} 条路径，选用第一条，实际距离: {main_distance:.2f}米, 新节点: {current_node}, 方向: {current_direction:.2f}°")
                    else:
                        print(f"      警告: 未生成有效路径，保持节点 {current_node}")
                else:
                    segment_path, segment_distance = walk_distance(
                        G, current_node, target_distance, current_direction=current_direction,
                        roads_gdf=roads_gdf, transport=transportation
                    )
                    if segment_path and len(segment_path) > 1:
                        path.extend(segment_path[1:])
                        total_distance = calculate_path_distance(G, path)
                        current_node = segment_path[-1]
                        current_direction = calculate_direction(segment_path[-2], segment_path[-1]) if len(segment_path) >= 2 else current_direction
                        print(f"      交通方式 + 时间: 非第一步，转换为距离 {target_distance:.2f}米，沿方向 {current_direction:.2f}°，实际距离: {segment_distance:.2f}米, 新节点: {current_node}")
                continue

            # 第二层：交通方式 + 参照性方向 + 距离
            elif reference_direction and step.get("distance") and not target_time and not action and not direction:
                target_distance = parse_distance(step, transportation)
                circle_needed = False
                circle_radius = None
                print(f"处理参照性方向 + 距离指令: 参照方向 {reference_direction}, 距离 {target_distance}米, 交通方式: {transportation}")

                ref_coords = get_amap_coordinates(reference_direction) if isinstance(reference_direction, str) else reference_direction
                if ref_coords:
                    ref_lng, ref_lat = ref_coords
                    ref_lng, ref_lat = gcj02_to_wgs84(ref_lng, ref_lat)
                    ref_coords = (round(ref_lng, 6), round(ref_lat, 6))
                    print(f"将参照方向名称 '{reference_direction}' 转换为坐标: {ref_coords}")
                else:
                    print(f"警告: 无法将参照方向名称 '{reference_direction}' 转换为坐标，跳过")
                    continue

                current_node = (round(current_node[0], 6), round(current_node[1], 6))
                nearest_start_node = insert_point_on_edge(G, current_node, tolerance=1e-4)
                nearest_ref_node = insert_point_on_edge(G, ref_coords, tolerance=1e-4)
                print(f"当前位置 {current_node} 映射到路网节点: {nearest_start_node}")
                print(f"参照点 {ref_coords} 映射到路网节点: {nearest_ref_node}")

                path = [nearest_start_node]
                current_node = nearest_start_node
                visited = {current_node}
                remaining_distance = target_distance
                total_distance = 0.0
                max_steps = 1000
                max_backtracks = 50
                backtrack_count = 0
                current_road_id = G.nodes.get(nearest_start_node, {}).get("road_id")
                current_path_direction = None

                target_direction = calculate_direction(nearest_start_node, nearest_ref_node)
                print(f"目标方向: {target_direction:.2f}°")

                while remaining_distance > 0 and max_steps > 0 and backtrack_count < max_backtracks:
                    neighbors = list(G.neighbors(current_node))
                    if not neighbors:
                        print(f"节点 {current_node} 无邻居，尝试回溯")
                        if len(path) <= 1:
                            print("已到起点，无法回溯")
                            break
                        path.pop()
                        current_node = path[-1]
                        visited.remove(current_node)
                        if len(path) >= 2:
                            total_distance -= G[path[-2]][current_node].get('weight', 0.0)
                            remaining_distance += G[path[-2]][current_node].get('weight', 0.0)
                            current_path_direction = calculate_direction(path[-2], path[-1])
                        backtrack_count += 1
                        print(f"回溯到 {current_node}，剩余距离: {remaining_distance:.2f}米，当前方向: {current_path_direction:.2f}°")
                        continue

                    is_intersection = False
                    road_ids = set()
                    for neighbor in neighbors:
                        road_ids.add(G.nodes[neighbor].get("road_id", current_road_id))
                    if len(road_ids) > 1:
                        is_intersection = True
                        print(f"当前节点 {current_node} 是路口，检测到不同道路 ID: {road_ids}")

                    candidates = []
                    for neighbor in neighbors:
                        if neighbor in visited:
                            print(f"邻居 {neighbor} 已被访问，跳过")
                            continue
                        neighbor_direction = calculate_direction(current_node, neighbor)
                        raw_diff = abs(neighbor_direction - target_direction)
                        angle_diff_to_target = raw_diff if raw_diff <= 180 else 360 - raw_diff
                        dist_to_ref = haversine(neighbor[0], neighbor[1], nearest_ref_node[0], nearest_ref_node[1])
                        candidates.append((neighbor, angle_diff_to_target, dist_to_ref))
                        print(f"邻居 {neighbor}，方向: {neighbor_direction:.2f}°，与目标方向绝对夹角: {angle_diff_to_target:.2f}°，到参照点距离: {dist_to_ref:.2f}米")

                    if not candidates:
                        print(f"节点 {current_node} 无有效邻居，尝试回溯")
                        if len(path) <= 1:
                            print("已到起点，无法回溯")
                            break
                        path.pop()
                        current_node = path[-1]
                        visited.remove(current_node)
                        if len(path) >= 2:
                            total_distance -= G[path[-2]][current_node].get('weight', 0.0)
                            remaining_distance += G[path[-2]][current_node].get('weight', 0.0)
                            current_path_direction = calculate_direction(path[-2], path[-1])
                        backtrack_count += 1
                        print(f"回溯到 {current_node}，剩余距离: {remaining_distance:.2f}米，当前方向: {current_path_direction:.2f}°")
                        continue

                    best_neighbor, best_angle_diff_to_target, best_dist_to_ref = min(candidates, key=lambda x: x[1])
                    print(f"选择最佳邻居 {best_neighbor}，绝对夹角: {best_angle_diff_to_target:.2f}°，到参照点距离: {best_dist_to_ref:.2f}米")
                    edge_weight = G[current_node][best_neighbor].get('weight', haversine(current_node[0], current_node[1], best_neighbor[0], best_neighbor[1]))
                    if edge_weight <= 0:
                        print(f"警告: 边 {current_node} -> {best_neighbor} 权重为 {edge_weight}，跳过")
                        visited.add(best_neighbor)
                        continue

                    if edge_weight <= remaining_distance:
                        path.append(best_neighbor)
                        total_distance += edge_weight
                        remaining_distance -= edge_weight
                        current_node = best_neighbor
                        visited.add(current_node)
                        current_path_direction = calculate_direction(path[-2], path[-1]) if len(path) >= 2 else current_path_direction
                        print(f"移动到 {current_node}，累计距离: {total_distance:.2f}米")
                    else:
                        ratio = remaining_distance / edge_weight
                        new_lon = current_node[0] + (best_neighbor[0] - current_node[0]) * ratio
                        new_lat = current_node[1] + (best_neighbor[1] - current_node[1]) * ratio
                        interpolated_node = (round(new_lon, 6), round(new_lat, 6))
                        G.add_node(interpolated_node, road_id=current_road_id)
                        dist_to_current = haversine(current_node[0], current_node[1], interpolated_node[0], interpolated_node[1])
                        dist_to_neighbor = haversine(interpolated_node[0], interpolated_node[1], best_neighbor[0], best_neighbor[1])
                        G.add_edge(current_node, interpolated_node, weight=dist_to_current, road_id=current_road_id)
                        G.add_edge(interpolated_node, best_neighbor, weight=dist_to_neighbor, road_id=current_road_id)
                        path.append(interpolated_node)
                        total_distance += remaining_distance
                        current_node = interpolated_node
                        remaining_distance = 0
                        print(f"截断边到 {interpolated_node}，剩余距离: {remaining_distance:.2f}米，累计距离: {total_distance:.2f}米，添加到路网")

                    max_steps -= 1

                if total_distance > 0:
                    current_road_id = G.nodes.get(current_node, {}).get("road_id", current_road_id)
                    current_direction = calculate_direction(path[-2], path[-1]) if len(path) >= 2 else current_direction
                    all_paths = [(path, total_distance, current_direction)]
                    print(f"朝参照方向移动完成，距离: {total_distance:.2f}米, 新节点: {current_node}, 方向: {current_direction:.2f}°")
                    print(f"返回 all_paths: {all_paths}")
                    return path, total_distance, current_direction, all_paths, circle_needed, circle_radius
                else:
                    print("警告: 无法移动任何距离")
                    raise ValueError("无法生成有效路径")
                continue

            # 第二层：交通方式 + 参照性方向 + 时间
            elif reference_direction and target_time and not step.get("distance") and not action and not direction:
                circle_needed = False
                circle_radius = None
                print(f"处理参照性方向 + 时间 + 交通方式指令: 参照方向 {reference_direction}, 时间 {target_time}, 交通方式: {transportation}")

                time_value = parse_time_to_minutes(target_time)
                if time_value is None:
                    print(f"警告: 无法解析时间 {target_time}，跳过")
                    continue
                speed = get_speed_by_transportation(transportation)
                target_distance = time_value * speed
                print(f"根据时间 {target_time}，转换为距离: {target_distance:.2f}米")

                ref_coords = get_amap_coordinates(reference_direction) if isinstance(reference_direction, str) else reference_direction
                if ref_coords:
                    ref_lng, ref_lat = ref_coords
                    ref_lng, ref_lat = gcj02_to_wgs84(ref_lng, ref_lat)
                    ref_coords = (round(ref_lng, 6), round(ref_lat, 6))
                    print(f"将参照方向名称 '{reference_direction}' 转换为坐标: {ref_coords}")
                else:
                    print(f"警告: 无法将参照方向名称 '{reference_direction}' 转换为坐标，跳过")
                    continue

                current_node = (round(current_node[0], 6), round(current_node[1], 6))
                nearest_start_node = insert_point_on_edge(G, current_node, tolerance=1e-4)
                nearest_ref_node = insert_point_on_edge(G, ref_coords, tolerance=1e-4)
                print(f"当前位置 {current_node} 映射到路网节点: {nearest_start_node}")
                print(f"参照点 {ref_coords} 映射到路网节点: {nearest_ref_node}")

                path = [nearest_start_node]
                current_node = nearest_start_node
                visited = {current_node}
                remaining_distance = target_distance
                total_distance = 0.0
                max_steps = 1000
                max_backtracks = 50
                backtrack_count = 0
                current_road_id = G.nodes.get(nearest_start_node, {}).get("road_id")
                current_path_direction = None

                target_direction = calculate_direction(nearest_start_node, nearest_ref_node)
                print(f"目标方向: {target_direction:.2f}°")

                while remaining_distance > 0 and max_steps > 0 and backtrack_count < max_backtracks:
                    neighbors = list(G.neighbors(current_node))
                    if not neighbors:
                        print(f"节点 {current_node} 无邻居，尝试回溯")
                        if len(path) <= 1:
                            print("已到起点，无法回溯")
                            break
                        path.pop()
                        current_node = path[-1]
                        visited.remove(current_node)
                        if len(path) >= 2:
                            total_distance -= G[path[-2]][current_node].get('weight', 0.0)
                            remaining_distance += G[path[-2]][current_node].get('weight', 0.0)
                            current_path_direction = calculate_direction(path[-2], path[-1])
                        backtrack_count += 1
                        print(f"回溯到 {current_node}，剩余距离: {remaining_distance:.2f}米，当前方向: {current_path_direction:.2f}°")
                        continue

                    is_intersection = False
                    road_ids = set()
                    for neighbor in neighbors:
                        road_ids.add(G.nodes[neighbor].get("road_id", current_road_id))
                    if len(road_ids) > 1:
                        is_intersection = True
                        print(f"当前节点 {current_node} 是路口，检测到不同道路 ID: {road_ids}")

                    candidates = []
                    for neighbor in neighbors:
                        if neighbor in visited:
                            print(f"邻居 {neighbor} 已被访问，跳过")
                            continue
                        neighbor_direction = calculate_direction(current_node, neighbor)
                        raw_diff = abs(neighbor_direction - target_direction)
                        angle_diff_to_target = raw_diff if raw_diff <= 180 else 360 - raw_diff
                        dist_to_ref = haversine(neighbor[0], neighbor[1], nearest_ref_node[0], nearest_ref_node[1])
                        candidates.append((neighbor, angle_diff_to_target, dist_to_ref))
                        print(f"邻居 {neighbor}，方向: {neighbor_direction:.2f}°，与目标方向绝对夹角: {angle_diff_to_target:.2f}°，到参照点距离: {dist_to_ref:.2f}米")

                    if not candidates:
                        print(f"节点 {current_node} 无有效邻居，尝试回溯")
                        if len(path) <= 1:
                            print("已到起点，无法回溯")
                            break
                        path.pop()
                        current_node = path[-1]
                        visited.remove(current_node)
                        if len(path) >= 2:
                            total_distance -= G[path[-2]][current_node].get('weight', 0.0)
                            remaining_distance += G[path[-2]][current_node].get('weight', 0.0)
                            current_path_direction = calculate_direction(path[-2], path[-1])
                        backtrack_count += 1
                        print(f"回溯到 {current_node}，剩余距离: {remaining_distance:.2f}米，当前方向: {current_path_direction:.2f}°")
                        continue

                    best_neighbor, best_angle_diff_to_target, best_dist_to_ref = min(candidates, key=lambda x: x[1])
                    print(f"选择最佳邻居 {best_neighbor}，绝对夹角: {best_angle_diff_to_target:.2f}°，到参照点距离: {best_dist_to_ref:.2f}米")
                    edge_weight = G[current_node][best_neighbor].get('weight', haversine(current_node[0], current_node[1], best_neighbor[0], best_neighbor[1]))
                    if edge_weight <= 0:
                        print(f"警告: 边 {current_node} -> {best_neighbor} 权重为 {edge_weight}，跳过")
                        visited.add(best_neighbor)
                        continue

                    if edge_weight <= remaining_distance:
                        path.append(best_neighbor)
                        total_distance += edge_weight
                        remaining_distance -= edge_weight
                        current_node = best_neighbor
                        visited.add(current_node)
                        current_path_direction = calculate_direction(path[-2], path[-1]) if len(path) >= 2 else current_path_direction
                        print(f"移动到 {current_node}，累计距离: {total_distance:.2f}米")
                    else:
                        ratio = remaining_distance / edge_weight
                        new_lon = current_node[0] + (best_neighbor[0] - current_node[0]) * ratio
                        new_lat = current_node[1] + (best_neighbor[1] - current_node[1]) * ratio
                        interpolated_node = (round(new_lon, 6), round(new_lat, 6))
                        G.add_node(interpolated_node, road_id=current_road_id)
                        dist_to_current = haversine(current_node[0], current_node[1], interpolated_node[0], interpolated_node[1])
                        dist_to_neighbor = haversine(interpolated_node[0], interpolated_node[1], best_neighbor[0], best_neighbor[1])
                        G.add_edge(current_node, interpolated_node, weight=dist_to_current, road_id=current_road_id)
                        G.add_edge(interpolated_node, best_neighbor, weight=dist_to_neighbor, road_id=current_road_id)
                        path.append(interpolated_node)
                        total_distance += remaining_distance
                        current_node = interpolated_node
                        remaining_distance = 0
                        print(f"截断边到 {interpolated_node}，剩余距离: {remaining_distance:.2f}米，累计距离: {total_distance:.2f}米，添加到路网")

                    max_steps -= 1

                if total_distance > 0:
                    current_road_id = G.nodes.get(current_node, {}).get("road_id", current_road_id)
                    current_direction = calculate_direction(path[-2], path[-1]) if len(path) >= 2 else current_direction
                    all_paths = [(path, total_distance, current_direction)]
                    print(f"朝参照方向移动完成，距离: {total_distance:.2f}米, 新节点: {current_node}, 方向: {current_direction:.2f}°")
                    print(f"返回 all_paths: {all_paths}")
                    return path, total_distance, current_direction, all_paths, circle_needed, circle_radius
                else:
                    print("警告: 无法移动任何距离")
                    raise ValueError("无法生成有效路径")
                continue

            else:
                print(f"      交通方式指令未匹配任何子条件，跳过")
                continue

        # 第一层：只有绝对方向
        elif direction and not transportation and not reference_direction and not action:
            # 第二层：绝对方向 + 距离
            if step.get("distance") and not target_time:
                target_distance = parse_distance(step, transportation)
                abs_direction = None
                if isinstance(direction, str):
                    direction = direction.lower().replace("°", "")
                    abs_direction = direction_mapping.get(direction, None)
                    if abs_direction is None:
                        try:
                            abs_direction = float(direction)
                        except ValueError:
                            print(f"      警告: 方向 {direction} 格式错误，跳过")
                            continue
                elif isinstance(direction, (int, float)):
                    abs_direction = float(direction)
                else:
                    print(f"      警告: 无效方向格式 {direction}，跳过")
                    continue

                segment_path, segment_distance = walk_distance(
                    G, current_node, target_distance, direction=abs_direction,
                    current_direction=current_direction, roads_gdf=roads_gdf, transport=transportation
                )
                if segment_path and len(segment_path) > 1:
                    path.extend(segment_path[1:])
                    total_distance = calculate_path_distance(G, path)
                    current_node = segment_path[-1]
                    current_road_id = G.nodes[current_node].get("road_id") if current_node in G.nodes else current_road_id
                    current_direction = calculate_direction(segment_path[-2], current_node) if len(segment_path) >= 2 else current_direction
                    print(f"      绝对方向 + 距离: 按绝对方向 {abs_direction}° 移动，距离: {segment_distance:.2f}米, 新节点: {current_node}, 方向: {current_direction:.2f}°")
                continue

            # 第二层：绝对方向 + 时间
            elif target_time and not step.get("distance"):
                abs_direction = None
                if isinstance(direction, str):
                    direction = direction.lower().replace("°", "")
                    abs_direction = direction_mapping.get(direction, None)
                    if abs_direction is None:
                        print(f"      警告: 方向 {direction} 未在映射表中，跳过")
                        continue
                elif isinstance(direction, (int, float)):
                    abs_direction = float(direction)
                else:
                    print(f"      警告: 无效方向格式 {direction}，跳过")
                    continue

                time_value = parse_time_to_minutes(target_time)
                if time_value is None:
                    print(f"      警告: 无法解析时间 {target_time}，跳过")
                    continue
                speed = get_speed_by_transportation(transportation)
                target_distance = time_value * speed

                segment_path, segment_distance = walk_distance(
                    G, current_node, target_distance, direction=abs_direction,
                    current_direction=current_direction, roads_gdf=roads_gdf, transport=transportation
                )
                if segment_path and len(segment_path) > 1:
                    path.extend(segment_path[1:])
                    total_distance = calculate_path_distance(G, path)
                    current_node = segment_path[-1]
                    current_road_id = G.nodes[current_node].get("road_id") if current_node in G.nodes else current_road_id
                    current_direction = calculate_direction(segment_path[-2], current_node) if len(segment_path) >= 2 else current_direction
                    print(f"      绝对方向 + 时间: 按绝对方向 {abs_direction}° 移动，距离: {segment_distance:.2f}米, 新节点: {current_node}, 方向: {current_direction:.2f}°")
                continue

            else:
                print(f"      只有绝对方向，但未匹配任何子条件，跳过")
                continue

        # 如果没有匹配任何条件
        else:
            print(f"      警告: 未知指令 {action}，跳过")
            continue

    for point in path:
        if not isinstance(point, tuple) or len(point) != 2:
            print(f"      错误: 最终路径包含无效点 {point}")
            raise ValueError(f"无效路径点: {point}")

    print(f"      指令序列执行完毕，路径长度: {len(path)}, 总距离: {total_distance:.2f}米，最终方向: {current_direction:.2f}°")
    if total_distance == 0 and len(path) <= 1 and not circle_needed:
        print(f"      警告: 仅生成圆形范围，未生成路径点")
    return path, total_distance, current_direction, all_paths, circle_needed, circle_radius

def handle_turn(G, current_node, current_direction, direction, roads_gdf=None, transport="walk"):
    print(f"      执行转弯: {direction}，当前方向: {current_direction:.2f}°，交通方式: {transport}，开始大胆探索，时间: 09:02 PM +08, Saturday, May 17, 2025")

    # 处理相对方向，计算新方向
    new_direction = process_relative_direction(current_direction, direction)
    if new_direction is None:
        print(f"      警告: 无法处理转弯方向 {direction}")
        return current_node, current_direction, [current_node], None

    # 获取当前节点的道路 ID
    current_road_id = G.nodes[current_node].get("road_id") if current_node in G.nodes else None
    if current_road_id is None:
        print(f"      错误: 当前节点 {current_node} 没有 road_id")
        return current_node, current_direction, [current_node], None

    # 初始路径
    path = [current_node]
    total_distance = 0.0
    current_node = tuple(current_node)  # 确保是 tuple
    start_node = current_node  # 记录起点
    visited = {current_node}  # 记录访问过的节点

    # 步骤 1: 检查当前点是否能转弯
    angle_tolerance = 60  # 角度容差 ±60°
    neighbors = list(G.neighbors(current_node))
    if neighbors:
        filtered_neighbors = []
        target_angle = new_direction
        for neighbor in neighbors:
            neighbor_road_id = G.nodes[neighbor].get("road_id")
            if neighbor_road_id == current_road_id:
                continue
            angle = calculate_direction(current_node, neighbor)
            angle_diff = min((angle - target_angle) % 360, (target_angle - angle) % 360)
            if angle_diff <= angle_tolerance:
                filtered_neighbors.append((neighbor, angle_diff, neighbor_road_id))
        
        if filtered_neighbors:
            filtered_neighbors.sort(key=lambda x: x[1])
            next_node, angle_diff, new_road_id = filtered_neighbors[0]
            edge_weight = G[current_node][next_node]['weight']
            path.append(next_node)
            total_distance += edge_weight
            current_node = next_node
            print(f"      在当前点转弯成功，选择新道路 ID: {new_road_id}，方向偏差: {angle_diff:.2f}°，新方向: {new_direction:.2f}°，距离: {total_distance:.2f}米")
            return current_node, new_direction, path, new_road_id

    # 步骤 2: 找到最近的节点（路口）
    path_to_intersection, distance_to_intersection = find_intersection(
        G, current_node, current_road_id, max_distance=500
    )
    if path_to_intersection and len(path_to_intersection) > 1:
        path.extend(path_to_intersection[1:])
        nearest_node = path_to_intersection[-1]
        total_distance += distance_to_intersection
        current_node = nearest_node
        visited.add(current_node)
        print(f"      移动到最近节点 {nearest_node}，距离: {distance_to_intersection:.2f}米")
    else:
        print(f"      当前节点 {current_node} 已在路口或未找到最近节点")

    # 步骤 3: 在最近节点尝试转弯
    neighbors = list(G.neighbors(current_node))
    if neighbors:
        filtered_neighbors = []
        target_angle = new_direction
        for neighbor in neighbors:
            neighbor_road_id = G.nodes[neighbor].get("road_id")
            if neighbor_road_id == current_road_id:
                continue
            angle = calculate_direction(current_node, neighbor)
            angle_diff = min((angle - target_angle) % 360, (target_angle - angle) % 360)
            if angle_diff <= angle_tolerance:
                filtered_neighbors.append((neighbor, angle_diff, neighbor_road_id))
        
        if filtered_neighbors:
            filtered_neighbors.sort(key=lambda x: x[1])
            next_node, angle_diff, new_road_id = filtered_neighbors[0]
            edge_weight = G[current_node][next_node]['weight']
            path.append(next_node)
            total_distance += edge_weight
            current_node = next_node
            print(f"      在最近节点转弯成功，选择新道路 ID: {new_road_id}，方向偏差: {angle_diff:.2f}°，新方向: {new_direction:.2f}°，距离: {total_distance:.2f}米")
            return current_node, new_direction, path, new_road_id

    # 步骤 4: 沿起点到最近节点的方向持续探索
    explore_direction = calculate_direction(start_node, current_node)
    print(f"      最近节点无法转弯，沿起点到最近节点方向 {explore_direction:.2f}° 持续探索")
    while True:
        neighbors = list(G.neighbors(current_node))
        if not neighbors:
            print(f"      警告: 节点 {current_node} 没有邻居，停止搜索")
            break

        # 选择方向最接近 explore_direction 的邻居
        best_neighbor = None
        best_angle_diff = float('inf')
        best_edge_weight = 0
        for neighbor in neighbors:
            angle = calculate_direction(current_node, neighbor)
            angle_diff = min((angle - explore_direction) % 360, (explore_direction - angle) % 360)
            if angle_diff < best_angle_diff:
                best_angle_diff = angle_diff
                best_neighbor = neighbor
                best_edge_weight = G[current_node][neighbor]['weight']

        if best_neighbor is None:
            print(f"      错误: 无法找到任何邻居，停止搜索")
            break

        # 移动到下一个节点
        path.append(best_neighbor)
        total_distance += best_edge_weight
        current_node = best_neighbor

        # 检查是否陷入循环
        if current_node in visited:
            print(f"      沿方向 {explore_direction:.2f}° 探索回到已访问节点 {current_node}，终止探索，累计距离: {total_distance:.2f}米")
            break
        visited.add(current_node)

        # 检查是否能转弯
        neighbors = list(G.neighbors(current_node))
        if neighbors:
            filtered_neighbors = []
            target_angle = new_direction
            for neighbor in neighbors:
                neighbor_road_id = G.nodes[neighbor].get("road_id")
                if neighbor_road_id == current_road_id:
                    continue
                angle = calculate_direction(current_node, neighbor)
                angle_diff = min((angle - target_angle) % 360, (target_angle - angle) % 360)
                if angle_diff <= angle_tolerance:
                    filtered_neighbors.append((neighbor, angle_diff, neighbor_road_id))
            if filtered_neighbors:
                filtered_neighbors.sort(key=lambda x: x[1])
                next_node, angle_diff, new_road_id = filtered_neighbors[0]
                edge_weight = G[current_node][next_node]['weight']
                path.append(next_node)
                total_distance += edge_weight
                current_node = next_node
                print(f"      沿方向 {explore_direction:.2f}° 探索到可转弯点，选择新道路 ID: {new_road_id}，方向偏差: {angle_diff:.2f}°，总距离: {total_distance:.2f}米")
                return current_node, new_direction, path, new_road_id

    # 步骤 5: 沿反方向（起点到最近节点的反向）探索
    print(f"      沿方向 {explore_direction:.2f}° 探索失败，沿反方向 {((explore_direction + 180) % 360):.2f}° 持续探索")
    reverse_direction = (explore_direction + 180) % 360
    current_node = start_node  # 回到起点
    path = [current_node]
    total_distance = 0.0
    visited = {current_node}

    while True:
        neighbors = list(G.neighbors(current_node))
        if not neighbors:
            print(f"      警告: 节点 {current_node} 没有邻居，停止搜索")
            break

        # 选择方向最接近 reverse_direction 的邻居
        best_neighbor = None
        best_angle_diff = float('inf')
        best_edge_weight = 0
        for neighbor in neighbors:
            angle = calculate_direction(current_node, neighbor)
            angle_diff = min((angle - reverse_direction) % 360, (reverse_direction - angle) % 360)
            if angle_diff < best_angle_diff:
                best_angle_diff = angle_diff
                best_neighbor = neighbor
                best_edge_weight = G[current_node][neighbor]['weight']

        if best_neighbor is None:
            print(f"      错误: 无法找到任何邻居，停止搜索")
            break

        # 移动到下一个节点
        path.append(best_neighbor)
        total_distance += best_edge_weight
        current_node = best_neighbor

        # 检查是否陷入循环
        if current_node in visited:
            print(f"      沿反方向 {reverse_direction:.2f}° 探索回到已访问节点 {current_node}，终止探索，累计距离: {total_distance:.2f}米")
            break
        visited.add(current_node)

        # 检查是否能转弯
        neighbors = list(G.neighbors(current_node))
        if neighbors:
            filtered_neighbors = []
            target_angle = new_direction
            for neighbor in neighbors:
                neighbor_road_id = G.nodes[neighbor].get("road_id")
                if neighbor_road_id == current_road_id:
                    continue
                angle = calculate_direction(current_node, neighbor)
                angle_diff = min((angle - target_angle) % 360, (target_angle - angle) % 360)
                if angle_diff <= angle_tolerance:
                    filtered_neighbors.append((neighbor, angle_diff, neighbor_road_id))
            if filtered_neighbors:
                filtered_neighbors.sort(key=lambda x: x[1])
                next_node, angle_diff, new_road_id = filtered_neighbors[0]
                edge_weight = G[current_node][next_node]['weight']
                path.append(next_node)
                total_distance += edge_weight
                current_node = next_node
                print(f"      沿反方向 {reverse_direction:.2f}° 探索到可转弯点，选择新道路 ID: {new_road_id}，方向偏差: {angle_diff:.2f}°，总距离: {total_distance:.2f}米")
                return current_node, new_direction, path, new_road_id

    # 如果所有方向都找不到，返回当前状态
    print(f"      所有方向探索失败，保持原点，方向: {new_direction:.2f}°")
    return current_node, new_direction, path, current_road_id

def calculate_angle(from_node, to_node):
    """计算两点之间的角度（北向为 0 度，顺时针增加）"""
    dx = to_node[0] - from_node[0]
    dy = to_node[1] - from_node[1]
    angle = (90 - math.degrees(math.atan2(dy, dx))) % 360
    return angle

def explore_path(G, current_path, direction, target_distance, visited, 
                last_direction, result_paths, current_distance, max_paths=5):
    """探索满足方向和距离的路径"""
    if current_distance >= target_distance * 0.8:
        result_paths.append(current_path.copy())
        if len(result_paths) >= max_paths:
            return
    
    if current_distance > target_distance * 1.5:
        return
    
    current_node = current_path[-1]
    neighbors = list(G.neighbors(current_node))
    filtered_neighbors = filter_by_direction(current_node, neighbors, direction, last_direction)
    filtered_neighbors = [n for n in filtered_neighbors if n not in visited]
    
    if not filtered_neighbors:
        filtered_neighbors = [n for n in neighbors if n not in visited]
    
    for next_node in filtered_neighbors:
        if next_node in visited:
            continue
            
        current_path.append(next_node)
        visited.add(next_node)
        
        if len(current_path) >= 2:
            new_direction = get_direction_from_angle(calculate_angle(current_path[-2], current_path[-1]))
        else:
            new_direction = last_direction
            
        edge_dist = G[current_path[-2]][current_path[-1]]['weight']
        new_distance = current_distance + edge_dist
        
        explore_path(G, current_path, direction, target_distance, 
                    visited, new_direction, result_paths, new_distance)
        
        visited.remove(next_node)
        current_path.pop()


def convert_and_print_coordinates(start_name):
    """转换并打印坐标（GCJ-02 到 WGS-84）"""
    start_coord_gcj = get_amap_coordinates(start_name)
    if not start_coord_gcj:
        print(f"无法获取起点 '{start_name}' 的坐标")
        return None
    print(f"高德地图返回的起点坐标 (GCJ-02): {start_coord_gcj}")
    
    start_wgs = gcj02_to_wgs84(*start_coord_gcj)
    print(f"转换后的起点坐标 (WGS-84): {start_wgs}")
    return start_wgs
def compute_full_path(data, shapefile_path):
    """计算完整路径"""
    try:
        G, roads_gdf = build_graph(shapefile_path)
        print(f"成功构建路网图: {len(G.nodes)} 节点, {len(G.edges)} 边")
        
        start_name = data["start"]["name"]
        print(f"查询起点: {start_name}")
        
        start_coord_gcj = get_amap_coordinates(start_name)
        if not start_coord_gcj:
            raise HTTPException(status_code=404, detail=f"无法获取起点 '{start_name}' 的坐标")
        
        start_wgs = gcj02_to_wgs84(*start_coord_gcj)
        print(f"起点坐标 (WGS84): {start_wgs}")
        current_node = insert_point_on_edge(G, start_wgs)

        full_coords = []
        instructions = []
        current_direction = None  # 初始化方向

        for seg_id, seg in data["segments"].items():
            print(f"\n开始处理 {seg_id}...")
            try:
                path_instr = seg.get("path_instructions", [])
                if path_instr:
                    print(f"  发现路径指令: {path_instr}")
                    # 调用 advance_segment，返回路径、距离和方向
                    segment_path, segment_distance, current_direction = advance_segment(
                        G, current_node, path_instr, prev_direction=current_direction, roads_gdf=roads_gdf
                    )
                    # 验证路径点
                    for point in segment_path:
                        if not isinstance(point, tuple) or len(point) != 2:
                            print(f"错误: 段 {seg_id} 的路径点 {point} 无效")
                            raise ValueError(f"无效路径点: {point}")
                    if not segment_path or len(segment_path) <= 1:
                        print(f"警告: {seg_id} 未生成有效路径")
                        continue
                    # 合并路径，排除重复的起点
                    full_coords.extend(segment_path if not full_coords else segment_path[1:])
                    current_node = segment_path[-1]
                    print(f"更新 current_node: {current_node}, 是否在图中: {current_node in G.nodes}")

                    # 生成指令描述
                    parts = []
                    for step in path_instr:
                        step_desc = []
                        if "action" in step: step_desc.append(step["action"])
                        if "direction" in step: step_desc.append(step["direction"])
                        if "distance" in step: step_desc.append(step["distance"])
                        if "time" in step: step_desc.append(f"{step['time']}分钟")
                        parts.append(" ".join(step_desc))
                    text = " → ".join(parts)
                    instructions.append({
                        "segment": seg_id,
                        "description": text,
                        "point_count": len(segment_path)
                    })
                    print(f"完成路径指令，路径点数: {len(segment_path)}")
                else:
                    print(f"{seg_id} 没有路径指令")
            except Exception as e:
                print(f"处理 {seg_id} 时出错: {e}")
                raise e
        
        if not full_coords or len(full_coords) < 2:
            raise Exception("未能生成有效的路径点")
            
        print(f"完成路径生成，总路径点数: {len(full_coords)}")
        return {
            "geojson": path_to_geojson(full_coords),
            "instructions": instructions,
            "path": full_coords  # 显式返回路径点列表
        }
    except Exception as e:
        print(f"路径计算过程中出错: {e}")
        raise e