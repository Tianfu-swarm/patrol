#include <vector>
#include <iostream>
#include "ros/ros.h"
#include <patrol/Patrol_area.h>
#include <patrol/Obstacle_area.h>
#include "patrol/Adjacency_matrix.h"
#include <visualization_msgs/Marker.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/Voronoi_diagram_2/Halfedge.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_2<K> DT;
typedef CGAL::Delaunay_triangulation_adaptation_traits_2<DT> AT;
typedef CGAL::Delaunay_triangulation_caching_degeneracy_removal_policy_2<DT> AP;
typedef CGAL::Voronoi_diagram_2<DT, AT, AP> VoronoiDiagram;
typedef VoronoiDiagram::Face_handle Face_handle;
typedef AT::Point_2 Point_2;

ros::Publisher mapEdge_pub;
ros::Publisher obstacleEdge_pub;
ros::Publisher voronoiEdge_pub;
ros::Publisher seedPoint_pub;

//创建种子点
std::vector<std::pair<double, double>> seedPoints;

// 计算叉积
long double cross_product(std::pair<double, double> p1, std::pair<double, double> p2)
{
    return p1.first * p2.second - p1.second * p2.first;
}

// 辅助函数，判断点q是否在线段pr上
bool on_segment(std::pair<double, double> p, std::pair<double, double> q, std::pair<double, double> r) {
    return (q.first <= std::max(p.first, r.first) && q.first >= std::min(p.first, r.first) &&
            q.second <= std::max(p.second, r.second) && q.second >= std::min(p.second, r.second));
}

// 判断两个线段是否相交
bool segments_intersect(std::pair<std::pair<double, double>, std::pair<double, double>> seg1,
                        std::pair<std::pair<double, double>, std::pair<double, double>> seg2) {
    const double EPSILON = 1e-15;
    auto p1 = seg1.first, p2 = seg1.second, p3 = seg2.first, p4 = seg2.second;

    // 计算四个叉积
    double d1 = cross_product({p3.first - p1.first, p3.second - p1.second},
                              {p4.first - p1.first, p4.second - p1.second});
    double d2 = cross_product({p3.first - p2.first, p3.second - p2.second},
                              {p4.first - p2.first, p4.second - p2.second});
    double d3 = cross_product({p1.first - p3.first, p1.second - p3.second},
                              {p2.first - p3.first, p2.second - p3.second});
    double d4 = cross_product({p1.first - p4.first, p1.second - p4.second},
                              {p2.first - p4.first, p2.second - p4.second});

    // 一般情况
    if (d1 * d2 < -EPSILON && d3 * d4 < -EPSILON) {
        return true;
    }

    // 检查共线情况下的点是否在线段上
    if (fabs(d1) < EPSILON && on_segment(p1, p3, p2)) return true;
    if (fabs(d2) < EPSILON && on_segment(p1, p4, p2)) return true;
    if (fabs(d3) < EPSILON && on_segment(p3, p1, p4)) return true;
    if (fabs(d4) < EPSILON && on_segment(p3, p2, p4)) return true;

    return false;
}

// 点到线段的距离
double distancePointToLineSegment(std::pair<double, double> &point,
                                  std::pair<std::pair<double, double>, std::pair<double, double>> lineSegment)
{
    // 提取点的坐标和线段端点的坐标
    auto &[px, py] = point;
    auto &[A, B] = lineSegment;
    auto &[ax, ay] = A;
    auto &[bx, by] = B;

    // 计算向量 AB
    double ABx = bx - ax;
    double ABy = by - ay;

    // 计算向量 AP
    double APx = px - ax;
    double APy = py - ay;

    // 计算向量 AP 和 AB 的点积
    double AB_AB = ABx * ABx + ABy * ABy;
    double AP_AB = APx * ABx + APy * ABy;

    // 计算点 P 在线段 AB 上的投影参数 t
    double t = AP_AB / AB_AB;

    if (t < 0.0)
    {
        // 如果 t < 0，最近点是 A
        return std::sqrt(APx * APx + APy * APy);
    }
    else if (t > 1.0)
    {
        // 如果 t > 1，最近点是 B
        double BPx = px - bx;
        double BPy = py - by;
        return std::sqrt(BPx * BPx + BPy * BPy);
    }
    else
    {
        // 如果 0 <= t <= 1，最近点在线段上
        double closestX = ax + t * ABx;
        double closestY = ay + t * ABy;
        double dx = px - closestX;
        double dy = py - closestY;
        return std::sqrt(dx * dx + dy * dy);
    }
}


// 计算两个线段的交点
std::pair<double, double> intersection_point(
    std::pair<std::pair<double, double>, std::pair<double, double>> seg1,
    std::pair<std::pair<double, double>, std::pair<double, double>> seg2)
{
    auto p1 = seg1.first, p2 = seg1.second, p3 = seg2.first, p4 = seg2.second;
    double A1 = p2.second - p1.second;
    double B1 = p1.first - p2.first;
    double C1 = A1 * p1.first + B1 * p1.second;

    double A2 = p4.second - p3.second;
    double B2 = p3.first - p4.first;
    double C2 = A2 * p3.first + B2 * p3.second;

    double det = A1 * B2 - A2 * B1;

    double x = (B2 * C1 - B1 * C2) / det;
    double y = (A1 * C2 - A2 * C1) / det;
    return {x, y};
}

ros::Publisher adjacencyMatrix_pub;
// 计算adjacency_matrix
void generate_adjacencyMatrix(std::vector<std::pair<double, double>> &points,
                              std::vector<std::pair<std::pair<double, double>, std::pair<double, double>>> &obstacleEdges,
                              std::vector<std::pair<std::pair<double, double>, std::pair<double, double>>> voronoiEdges)
{

    patrol::Adjacency_matrix Adjacency_matrix;
    
    int matrixDimension = seedPoints.size();
    std::vector<std::vector<int>> adjacency_matrix(matrixDimension, std::vector<int>(matrixDimension, 0));
    // 将相邻的point的紧邻矩阵设置为1，不相邻设置为0
    for (int i = 0; i < points.size() - 1; i++)
    {
        std::pair<std::pair<double, double>, std::pair<double, double>> seedSeg;
        for (int j = 1; j < points.size(); j++)
        {
            seedSeg = {points[i], points[j]};
            int num_intersect = 0;
            int num_intersect_obstacle = 0;
            for (const auto &seg : voronoiEdges)
            {
                if (segments_intersect(seedSeg, seg))
                {
                    num_intersect++;
                }
            }
            for (const auto &seg2 : obstacleEdges)
            {
                if (segments_intersect(seedSeg, seg2))
                {
                    num_intersect_obstacle++;
                }
            }
            if (num_intersect == 1 && num_intersect_obstacle == 0)
            {
                adjacency_matrix[i][j] = 1;
                adjacency_matrix[j][i] = 1;
            }
        }
    }

    //将points与adjacency_matrix转为patrol：：Adjacency_matrix格式并发布
    for(auto &point : points){
        Adjacency_matrix.x.push_back(point.first);
        Adjacency_matrix.y.push_back(point.second);
    }

    for (int i = 0; i < matrixDimension; ++i)
    {
        for (int j = 0; j < matrixDimension; ++j)
        {
            Adjacency_matrix.matrix.push_back(adjacency_matrix[i][j]);
        }
    }
    adjacencyMatrix_pub.publish(Adjacency_matrix);
}

//找到对应点被哪个最小的封闭多边形包围
std::vector<std::pair<std::pair<double, double>, std::pair<double, double>>> closedPolygon(std::pair<double, double> &points,
                                                                                           std::pair<std::pair<double, double>, std::pair<double, double>> &lineSegment,
                                                                                           std::vector<std::pair<std::pair<double, double>, std::pair<double, double>>> &voronoiEdges)
{

     double endX;
     double endY;
     std::vector<std::pair<std::pair<double, double>, std::pair<double, double>>> generateRays;
     double radius = 1000.0;
     int numRays = static_cast<int>(360 / 0.5);
     for (int i = 0; i < numRays; ++i)
     {
         double angleDeg = i * 0.5;
         double angleRad = angleDeg * M_PI / 180.0;
         endX= points.first + radius * cos(angleRad);
         endY= points.second + radius * sin(angleRad);
         generateRays.push_back({points,{endX,endY}});
     }
    std::vector<std::pair<std::pair<double, double>, std::pair<double, double>>> polygons;
    std::pair<std::pair<double, double>, std::pair<double, double>> polygon;
    for(const auto&ray :generateRays){
        double minDistance = 1e9;
        for(const auto&edge:voronoiEdges){
            if(segments_intersect(ray,edge)){
                std::pair<double, double> inter = intersection_point(ray,edge);
                double dis = std::sqrt((inter.first-points.first)*(inter.first-points.first)+(inter.second-points.second)*(inter.second-points.second));
                if (dis < minDistance)
                {
                    minDistance = dis;
                    polygon = edge;
                }
            }
            
        }
        polygons.push_back(polygon);
    }
    //去除polygos中重复的边
    // 自定义比较器，用于 std::set
    auto compareEdges = [](const std::pair<std::pair<double, double>, std::pair<double, double>> &a,
                           const std::pair<std::pair<double, double>, std::pair<double, double>> &b)
    {
        auto min_a = std::min(a.first, a.second);
        auto max_a = std::max(a.first, a.second);
        auto min_b = std::min(b.first, b.second);
        auto max_b = std::max(b.first, b.second);
        return std::tie(min_a, max_a) < std::tie(min_b, max_b);
    };

    // 使用 std::set 去除重复边
    std::set<std::pair<std::pair<double, double>, std::pair<double, double>>, decltype(compareEdges)> uniqueEdges(compareEdges);

    for (const auto &edge : polygons)
    {
        auto min_edge = std::min(edge.first, edge.second);
        auto max_edge = std::max(edge.first, edge.second);
        uniqueEdges.emplace(min_edge, max_edge);
    }

    // 将集合中的唯一边转换回向量
    std::vector<std::pair<std::pair<double, double>, std::pair<double, double>>> closePolygons(uniqueEdges.begin(), uniqueEdges.end());

    return closePolygons;
}

// 计算每个seed对应的voronoi边
std::vector<std::pair<double, double>> lloyd(std::vector<std::pair<double, double>> &points,
                                             std::vector<std::pair<std::pair<double, double>, std::pair<double, double>>> voronoiEdges)
{
    std::vector<std::pair<double, double>> point = points;
    std::vector<std::pair<double, double>> newpoint;
    for (int i = 0; i < point.size(); i++)
    {
        //找到point距离最近的边
        
        double minDis = 1e9;
        int numOfEdge;
        for (int j = 0; j < voronoiEdges.size(); j++)
        {
            double dis = distancePointToLineSegment(point[i], voronoiEdges[j]);
            if (dis < minDis)
            {
                minDis = dis;
                numOfEdge = j;
            }
        }
        //找到points在哪个由voronoi边围成的多边形中
        std::vector<std::pair<std::pair<double, double>, std::pair<double, double>>> polygons;
        polygons = closedPolygon(point[i],voronoiEdges[numOfEdge],voronoiEdges);
        double newPointX = 0;
        double newPointY = 0;
        int num = 0;
        for (const auto &polygon : polygons)
        {
            newPointX = newPointX + polygon.first.first + polygon.second.first;
            newPointY = newPointY + polygon.first.second + polygon.second.second;
            num++;
            //ROS_WARN("(%f,%f),(%f,%f)",polygon.first.first,polygon.first.first,polygon.second.first,polygon.second.second);
        }
        newpoint.push_back({newPointX / num / 2,newPointY / num / 2});

    }
    return newpoint;
}

//计算voronoi的边
std::vector<std::pair<std::pair<double, double>, std::pair<double, double>>> voronoiEdges;
void generate_voronoi(std::vector<std::pair<double, double>> &points,
                      std::vector<std::pair<std::pair<double, double>, std::pair<double, double>>> &mapEdges,
                      std::vector<std::pair<std::pair<double, double>, std::pair<double, double>>> &obstacleEdges,
                      std::vector<std::vector<std::pair<std::pair<double, double>, std::pair<double, double>>>> &obstacleGroup)
{

    // 将输入点转换为CGAL点类型
    std::vector<K::Point_2> seeds;
    for (const auto &p : points)
    {
        seeds.push_back(K::Point_2(p.first, p.second));
    }

    // 使用seeds构造Delaunay三角剖分和Voronoi图
    DT dt(seeds.begin(), seeds.end());
    VoronoiDiagram vd(dt);

    //提取所有voronoi边
    
    VoronoiDiagram::Edge_iterator it;
    for (it = vd.edges_begin(); it != vd.edges_end(); ++it)
    {
        // 有界边
        if (it->has_source() && it->has_target())
        {
            auto source = it->source()->point();
            auto target = it->target()->point();
            voronoiEdges.push_back({{source.x(), source.y()},
                                    {target.x(), target.y()}});
        }
        //无界边
        if (!it->has_source() || !it->has_target())
        {
            // 无界边的源点
            auto source = it->has_source()
                              ? it->source()->point()
                              : it->target()->point();

            // 原始的Delaunay边
            auto dual_edge = it->dual();

            // 对应Delaunay边的两个顶点
            auto dual_source = dual_edge.first->vertex((dual_edge.second + 1) % 3)->point();
            auto dual_target = dual_edge.first->vertex((dual_edge.second + 2) % 3)->point();

            // 使用map的质心代替内部
            double mx = 0;
            double my = 0;
            int numMapPoint = 0;
            for (const auto &map : mapEdges)
            {
                mx = mx + map.first.first + map.second.first;
                my = my + map.first.second + map.second.second;
                numMapPoint++;
            }
            mx = mx / numMapPoint / 2;
            my = my / numMapPoint / 2;

            // 计算Delaunay边的中点
            double midpoint_x = (dual_source.x() + dual_target.x()) / 2.0;
            double midpoint_y = (dual_source.y() + dual_target.y()) / 2.0;
            // 计算从Delaunay边的中点指向质心的方向
            auto midDirection = std::make_pair(midpoint_x - mx, midpoint_y - my);
            // 计算Delaunay边的方向，垂直于Delaunay边，并从midDirection获取起始方向
            auto direction = std::make_pair(-(dual_target.y() - dual_source.y()),
                                            dual_target.x() - dual_source.x());

            // 计算direction向量和midDirection向量的点积
            double dot_product = direction.first * midDirection.first + direction.second * midDirection.second;

            // 如果点积小于零，说明两向量夹角小于90度，因此需要调整direction向量的方向
            // 注意：我们无需判断点积是否大于零，因为如果两个向量夹角大于90度，我们希望保持当前的direction向量方向
            if (dot_product < 0)
            {
                direction.first = -direction.first;
                direction.second = -direction.second;
            }

            // 计算无界边的目标点
            double fixed_length = 1000.0; // 设定一个最远距离
            auto target_point = std::make_pair(source.x() + fixed_length * direction.first,
                                               source.y() + fixed_length * direction.second);

            voronoiEdges.push_back({{source.x(), source.y()}, {target_point.first, target_point.second}});
        }
    }

    //存储在地图与障碍物边界上的点
    std::vector<std::pair<double, double>> newPoints;
    //新生成的边
    std::vector<std::pair<std::pair<double, double>, std::pair<double, double>>> shouldAddEdge;
    std::vector<int> shouldRemove;
    int numofinster = 0;
    //处理相交
    for (size_t i = 0; i < voronoiEdges.size(); ++i)
    {
        auto &edge = voronoiEdges[i];

        int flag1 = 0;
        int flag2 = 0;
        double minDistance = 1e9;
        // 存储map最近的交点
        std::pair<double, double> minIntersectPoint;
        //存储obstacle最近的交点
        std::pair<double, double> minIntersectPointObstacle;

        // 删除在map范围外的边
        // map的质心
        double mx = 0;
        double my = 0;
        int numMapPoint = 0;
        for (const auto &map : mapEdges)
        {
            mx = mx + map.first.first + map.second.first;
            my = my + map.first.second + map.second.second;
            numMapPoint++;
        }
        mx = mx / numMapPoint / 2;
        my = my / numMapPoint / 2;
        std::pair<std::pair<double, double>, std::pair<double, double>> seg1 = {edge.first, {mx, my}};
        std::pair<std::pair<double, double>, std::pair<double, double>> seg2 = {edge.second, {mx, my}};
        bool inter1 = false;
        bool inter2 = false;
        for (const auto &mapedge : mapEdges)
        {
            if (segments_intersect(seg1, mapedge))
            {
                inter1 = true;
            }
            if (segments_intersect(seg2, mapedge))
            {
                inter2 = true;
            }
        }
        if (inter1 == true && inter2 == true)
        {
            shouldRemove.push_back(i);
            continue;
        }
        //删除在obstacle内的边
        for (const auto &group : obstacleGroup)
        {
            double mx = 0;
            double my = 0;
            int numMapPoint = 0;
            for (const auto &map : group)
            {
                mx = mx + map.first.first + map.second.first;
                my = my + map.first.second + map.second.second;
                numMapPoint++;
            }
            mx = mx / numMapPoint / 2;
            my = my / numMapPoint / 2;
            std::pair<std::pair<double, double>, std::pair<double, double>> seg1 = {edge.first, {mx, my}};
            std::pair<std::pair<double, double>, std::pair<double, double>> seg2 = {edge.second, {mx, my}};
            bool inter1 = false;
            bool inter2 = false;
            for (const auto &edge : group)
            {
                if (segments_intersect(seg1, edge))
                {
                    inter1 = true;
                    continue;
                }
                if (segments_intersect(seg2, edge))
                {
                    inter2 = true;
                    continue;
                }
            }
            if (inter1 == false && inter2 == false)
            {
                shouldRemove.push_back(i);
                continue;
            }
        }

        //判断与地图边界的交点
        for (const auto &mapEdge : mapEdges)
        {
            if (segments_intersect(edge, mapEdge))
            {
                std::pair<double, double> intersectPoint = intersection_point(edge, mapEdge);
                double dx = edge.first.first - intersectPoint.first;
                double dy = edge.first.second - intersectPoint.second;
                double distance = std::sqrt(dx * dx + dy * dy);
                if (distance < minDistance)
                {
                    minIntersectPoint = intersectPoint;
                }
                flag1 = 1;
            }
        }

        // 判断与障碍物边界的交点
        int numObstacleEdge;
        int numObstacleGroup;
        // 存储所有交点
        std::vector<std::pair<double, double>> intersectPointObstacle;
        for (int k = 0; k < obstacleGroup.size(); k++)
        {
            std::vector<std::pair<std::pair<double, double>, std::pair<double, double>>> obstacle = obstacleGroup[k];
            for (int j = 0; j < obstacle.size(); j++)
            {
                auto &obstacleEdge = obstacle[j];
                if (segments_intersect(edge, obstacleEdge))
                {
                    intersectPointObstacle.push_back(intersection_point(edge, obstacleEdge));
                    numObstacleEdge = j;
                    numObstacleGroup = k;
                }
            }
        }
        if (intersectPointObstacle.size() > 1)
        {
            flag2 = 1;
        }
        else if (intersectPointObstacle.size() == 1)
        {
            flag2 = 2;
        }

        // 如果出现交点在map边
        if (flag1 == 1)
        {
            newPoints.push_back(minIntersectPoint);
            
            // 将边的终点修改为交点
            bool inter1 = false;
            bool inter2 = false;
            for (const auto &edge : mapEdges)
            {
                if (segments_intersect(seg1, edge))
                {
                    inter1 = true;
                }
                if (segments_intersect(seg2, edge))
                {
                    inter2 = true;
                }
            }
            if (inter1 == false)
            {
                edge.second.first = minIntersectPoint.first;
                edge.second.second = minIntersectPoint.second;
            }
            if (inter2 == false)
            {
                edge.first.first = minIntersectPoint.first;
                edge.first.second = minIntersectPoint.second;
            }
        }

        //如果出现交点在obstacle边且有多个交点
        if (flag2 ==1){

            for (int x; x < intersectPointObstacle.size() - 1; x++)
            {
                // 将与障碍物的所有交点依次连接为边，并加入voronoi
                shouldAddEdge.push_back({intersectPointObstacle[x], intersectPointObstacle[x + 1]});
            }

            for (const auto &point : intersectPointObstacle)
            {
                newPoints.push_back(point);
            }

            //删除存放的交点
            intersectPointObstacle.clear();
            // 删除当前边
            shouldRemove.push_back(i);
            
        }
        //如果出现交点在obstacle边且只有一个交点
        if(flag2 ==2){

            minIntersectPoint.first = intersectPointObstacle[0].first;
            minIntersectPoint.second = intersectPointObstacle[0].second;
            newPoints.push_back(minIntersectPoint);

            double mx = 0;
            double my = 0;
            int numObstaclePoint = 0;
    
            for (const auto &obstacle : obstacleGroup[numObstacleGroup])
            {
                mx = mx + obstacle.first.first + obstacle.second.first;
                my = my + obstacle.first.second + obstacle.second.second;
                numObstaclePoint++;
            }
            mx = mx / numObstaclePoint / 2;
            my = my / numObstaclePoint / 2;
            std::pair<std::pair<double, double>, std::pair<double, double>> seg3 = {edge.first, {mx, my}};
            std::pair<std::pair<double, double>, std::pair<double, double>> seg4 = {edge.second, {mx, my}};

            // 将边的终点修改为交点
            bool inter1 = false;
            bool inter2 = false;
            for (const auto &mapedge : obstacleGroup[numObstacleGroup])
            {
                if (segments_intersect(seg3, mapedge))
                {
                    inter1 = true;
                }
                if (segments_intersect(seg4, mapedge))
                {
                    inter2 = true;
                }
            }

            if (inter2 == false)
            {
                edge.second.first = minIntersectPoint.first;
                edge.second.second = minIntersectPoint.second;
            }
            if (inter1 == false)
            {
                edge.first.first = minIntersectPoint.first;
                edge.first.second = minIntersectPoint.second;
            }
        }
    }

    //从voronoiEdges中移除shouldRemove对应下标的点从最大的索引开始删除
    std::sort(shouldRemove.rbegin(), shouldRemove.rend());

    // 从后向前遍历 shouldRemove，并删除对应的元素
    for (auto i : shouldRemove)
    {
        voronoiEdges.erase(voronoiEdges.begin() + i);
    }

    // 添加新增的边
    for (const auto &edge : shouldAddEdge)
    {
        voronoiEdges.push_back(edge);
    }

    // 处理地图/障碍物边界的交点生成的新边
    //与地图边界的新边
    for (const auto &mapEdge : mapEdges)
    {
        // 定义在同一条边的所有点
        std::vector<std::pair<double, double>> pointOnEdge;
        // 找到在同一条边界边的所有点
        for (const auto &point : newPoints)
        {
            // 计算点到线段两个端点的距离，如果距离和与线段长度相等则表示该点在线段上
            //  计算点到两个端点的距离
            double distance1 = std::sqrt((point.first - mapEdge.first.first) * (point.first - mapEdge.first.first) + (point.second - mapEdge.first.second) * (point.second - mapEdge.first.second));
            double distance2 = std::sqrt((point.first - mapEdge.second.first) * (point.first - mapEdge.second.first) + (point.second - mapEdge.second.second) * (point.second - mapEdge.second.second));

            // 计算线段长度
            double segment_length = std::sqrt((mapEdge.first.first - mapEdge.second.first) * (mapEdge.first.first - mapEdge.second.first) + (mapEdge.first.second - mapEdge.second.second) * (mapEdge.first.second - mapEdge.second.second));

            // 如果点到两个端点的距离之和等于线段长度，则点在线段上
            if (std::fabs(distance1 + distance2 - segment_length) < 1e-9)
            {
                pointOnEdge.push_back(point);
            }
        }
        //判断是否有点
        if (!pointOnEdge.empty())
        {
            // 将pointOnEdge中的点按照距离mapEdge.first的距离重新排列
            std::sort(pointOnEdge.begin(), pointOnEdge.end(), [&](const std::pair<double, double> &p1, const std::pair<double, double> &p2)
                      {    
                        double distance1 = std::sqrt((p1.first - mapEdge.first.first) * (p1.first - mapEdge.first.first) + (p1.second - mapEdge.first.second) * (p1.second - mapEdge.first.second));
                        double distance2 = std::sqrt((p2.first - mapEdge.first.first) * (p2.first - mapEdge.first.first) + (p2.second - mapEdge.first.second) * (p2.second - mapEdge.first.second));
                        return distance1 < distance2; });

            // 将同一条线上的点组合voronoi边
            if (pointOnEdge.size() > 1)
            {
                voronoiEdges.push_back({mapEdge.first, pointOnEdge[0]});
                for (int i = 0; i < pointOnEdge.size() - 1; i++)
                {
                    voronoiEdges.push_back({pointOnEdge[i], pointOnEdge[i + 1]});
                }
                int i = pointOnEdge.size() - 1;
                voronoiEdges.push_back({pointOnEdge[i], mapEdge.second});
            }
            else
            {
                voronoiEdges.push_back({mapEdge.first, pointOnEdge[0]});
                voronoiEdges.push_back({pointOnEdge[0], mapEdge.second});
            }
        }
        //清空pointOnEdge
        pointOnEdge.clear();
    }
    // 与障碍物边界的新边
    for (const auto &obstacleEdge : obstacleEdges)
    {
        // 定义在同一条边的所有点
        std::vector<std::pair<double, double>> pointOnEdge;
        // 找到在同一条边界边的所有点
        for (const auto &point : newPoints)
        {
            // 计算点到线段两个端点的距离，如果距离和与线段长度相等则表示该点在线段上
            //  计算点到两个端点的距离
            double distance1 = std::sqrt((point.first - obstacleEdge.first.first) * (point.first - obstacleEdge.first.first) + (point.second - obstacleEdge.first.second) * (point.second - obstacleEdge.first.second));
            double distance2 = std::sqrt((point.first - obstacleEdge.second.first) * (point.first - obstacleEdge.second.first) + (point.second - obstacleEdge.second.second) * (point.second - obstacleEdge.second.second));

            // 计算线段长度
            double segment_length = std::sqrt((obstacleEdge.first.first - obstacleEdge.second.first) * (obstacleEdge.first.first - obstacleEdge.second.first) + (obstacleEdge.first.second - obstacleEdge.second.second) * (obstacleEdge.first.second - obstacleEdge.second.second));

            // 如果点到两个端点的距离之和等于线段长度，则点在线段上
            if (std::fabs(distance1 + distance2 - segment_length) < 1e-9)
            {
                pointOnEdge.push_back(point);
            }
        }
        //判断是否有点
        if (!pointOnEdge.empty())
        {
            // 将pointOnEdge中的点按照距离mapEdge.first的距离重新排列
            std::sort(pointOnEdge.begin(), pointOnEdge.end(), [&](const std::pair<double, double> &p1, const std::pair<double, double> &p2)
                      {    
                        double distance1 = std::sqrt((p1.first - obstacleEdge.first.first) * (p1.first - obstacleEdge.first.first) + (p1.second - obstacleEdge.first.second) * (p1.second - obstacleEdge.first.second));
                        double distance2 = std::sqrt((p2.first - obstacleEdge.first.first) * (p2.first - obstacleEdge.first.first) + (p2.second - obstacleEdge.first.second) * (p2.second - obstacleEdge.first.second));
                        return distance1 < distance2; });

            // 将同一条线上的点组合voronoi边
            if (pointOnEdge.size() > 1)
            {
                voronoiEdges.push_back({obstacleEdge.first, pointOnEdge[0]});
                for (int i = 0; i < pointOnEdge.size() - 1; i++)
                {
                    voronoiEdges.push_back({pointOnEdge[i], pointOnEdge[i + 1]});
                }
                int i = pointOnEdge.size() - 1;
                voronoiEdges.push_back({pointOnEdge[i], obstacleEdge.second});
            }
            else
            {
                voronoiEdges.push_back({obstacleEdge.first, pointOnEdge[0]});
                voronoiEdges.push_back({pointOnEdge[0], obstacleEdge.second});
            }
        }
        //清空pointOnEdge
        pointOnEdge.clear();
    }
}

double calculateTriangleArea(const std::pair<double, double> &p1, const std::pair<double, double> &p2, const std::pair<double, double> &p3)
{
    // 计算三角形的面积
    return 0.5 * std::abs((p1.first * (p2.second - p3.second) + p2.first * (p3.second - p1.second) + p3.first * (p1.second - p2.second)));
}
// 计算区域面积
double calculatePolygonArea(const std::vector<std::pair<double, double>> &points)
{
    int numPoints = points.size();
    double area = 0.0;

    // 判断是否有足够的点构成多边形
    if (numPoints < 3)
    {
        return 0.0; // 无法构成多边形，面积为0
    }

    // 遍历顶点，划分成三角形并计算面积
    for (int i = 1; i < numPoints - 1; ++i)
    {
        area += calculateTriangleArea(points[0], points[i], points[i + 1]);
    }

    return area;
}

double patrolPointsDensity;
// 生成种子点
std::vector<std::pair<double, double>> generatePoint(double &patrolArea,
                                                     std::vector<std::pair<std::pair<double, double>, std::pair<double, double>>> &mapEdges,
                                                     std::vector<std::vector<std::pair<std::pair<double, double>, std::pair<double, double>>>> obstacle)
{
    std::vector<std::pair<double, double>> points;
    // 计算应该生成的点的数量,每10m^2 设置一个巡逻点
    int numPointsToGenerate = static_cast<int>(patrolArea / 10);
    // 生成在map内接圆内，且不在obstacle外接园内的点
    std::random_device rd;
    std::mt19937 gen(rd());
    while (!(points.size() >= numPointsToGenerate))
    {
        // 计算内接圆的圆心与直径
        double cx = 0.0;
        double cy = 0.0;
        double area = 0.0;

        for (size_t i = 0; i < mapEdges.size(); ++i)
        {
            double x1 = mapEdges[i].first.first;
            double y1 = mapEdges[i].first.second;
            double x2 = mapEdges[(i + 1) % mapEdges.size()].first.first;
            double y2 = mapEdges[(i + 1) % mapEdges.size()].first.second;

            double crossProduct = (x1 * y2 - x2 * y1);
            area += crossProduct;

            cx += (x1 + x2) * crossProduct;
            cy += (y1 + y2) * crossProduct;
        }

        area *= 0.5;
        cx /= (6.0 * area);
        cy /= (6.0 * area);
        // 圆心用重心代替
        std::pair<double, double> centroid = {cx, cy};
        
        ROS_WARN("%f,%f",cx,cy);
        // 计算多边形的最近顶点到重心的距离
        double minDistanceSquared = 1000000.0;
        for (const auto &edge : mapEdges)
        {
            double dx = centroid.first - edge.first.first;
            double dy = centroid.second - edge.first.second;
            double distanceSquared = dx * dx + dy * dy;
            minDistanceSquared = std::min(minDistanceSquared, distanceSquared);
        }

        // 计算半径（最近顶点到重心的距离）
        double radius = std::sqrt(minDistanceSquared)/1.5;
        
        // 生成x坐标
        double min = cx - radius;
        double max = cx + radius;
        std::uniform_real_distribution<double> dis1(min, max);
        double dx = dis1(gen);
        // 生成y坐标
        min = cy - radius;
        max = cy + radius;
        std::uniform_real_distribution<double> dis2(min, max);
        double dy = dis2(gen);

        bool isInobstacle = false;
        // 计算外接接圆的圆心与直径
        for (const auto &obstacle : obstacle)
        {
            double cx = 0.0;
            double cy = 0.0;
            double area = 0.0;

            for (size_t i = 0; i < obstacle.size(); ++i)
            {
                double x1 = obstacle[i].first.first;
                double y1 = obstacle[i].first.second;
                double x2 = obstacle[(i + 1) % obstacle.size()].first.first;
                double y2 = obstacle[(i + 1) % obstacle.size()].first.second;

                double crossProduct = (x1 * y2 - x2 * y1);
                area += crossProduct;

                cx += (x1 + x2) * crossProduct;
                cy += (y1 + y2) * crossProduct;
            }

            area *= 0.5;
            cx /= (6.0 * area);
            cy /= (6.0 * area);
            // 圆心用重心代替
            std::pair<double, double> centroid = {cx, cy};
            

            // 计算多边形的最远顶点到重心的距离
            double maxDistanceSquared = 0.0;
            for (const auto &edge : obstacle)
            {
                double dx = centroid.first - edge.first.first;
                double dy = centroid.second - edge.first.second;
                double distanceSquared = dx * dx + dy * dy;
                maxDistanceSquared = std::max(maxDistanceSquared, distanceSquared);
            }

            // 半径 （最远顶点到重心的距离）
            double radius = std::sqrt(maxDistanceSquared);

            // 判断（dx,dy）是否在外接圆内
            double distance = std::sqrt((dx - centroid.first) * (dx - centroid.first) + (dy - centroid.second) * (dy - centroid.second));
            if (distance < radius)
            {
                isInobstacle = true;
            }
        }
        bool isTheSame = false;
        for (const auto &point : points)
        {
            if (point.first == dx && point.second == dy)
            {
                isTheSame = true;
            }
        }
        if (isInobstacle == false && isTheSame == false)
        {
            points.push_back({dx, dy});
        }
    }

    return points;
}

int numOfIteration = 100;

double patrolArea = 0.0;
// 创建障碍物边界线段
std::vector<std::pair<std::pair<double, double>, std::pair<double, double>>> obstacleEdges;
// 多个障碍物放进一个vector
std::vector<std::vector<std::pair<std::pair<double, double>, std::pair<double, double>>>> obstacleGroup;

// 创建地图边界线段
std::vector<std::pair<std::pair<double, double>, std::pair<double, double>>> mapEdges;

void patrolAreaCallback(const patrol::Patrol_area::ConstPtr &msg)
{
    // 初始化 mapEdge
    visualization_msgs::Marker mapEdge;
    mapEdge.header.frame_id = "world";
    mapEdge.ns = "edges";
    mapEdge.id = 0;
    mapEdge.type = visualization_msgs::Marker::LINE_LIST;
    mapEdge.action = visualization_msgs::Marker::ADD;
    mapEdge.pose.orientation.w = 1.0;
    mapEdge.scale.x = 0.05;
    mapEdge.color.a = 1.0;
    mapEdge.color.r = 0.0;
    mapEdge.color.g = 0.0;
    mapEdge.color.b = 1.0;

    std::vector<std::pair<double, double>> mapPoints;
    size_t num_mapPoints = std::min(msg->mapVertex_x.size(), msg->mapVertex_y.size());
    for (size_t i = 0; i < num_mapPoints; ++i)
    {
        double x = msg->mapVertex_x[i];
        double y = msg->mapVertex_y[i];
        mapPoints.emplace_back(x, y);
        // 创建一个Point实例并设置其x, y, z值
        if (i < num_mapPoints - 1)
        {
            geometry_msgs::Point first_point;
            geometry_msgs::Point second_point;
            first_point.x = msg->mapVertex_x[i];
            first_point.y = msg->mapVertex_y[i];
            first_point.z = 0; // 设z值为0

            second_point.x = msg->mapVertex_x[i + 1];
            second_point.y = msg->mapVertex_y[i + 1];
            second_point.z = 0;
            // 现在将这个点添加到mapEdge.points中
            mapEdge.points.push_back(first_point);
            mapEdge.points.push_back(second_point);
        }
    }
    // 确保图形是封闭的，将第一个点再次添加到末尾
    geometry_msgs::Point last_point;
    geometry_msgs::Point first_point;
    last_point.x = msg->mapVertex_x[num_mapPoints - 1];
    last_point.y = msg->mapVertex_y[num_mapPoints - 1];
    last_point.z = 0; 

    first_point.x = msg->mapVertex_x[0];
    first_point.y = msg->mapVertex_y[0];
    first_point.z = 0; 

    mapEdge.points.push_back(last_point);
    mapEdge.points.push_back(first_point);

    mapEdge_pub.publish(mapEdge);
    mapEdge.points.clear();

    double mapArea = calculatePolygonArea(mapPoints);
    patrolArea = mapArea;
    std::cout << "Area of the mapArea: " << mapArea << std::endl;

    // 将短点依次连成线段
    for (size_t i = 0; i < mapPoints.size() - 1; ++i)
    {
        mapEdges.emplace_back(mapPoints[i], mapPoints[i + 1]);
    }
    // 首尾相连
    mapEdges.emplace_back(mapPoints.back(), mapPoints.front());

    //生成种子点然后进行voronoi
    seedPoints = generatePoint(patrolArea, mapEdges, obstacleGroup);
    int i = 0;
    while (1)
    {
        generate_voronoi(seedPoints, mapEdges, obstacleEdges, obstacleGroup);

        // 初始化 seedPoint_pub
        visualization_msgs::Marker seedPoint;
        seedPoint.header.frame_id = "world";
        seedPoint.ns = "goals";
        seedPoint.id = 0;
        seedPoint.type = visualization_msgs::Marker::POINTS;
        seedPoint.action = visualization_msgs::Marker::ADD;
        seedPoint.pose.orientation.w = 1.0;
        seedPoint.scale.x = 0.1; // 在 x 方向上的大小
        seedPoint.scale.y = 0.1; // 在 y 方向上的大小
        seedPoint.scale.z = 0.1; // 在 z 方向上的大小
        seedPoint.color.a = 1.0; // 透明度
        seedPoint.color.r = 1.0; // 红色
        seedPoint.color.g = 0.0; // 绿色
        seedPoint.color.b = 0.0; // 蓝色

        for (const auto &point : seedPoints)
        {
            geometry_msgs::Point seedpoint;
            seedpoint.x = point.first;
            seedpoint.y = point.second;
            seedpoint.z = 0.0;
            seedPoint.points.push_back(seedpoint);
        }
        seedPoint_pub.publish(seedPoint);
        seedPoint.points.clear();

        // 初始化 voronoiEdge
        visualization_msgs::Marker voronoiEdge;
        voronoiEdge.header.frame_id = "world";
        voronoiEdge.ns = "edges";
        voronoiEdge.id = 0;
        voronoiEdge.type = visualization_msgs::Marker::LINE_LIST;
        voronoiEdge.action = visualization_msgs::Marker::ADD;
        voronoiEdge.pose.orientation.w = 1.0;
        voronoiEdge.scale.x = 0.05;
        voronoiEdge.color.a = 1.0;
        voronoiEdge.color.r = 1.0;
        voronoiEdge.color.g = 0.0;
        voronoiEdge.color.b = 1.0;

        for (const auto &edge : voronoiEdges)
        {
            geometry_msgs::Point first_point;
            geometry_msgs::Point second_point;
            first_point.x = edge.first.first;
            first_point.y = edge.first.second;
            first_point.z = 0;
            second_point.x = edge.second.first;
            second_point.y = edge.second.second;
            second_point.z = 0;
            voronoiEdge.points.push_back(first_point);
            voronoiEdge.points.push_back(second_point);
        }

        voronoiEdge_pub.publish(voronoiEdge);

        // lloyd
        i++;
        ROS_WARN("___________________________________________LLOYD %d", i);
        // for (const auto &point : seedPoints)
        // {
        //     ROS_WARN("%f,%f", point.first, point.second);
        // }

        seedPoints = lloyd(seedPoints, voronoiEdges);
        if(i>numOfIteration){break;}
        voronoiEdges.clear();
        
    }
    generate_adjacencyMatrix(seedPoints,obstacleEdges,voronoiEdges);
}


// 初始化 obstacleEdge 为全局
visualization_msgs::Marker obstacleEdge;
void obstacleAreaCallback(const patrol::Obstacle_area::ConstPtr &msg)
{

    obstacleEdge.header.frame_id = "world";
    obstacleEdge.ns = "edges";
    obstacleEdge.id = 0;
    obstacleEdge.type = visualization_msgs::Marker::LINE_LIST;
    obstacleEdge.action = visualization_msgs::Marker::ADD;
    obstacleEdge.pose.orientation.w = 1.0;
    obstacleEdge.scale.x = 0.05;
    obstacleEdge.color.a = 1.0;
    obstacleEdge.color.r = 0.0;
    obstacleEdge.color.g = 0.0;
    obstacleEdge.color.b = 1.0;

    std::vector<std::pair<double, double>> obstaclePoints;
    size_t num_obstaclePoints = std::min(msg->obstacleVertex_x.size(), msg->obstacleVertex_y.size());
    for (size_t i = 0; i < num_obstaclePoints; ++i)
    {
        double x = msg->obstacleVertex_x[i];
        double y = msg->obstacleVertex_y[i];
        obstaclePoints.emplace_back(x, y);

        // 创建一个Point实例并设置其x, y, z值
        if (i < num_obstaclePoints - 1)
        {
            geometry_msgs::Point first_point;
            geometry_msgs::Point second_point;
            first_point.x = msg->obstacleVertex_x[i];
            first_point.y = msg->obstacleVertex_y[i];
            first_point.z = 0; // 设z值为0

            second_point.x = msg->obstacleVertex_x[i + 1];
            second_point.y = msg->obstacleVertex_y[i + 1];
            second_point.z = 0;
            // 现在将这个点添加到mapEdge.points中
            obstacleEdge.points.push_back(first_point);
            obstacleEdge.points.push_back(second_point);
        }
    }

    // 确保图形是封闭的，将第一个点再次添加到末尾
    geometry_msgs::Point last_point;
    geometry_msgs::Point first_point;
    last_point.x = msg->obstacleVertex_x[num_obstaclePoints - 1];
    last_point.y = msg->obstacleVertex_y[num_obstaclePoints - 1];
    last_point.z = 0;

    first_point.x = msg->obstacleVertex_x[0];
    first_point.y = msg->obstacleVertex_y[0];
    first_point.z = 0;

    obstacleEdge.points.push_back(last_point);
    obstacleEdge.points.push_back(first_point);

    obstacleEdge_pub.publish(obstacleEdge);

    double obstacleArea = calculatePolygonArea(obstaclePoints);
    patrolArea = patrolArea - obstacleArea;
    std::cout << "Area of the obstacleArea: " << obstacleArea << std::endl;
    std::cout << "Area of the patrolArea: " << patrolArea << std::endl;

    std::vector<std::pair<std::pair<double, double>, std::pair<double, double>>> obstacle;
    // 将短点依次连成线段
    for (size_t i = 0; i < obstaclePoints.size() - 1; ++i)
    {
        obstacleEdges.emplace_back(obstaclePoints[i], obstaclePoints[i + 1]);
        obstacle.emplace_back(obstaclePoints[i], obstaclePoints[i + 1]);
    }
    // 首尾相连
    obstacleEdges.emplace_back(obstaclePoints.back(), obstaclePoints.front());
    obstacle.emplace_back(obstaclePoints.back(), obstaclePoints.front());
    
    //将障碍物边分组装进vector
    obstacleGroup.push_back(obstacle);
    obstacle.clear();

    // 生成种子点然后进行voronoi
    seedPoints = generatePoint(patrolArea, mapEdges, obstacleGroup);

    int i = 0;
    while (1)
    {
        generate_voronoi(seedPoints, mapEdges, obstacleEdges, obstacleGroup);

        // 初始化 seedPoint_pub
        visualization_msgs::Marker seedPoint;
        seedPoint.header.frame_id = "world";
        seedPoint.ns = "goals";
        seedPoint.id = 0;
        seedPoint.type = visualization_msgs::Marker::POINTS;
        seedPoint.action = visualization_msgs::Marker::ADD;
        seedPoint.pose.orientation.w = 1.0;
        seedPoint.scale.x = 0.1; // 在 x 方向上的大小
        seedPoint.scale.y = 0.1; // 在 y 方向上的大小
        seedPoint.scale.z = 0.1; // 在 z 方向上的大小
        seedPoint.color.a = 1.0; // 透明度
        seedPoint.color.r = 1.0; // 红色
        seedPoint.color.g = 0.0; // 绿色
        seedPoint.color.b = 0.0; // 蓝色

        for (const auto &point : seedPoints)
        {
            geometry_msgs::Point seedpoint;
            seedpoint.x = point.first;
            seedpoint.y = point.second;
            seedpoint.z = 0.0;
            seedPoint.points.push_back(seedpoint);
        }
        seedPoint_pub.publish(seedPoint);
        seedPoint.points.clear();

        // 初始化 voronoiEdge
        visualization_msgs::Marker voronoiEdge;
        voronoiEdge.header.frame_id = "world";
        voronoiEdge.ns = "edges";
        voronoiEdge.id = 0;
        voronoiEdge.type = visualization_msgs::Marker::LINE_LIST;
        voronoiEdge.action = visualization_msgs::Marker::ADD;
        voronoiEdge.pose.orientation.w = 1.0;
        voronoiEdge.scale.x = 0.05;
        voronoiEdge.color.a = 1.0;
        voronoiEdge.color.r = 1.0;
        voronoiEdge.color.g = 0.0;
        voronoiEdge.color.b = 1.0;

        for (const auto &edge : voronoiEdges)
        {
            geometry_msgs::Point first_point;
            geometry_msgs::Point second_point;
            first_point.x = edge.first.first;
            first_point.y = edge.first.second;
            first_point.z = 0;
            second_point.x = edge.second.first;
            second_point.y = edge.second.second;
            second_point.z = 0;
            voronoiEdge.points.push_back(first_point);
            voronoiEdge.points.push_back(second_point);
        }

        voronoiEdge_pub.publish(voronoiEdge);
        

        //lloyd
        i++;
        ROS_WARN("___________________________________________LLOYD %d",i);
        // for (const auto &point : seedPoints)
        // {
        //     ROS_WARN("%f,%f", point.first, point.second);
        // }

        seedPoints = lloyd(seedPoints,voronoiEdges);

        if (i>numOfIteration){break;}
        voronoiEdges.clear();
        
    }
    generate_adjacencyMatrix(seedPoints,obstacleEdges,voronoiEdges);
}
