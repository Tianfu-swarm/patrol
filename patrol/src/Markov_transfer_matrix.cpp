#include <cmath> 
#include <vector>
#include <random>
#include "ros/ros.h"
#include "patrol.hpp"
#include <nav_msgs/Odometry.h>
#include <geometry_msgs/Point.h>
#include <geometry_msgs/PoseStamped.h>
#include "patrol/Adjacency_matrix.h"

// 定义全局变量 markov_matrix
std::vector<std::vector<double>> markov_matrix;
void generateMarkovTransitionMatrix(const std::vector<std::vector<int>> &input_matrix)
{
    // 生成马尔可夫转移矩阵
    // 获取输入矩阵的维度
    int n = input_matrix.size();

    // 初始化马尔可夫转移矩阵为和输入矩阵相同大小，并填充为0
    markov_matrix = std::vector<std::vector<double>>(n, std::vector<double>(n, 0.0));

    // 使用随机数生成器
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    // 遍历输入矩阵
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            // 如果输入矩阵中的元素为1，生成一个随机小数
            if (input_matrix[i][j] == 1)
            {
                markov_matrix[i][j] = dis(gen);
            }
            // 如果输入矩阵中的元素为0，保持为0
            else
            {
                markov_matrix[i][j] = 0.0;
            }
        }
    }

    // 打印邻接矩阵
    ROS_INFO("Received markov matrix:");
    for (int i = 0; i < n; ++i)
    {
        std::string row_string;
        for (int j = 0; j < n; ++j)
        {
            row_string += std::to_string(markov_matrix[i][j]) + " ";
        }
        ROS_INFO("%s", row_string.c_str());
    }

}

// 定义全局变量 patrol_point
std::vector<std::pair<double, double>> patrol_points;
void matrixCallback(const patrol::Adjacency_matrix::ConstPtr &msg)
{
    // 获取 x 和 y 坐标
    patrol_points.clear(); // 清空之前的数据
    for (size_t i = 0; i < msg->x.size(); ++i)
    {
        patrol_points.push_back(std::make_pair(msg->x[i], msg->y[i]));
    }

    // 获取 adjacency_matrix 数据并转换成二维数组
    std::vector<int> adjacency_matrix_data = msg->matrix;
    int n = std::sqrt(adjacency_matrix_data.size()); // 计算矩阵维度
    std::vector<std::vector<int>> adjacency_matrix(n, std::vector<int>(n));

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            int index = i * n + j;
            adjacency_matrix[i][j] = adjacency_matrix_data[index];
        }
    }

    // 打印坐标
    ROS_INFO("Received patrol points:");
    for (const auto &point : patrol_points)
    {
        ROS_INFO("(%f, %f)", point.first, point.second);
    }

    // 打印矩阵
    ROS_INFO("Received Adjacency Matrix:");
    for (int i = 0; i < n; ++i)
    {
        std::string row_string;
        for (int j = 0; j < n; ++j)
        {
            row_string += std::to_string(adjacency_matrix[i][j]) + " ";
        }
        ROS_INFO("%s", row_string.c_str());
    }

    generateMarkovTransitionMatrix(adjacency_matrix);
}

//计算两个点之间的距离
double distance(const std::pair<double, double> &p1, const std::pair<double, double> &p2)
{
    return std::sqrt(std::pow(p1.first - p2.first, 2) + std::pow(p1.second - p2.second, 2));
}

ros::Publisher goal_pub;
geometry_msgs::PoseStamped goal_point;

int current_goal_index = -1; // 初始值为 -1，表示尚未选择目标点
void pubPatrolGoalPoint(const std::pair<double, double> &current_position,
                        const std::vector<std::pair<double, double>> &patrol_point,
                        const std::vector<std::vector<double>> &transfor_matrix)
{
    
    // 如果还没有发布目标点
    if (goal_point.pose.position.x == 0 && goal_point.pose.position.y == 0)
    {
        double min_distance = std::numeric_limits<double>::max();
        std::pair<double, double> closest_point;
        for (int i = 0; i < patrol_point.size(); ++i)
        {
            const auto &point = patrol_point[i];

            double dist = distance(current_position, point);
            if (dist < min_distance)
            {
                min_distance = dist;
                closest_point = point;
                current_goal_index = i; // 更新当前目标点的索引
            }
        }
	ROS_INFO("Pub the first goal");
        goal_point.pose.position.x = closest_point.first;
        goal_point.pose.position.y = closest_point.second;
        goal_pub.publish(goal_point);
    }

    std::pair<double, double> goal;
    goal.first = goal_point.pose.position.x;
    goal.second = goal_point.pose.position.y;
    double distanceToGoal = distance(current_position, goal);

    double distanceThreshold;
    ros::param::get("/parameters/distance_threshold", distanceThreshold);
    double Threshold = distanceThreshold;
    if (distanceToGoal < Threshold)
    {
        // 到达目标点时的操作
        // 检查转换矩阵的第i行
        const std::vector<double> &probability_row = transfor_matrix[current_goal_index];

        // 寻找所有概率不为零的元素，即非零概率对应的巡逻点索引
        std::vector<int> nonzero_indices;
        for (int j = 0; j < probability_row.size(); ++j)
        {
            if (probability_row[j] != 0)
            {
                nonzero_indices.push_back(j);
            }
        }

        std::random_device rd;
        std::mt19937 gen(rd());
        // 根据概率选择下一个巡逻点
        if (!nonzero_indices.empty())
        {
            // 创建一个离散概率分布，概率由转换矩阵的非零元素提供
            std::vector<double> probabilities;
            for (int index : nonzero_indices)
            {
                probabilities.push_back(transfor_matrix[current_goal_index][index]);
            }
            std::discrete_distribution<int> distribution(probabilities.begin(), probabilities.end());

            // 使用离散分布进行加权随机选择
            int next_patrol_index = nonzero_indices[distribution(gen)];

            // 发布下一个巡逻点
            goal_point.pose.position.x = patrol_point[next_patrol_index].first;
            goal_point.pose.position.y = patrol_point[next_patrol_index].second;
            // 更新当前目标点的索引
            current_goal_index = next_patrol_index;
        }
	ROS_INFO("patrlling~");
        goal_pub.publish(goal_point);
    }
}

// 回调函数，处理接收到的odometry消息
void odometryCallback(const nav_msgs::Odometry::ConstPtr &msg)
{

    // 机器人的姿态
    std::pair<double, double> robot_pose;

    // 提取位置信息（x 和 y 坐标）
    robot_pose.first = msg->pose.pose.position.x;
    robot_pose.second = msg->pose.pose.position.y;

    if(!patrol_points.empty()){
    pubPatrolGoalPoint(robot_pose, patrol_points, markov_matrix);}
}


int main(int argc, char **argv)
{

    // 初始化ROS节点
    ros::init(argc, argv, "markov_transfer_matrix");

    // 创建节点句柄
    ros::NodeHandle nh;

    std::string patrol_goal_topic;
    std::string Adjacency_matrix_topic;
    std::string Odometry_topic;

    goal_point.header.frame_id = "world";
    goal_point.pose.position.x = 0;
    goal_point.pose.position.y = 0;

    nh.param<std::string>("common/patrol_goal_topic", patrol_goal_topic, "/patrol_goal");
    nh.param<std::string>("common/Adjacency_matrix_topic", Adjacency_matrix_topic, "/Adjacency_matrix");
    nh.param<std::string>("common/Odometry_topic", Odometry_topic, "/Odometry");
    nh.param("common/numOfIteration", numOfIteration);
    nh.param("common/patrolPointsDensity", patrolPointsDensity);

    goal_pub = nh.advertise<geometry_msgs::PoseStamped>(patrol_goal_topic, 10);
    ros::Subscriber Adjacency_matrix_sub = nh.subscribe(Adjacency_matrix_topic, 10, matrixCallback);
    ros::Subscriber Odometry_sub = nh.subscribe(Odometry_topic, 10, odometryCallback);
    
    ros::Subscriber patrol_sub = nh.subscribe("/patrol_area", 10, patrolAreaCallback);
    ros::Subscriber obstacle_sub = nh.subscribe("/obstacle_area", 10, obstacleAreaCallback);
    adjacencyMatrix_pub = nh.advertise<patrol::Adjacency_matrix>(Adjacency_matrix_topic, 10);

    mapEdge_pub = nh.advertise<visualization_msgs::Marker>("/mapEdge", 1);
    obstacleEdge_pub = nh.advertise<visualization_msgs::Marker>("/obstacleEdge", 1);
    voronoiEdge_pub = nh.advertise<visualization_msgs::Marker>("/voronoiEdge", 1);
    seedPoint_pub = nh.advertise<visualization_msgs::Marker>("/seedPoint", 1);

    // 循环等待回调函数
    ros::spin();

    return 0;
}