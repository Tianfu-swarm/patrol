cmake_minimum_required(VERSION 2.8.3)
project(patrol)  # 替换为你的包名

## 找到catkin宏，这个宏用来编译catkin软件包
## 和在ROS中使用catkin软件包
find_package(catkin  REQUIRED COMPONENTS
  roscpp
  std_msgs
  message_generation  # 添加message_generation
  
)
find_package(CGAL)

## 添加消息文件
add_message_files(
  FILES
  Adjacency_matrix.msg
  Patrol_area.msg
  Obstacle_area.msg
)

## 在运行编译器时为C++代码生成消息头文件
generate_messages(
  DEPENDENCIES
  std_msgs
)

## 声明一个catkin包
## 填写在构建过程中需要的其它CMake包
catkin_package(
  CATKIN_DEPENDS roscpp std_msgs message_runtime  # 添加message_runtime
)

## 指定构建目录
include_directories(
  ${catkin_INCLUDE_DIRS}
)



## 添加可执行文件
add_executable(Markov_transfer_matrix src/Markov_transfer_matrix.cpp)


## 链接catkin库
target_link_libraries(Markov_transfer_matrix
  ${catkin_LIBRARIES}
  ${CGAL_LIBRARIES}
)

## 声明消息依赖
add_dependencies(Markov_transfer_matrix ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})

