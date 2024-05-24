
example 1:
手动发布拓扑图的紧邻矩阵
note：在发布话题时的终端需要source
3--------4
|\      /|
|  \  /  |
|   / \  |
|  /   \ |
|/      \|
1        2
rostopic pub /Adjacency_matrix patrol/Adjacency_matrix '{x: [0.1, 1.0, 0.1, 1.0], y: [0.1, 0.1, 1.0, 1.0], matrix: [0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0]}'
                                                       {x: [0.1, 1.0, 0.1, 1.0], 
                                                        y: [0.1, 0.1, 1.0, 1.0], matrix: [0, 0, 1, 1, 
                                                                                          0, 0, 1, 1, 
                                                                                          1, 1, 0, 1, 
                                                                                          1, 1, 1, 0]}


发布map 与 obstacle
rostopic pub /patrol_area patrol/Patrol_area "{mapVertex_x: [0.0, 0.0, 20.0, 20.0], mapVertex_y: [0.0, 20.0, 20.0, 0.0]}"

rostopic pub /obstacle_area patrol/Obstacle_area "{ obstacleVertex_x: [5, 5, 10.0, 10.0], obstacleVertex_y:[5, 10.0, 10.0, 5]}"

rostopic pub /obstacle_area patrol/Obstacle_area "{ obstacleVertex_x: [11.0, 11.0, 15.0, 15.0], obstacleVertex_y:[11.0, 15.0, 15.0, 11.0]}"

rostopic pub /patrol_area patrol/Patrol_area "{mapVertex_x: [6.0, 0.0, 0.0, 3.0, 6.0], 
                                               mapVertex_y: [4.0, 4.0, 1.0, -2.0, -2.0]}"

##ROBOT1
rostopic pub /patrol_area patrol/Patrol_area "{mapVertex_x: [6.0, 5.0, 5.0, 0.0, 0.0, 3.0, 6.0], 
                                               mapVertex_y: [3.0, 3.0, 4.0, 4.0, 1.0, -2.0, -2.0]}"
##ROBOT2
rostopic pub /patrol_area patrol/Patrol_area "{mapVertex_x: [6.0, 5.0, 5.0, 0.0, 0.0, 3.0, 6.0], 
                                               mapVertex_y: [3.0, 3.0, 4.0, 4.0, 1.0, -2.0, -2.0]}"