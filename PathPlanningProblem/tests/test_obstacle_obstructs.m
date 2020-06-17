function success = test_obstacle_obstructs(obstacle, p1, p2, expected)
    success = (obstacle.obstructs(p1, p2) == expected) && (obstacle.obstructs(p2, p1)==expected);
end

