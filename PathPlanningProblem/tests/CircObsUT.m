%%
%%%%%%%%%%%%%%%%%%
% CircObs
%%%%%%%%%%%%%%%%%%
test_num = 0;
failures = 0;

fprintf('Testing CircObs ''obstructs'' method...\n')
testc1 = CircObs(1,[0,0]);

%test interior and interior (always obstructs)
test_num = test_num + 1;
if ~test_obstacle_obstructs(testc1, [0.75,0.25], [0.25, -0.25], 1)
    failures = failures + 1;
    fprintf('Failed Test %d\n', test_num);
end

%Test interior and exterior (always obstructs)
test_num = test_num + 1;
if ~test_obstacle_obstructs(testc1, [2,2], [0.5, 0.5], 1)
    failures = failures + 1;
    fprintf('Failed Test %d\n', test_num);
end

test_num = test_num + 1;
if ~test_obstacle_obstructs(testc1, [-2,2], [-0.5, 0.5], 1)
    failures = failures + 1;
    fprintf('Failed Test %d\n', test_num);
end

test_num = test_num + 1;
if ~test_obstacle_obstructs(testc1, [-2,-2], [-0.5, -0.5], 1)
    failures = failures + 1;
    fprintf('Failed Test %d\n', test_num);
end

test_num = test_num + 1;
if ~test_obstacle_obstructs(testc1, [2,-2], [0.5, -0.5], 1)
    failures = failures + 1;
    fprintf('Failed Test %d\n', test_num);
end

%Test exterior, does obstrcut
test_num = test_num + 1;
if ~test_obstacle_obstructs(testc1, [-1.1,0], [0, 1.1], 1)
    failures = failures + 1;
    fprintf('Failed Test %d\n', test_num);
end

test_num = test_num + 1;
if ~test_obstacle_obstructs(testc1, [-1.1,0], [0, -1.1], 1)
    failures = failures + 1;
    fprintf('Failed Test %d\n', test_num);
end

test_num = test_num + 1;
if ~test_obstacle_obstructs(testc1, [1.1,0], [0, -1.1], 1)
    failures = failures + 1;
    fprintf('Failed Test %d\n', test_num);
end

test_num = test_num + 1;
if ~test_obstacle_obstructs(testc1, [1.1,0], [0, 1.1], 1)
    failures = failures + 1;
    fprintf('Failed Test %d\n', test_num);
end

%Test exterior, does not obstruct
test_num = test_num + 1;
if ~test_obstacle_obstructs(testc1, [2, 0], [0, 2], 0)
    failures = failures + 1;
    fprintf('Failed Test %d\n', test_num);
end

test_num = test_num + 1;
if ~test_obstacle_obstructs(testc1, [-2, 0], [0, 2], 0)
    failures = failures + 1;
    fprintf('Failed Test %d\n', test_num);
end

test_num = test_num + 1;
if ~test_obstacle_obstructs(testc1, [-2, 0], [0, -2], 0)
    failures = failures + 1;
    fprintf('Failed Test %d\n', test_num);
end

test_num = test_num + 1;
if ~test_obstacle_obstructs(testc1, [2, 0], [0, -2], 0)
    failures = failures + 1;
    fprintf('Failed Test %d\n', test_num);
end


fprintf('Ran %d tests.\n',test_num);
if (test_num - failures == test_num)
    fprintf('All tests passed!\n');
else
    fprintf('Failed %d of %d tests\n',failures, test_num);
end
