%%
%%%%%%%%%%%%%%%%%%
% RectObs
%%%%%%%%%%%%%%%%%%
test_num = 0;
failures = 0;

fprintf('Testing RectObs ''obstructs'' method...\n')
test_rect = RectObs([-2, 2, -3, 3]);

%both points instide
test_num = test_num + 1;
if ~test_obstacle_obstructs(test_rect, [1, 0], [0, -1], 1)
    failures = failures + 1;
    fprintf('Failed Test %d\n', test_num);
end

%test one point in , one out
test_num = test_num + 1;
if ~test_obstacle_obstructs(test_rect, [1, 0], [-4, 0], 1)
    failures = failures + 1;
    fprintf('Failed Test %d\n', test_num);
end

test_num = test_num + 1;
if ~test_obstacle_obstructs(test_rect, [1, 0], [1, 10], 1)
    failures = failures + 1;
    fprintf('Failed Test %d\n', test_num);
end

test_num = test_num + 1;
if ~test_obstacle_obstructs(test_rect, [1, 0], [4, 2], 1)
    failures = failures + 1;
    fprintf('Failed Test %d\n', test_num);
end

test_num = test_num + 1;
if ~test_obstacle_obstructs(test_rect, [1, 0], [1.5, -10], 1)
    failures = failures + 1;
    fprintf('Failed Test %d\n', test_num);
end

%Both out, rectangle obstructs

%vertical line
test_num = test_num + 1;
if ~test_obstacle_obstructs(test_rect, [1, -10], [1, 10], 1)
    failures = failures + 1;
    fprintf('Failed Test %d\n', test_num);
end

%horizontal line
test_num = test_num + 1;
if ~test_obstacle_obstructs(test_rect, [-10, 0], [10, 0], 1)
    failures = failures + 1;
    fprintf('Failed Test %d\n', test_num);
end

%diagonal
test_num = test_num + 1;
if ~test_obstacle_obstructs(test_rect, [-3, -1], [3.25, 1.75], 1)
    failures = failures + 1;
    fprintf('Failed Test %d\n', test_num);
end

%single point on the edge
test_num = test_num + 1;
if ~test_obstacle_obstructs(test_rect, [-2, 3], [-2, 3], 1)
    failures = failures + 1;
    fprintf('Failed Test %d\n', test_num);
end

% no obstruction
test_num = test_num + 1;
if ~test_obstacle_obstructs(test_rect, [-5, -10], [7, 100], 0)
    failures = failures + 1;
    fprintf('Failed Test %d\n', test_num);
end

%vertical line, no obstruction
test_num = test_num + 1;
if ~test_obstacle_obstructs(test_rect, [-5, -10], [-5, 20], 0)
    failures = failures + 1;
    fprintf('Failed Test %d\n', test_num);
end

%Horizontal line, no obstruction
test_num = test_num + 1;
if ~test_obstacle_obstructs(test_rect, [-1, 10], [1, 10], 0)
    failures = failures + 1;
    fprintf('Failed Test %d\n', test_num);
end


fprintf('Ran %d tests.\n',test_num);
if (test_num - failures == test_num)
    fprintf('All tests passed!\n');
else
    fprintf('Failed %d of %d tests\n',failures, test_num);
end

