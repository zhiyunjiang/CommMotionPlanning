function dist = ManhattanDistance(n1, n2, path)
    if nargin == 3
        [dist, ~] = size(path);
        dist = dist - 1;
    else
        p1 = n1.getPos();
        p2 = n2.getPos();
        dist = norm(p1-p2, 1);
    end
end

