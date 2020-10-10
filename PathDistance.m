function dist = PathDistance(path)
    [path_length, ~] = size(path);
    dist = 0;
    for i=2:path_length
        dist = dist + norm(path(i-1,:) - path(i,:)); 
    end
end

