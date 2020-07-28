function result = IntegrateWSimpsonsIrregular(f, h)
    %number of intervals
    n = length(h);

    w = SimpsonWeights(n, h);
    
    %now just sum those bad bois
    result = w*f;
end

