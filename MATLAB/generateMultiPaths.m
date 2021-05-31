function leaves = generateMultiPaths(RDT, source, dest, star, dest_freq, n_paths)
    leaves(1:n_paths) = RDTNode;

    for i=1:n_paths
        [~, leaf, ~] = RDT.FindPath(source, dest, star, [], dest_freq);
        leaves(i) = leaf;
    end
end

