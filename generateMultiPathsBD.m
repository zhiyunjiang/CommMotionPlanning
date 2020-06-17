function [s_roots, d_roots, s_leaves, d_leaves] = generateMultiPathsBD(RDT, source, dest, star, dest_freq, n_paths)
    s_roots(1:n_paths) = RDTNode;
    d_roots(1:n_paths) = RDTNode;
    s_leaves(1:n_paths) = RDTNode;
    d_leaves(1:n_paths) = RDTNode;

    for i=1:n_paths
        [Vs, Vd, ~] = RDT.FindPathBD(source, dest, star, [], dest_freq);
        s_roots(i) = Vs(1);
        d_roots(i) = Vd(1);
        s_leaves(i) = Vs(end);
        d_leaves(i) = Vd(end);
    end
end

