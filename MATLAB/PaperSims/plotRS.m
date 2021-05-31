function scatter_handle = plotRS(q, color)
    if nargin == 2
        marker_face_color = color;
    else
        marker_face_color = [1,0,1];
    end
    marker_edge_color = 'k';
    marker_size = 90;
    hold on
    scatter_handle = scatter(q(1), q(2), marker_size, '^', 'filled', ...
                            'MarkerEdgeColor', marker_edge_color,...
                            'MarkerFaceColor', marker_face_color ,...
                            'DisplayName', 'Remote Station');
end