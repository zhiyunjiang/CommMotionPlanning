function plotComAwareRRT(pppi,q1, q2)
    hold on
    pppi.plotProb(0)
    plotBS(q1);
    
    if nargin == 3
        plotBS(q2);
    end
end

function plotBS(q)
    marker_edge_color = 'k';
    hold on
    scatter(q(1), q(2), 90, '^', 'filled', 'MarkerEdgeColor', marker_edge_color);
end

