%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getRRTSolver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots the problem instance as well as the remote station markers. Assumes
% the Communication Power field has already been plotted
%
% Inputs:
% pppi - the path planning problem instance to be solved
% q1 - [x,y] of the first remote station
% q2 - [x,y] of the second remote station
%
% Output:
% rshandle - handle to the first remote station series. May be useful when
%               creating legends


function rshandle = plotComAwareRRT(pppi,q1, q2)
    hold on
    pppi.plotProb()
    rshandle = plotRS(q1);
    
    if nargin == 3
        plotRS(q2);
    end
end

function scatter_handle = plotRS(q)
    marker_edge_color = 'k';
    marker_face_color = [1,0,1];
    marker_size = 90;
    hold on
    scatter_handle = scatter(q(1), q(2), marker_size, '^', 'filled', ...
                            'MarkerEdgeColor', marker_edge_color,...
                            'MarkerFaceColor', marker_face_color ,...
                            'DisplayName', 'Remote Station');
end

