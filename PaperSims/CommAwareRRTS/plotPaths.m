%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotPaths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot paths over the workspace. Assumes the Communication Power field and 
% problem instance + remote stations have already been plotted. Removes any
% title
% 
% Inputs:
% path1 - nx2 array of waypoints for path 1
% label1 - legend label for path 1
% path2 - nx2 array of waypoints for path 2 (optional)
% label2 - legend label for path 2 (optional)
% 
% Output:
% path_handles - either the handle to plotted path1 if only 1 path is given
%                   or an array of handle objects to the paths if two paths
%                   given

function path_handles = plotPaths( path1, label1,  path2, label2)
    line_width = 2.5;
    
    hold on
    if nargin < 3
        path_handles = plot(path1(:,1), path1(:,2), 'LineWidth', line_width, ...
                            'DisplayName', label1);
    else
        plot_handle1 = plot(path1(:,1), path1(:,2), 'LineWidth', line_width, ...
                            'DisplayName', label1);

        plot_handle2 = plot(path2(:,1), path2(:,2), 'LineWidth', line_width, ...
                            'DisplayName', label2);

        path_handles = [plot_handle1, plot_handle2];
    end
    
    title('')
    
end

