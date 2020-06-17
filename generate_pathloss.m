function gamma_PL_dB = generate_pathloss(region, q_b, res, K_PL, n_PL)
%generate_pathloss - Modulated from ChannelSim/channel_simulator

    % the rectangular region specifying the environment
    x_max = region(1);
    x_min = region(2);
    y_max = region(3);
    y_min = region(4);

    % the position of the base station
    x_b = q_b(1);
    y_b = q_b(2);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % generating the grid
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [g_x, g_y] = meshgrid(x_min:1/res:x_max,y_min:1/res:y_max);
    [M, N] = size(g_x);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % path loss
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\ngrid size = %d pixels\n', M*N*res)
    fprintf('generating path loss...\n')

    g_d = sqrt((g_x - x_b).^2 + (g_y - y_b).^2);

    % prevent the path loss to become very large if the samples are very close
    % to the base station
    g_d(g_d < 1/res) = 1/res;

    % generating the path loss
    gamma_PL_dB = K_PL - 10*n_PL*log10(g_d);
end

