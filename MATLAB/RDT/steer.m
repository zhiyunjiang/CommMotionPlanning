function path = steer(steer_rad, dest, start, pppi)

    diff = dest - start;
    %don't extend beyond the newly sampled point
    if norm(diff) > steer_rad

        %find the direction
        phi = atan2(diff(2), diff(1));

        new_diff = [steer_rad*cos(phi), steer_rad*sin(phi)];
        new = start + new_diff;
        %snap to grid
        new = pppi.toRawCoordinate(pppi.toGridCoordinate(new));
    else
       new =  dest;
    end

   path = pppi.getSLPath(start, new);

end

