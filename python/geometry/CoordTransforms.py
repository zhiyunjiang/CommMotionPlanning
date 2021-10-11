import numpy as np

import warnings
import math

def toGridFromRaw(region, res, raw_pt):
    """
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % toGridFromRaw
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Translates arbitrary point to a point on discretized version of the
    % space. The discretized space is 0-indexed.
    %
    % Inputs:
    % region - [x_max, x_min, y_max, y_min]
    % res - resolution of the discretized region, i.e. 1/res = step size
    % raw_pt - arbitrary point x,y
    %
    % Output:
    % grid_pt - x,y in grid (index) coordiantes
    """
    raw_pt = np.array(raw_pt)
    if raw_pt.ndim == 1:
        raw_pt = np.reshape(raw_pt, (1,2))
    _check_region(region, raw_pt)

    if math.isinf(res): # just offset so that min is 0, 0
        grid_pt = raw_pt - np.array([region[1], region[3]])
    else:
        xGMax = np.ceil( (region[0] - region[1])*res)-2
        yGMax = np.ceil( (region[2] - region[3])*res)-2
        xG = np.minimum(np.rint((raw_pt[:,0] - region[1])*res), xGMax)
        yG = np.minimum(np.rint((raw_pt[:,1] - region[3])*res), yGMax)
        grid_pt = np.array([xG,yG], dtype=int).T
    return grid_pt


def toRawFromGrid(region, res, grid_pt):
    """
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % toRawFromGrid
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Translates point in discretized, normalized space to a point in
    % continuous space
    %
    % Inputs:
    % region - [x_max, x_min, y_max, y_min]
    % res - resolution of the discretized region, i.e. 1/res = step size
    % grid_pt - point in discretized, normalized space
    %
    % Output:
    % raw_pt - x,y in continuous, non-shifted space
    """
    if math.isinf(res):
       res = 1 #we just need to shift 

    grid_pt = np.array(grid_pt)
    if grid_pt.ndim == 1:#it's a single point
        grid_pt = np.reshape(grid_pt, (1,2))
    raw_pt = (grid_pt)/res + [region[1], region[3]]
    
    _check_region(region, raw_pt)    

    return raw_pt


def changes_res(pts, region, res):
    grid = toGridFromRaw(region, res, pts)
    new_raw = toRawFromGrid(region, res, grid)
    return np.unique(new_raw, axis = 0)


def _check_region(region, raw_pt):
    if any( raw_pt[:,0] > region[0] ) or any( raw_pt[:,0] < region[1] )\
    or any( raw_pt[:,1] > region[2] ) or any( raw_pt[:,1] < region[3] ):
        warnings.warn('Some points are outside of the region')


