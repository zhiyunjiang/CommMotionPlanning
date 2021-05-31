%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% runOneSim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runs a single instance of RRT on a path planning problem instance
% 
% Inputs:
% rrt_solver - the solver to use to path planning problem
% pppi - the path planning problem instance to be solved
%
% Outputs:
% rrt_path - the solution (path) found by rrt_solver, represented by a
%               series of waypoints
function rrt_path = runOneSim(rrt_solver, pppi)
    rrt_solver.solve(pppi);
    bsf = rrt_solver.getBSF();
    rrt_path = bsf.pathToRoot(0);
end