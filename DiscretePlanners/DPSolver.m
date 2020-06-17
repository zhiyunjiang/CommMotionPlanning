classdef DPSolver < handle
    
    properties
        MUstar;
        JStar;
        BSF;
        pppi;
    end
    
    methods (Access = public)
        
        function solve(o, pppi)
            o.pppi = pppi;
            loc_0 = pppi.getSourceGrid();
            visited_0 = containers.Map('KeyType','double','ValueType','any');
            o.BSF = Inf;
            
            x_0 = PPDPState(loc_0, visited_0);
            
            o.J(x_0,0);
        end
        
        function setBSF(o, val)
           o.BSF =  val;
        end
    end
    
    methods (Access = private)
        function min =  J(o, state, cost2here)
            if o.pppi.pointInDestGrid(state.loc)
                min = 0;
                if cost2here<o.BSF
                   o.setBSF(cost2here); 
                end
            else
                min = Inf;
                controls = o.getAdmissibleControls(state);
                for i = 1:length(controls)
                    stage_cost = 1;
                    next_loc = state.loc + controls(i);
                    optimistic_c2g = sum(abs(next_loc - o.pppi.getSourceGrid()));

                    if optimistic_c2g + cost2here + stage_cost < o.BSF
                       %only bother going down this branch if we have a shot 
                       next_visited = state.calcNextVisited();
                       next_state = PPDPState(next_loc, next_visited);
                       
                       minc2g = o.J(next_state, cost2here + 1);
                       total = stage_cost + minc2g;
                       if min < total
                            min = total;
                       end
                       
                    end

                end
            end
            
        end
        
        function controls = getAdmissibleControls(o, state)
            %can go up, down, left or, right, now see which are disasslowed
            init = state.loc;
            controls = [];
            for i = 1:4
               xp = sign(i-3)*mod(i+1,2);
               yp = sign(i-2)*mod(i,2);
               control = [xp, yp];
               next_loc = init+control; 
               
               if o.pppi.gridPtInRegion(next_loc) && (~state.locIsVisited(next_loc))
                   free_sub = o.pppi.collisionFree([init, next_loc]);
                   if sum(size(free_sub)) == 4 %i.e. it's the same full path
                        controls = [control, controls]; 
                   end
               end
            end
        end
    end
end

