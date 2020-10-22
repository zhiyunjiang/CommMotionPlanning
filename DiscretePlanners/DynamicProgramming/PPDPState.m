classdef PPDPState < handle
    
    properties (Access = public)
        loc;
        visited;
    end
    
    methods
        function obj = PPDPState(loc, visited)
            obj.loc = loc;
            obj.visited = visited;
        end
        
        function next_visited = calcNextVisited(o)
           %create a deep copy of the map
           xkeys = cell2mat(o.visited.keys());
           
           next_visited = containers.Map('KeyType','double','ValueType','any');
           for i = 1:length(xkeys)
               key = xkeys(i);
               submap = o.visited(key);
               next_visited(key) = containers.Map(submap.keys, submap.values);
           end
           x = o.loc(1);
           y = o.loc(2);
           if next_visited.isKey(x)
                submap = next_visited(x);
                submap(y) = 1;
           else
                submap = containers.Map(y,1);
                next_visited([x]) = submap;
           end
        end
        
        function tf = locIsVisited(o, loc)
            tf = 0;
             x = loc(1);
             if o.visited.isKey(x)
                y = loc(2);
                submap = o.visited(x);
                if submap.isKey(y)
                   tf = 1; 
                end
             end
        end
        
    end
end

