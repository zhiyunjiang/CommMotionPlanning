classdef VisitedTable
    %VISITEDTABLE - keeps tracks of points interested or
    
    properties
        map;
    end
    
    methods
        function this = VisitedTable()
            this.map = containers.Map('KeyType','double','ValueType','any');
        end
        
        function add(this, pos, do_add)
            if nargin == 2
                do_add = 0;
            end
            
            x = pos(1);
            y = pos(2);
            if this.map.isKey(x)
                submap = this.map(x);
                %submap(y) is a handle object - modification persists 
                %outside this function call
                if submap.isKey(y) && do_add
                    submap(y) = submap(y) + 1;
                else
                    submap(y) = 1;
                end
            else
                submap = containers.Map(y,1);
                this.map(x) = submap;
            end
        end
        
        function [tf, val] = pointInTable(this, pos)
            tf = 0; val = 0;
            x = pos(1);

            if this.map.isKey(x)
                y = pos(2);
                submap = this.map(x);
                if submap.isKey(y)
                   tf = 1;
                   val = submap(y);
                end
            end
        end
        
        function new_table = deepCopy(this)
            new_table = VisitedTable();
            new_map = new_table.map;
            keys = cell2mat(this.map.keys);
            for i=1:length(keys)
                submap = this.map(keys(i));
                new_submap = containers.Map(submap.keys, submap.values);
                new_map(keys(i)) = new_submap;
            end
        end
    end
end

