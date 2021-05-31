classdef TreeNode < handle
    %TreeNode - represents a node in the tree created by a number of related
    % Random Dense Tree algorithms. Also includes convenient printing
    % methods
    %
    %PROPERTIES
    % wrkspcpos - x and y coordinates of the node in the workspace
    % parent - parent node. If this node is the root, then parent is itself
    % children - array of child nodes
    % isRoot - Flag set high for root node, low otherwise
    % distToHere - the distance to this node through the tree
    
    properties
        wrkspcpos;
        parent;
        pathFromParent;
        children = [];
        isRoot = 0;
        distToHere = Inf;
        isPropagating = 0;
        
        problemData;
    end

    
    methods (Access = public)
        %TreeNode Constructor
        function this = TreeNode(pos)
            if nargin>0
                this.wrkspcpos = pos;
            else
                this.wrkspcpos = [0, 0];
            end
            
            this.problemData = containers.Map;
        end
        
        function pos = getPos(this)
           pos = this.wrkspcpos; 
        end
        
        function dist = distTo(this, other, norm_num)
           if nargin == 2
               %make the default norm the 2-norm (euclidean distance)
              norm_num = 2; 
           end
           dist = norm(this.wrkspcpos - other.wrkspcpos, norm_num); 
        end

        function setParent(this, parent, dist, path_from_parent)
            %need to update cost to other nodes down the line, too
            this.updateDist(parent.distToHere, dist);
            if ~isempty(this.parent)
                this.parent.removeChild(this);
            end
            
            this.parent = parent;
            parent.children = [parent.children, this];
            
            %if the path_to_parent is actually passed
            if nargin == 4
                this.pathFromParent = path_from_parent;
            else
                %otherwise assume it's just the line between them
                this.pathFromParent = [parent.getPos(); this.getPos()];
            end
        end
        
        function setDist(this, new_dist)
            
            if new_dist < 0
                error("Distance must always be non-negative");
            else
                this.distToHere = new_dist;
            end
        end
          
        function root = getRootNode(this)
            node = this;
            while ~node.isRoot
               node = node.parent; 
            end
            
            root = node;
        end
         
        function tf = eq(this, that)
           tf = (norm(this.wrkspcpos - that.wrkspcpos) == 0);
        end
        
        function removeChild(this, thatChild)
           for i=1:length(this.children)
                thisChild = this.children(i);
               if  thisChild == thatChild
                    this.children(i) = [];
                   break;
               end
           end
        end
        
        function path = pathToRoot(this, do_plot, scale)
            if nargin == 2
               scale = 1; 
            end
            current = this;
            path = [];
            while ~current.isRoot
               %only take the path up to but not including the parent
               %parent is first, so don't add it (will be added next
               %iteration
               path = [ current.pathFromParent(2:end,:); path];
               current = current.parent;
            end
            %tack on the root
            path = [current.getPos(); path];
            %now scale
            path = scale*path;
            if do_plot
                plot(path(:,1), path(:,2));
            end
        end
        
        function plotTree(obj, style)
           queue = [obj];
           while ~isempty(queue)
               %pop the first off the queue
               current = queue(1);
               queue(1)=[];
               current_children = current.children;
               if ~isempty(current_children)
                   for i=1:length(current_children)
                       child = current_children(i);
                       path_to_parent = child.pathFromParent;
                       plot(path_to_parent(:,1), path_to_parent(:,2), style);
                       hold on
                       queue = [queue, child];
                   end
               end
           end
           hold off
        end
        
    end
    
    methods (Access = private)
        
        function updateDist(this, parent_dist, dist_btwn)
            if isnan(parent_dist) || parent_dist == Inf || parent_dist < 0 || ...
               isnan(dist_btwn) || dist_btwn == Inf || dist_btwn < 0
                error("Distance must be non-negative, finite number.");
            end
            
            if this.distToHere == Inf
                %we've just connected, we have no children!
                this.setDist(parent_dist + dist_btwn);
            else
                diff = this.distToHere - (parent_dist + dist_btwn);
                this.propogateDiff(diff);
            end
            
        end
        %if the distance to me gets updated, need to update distance to my
        %children
        %TODO - update to take into account the problem instance's distance
        %metric. Should only need to be calculated for rewire
        function propogateDiff(this, diff)
            if this.isPropagating
                error('Cycle detected in tree.');
            end
            this.isPropagating = 1;
            %avoid some very smoll rounding issues
            this.setDist( max(0,this.distToHere - diff) );

            for i = 1:length(this.children)
                this.children(i).propogateDiff(diff);
            end
            this.isPropagating = 0;
        end
    end
end

