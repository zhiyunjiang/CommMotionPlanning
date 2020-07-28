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
        
        PL;
    end

    
    methods (Access = public)
        %TreeNode Constructor
        function this = TreeNode(pos)
            if nargin>0
                this.wrkspcpos = pos;
            else
                this.wrkspcpos = [0, 0];
            end
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
        
        function PL_vals = getPLValsToRoot(this)
            current = this;
            PL_vals = [];
            while ~current.isRoot
               %only take the path up to but not including the parent
               %parent is first, so don't add it (will be added next
               %iteration
               PL_vals = [ current.PL; PL_vals];
               current = current.parent;
            end
            %tack on the root
            PL_vals = [current.PL; PL_vals];
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
            if this.distToHere == Inf
                %we've just connected, we have no children!
               this.distToHere = parent_dist + dist_btwn; 
            else
                diff = this.distToHere - (parent_dist + dist_btwn);
                this.propogateDiff(diff);
            end
            
        end
        %if the distance to me gets updated, need to update distance to my
        %children
        function propogateDiff(this, diff)
            this.distToHere = this.distToHere - diff;
            for i = 1:length(this.children)
                child = this.children(i);
                child.propogateDiff(diff);
            end
        end
    end
end

