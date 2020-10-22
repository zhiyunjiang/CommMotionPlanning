%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ExpectedDistance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the expected distance adding n2 to a tree path will add to the
% path. Not a valid RRT* cost function since it is not additive. The order
% paths are visited affects the total cost of hte path
%
% Inputs:
% TreeNode n1 - The node to be added to the graph
% TreeNode n2 - The node already in the graph. 
% path - points from n1 to n2, inclusive
% ChannelAnalyzer ca - channel analyzer with relevant channel data. Used to
% calculate probabilities.
% mode - If mode == 1, just compute the distance. If mode == 2, also update
% the 'J' and 'p_no_conn_to_here' attributes of n1
%
% OUTPUT:
% double exp_dist - expected distance added to the path 

function exp_dist = ExpectedDistance(n1, n2, path, ca, mode)


    if nargin == 3
        mode = 1;
    end

    no_conn_key = 'p_no_conn_to_here';
    j_key = 'J';
    n2_data = n2.problemData;
    if n2.isRoot && ~n2_data.isKey(no_conn_key)
        %root coming through for the first time, need to intialize J
        %functions
        n2_data(j_key) = ca.J0(n2.getPos());
        %n2_data is a handle object, so we can just modify here
        n2_data(no_conn_key) = ca.IntegrateJ(n2_data(j_key));
    end
    
    if n2_data.isKey(no_conn_key)
        p_no_conn_to_n2 = n2_data(no_conn_key);
    else
       error('p_no_conn_to_here has not been computed for this node'); 
    end
    
    nroot = n2.getRootNode();

    %find the conditional probability
    p_no_conn_cond = p_no_conn_to_n2/nroot.problemData(no_conn_key);
   
    %use the probability to scale the distance
    dist = n1.distTo(n2);
    exp_dist = p_no_conn_cond*dist;
    
    if mode == 2
        %also setup J, prob of no connectivity for new node, n1
        n1_data = n1.problemData;
        j_prev = n2_data(j_key);
        for i = 2:length(path)
           j_prev = ca.ItterativeJNextFromPoint(j_prev, path(i,:), dist);  
        end
        n1_data(j_key) = j_prev;
        n1_data(no_conn_key) = ca.IntegrateJ(n1_data(j_key));
    end
    
end

