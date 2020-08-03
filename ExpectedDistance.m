function exp_dist = ExpectedDistance(n1, n2, ca, method)
%ExpectedDistance - TreeNode X TreeNode -> positive real
% Calculate the expected distance adding n2 to a tree path will add to the
% path. Metric function to be used for RDT

%INPUT
% TreeNode n1 - The node to be added to the graph
% TreeNode n2 - The node already in the graph. 
% ChannelAnalyzer ca - channel analyzer with relevant channel data. Used to
% calculate probabilities.
% method - Flag indicacting which method to use to calculate the
% probability of no connection along the path up to that point.
%               1 - Use Botev's monte carlo based method. Can handle > 25
%               dimensions (path points), and may be more accurate for
%               higher dimensions, but will be slower
%               2 - Use matlab's built in multivarirate CDF function, which
%               will can only handle up to 25 dimensions. Path will be
%               truncated to the 25 most recent points.
%               3 - Assume the path is approximately Markovian, use Arjun's
%               characterization

%OUTPUT
% double exp_dist - expected distance added to the path 

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
    
    %also setup J, prob of no connectivity for new node, n1
    n1_data = n1.problemData;
    j_prev = n2_data(j_key);
    n1_data(j_key) = ca.ItterativeJNextFromPoint(j_prev, n1.getPos(), dist);
    n1_data(no_conn_key) = ca.IntegrateJ(n1_data(j_key));
    
end

