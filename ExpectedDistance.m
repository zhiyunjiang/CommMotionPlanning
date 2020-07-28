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

    %Retrieve the path all the way back to root
    path = n2.pathToRoot(0);
    %find the probability of no connection thusfar
    p_no_conn_2_here = ca.NoConnectionOnPath(path, method);
    
    %find the conditional probability
    p_no_conn_cond = p_no_conn_2_here/ca.NoConnectionPrior(path(1,:));
   
    %use the probability to scale the distance
    exp_dist = p_no_conn_cond*n1.distTo(n2);
end

