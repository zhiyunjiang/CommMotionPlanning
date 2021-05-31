%test the matrix weight sampler
%%
data_matrix = [1, 0; 0, 1];
mws = MatrixWeightSampler(data_matrix);
weights = mws.weights;
empirical_weights = zeros([4,1]);
for i = 1: 1000
    sub =  mws.sample() + 1;
    ind = sub2ind(size(data_matrix),sub(1), sub(2));
    empirical_weights(ind) = empirical_weights(ind) + 1;
end
(empirical_weights/sum(empirical_weights)) - weights