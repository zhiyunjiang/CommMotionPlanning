classdef MatrixWeightSampler < handle
    %MatrixWeightSampler 
    properties
        weights;
        vals;
        sz;
    end
    
    methods
        function this = MatrixWeightSampler(data_matrix)
            this.sz = size(data_matrix);
            
            min_data = min(data_matrix, [], 'all');
            max_data = max(data_matrix, [], 'all');
            
            %avoid 0 probability on min
            eps = (max_data - min_data)/sum(this.sz);
            weight_matrix = data_matrix - min_data + eps;
            this.weights = reshape(weight_matrix/(sum(weight_matrix, 'all')), this.sz(1)*this.sz(2), 1);
            this.vals = 1:length(this.weights);
        end
        
        function pos = sample(this)
            idx = randsample(this.vals, 1, 1, this.weights);
            [x, y] = ind2sub(this.sz',idx);
            pos = [x,y] - 1;%subtract 1 to keep this 0 indexed
        end
    end
end

