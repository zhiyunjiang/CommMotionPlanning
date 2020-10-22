%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IntegrateWSimpsonsIrregular
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uses Simpsons method to estimate an integral with function samples f
% irregularlay spaced
%
% Input:
% f - array of n+1 function values,
% h - array of n interval lengths
% 
% Output: 
% result - the aproximate integral value
 
function result = IntegrateWSimpsonsIrregular(f, h)
    %number intervals
    n = length(h);

    w = SimpsonWeights(n, h);
    
    %now just sum those bad bois
    result = w*f;
end

