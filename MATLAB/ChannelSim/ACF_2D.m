%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D ACF generator
%
% usage 
%-----------------------
% ACF = ACF_2D(s, x, y, d, res, perc)
%
% inputs
%--------------------------------
% s : 2D input signal 
% [x , y] : 2D grid of the positions
% d : the desired vector of distances at which the ACF needs to be calculated
% res : the resolution of the grid in (samples/m)
% perc : the percentage of the samples used for averaging in Monte Carlo method
%
% outputs
%--------------------------------
% ACF : the calculated ACF (the same size as input d vector)
%
% Original codes courtesy of Alireza Ghaffarkhah and Alejandro Gonzalez-Ruiz
% Updated by Herbert Cai (April 2016)
%
% If you have any questions or comments regarding the codes,
% please contact Herbert Cai at hcai@ece.ucsb.edu
%
% Dept. of Electrical and Computer Engineering
% University of California, Santa Barbara
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ACF = ACF_2D(s, x, y, d, res, perc)


% vectorize the input matrices
x = x(:);
y = y(:);
s = s(:)/std(s(:));

% choose the random samples used for averaging
M = fix(perc*length(x)/100);
I = randsample(length(x), M);       %indices of the samples we will pick

% calculate the ACF
ACF = zeros(size(d));
for i = 1:length(d)         %go through every possible distance
    
    b = [];
    for k = 1:M                                                 %go through each randomly picked point
        D = sqrt((x(I(k)) - x).^2 + (y(I(k)) - y).^2);          %calculate the distance between the samples and any other point in the workspace
        p = s(I(k)).*s;                                                 
        idx = find( abs(D - d(i)) <= sqrt(2)/(2*res) );         
        
        b = [b; p(idx)];
    end;
    
    ACF(i) = mean(b);
    
    fprintf('%d out of %d\n', i, length(d))

end;