function sobolevNorm = sobolevNorm(param, sobolevParam)
% Compute the sobolev norm of a set of parameterizations
%
% Define variables:
% param       - cell array of Fourier parameterizations (input)
% s           - H^s (input)
% sobolevNorm - sobolev norm array (output)
% x           - temp variable for matrices
% K           - number of Fourier modes
% sum         - temp variable to hold norm
% c           - holds both values for the parameterizations
% 
% Dependencies: 
% Fourier.m
    if nargin == 1
        sobolevParam = 1;
    end % end if
    sobolevNorm = zeros(size(param,2),1);
    for k = 1:size(param,2)
        x = [mat(param{1,k})';mat(param{2,k})'];
        K = length(param{1,k});
        sum = 0;
        for j = -K:K
            c = x(:,j+K+1);
            sum = sum + (1+abs(j)^2)^sobolevParam*max(abs(c))^2;
        end % end for loop
        sobolevNorm(k) = sqrt(sum);
    end % end for loop
end % end sobolevNorm

