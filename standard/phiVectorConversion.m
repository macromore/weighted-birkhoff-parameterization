function phi = phiVectorConversion(beta, scalars, param)
% Convert the cell arrays to a vector
%
% Define variables:
% beta          - beta from the map (input)
% scalars       - scalars for the map (intput)
% param         - cell array containing pairs of Fourier series (input)
% phi           - vector of input (output)
% K             - number of modes
% 
% Dependencies: 
% Fourier.m
    phi = beta;
    K = size(param,2);
    for k = 1:K
       phi = [phi; scalars(:,k); mat(param{1,k}); mat(param{2,k});...
           mat(param{3,k}); mat(param{4,k})]; 
    end % end for loop
end % end phiVectorConversion