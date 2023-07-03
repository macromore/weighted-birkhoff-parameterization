 function cosseries = fourierCosineApprox(series,N)
% Compute a truncated cosine series using a truncated fourier series
%
% Define variables:
% N           - order of approximation
% series      - Fourier series to be cosined
% 
% Dependencies: 
% Fourier.m
    cosseries = 0;
    for k = 0:N
        cosseries = cosseries + ((-1)^k)*(series^(2*k))/factorial(2*k);
    end % end for loop
 end % end fourierCosineApprox