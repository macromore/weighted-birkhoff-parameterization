function sinseries = fourierSineApprox(series, N)
% Compute a truncated cosine series using a truncated fourier series
%
% Define variables:
% N           - order of approximation
% series      - Fourier series to be sined
% 
% Dependencies: 
% Fourier.m
    sinseries = 0;
    for k = 0:N
       sinseries = sinseries + (-1)^k*series^(2*k+1)/factorial(2*k+1); 
    end % end for loop
end % fourierSineApprox