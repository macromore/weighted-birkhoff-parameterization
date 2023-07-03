function [ P_n ] = standardMap( P, alpha )
% Compute the image under the standard map
%
% Define variables:
% P         - base point (input)
% alpha     - standard map parameter (input)
% P_n       - image (output)
% 
% Dependencies: 
% fourierSineApprox.m
    %Standard Map of a point P with parameter alpha
    if nargin == 1
        alpha = 1/4;
    end % end if
    if isa(P, 'cell')
        P_n = cell(2,1);
        P_n{2} = P{2} + alpha*fourierSineApprox(P{1}, 50);
        P_n{1} = P{1} + P_n{2};
    else % not a cell array
        P_n = zeros(2,1);
        P_n(2) = P(2) + alpha*sin(P(1));%*fourierSineApprox(P(1),50);
        P_n(1) = P(1) + P_n(2);
%         P_n(1) = real(P_n(1)) + 1i*imag(P_n(1));
        P_n(1) = mod(real(P_n(1)), 2*pi) + 1i*imag(P_n(1));
    end % end if
end % end standardMap