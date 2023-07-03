function endVal = maxEndVal(param)
% Compute the maximum end value of a cell array
%
% Define variables:
% param  - cell array of Fourier series (input)
% endVal - maximum end value
% 
% Dependencies: 
% N/A
    endVal = max(max(cellfun(@endVals2, param(1:2,:))));
end % end maxEndVal