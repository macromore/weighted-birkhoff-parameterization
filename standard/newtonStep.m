function [betanext, scalarsnext, paramnext] = newtonStep...
        (beta, scalars, param, alpha, rho, phase)
% Preform a newton step on phi 
%
% Define variables:
% beta          - beta from previous step (input)
% scalars       - scalars for map from previous step (input)
% param         - previous parameterization (input)
% alpha         - standard map parameter (input)
% rho           - approximate rotation number (input)
% phase         - phase values for map (input)
% betanext      - new beta (output)
% scalarsnext   - new scalars (output)
% paramnext     - new parameterization (output)
% x             - vector representation of phi 
% xnext         - new vector repesentation of phi 
% N             - number of modes
% 
% Dependencies: 
% Fourier.m
% phiVectorConversion.m
% phi.m
% phiCellConversion
    x = phiVectorConversion(beta, scalars, param);
    [betatemp, scalarstemp, paramtemp] = ...
        phi(beta, scalars, param, alpha, rho, phase);
    xnext = x - dPhi(beta, scalars, param, alpha, rho)\...
          phiVectorConversion(betatemp, scalarstemp, paramtemp);
    N = length(param{1,1});
    [betanext, scalarsnext, paramnext] = phiCellConversion(xnext, N);
end % end newtonStep