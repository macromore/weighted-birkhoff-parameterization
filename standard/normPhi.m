function nor = normPhi(beta, scalars, param, alpha, rho, phase, nu)
% Compute the norm of the input for a phi... bad name
% Define variables:
% beta          - beta from the map (input)
% scalars       - scalars for the map (intput)
% param         - cell array containing pairs of Fourier series (input)
% nu            - function space parameter (input)
% nor           - norm of the input (output)
% 
% Dependencies: 
% Fourier.m
    [betatemp, scalarstemp, paramtemp] = ...
    phi(beta, scalars, param, alpha, rho, phase);
    nor = [abs(betatemp), abs(scalarstemp(:)')];
    for k = 1:size(paramtemp,2)
       nor = [nor, l1Norm(paramtemp{1,k},nu), l1Norm(paramtemp{2,k},nu)];
    end % end for loop
    nor = max(nor);
end % end normPhi