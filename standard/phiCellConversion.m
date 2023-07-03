function [beta, scalars, param] = phiCellConversion(vect, N)
% Convert the matrix to a cell array
%
% Define variables:
% vect          - vector of phi map (input)
% N             - number of circles (input)
% beta          - beta from the map (output)
% scalars       - scalars for the map (output)
% param         - cell array containing pairs of Fourier series (output)
% vec           - tail of the vect without beta
% K             - number of modes in Fourier series
% vectemp       - temp storange for vec
% 
% Dependencies: 
% N/A
    beta = vect(1);
    vec = vect(2:end);
    K = length(vec)/(2+4*(2*N+1));
    vec = reshape(vec,2+4*(2*N+1),K);
    param = cell(4,K);
    for k = 1:K
        vectemp = reshape(vec(3:end,k),2*N+1,4);
        param{1,k} = Fourier(vectemp(:,1));
        param{2,k} = Fourier(vectemp(:,2));
        param{3,k} = Fourier(vectemp(:,3));
        param{4,k} = Fourier(vectemp(:,4));
    end % end for loop
    scalars = vec(1:2,:);
end % end phiCellConversion