function dPhiMatrix = dPhi(beta, scalars, param, alpha, rho)
% Compute the derivative of the Phi map when the tori has arbitrary period.
%
% Define variables:
% N           - number of Fourier modes
% K           - period of Tori
% n           - total number of coeffients
% rho         - computed rotation number (input)
% alpha       - parameter from mapping
% beta        - parameter from system of equations
% vfill       - vertical zeros
% hfill       - horizontal zerso
% fill        - square zeros
% dPhiCell    - cell array to hold pieces
% rotMatrix   - matrix representing a rotation
% diffMatrix  - matrix representing a derivative
% dPhiMatrix  - assembled derivative
% 
% Dependencies: 
% Fourier.m
% FouerierOperator.m 
    % initialize variables
    N = length(param{1,1});
    K = size(param,2);
    n = 2*N+1;
    vfill = sparse(zeros(n,1));
    hfill = sparse(zeros(1,n));
    fill = sparse(zeros(n));
    dPhiCell = cell(6*K+1,6*K+1);
    % place filler matricies in each spot of the cell array
    for j = 1:size(dPhiCell,1)
        for k = 1:size(dPhiCell,2)
            dPhiCell{j,k} = fill;
        end % end for loop
    end % end for loop
    % Place 1x1 zeros where necessary
    for j = 1:3
        for k = 1:3
            dPhiCell{j,k} = 0;
        end % end for loop
        for k = 1:K-1
            dPhiCell{j,6*k+2} = 0;
            dPhiCell{j,6*k+3} = 0;
            dPhiCell{6*k+2,j} = 0;
            dPhiCell{6*k+3,j} = 0;
        end % end for loop
    end % end for loop
    for j = 1:K-1
        for k = 1:K-1
            dPhiCell{6*j+2,6*k+2} = 0;
            dPhiCell{6*j+3,6*k+2} = 0;
            dPhiCell{6*j+2,6*k+3} = 0;
            dPhiCell{6*j+3,6*k+3} = 0; 
        end % end for loop
    end % end for loop
    % Place columns and rows of zeros to match the filled in zeros
    for k1 = 0:K-1
        for k2 = 0:K-1
            for k = 1:4
                dPhiCell{1,6*k2+k+3} = hfill;
                dPhiCell{6*k2+k+3,1} = vfill;
                for j = 1:2
                    dPhiCell{6*k1+j+1,6*k2+k+3} = hfill;
                    dPhiCell{6*k2+k+3,6*k1+j+1} = vfill;
                end % end for loop
            end % end for loop
        end % end for loop
    end % end for loop
% Create the derivative and rotation matrices
    diffMat = sparse(zeros(n));
    rotMat = sparse(zeros(n));
    for j = -N:N
        rotMat(N+j+1,N+j+1) = exp(2*pi*1i*(K*rho)*j);
        diffMat(N+j+1,N+j+1) = 2*pi*1i*j;
    end
    % Place remaining derivatives for the system
    % Beta wrt beta
%     dPhiCell{1,1} = -1;
    % Beta wrt param{1,1}
    dPhiCell{1,4} = ones(1,n);
    % K_K^1 wrt Beta
    dPhiCell{6*K-2,1} = spMat(rotation(param{1,1},K*rho));
    % K_K^2 wrt Beta
    dPhiCell{6*K-1,1} = spMat(rotation(param{2,1},K*rho));
    for k = 1:K
        % l^1 wrt l^1
%         dPhiCell{6*k-4, 6*k-4} = -1;
%         dPhiCell{6*k-3, 6*k-3} = -1;
       % l^1 wrt S
       dPhiCell{6*k-4,6*k} = sparse(ones(1,n));
       % l^2 wrt C
       dPhiCell{6*k-3,6*k+1} = sparse(ones(1,n));
       % K^1 wrt K^1, K^2, S
       dPhiCell{6*k-2,6*k-2} = - speye(n);
       dPhiCell{6*k-2,6*k-1} = - speye(n);
       dPhiCell{6*k-2,6*k} = - alpha*speye(n);
       % K^2 wrt K^1, S
       dPhiCell{6*k-1,6*k-1} = - speye(n);
       dPhiCell{6*k-1,6*k} = - alpha*speye(n);
       % S wrt l^1 l^2 K^1 S C
       dPhiCell{6*k,6*k-4} = - spMat(param{3,k});
       dPhiCell{6*k,6*k-3} = spMat(param{4,k});
       dPhiCell{6*k,6*k-2} = spMat(FourierOperator(param{4,k})*diffMat);
       dPhiCell{6*k,6*k} = - diffMat - scalars(1,k)*speye(n);
       dPhiCell{6*k,6*k+1} = spMat(FourierOperator(diff(param{1,k}))) + scalars(2,k)*speye(n);
       % C wrt l^1 l^2 K^1 
       dPhiCell{6*k+1,6*k-4} = - spMat(param{4,k});
       dPhiCell{6*k+1,6*k-3} = - spMat(param{3,k});
       dPhiCell{6*k+1,6*k-2} = - spMat(FourierOperator(param{3,k})*diffMat);
       dPhiCell{6*k+1,6*k} = - spMat(FourierOperator(diff(param{1,k}))) - scalars(2,k)*speye(n);
       dPhiCell{6*k+1,6*k+1} = - diffMat - scalars(1,k)*speye(n);
    end % end for loop
    for k = 1:K-1
        % K^1 wrt 
        dPhiCell{6*k-2,6*(k+1)-2} = dPhiCell{6*k-2,6*(k+1)-2} + speye(n);%
        % K^2 wrt 
        dPhiCell{6*k-1,6*(k+1)-1} = dPhiCell{6*k-1,6*(k+1)-1} + speye(n);%
    end % end for loop
    % K^1_K wrt 
    dPhiCell{6*K-2,4} = dPhiCell{6*K-2,4} + (1+beta)*rotMat;
    % K^2_K wrt 
    dPhiCell{6*K-1,5} = dPhiCell{6*K-1,5} + (1+beta)*rotMat;
    % Convert the cell array into a matrix
    dPhiMatrix = cell2mat(dPhiCell);
%     s1 = whos('dPhiMatrix').bytes;
%     fdPM = full(dPhiMatrix);
%     s2 = whos('fdPM').bytes;
%     ratio = s1/s2
    figure(5)
    title('Sparsity Pattern of derivative')
    set(gcf,'Units','Normalized','OuterPosition',[.15 .1 .25 .4]);
    spy(dPhiMatrix);
end % end dPhi