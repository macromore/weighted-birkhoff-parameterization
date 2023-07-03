function conjError = conjugacyError(rho, param, alpha, numberOfPoints)
% Compute the conjugacy error of an estimated solution across the
% cohomology equations.
    K = size(param,2);
    conjError = zeros(K,numberOfPoints-1);
    for j = 1:K
        for p = 1:numberOfPoints-1
            phase = p/(numberOfPoints-1);
            pMap = [evaluate(param{1,j},phase); evaluate(param{2,j},phase)];
            for k = 1:K
                pMap = standardMap(pMap, alpha);
            end
            pConj = [evaluate(param{1,j},K*rho+phase); evaluate(param{2,j},K*rho+phase)];
            conjError(j,p) = norm(pMap-pConj);
        end
    end
    conjError = max(conjError,[],'all');
end