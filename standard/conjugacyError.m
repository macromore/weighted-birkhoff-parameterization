function conjError = conjugacyError(rho, param, alpha, numberOfPoints)
% Compute the conjugacy error of an estimated solution across the
% cohomology equations.
    K = size(param,2);
    conjError = zeros(numberOfPoints,K);
    figure(1)
    for phase = 1:numberOfPoints
        for k = 1:K
        pMap = [evaluate(param{1,k},phase/numberOfPoints); evaluate(param{2,k},phase/numberOfPoints)];
        pMap = standardMap(pMap, alpha);
%         plot(pMap(1),pMap(2),'*g')
        if k < K
            pConj = [evaluate(param{1,k+1},phase/numberOfPoints); evaluate(param{2,k+1},phase/numberOfPoints)];
        elseif k == K
            pConj = [evaluate(param{1,1},phase/numberOfPoints+K*rho); evaluate(param{2,1},phase/numberOfPoints+K*rho)];
        else
            continue
        end
%         plot(pConj(1),pConj(2),'*b')
        conjError(phase,k) = norm(pMap-pConj);
        end
    end
    conjError = max(conjError,[],'all');
end