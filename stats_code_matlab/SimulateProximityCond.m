function [CondX, CondY, muPost, SigmaPost] = SimulateProximityCond(tGrid, rm, winSize, alpha, N, dim)

M = 10000;
sigmaTot = sqrt(rm.acfRouse(0)) + rm.sigma_loc;
xGrid = linspace(-2*sigmaTot, 2*sigmaTot, M);

thrFun = @(x) abs(x)<alpha;

K = length(tGrid);
dx = xGrid(2) - xGrid(1);

sigma = sqrt(rm.acfRouse(0));




normTrunc = ~thrFun(xGrid).*normpdf(xGrid, 0, sigma);
normTrunc = normTrunc / sum(normTrunc*dx);

% Set statistics of the propagator
CondX = zeros(N, K, dim);
CondY = zeros(N, K, dim);
muPost = zeros(N, K, dim);

alphaMargin = 0.1;

parfor n=1:N

    contactInterval = zeros(size(tGrid));

    kStart = floor(K/2-winSize/2);
    kStart = randi([2, K-winSize-1]);
    kEnd = kStart + winSize-1;

    contactInterval(kStart:kEnd)=1;

    contactIdx = find(contactInterval);
    idx = [1, contactIdx, K];

    SigmaCond = rm.Sigma(idx, idx);
    idxComp = setdiff(1:K, idx);
    SigmaX = rm.Sigma(idxComp, idxComp);
    SigmaXY1 = rm.Sigma(idxComp, idx);
    SigmaXY2 = rm.Sigma(idx, idxComp);
    covRatio = SigmaXY1/SigmaCond;
    SigmaXCond = SigmaX - covRatio*SigmaXY2;

    for u=1:dim

        XEnd = randsample(xGrid, 1, true, normTrunc);
        XStart = randsample(xGrid, 1, true, normTrunc);

        lb = -alpha*ones(size(contactIdx))*alphaMargin;
        ub = -lb;

        yContact = mvrandn(lb, ub, rm.Sigma(contactIdx, contactIdx), 1);
        yCond = [XStart; yContact; XEnd];

        muXCond = covRatio*yCond;

        lb = zeros(1, length(muXCond));
        ub = zeros(1, length(muXCond));
        idxStart = idxComp<contactIdx(1);
        if (XStart<0)
            lb(idxStart) = -inf;
            ub(idxStart) = -alpha*alphaMargin;
        else
            lb(idxStart) = alpha*alphaMargin;
            ub(idxStart) = inf;
        end

        idxEnd = idxComp>contactIdx(end);
        if (XEnd<0)
            lb(idxEnd) = -inf;
            ub(idxEnd) = -alpha*alphaMargin;
        else
            lb(idxEnd) = alpha*alphaMargin;
            ub(idxEnd) = inf;
        end

        lb = lb - muXCond';
        ub = ub - muXCond';

        XCond = mvrandn(lb, ub, SigmaXCond, 1) + muXCond;

        X = zeros(1, K);
        X(idx) = yCond;
        X(idxComp) = XCond;

        CondX(n, :, u) = X;
        CondY(n, :, u) = CondX(n, :, u) + rm.sigma_loc*randn(1, K);

      
    end

    if (mod(n, 50) == 1)
        fprintf('Finished realization %d\n', n)
    end
end

% Calculate posterior covariance
SigmaXEst = cov([CondX(:, :, 1); CondX(:, :, 2); CondX(:, :, 3)]);
SigmaPost = inv(inv(SigmaXEst) + 1/rm.sigma_loc^2 * eye(K));

for n=1:N

    for u=1:dim
          % Calculate posterior mean
        muPost(n, :, u) = (1/rm.sigma_loc^2*CondY(n, :, u)) * SigmaPost;
    end

end