%% RunBayesianThresholding.m
% This script runs Bayesian thresholing on simulated 3D trajectories and
% compares it to simple thresholding. The trajectories are simulated from a
% Rouse model, but conditionally on a defined proximity interval. That
% means, the true distance is larger than the proximity threshold alpha
% everywhere except inside the considered proximity window, where it is
% smaller than alpha. We then use Bayesian and simple thresholding to
% analyze how the inferred proximity durations compare to the true one. The
% script uses parallelization to simulate 1000 trajectories, which should
% finish within a few minutes on a laptop / PC.

clear;
close all;

%rng(100);

addpath('mvrandn');

cols = GetDefaultColors;

dim = 3;

% Parameters for the Rouse dynamics.
deltaS = 16;
D = 291;
gamma = 1;
kappa = 1/17;
T = 60*60;
dt = 0.2*60;
tGrid = 0:dt:T;
K = length(tGrid);

% Localization error.
sigma_loc = sqrt(2)*20;

% Creates the Rouse model and associated statistics.
rm = CreateRouseModel(tGrid, deltaS, D, gamma, kappa, sigma_loc);

% Defines the vector of proximity intervals in frames.
winSizeVec = [2:2:20];

% Select which proximity intervals should be plotted as examples.
selWin = [10, 20];

% Define proximity threshold
alpha = 50;
M = 1000;

for u=1:length(winSizeVec)

    winSize = winSizeVec(u);


    % Simulates K trajectories with defined proximity duration and
    % calculate posterior distances using Bayesian inference (muPost). X is
    % the true distance and Y is the distance including localization error.
    [X, Y, muPost, SigmaPost] = SimulateProximityCond(tGrid, rm, winSize, alpha, M, dim);

    XMat = reshape(X, [3*M, K]);
    SigmaXEst = cov([X(:, :, 1); X(:, :, 2); X(:, :, 3)]);

    % Ground truth (must be the same as the prescribed contact interval by
    % construction).
    dX = max(abs(X), [], 3);
    dX = sqrt(X(:, :, 1).^2 + X(:, :, 2).^2 + X(:, :, 3).^2);

    % Thresholded distance Y.
    dY = max(abs(Y), [], 3);
    dY = sqrt(Y(:, :, 1).^2 + Y(:, :, 2).^2 + Y(:, :, 3).^2);

    % Thresholded posterior distance muPost.
    dmu = max(abs(muPost), [], 3);
    dmu = sqrt(muPost(:, :, 1).^2 + muPost(:, :, 2).^2 + muPost(:, :, 3).^2);

    %pProximity = chi2cdf(alpha^2, )
    sigma2Post = repmat(diag(SigmaPost)', M, 1);
    lambda = dmu.^2 ./ sigma2Post;
    prob = ncx2cdf(alpha^2 ./ sigma2Post, dim, lambda);

    % Extracts inferred proximity intervals.
    [timeI_True, EventTimesTotal_True, ~] = DetectEventTimes(dX<alpha);
    [timeI_Bayes, ~, EventTimesTotal_Bayes] = DetectEventTimes(dmu<alpha);
    [timeI_Threshold, ~, EventTimesTotal_Threshold] = DetectEventTimes(dY<alpha);

    % Calculate summary statistics
    dProbBayes(u) =  mean(sum(and(dX<alpha, dmu<alpha), 2))/winSize;
    meanTau_Bayes(u) = nanmean(EventTimesTotal_Bayes);
    quantiles_Bayes([1,2], u) = 2.5*nanstd(EventTimesTotal_Bayes)/sqrt(length(EventTimesTotal_Bayes));
    dProbThreshold(u) =  mean(sum(and(dX<alpha, dY<alpha), 2))/winSize;
    meanTau_Threshold(u) = nanmean(EventTimesTotal_Threshold);
    quantiles_Threshold([1, 2], u) = 2.5*nanstd(EventTimesTotal_Threshold)/sqrt(length(EventTimesTotal_Threshold));

    winSizeIdx = find(winSize == selWin);

    % Plotting
    if (~isempty(winSizeIdx))
        figure(1)
        subplot(3,length(selWin),winSizeIdx);
        dt = tGrid(2) - tGrid(1);
        plot(tGrid/60, dX(1, :), '-', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5); hold on;
        plot(tGrid/60, dY(1, :), '.', 'Color', cols(2, :));
        plot(tGrid/60, dmu(1, :), '--', 'Color', cols(1, :), 'LineWidth', 1.5);
        plot(tGrid/60, alpha*ones(size(tGrid)), '--', 'Color',[0.2, 0.2, 0.2], 'LineWidth', 1.5);
        ylim([0, 500]);
        ylabel('3D Distance');
        subplot(3,length(selWin),1*length(selWin)+winSizeIdx);
        bar((tGrid+dt/2)/60, dY(1, :)<alpha, 'BarWidth', 1, 'FaceColor', cols(2, :), 'LineStyle','none'); hold on;
        stairs(tGrid/60, dX(1, :)<alpha, 'LineWidth', 1.5, 'Color', [0.5, 0.5, 0.5]);
        legend('Inferred Interval', 'True Interval');
        ylabel('Inferred proximity (Threshold)');
        ylim([0, 1.5]);
        subplot(3,length(selWin),2*length(selWin)+winSizeIdx);
        bar((tGrid+dt/2)/60, dmu(1, :)<alpha, 'BarWidth', 1, 'FaceColor', cols(1, :), 'LineStyle','none'); hold on;
        stairs(tGrid/60, dX(1, :)<alpha, 'LineWidth', 1.5, 'Color', [0.5, 0.5, 0.5]);
        ylabel('Inferred proximity (Bayes)');
        ylim([0, 1.5]);
        xlabel('Time (min)');
        legend('Inferred Interval', 'True Interval');

    end

end


dim = 3;

r = alpha / sigma_loc;
% Calculate the theoretical detection probability for an ideal contact.
pEvent = ncx2cdf(r.^2, dim, 0);

% Calculate the mean inferred lifetime for an ideal contact.
tau_mri_L2 = 1./(1 - pEvent);
relDurationLimit = 1 - abs(tau_mri_L2 - winSizeVec)./winSizeVec;
dProbLimit = repmat(pEvent, 1, length(winSizeVec));

% Plot summary statistics / results
figure;
subplot(1,2,1);
plot(winSizeVec*dt/60, meanTau_Bayes*dt/60, '.-', 'Color', cols(1, :), 'LineWidth', 1.5, 'MarkerSize', 12); hold on;
plot(winSizeVec*dt/60, meanTau_Threshold*dt/60, '.-', 'Color', cols(2, :), 'LineWidth', 1.5, 'MarkerSize', 12); hold on;
plot([0, winSize*dt/60], [0, winSize*dt/60], '--', 'LineWidth', 1.5, 'Color', [0.5, 0.5, 0.5]); hold on;
plot([0, max(winSizeVec)*dt]/60, [1, 1]*tau_mri_L2*dt/60, 'k--', 'LineWidth',1.5);
xlim([0, max(winSizeVec)*dt]/60);
ylim([0, max(winSizeVec)*dt]/60);
xlabel('True contact lifetime [min]');
ylabel('Inferred contact lifetime [min]');

subplot(1,2,2);
relDurationBayes = 1-abs(meanTau_Bayes-winSizeVec)./winSizeVec;
relDurationThreshold = 1-abs(meanTau_Threshold-winSizeVec)./winSizeVec;
plot(dProbBayes, relDurationBayes, '.-', 'Color', cols(1, :), 'LineWidth', 1.5); hold on;
plot(dProbThreshold, relDurationThreshold, '-', 'Color', cols(2, :), 'LineWidth', 1.5); hold on;

for k=1:length(dProbBayes)
    plot(dProbBayes(k), relDurationBayes(k), '.', 'MarkerSize', 8 + 12*(winSizeVec(k))/30, 'Color', cols(1, :));
    plot(dProbThreshold(k), relDurationThreshold(k), '.', 'MarkerSize', 8 + 12*(winSizeVec(k))/30, 'Color', cols(2, :));
end

xlabel('True positive rate');
ylabel('Interval agreement');
