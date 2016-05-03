function [stats, fpr, tpr, mcr, fdr, f1s, mcc, val] = classifierperformance(Y, P, Yp, fignum)

% sort scores and hits
[P, si] = sort(P, 'descend');
val = P;
hit = Y(si) == 1;

% basic calculations
population = numel(hit);

actualpositives = sum(hit);
actualnegatives = sum(~hit);

predictedpositives = (1:1:population)';
% predictednegatives = population - predictedpositives;

truepositives = cumsum(hit);
falsepositives = cumsum(~hit);

truenegatives = actualnegatives - falsepositives;
falsenegatives = actualpositives - truepositives;

% advanced calculations
tpr = truepositives/actualpositives; % sensitivity, recall
fpr = falsepositives/actualnegatives; % fall-out

% fnr = falsenegatives/actualpositives; % miss rate
% tnr = truenegatives/actualnegatives; % specificity

% accuracy = (truepositives + truenegatives)/population;
mcr = (falsepositives + falsenegatives)/population;

% positivepredictivevalue = truepositives./predictedpositives;
% falseommisionrate = falsenegatives./predictednegatives;
fdr = falsepositives./predictedpositives;
% negativepredictivevalue = truenegatives./predictednegatives;

% positivelikelihoodratio = tpr./fpr;
% negativelikelihoodratio = fnr./tnr;
% diagnosticoddsratio = positivelikelihoodratio./negativelikelihoodratio;

f1s = 2*truepositives./(2*truepositives + falsepositives + falsenegatives);
mcc = (truepositives.*truenegatives - falsepositives.*falsenegatives)./sqrt((truepositives + falsepositives).*(truepositives + falsenegatives).*(truenegatives + falsepositives).*(truenegatives + falsenegatives));

stats.idx = find(f1s+mcc >= 0.99*max(f1s+mcc), 1, 'first');
stats.val = P(stats.idx);
stats.f1s = f1s(stats.idx);
stats.mcc = mcc(stats.idx);
stats.fdr = fdr(stats.idx);
stats.mcr = mcr(stats.idx);
stats.tpr = tpr(stats.idx);
stats.fpr = fpr(stats.idx);

stats.idx50 = find(fdr <= 0.5, 1, 'last');
if isempty(stats.idx50)
    stats.idx50 = NaN;
    stats.val50 = NaN;
    stats.f1s50 = NaN;
    stats.mcc50 = NaN;
    stats.fdr50 = NaN;
    stats.mcr50 = NaN;
    stats.tpr50 = NaN;
    stats.fpr50 = NaN;
else
    stats.val50 = P(stats.idx50);
    stats.f1s50 = f1s(stats.idx50);
    stats.mcc50 = mcc(stats.idx50);
    stats.fdr50 = fdr(stats.idx50);
    stats.mcr50 = mcr(stats.idx50);
    stats.tpr50 = tpr(stats.idx50);
    stats.fpr50 = fpr(stats.idx50);
end

stats.fdrmin = min(fdr);
stats.idxfdrmin = find(fdr == stats.fdrmin, 1, 'last');
stats.valfdrmin = P(stats.idxfdrmin);

stats.valmax = max(P);

stats.auc = trapz(fpr, tpr);
stats.dev = -2*mean(hit.*log10(P) + (1 - hit).*log10(1 - P));

if ~isempty(Yp)
    stats.MCR = sum(Y~=Yp)/numel(Y);
end

if ~isempty(fignum)
    if population > 1e8
        pltidx = (1:1e5:population)';
    elseif population > 1e7
        pltidx = (1:1e4:population)';
    elseif population > 1e6
        pltidx = (1:1e3:population)';
    elseif population > 1e5
        pltidx = (1:1e2:population)';
    elseif population > 1e4
        pltidx = (1:1e1:population)';
    else
        pltidx = (1:1:population)';
    end
    figure(fignum);
    clf;
    subplot(2,2,1);
    plot(fpr(pltidx), tpr(pltidx), '-k');
    xlabel('fpr');
    ylabel('tpr');
    axis([0 1 0 1]);
    subplot(2,2,2);
    plot(predictedpositives(pltidx)/population, mcr(pltidx), '-k');
    xlabel('predicted positives');
    ylabel('mcr');
    xlim([0 1]);
    subplot(2,2,3);
    plot(predictedpositives(pltidx)/population, fdr(pltidx), '-k');
    xlabel('predicted positives');
    ylabel('fdr');
    xlim([0 1]);
    subplot(2,2,4);
    plot(predictedpositives(pltidx)/population, f1s(pltidx), '-k', predictedpositives(pltidx)/population, mcc(pltidx), '--r');
    xlabel('predicted positives');
    ylabel('f1s (black) and mcc (red)');
    xlim([0 1]);
end


