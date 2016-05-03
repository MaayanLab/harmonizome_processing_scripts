function cm = cm_ksdensity_standardization_sparsematrix(cm, support, fignum)

if exist('support', 'var') == 0 || isempty(support)
    support = 'unbounded';
end
        
X = cm.matrix(cm.matrix ~= 0);
if numel(X) <= 3000
    P = ksdensity(X, X, 'function', 'cdf', 'support', support);
else
    [p, x] = ksdensity(X, 'npoints', 3000, 'function', 'cdf', 'support', support);
    P = interp1(x, p, X, 'spline');
end
cm.matrix(cm.matrix ~= 0) = P;

if exist('fignum', 'var') == 1 || ~isempty(fignum)
    [f, x] = ecdf(X, 'function', 'cdf');
    [xfit, ui] = unique(X);
    ffit = P(ui);
    figure(fignum);
    clf;
    plot(x, f, 'ok', xfit, ffit, '-r', 'linewidth', 2);
    ylim([0 1]);
    legend('empirical', 'fit');
end


