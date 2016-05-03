function P = kscumulativeprobability_matrix_ignorezeros(X, support, fignum)

if exist('support', 'var') == 0 || isempty(support)
    support = 'unbounded';
end
        
Xnz = X(X ~= 0);

if numel(Xnz) <= 3000
    Pnz = ksdensity(Xnz, Xnz, 'function', 'cdf', 'support', support);
else
    [p, x] = ksdensity(Xnz, 'npoints', 3000, 'function', 'cdf', 'support', support);
    Pnz = interp1(x, p, Xnz);
end

P = zeros(size(X));
P(X ~= 0) = Pnz;

if exist('fignum', 'var') == 1 || ~isempty(fignum)
    [f, x] = ecdf(Xnz, 'function', 'cdf');
    [xfit, ui] = unique(Xnz);
    ffit = Pnz(ui);
    figure(fignum);
    clf;
    plot(x, f, 'ok', xfit, ffit, '-.r');
    ylim([0 1]);
    legend('empirical', 'fit');
end


