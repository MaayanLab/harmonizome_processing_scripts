function sm = smcluster(sm)

if sm.numterms < 2000
    
    try
        
        d = 1 - sm.matrix(tril(ones(size(sm.matrix)), -1) == 1)';

        lnk = linkage(d, 'average');

        idxorder = optimalleaforder(lnk, d, 'Criteria', 'adjacent', 'Transformation', 'linear');

        sm = smorder(sm, idxorder);

        disp('computed optimal leaf order');
        
    catch e
        
        disp(getReport(e));

        figure(666);
        clf;
        [~, ~, idxorder] = dendrogram(lnk,0);
        close(gcf);

        sm = smorder(sm, idxorder);

        disp('computed dendrogram default order');
        
    end
	
elseif sm.numterms < 25000
    
    d = 1 - sm.matrix(tril(ones(size(sm.matrix)), -1) == 1)';
    
    lnk = linkage(d, 'average');
	
	figure(666);
	clf;
	[~, ~, idxorder] = dendrogram(lnk,0);
	close(gcf);

	sm = smorder(sm, idxorder);

	disp('computed dendrogram default order');

else

%     sm.vector = sm.matrix(tril(ones(size(sm.matrix)), -1) == 1);
    
%     sm.scaledmatrix = squarematrixscale(sm.matrix);
    
	disp('matrix too large to order');

end


